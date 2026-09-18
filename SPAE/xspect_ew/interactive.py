"""
Interactive matplotlib widget for xspect-ew: step through Load -> Normalize
-> RV Shift -> Measure EW, with live-adjustable parameters at each stage.

Same architecture as moog/interactive.py's SynthWidget -- pure
matplotlib.widgets (Slider/Button/TextBox), no separate GUI framework --
so the exact same code opens as a standalone window from a terminal
script, or embeds in a Jupyter cell via `%matplotlib widget`.

Usage (terminal), single exposure:
    from SPAE.xspect_ew.interactive import ew_interactive
    ew_interactive('spectrum.fits', linelist='linelist_star.txt')

Usage, multiple exposures of the same star (Stage 1 lets you inspect each
one individually or overlaid, and combine them, before continuing):
    ew_interactive(spectra=['exp1.fits', 'exp2.fits'], linelist='linelist_star.txt')

Usage (Jupyter notebook):
    %matplotlib widget
    from SPAE.xspect_ew.interactive import ew_interactive
    ew_interactive('spectrum.fits', linelist='linelist_star.txt')

Every stage calls into the SAME Spectrum_Data methods an automated script
would call (normalize()/normalize_all(), apply_rv_shift(), measure_ew()) --
this widget is a thin, optional interactive layer on top of that one
shared API, not a second code path, so it can never drift out of sync
with the automated end-to-end mode. In particular, Stage 4 (EW
measurement) treats measure_ew() as a black box: it only relies on the
axes= parameter added alongside this widget (see spectrum_data.py), not
on any of that routine's internal fitting logic -- which is expected to
keep changing as the continuum-fitting-within-EW-measurement work
continues on its own branch.
"""

import copy
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox

from .spectrum_data import Spectrum_Data

_STAGES = ['Load', 'Normalization', 'RV Shift', 'Measure EW']
# Reject Back/Next/Combine clicks that land within this long after any
# stage rebuild -- see _on_nav_clicked()/_on_combine_clicked()'s
# docstring note for why this exists (it's a different failure mode from
# the self._busy reentrancy guard, and neither alone is sufficient).
_TRANSITION_DEBOUNCE_SECONDS = 0.5
_DEFAULT_BTN_COLOR = '0.85'   # matplotlib Button's own default
_COMBINE_COLOR = '#c9e8b0'
_BUSY_COLOR = '#ffd966'
# 'normalized' reuses _COMBINE_COLOR -- both mark "showing the processed
# result", same idea as Combine Spectra -> the Combined view.
_NORM_VIEW_COLORS = {'raw': ('#a8d0e8', '#7fb8dc'), 'normalized': (_COMBINE_COLOR, '#a8d98a')}
# Same two colors, same reasoning, for the RV stage's Shifted/Observed
# toggle -- Shifted (rest-frame, the corrected result) = green like
# Normalized; Observed (uncorrected) = the same blue as Raw+Fit.
_RV_VIEW_COLORS = {True: _NORM_VIEW_COLORS['normalized'], False: _NORM_VIEW_COLORS['raw']}


def _read_instrument_label(filename):
    """Best-effort 'Telescope / Instrument' string from a FITS primary
    header, for the window's title bar -- purely cosmetic, so any failure
    (not a real path, not FITS, missing keywords) just yields None rather
    than raising."""
    try:
        from astropy.io import fits
        with fits.open(filename) as hdul:
            header = hdul[0].header
        parts = [p for p in (header.get('TELESCOP'), header.get('INSTRUME')) if p]
        return ' / '.join(parts) if parts else None
    except Exception:
        return None


class EWWidget:
    """
    Interactive matplotlib GUI for stepping through the xspect-ew pipeline.

    Keep a reference to the returned object -- garbage-collecting it will
    destroy the figure's callbacks.
    """

    def __init__(self, filename=None, spectrum=None, spectra=None, linelist=None,
                 KECK_file=True, **reader_kwargs):
        """
        filename/spectrum : a single exposure (path, or an already-
            constructed Spectrum_Data) -- the common single-exposure case.
        spectra : a list of exposures of the SAME star (paths and/or
            already-constructed Spectrum_Data, freely mixed) when more
            than one is available. The FIRST is the "primary" spectrum
            (self.spec) -- it receives the linelist and, after Stage 1's
            Combine step, the co-added result. Stage 1 (Load) then lets
            you inspect each exposure individually or overlaid, and
            combine them once you've visually confirmed they line up,
            before continuing to Normalize/RV/EW.
        """
        def _as_spectrum(x):
            if isinstance(x, Spectrum_Data):
                return x
            return Spectrum_Data(x, KECK_file=KECK_file, **reader_kwargs)

        if spectra is not None:
            self._exposures = [_as_spectrum(x) for x in spectra]
        elif spectrum is not None:
            self._exposures = [_as_spectrum(spectrum)]
        elif filename is not None:
            self._exposures = [_as_spectrum(filename)]
        else:
            raise ValueError("Provide filename=/spectrum= (a single exposure), "
                              "or spectra=[...] (two or more exposures of the same star).")

        self.spec = self._exposures[0]
        self._n_exposures = len(self._exposures)
        self._instrument_label = _read_instrument_label(self.spec.filename)
        # Raw-flux snapshot per exposure, captured now (before any
        # normalization/combining touches self.spec.flux in place) so
        # Stage 1 can still show each exposure's own original data even
        # after combining has overwritten self.spec.flux with the
        # co-added result.
        self._raw_flux_snapshots = [copy.deepcopy(e.flux) for e in self._exposures]
        self._exposure_names = [os.path.basename(str(e.filename)) for e in self._exposures]
        # Computed unconditionally (not just when > 1 exposure) so the
        # single-exposure Load view can also color its line to match --
        # see _redraw_load().
        cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
        self._exposure_colors = [cycle[i % len(cycle)] for i in range(self._n_exposures)]
        self._combined = False
        # 'overlay' (only offered when there's >1 exposure) or an int
        # index into self._exposures.
        self._view_idx = 0 if self._n_exposures == 1 else 'overlay'

        self._linelist_path = linelist
        if linelist is not None:
            self.spec.load_lines(linelist)

        self._stage = 0
        self._order = 0
        self._line_idx = 0

        self._normalized = False
        self._rv_applied = False
        self._show_shifted = True
        self._norm_view = 'raw'  # or 'normalized' -- see _toggle_norm_view()
        # order -> (lam, p) actually used for THAT order's own pred_all/
        # normalized_flux -- since the lam/p sliders are shared, single
        # widget state (not per-order), stepping to a different order
        # would otherwise leave them showing whatever was last touched
        # for a DIFFERENT order, misrepresenting what actually produced
        # the fit currently on screen. Populated in _build_stage_normalize()
        # (normalize_all()/Apply to ALL orders) and _refit_current_order()
        # (a single order's own slider tweak); consulted by
        # _on_normalize_order_change() whenever the displayed order changes.
        self._order_norm_params = {}

        # True once the user has typed explicit axis limits into the zoom
        # strip -- while set, redraw functions restore those limits after
        # their ax.clear()+replot instead of falling back to autoscale
        # (mirrors moog/interactive.py's SynthWidget._zoom_locked). Reset
        # on every stage switch and on every order/line step, since a
        # different order/line has a different natural wavelength range.
        # Also set to True by _on_axes_view_changed() -- the matplotlib
        # navigation toolbar's own rectangle-zoom/pan tools change the
        # axes' xlim/ylim directly, bypassing the zoom strip entirely, so
        # without that hook a toolbar-drawn zoom silently reset itself on
        # the very next order-step or slider tweak.
        self._zoom_locked = False
        # True for the duration of each redraw function's own
        # ax.clear()+replot+_apply_zoom() sequence -- guards
        # _on_axes_view_changed() against mistaking OUR OWN programmatic
        # set_xlim()/autoscale() calls for a toolbar-driven user zoom.
        self._redrawing = False
        # True while a Back/Next/Combine/Apply-to-ALL click is being
        # handled. _flash_busy() calls canvas.flush_events() to force the
        # amber color to actually paint before the slow work underneath
        # it runs -- but flush_events() can ALSO reentrantly dispatch a
        # second click that arrived while the first was still being
        # handled (e.g. an impatient re-click during a multi-second
        # combine_spectra()/normalize_all() call), running the whole
        # handler a second time NESTED inside the first. That nested
        # call's own _goto_stage() completes before the outer call's own
        # _goto_stage() even runs. This flag makes any click that arrives
        # while one is already being handled (nested, on the same call
        # stack) a safe no-op instead.
        self._busy = False
        # Timestamp of the last _goto_stage() rebuild -- see
        # _TRANSITION_DEBOUNCE_SECONDS. Covers a DIFFERENT failure mode
        # than self._busy above: Combine Spectra and Next share the same
        # screen position (by design, so Next visually "takes over" once
        # combining is done), so an impatient real double-click -- first
        # click starts Combine, second click lands moments later, once
        # the rebuild has ALREADY made Next appear in that exact spot --
        # advances straight past the Combined view to Normalization
        # entirely correctly, from Next's own perspective. self._busy
        # alone doesn't catch this: the first click's handler has
        # already fully returned (busy back to False) by the time the
        # second, separate click event is dispatched.
        self._last_transition_time = 0.0

        # Per-order normalization knobs -- live-adjustable, seeded with
        # normalize()/normalize_all()'s own defaults.
        self._lam = 2e3
        self._p = 0.01

        # EW-measurement ex_params, see measure_ew()'s own docstring:
        # [continuum shift, left bound (A), right bound (A), center (A)]
        self._ex_params = [0.0, 0.0, 0.0, 0.0]

        self._main_axes = []       # main-area axes, rebuilt fresh per stage
        self._dynamic_axes = []    # stage-specific control axes
        self._dynamic_widgets = []  # Button/Slider/TextBox living on those axes
        self._dynamic_texts = []   # fig.text() labels living alongside them

        self._build_figure()
        self._goto_stage(0)

    # ------------------------------------------------------------------ #
    # Figure scaffold (persistent across stages)                          #
    # ------------------------------------------------------------------ #

    def _build_figure(self):
        self._fig = plt.figure('xspect-ew interactive', figsize=(13, 8))
        self._fig.subplots_adjust(left=0, bottom=0, right=1, top=1)

        title = 'XSpect-EW' + (f': {self._instrument_label}' if self._instrument_label else '')
        # Left-aligned, not centered -- a long instrument string centered
        # on the figure crowds into the Combine Spectra button on the right.
        self._stage_text = self._fig.text(
            0.08, 0.965, title, fontsize=13, ha='left', va='center', weight='bold')

        # Widened from the original generic "Back"/"Next" to fit the
        # stage-name labels _goto_stage() fills in (e.g. "◀ Load",
        # "Normalization ▶") -- see there for why the text isn't fixed here.
        # y=0.895, not 0.93: the header title's own rendered text (measured
        # via get_window_extent(), not just its anchor point) extends down
        # to y~0.954, and a button row at 0.93 genuinely overlapped it.
        ax_back = self._fig.add_axes([0.08, 0.895, 0.16, 0.04])
        ax_next = self._fig.add_axes([0.82, 0.895, 0.14, 0.04])
        self._btn_back = Button(ax_back, '← Back')
        # Green like Combine Spectra -- all forward-moving navigation
        # (Normalization/RV Shift/Measure EW ▶) reads as one consistent
        # "go" action; Back stays the neutral default to read as distinct.
        self._btn_next = Button(ax_next, 'Next →', color=_COMBINE_COLOR, hovercolor='#a8d98a')
        self._btn_back.on_clicked(lambda _e: self._on_nav_clicked(self._btn_back, self._stage - 1))
        self._btn_next.on_clicked(lambda _e: self._on_nav_clicked(self._btn_next, self._stage + 1))

        # Combine Spectra lives here (persistent, shown/hidden per stage in
        # _goto_stage()) rather than down among Stage 1's own per-stage
        # controls -- it's the one Load-stage action important enough to
        # want next to the main Back/Next navigation. Same x AND width as
        # Next -- the two are never visible simultaneously (Combine hides
        # once combined, Next only appears once combined), so this makes
        # Next take over the exact same slot Combine just vacated.
        ax_combine = self._fig.add_axes([0.82, 0.895, 0.14, 0.04])
        self._btn_combine = Button(ax_combine, 'Combine Spectra',
                                    color=_COMBINE_COLOR, hovercolor='#a8d98a')
        self._btn_combine.on_clicked(lambda _e: self._on_combine_clicked())

        # Tab-cycling and double-click-to-clear for the Order/zoom TextBoxes
        # -- populated fresh per stage (see _goto_stage()) by _order_pager()
        # and _build_zoom_strip(), in the order they should tab through.
        self._tab_order = []
        self._fig.canvas.mpl_connect('key_press_event', self._on_tab_key)
        self._fig.canvas.mpl_connect('button_press_event', self._on_textbox_dblclick)

        # Main plotting region reserved here; individual stages create
        # whatever axes they need inside it (one wide axes, or two
        # side-by-side, etc.) -- see _clear_main()/_main_rect(). Bottom is
        # 0.26, not right at the top of the control rows (0.135-0.175),
        # to leave room for the main axes' own x-axis tick/label text
        # (otherwise it overlaps the order/line-pager buttons directly
        # below it). Top is 0.85, not 0.88 -- the axes' own title (loc=
        # 'left') renders ABOVE the axes top, and at 0.88 it collided with
        # the Back/Next/Combine row at y=0.895 (confirmed via
        # get_window_extent() on the actual title text, not just the
        # axes patch -- same class of issue as the header/xlabel ones).
        self._main_rect = (0.08, 0.26, 0.88, 0.59)

    def _main_rect_split(self, n):
        """n side-by-side axes rects spanning the main region."""
        left, bottom, width, height = self._main_rect
        gap = 0.03
        w = (width - gap * (n - 1)) / n
        return [(left + i * (w + gap), bottom, w, height) for i in range(n)]

    def _clear_main(self):
        for ax in self._main_axes:
            ax.remove()
        self._main_axes = []

    def _clear_dynamic(self):
        # Button/Slider/TextBox stay connected to the figure's canvas-wide
        # mouse-motion/click events even after their axes is removed --
        # matplotlib.widgets' own hover-highlight handlers (e.g. Button.
        # _motion) then crash on the next mouse move over that dead screen
        # location (self.ax.figure is None -> AttributeError). Disconnect
        # each widget's own event handlers BEFORE removing its axes.
        for w in self._dynamic_widgets:
            w.disconnect_events()
        self._dynamic_widgets = []
        for t in self._dynamic_texts:
            t.remove()
        self._dynamic_texts = []
        for ax in self._dynamic_axes:
            ax.remove()
        self._dynamic_axes = []

    def _add_main_axes(self, n=1):
        rects = self._main_rect_split(n) if n > 1 else [self._main_rect]
        axes = [self._fig.add_axes(list(r)) for r in rects]
        for ax in axes:
            self._connect_view_watch(ax)
        self._main_axes.extend(axes)
        return axes[0] if n == 1 else axes

    def _connect_view_watch(self, ax):
        """(Re)connect the view-change watch on a main axes -- ax.clear()
        wipes matplotlib's ENTIRE callback registry for that axes
        (confirmed empirically), so this must be called again after
        every single clear(), not just once when the axes is created."""
        ax.callbacks.connect('xlim_changed', self._on_axes_view_changed)
        ax.callbacks.connect('ylim_changed', self._on_axes_view_changed)

    def _on_axes_view_changed(self, ax):
        """Fires on ANY xlim/ylim change to a main axes -- including the
        matplotlib navigation toolbar's own rectangle-zoom/pan tools,
        which set limits directly and have no other way to notify us.
        Ignored during our own redraws (self._redrawing) so a fresh
        autoscale doesn't get mistaken for a user-driven zoom."""
        if self._redrawing or not hasattr(self, '_zoom_sync'):
            return
        self._zoom_locked = True
        self._zoom_sync()

    def _add_control_axes(self, rect):
        ax = self._fig.add_axes(rect)
        self._dynamic_axes.append(ax)
        return ax

    def _track(self, widget):
        """Register a stage-local Button/Slider/TextBox for teardown."""
        self._dynamic_widgets.append(widget)
        return widget

    def _relabel_above(self, widget):
        """Move a Slider/TextBox's own label from matplotlib's default
        position (outside its axes, to the LEFT, in axes-fraction coords)
        to sit centered ABOVE the axes instead. The default position only
        looks fine in isolation -- pack two such widgets within roughly a
        label-width of each other (as several rows here do) and the
        left-bleeding label of the second one lands on top of the first
        one's axes (confirmed via get_window_extent(), not just visual
        guessing). Repositioning is a general fix for the whole class of
        bug, rather than hand-tuning gaps around each label's guessed width.

        1.3, not some larger offset -- far enough above the axes to
        clear it, but a p-slider whose only close neighbor is a row well
        above it (not directly touching) still needs its OWN label to
        read as clearly attached to it, not to whatever's in that row
        instead (confirmed as a real visual-clarity complaint, not just
        an overlap: 'p' floated distractingly close to a Next-Order
        button two full rows above the actual p-slider)."""
        widget.label.set_position((0.5, 1.3))
        widget.label.set_ha('center')
        widget.label.set_va('bottom')
        return widget

    def _add_text(self, *args, **kwargs):
        """fig.text() label that gets removed on the next stage switch."""
        t = self._fig.text(*args, **kwargs)
        self._dynamic_texts.append(t)
        return t

    def _tabbable(self, tb):
        """Register a TextBox for Tab-cycling and double-click-to-clear
        (see _on_tab_key()/_on_textbox_dblclick()) -- used for the Order
        and zoom boxes specifically, in the order they should tab through.
        Does NOT change teardown; pair with self._track() as usual."""
        self._tab_order.append(tb)
        return tb

    def _on_tab_key(self, event):
        """Tab/Shift+Tab cycles focus among the current stage's Order/zoom
        TextBoxes (self._tab_order) -- matplotlib's TextBox has no native
        notion of a tab order across independent widgets, so this walks
        the list ourselves: stop_typing() on whichever box currently has
        focus (this also submits its value, same as clicking away would)
        and begin_typing() on the next one."""
        # macOS's native backend reports Shift+Tab as 'shift+backtab', not
        # 'shift+tab' (confirmed empirically -- other backends may differ,
        # so both spellings are accepted).
        if event.key not in ('tab', 'shift+tab', 'shift+backtab', 'backtab') or not self._tab_order:
            return
        focused = [i for i, tb in enumerate(self._tab_order) if tb.capturekeystrokes]
        if not focused:
            return
        i = focused[0]
        self._tab_order[i].stop_typing()
        reverse = event.key in ('shift+tab', 'shift+backtab', 'backtab')
        nxt = self._tab_order[(i + (-1 if reverse else 1)) % len(self._tab_order)]
        nxt.begin_typing()
        nxt.cursor_index = len(nxt.text)
        nxt._rendercursor()

    def _on_textbox_dblclick(self, event):
        """Double-clicking one of the Order/zoom TextBoxes clears it and
        starts typing fresh -- matplotlib's TextBox has no concept of a
        text selection/highlight to click-and-replace, so this is a
        practical stand-in for the same "quickly overwrite the number"
        workflow rather than a literal highlight."""
        if not event.dblclick:
            return
        for tb in self._tab_order:
            if tb.ax.get_visible() and event.inaxes is tb.ax:
                # NOT tb.set_val('') -- that fires the 'submit' observer
                # synchronously, and a box whose own on_submit reverts
                # invalid/empty text (e.g. the Order box's int() parse)
                # would immediately undo the clear before the user ever
                # gets to type. Clear the display directly instead.
                tb.text_disp.set_text('')
                tb.cursor_index = 0
                tb.begin_typing()
                tb._rendercursor()
                break

    # ------------------------------------------------------------------ #
    # Stage navigation                                                     #
    # ------------------------------------------------------------------ #

    def _set_btn_visible(self, btn, visible):
        """set_visible() ALONE is not enough for Back/Next/Combine: Next
        and Combine deliberately share the same screen position (Next
        takes over Combine's slot once combining is done), and
        matplotlib's Button._click()/_release() only check `self.active`
        (via .ignore()) -- NOT axes visibility -- before processing a
        click. Confirmed via a real 'RuntimeError: Another Axes already
        grabs mouse input': the INVISIBLE button's own _click() still
        fired first (by construction/connection order), grabbed the
        mouse, and then WON the release too, firing its own action
        instead of the visible button's -- e.g. clicking "Combine
        Spectra" silently advancing straight to Normalization instead,
        because the hidden Next button underneath it grabbed the click.
        Setting .active = visible makes ignore() correctly skip the
        hidden widget's own event handling entirely."""
        btn.ax.set_visible(visible)
        btn.active = visible

    def _flash_busy(self, btn):
        """Immediate visual feedback that a click registered and slow work
        (normalize_all()/apply_rv_shift()/combine_spectra(), all called
        synchronously from _goto_stage()/_do_combine()) has started.
        matplotlib doesn't repaint mid-callback on its own -- without the
        explicit draw()+flush_events() here, a color/text change made
        just before a blocking call wouldn't actually appear on screen
        until AFTER that call returns, which defeats the whole point."""
        btn.color = _BUSY_COLOR
        btn.ax.set_facecolor(_BUSY_COLOR)
        self._fig.canvas.draw()
        self._fig.canvas.flush_events()

    def _on_nav_clicked(self, btn, target_stage):
        if self._busy or self._too_soon_after_transition():
            return
        self._busy = True
        self._flash_busy(btn)
        self._goto_stage(target_stage)
        self._busy = False

    def _on_combine_clicked(self):
        if self._busy or self._too_soon_after_transition():
            return
        self._busy = True
        self._flash_busy(self._btn_combine)
        self._do_combine()
        self._busy = False

    def _too_soon_after_transition(self):
        return (time.time() - self._last_transition_time) < _TRANSITION_DEBOUNCE_SECONDS

    def _goto_stage(self, stage):
        stage = max(0, min(len(_STAGES) - 1, stage))
        self._stage = stage
        self._last_transition_time = time.time()
        self._tab_order = []

        # Whichever button triggered this (if any) was just set to
        # _BUSY_COLOR by _flash_busy() -- restore normal colors now that
        # the slow work (if any) is done and the new stage is ready.
        self._btn_back.color = _DEFAULT_BTN_COLOR
        self._btn_back.ax.set_facecolor(_DEFAULT_BTN_COLOR)
        self._btn_next.color = _COMBINE_COLOR
        self._btn_next.ax.set_facecolor(_COMBINE_COLOR)
        self._btn_combine.color = _COMBINE_COLOR
        self._btn_combine.ax.set_facecolor(_COMBINE_COLOR)

        self._set_btn_visible(self._btn_back, stage > 0)
        if stage > 0:
            self._btn_back.label.set_text(f'◀ {_STAGES[stage - 1]}')

        # On the Load stage with multiple exposures, Next is pointless
        # (and misleading) until Combine Spectra has actually been run --
        # there's nothing to normalize yet.
        pending_combine = stage == 0 and self._n_exposures > 1 and not self._combined
        next_visible = stage < len(_STAGES) - 1 and not pending_combine
        self._set_btn_visible(self._btn_next, next_visible)
        if next_visible:
            self._btn_next.label.set_text(f'{_STAGES[stage + 1]} ▶')

        self._set_btn_visible(self._btn_combine, pending_combine)

        self._zoom_locked = False
        self._clear_dynamic()
        self._clear_main()

        if stage == 0:
            self._build_stage_load()
        elif stage == 1:
            self._build_stage_normalize()
        elif stage == 2:
            self._build_stage_rv()
        elif stage == 3:
            self._build_stage_ew()

        self._fig.canvas.draw_idle()

    def _order_pager(self, y, on_change):
        """Shared Previous/Next-order control plus a jump-to-order box,
        used by stages 1-3. The box always shows the current order (kept
        in sync from the Prev/Next buttons too), and typing a number + Enter
        jumps straight there."""
        ax_prev = self._add_control_axes([0.08, y, 0.16, 0.04])
        ax_next = self._add_control_axes([0.26, y, 0.16, 0.04])
        btn_prev = self._track(Button(ax_prev, '◀ Previous Order'))
        btn_next = self._track(Button(ax_next, 'Next Order ▶'))

        # Orders are always <=3 digits text-wise, but a too-narrow box is
        # an unreliable double-click target -- 0.045 balances the two.
        ax_tb = self._add_control_axes([0.485, y, 0.045, 0.04])
        # NOT _relabel_above() -- this row (y=0.135) sits close enough
        # below the main plot's own x-axis label/tick text that the
        # relabeled-above position collided with it (confirmed via
        # get_window_extent() on the actual xlabel, not just the axes
        # patch -- easy to miss since the patch itself has plenty of
        # clearance). The default left-bleeding position is safe here.
        tb_order = self._tabbable(self._track(
            TextBox(ax_tb, 'Order:  ', initial=str(self._order))))

        def _goto_order(new_order):
            n = len(self.spec.wavelength)
            new_order = max(0, min(n - 1, new_order))
            self._order = new_order
            self._zoom_locked = False  # a different order has a different natural range
            tb_order.eventson = False
            tb_order.set_val(str(self._order))
            tb_order.eventson = True
            on_change()

        def _step(delta):
            _goto_order(self._order + delta)

        def _submit_order(text):
            try:
                new_order = int(text)
            except ValueError:
                tb_order.set_val(str(self._order))
                return
            _goto_order(new_order)

        btn_prev.on_clicked(lambda _e: _step(-1))
        btn_next.on_clicked(lambda _e: _step(1))
        tb_order.on_submit(_submit_order)

    def _capture_zoom(self, axes):
        """Grab current view limits, if zoom is locked, to survive the
        upcoming ax.clear()+replot -- see self._zoom_locked's docstring."""
        if not self._zoom_locked:
            return None
        return [(ax.get_xlim(), ax.get_ylim()) for ax in axes]

    def _apply_zoom(self, axes, saved):
        if saved is None:
            # Force autoscale to compute NOW, synchronously, while the
            # caller's self._redrawing guard is still active. Left to
            # matplotlib's own default behavior, this computation happens
            # lazily at actual render/draw time -- by which point
            # self._redrawing has already been reset to False, so the
            # resulting xlim_changed event would be mistaken for a
            # toolbar-driven user zoom (confirmed: caused EVERY ordinary
            # redraw to spuriously set zoom_locked=True).
            for ax in axes:
                ax.relim()
                ax.autoscale_view()
            return
        for ax, (xl, yl) in zip(axes, saved):
            ax.set_xlim(xl)
            ax.set_ylim(yl)

    def _build_zoom_strip(self, axes_getter):
        """Compact wavelength/flux zoom strip (xmin-xmax, ymin-ymax + Reset),
        mirroring moog/interactive.py's SynthWidget zoom controls. Used by
        every stage -- `axes_getter()` returns whichever Axes are current
        for that stage (recomputed live, since main-area axes are torn
        down and rebuilt on every stage switch)."""
        # Box widths sized to their actual content: wavelengths never
        # exceed "xxxx.xx" (7 chars), orders/order-box elsewhere never
        # exceed 3 digits, flux can vary widely but never exceeds 7 digits
        # (plus room for a sign/decimal).
        y, h = 0.01, 0.03
        self._add_text(0.08, y + h / 2, 'λ:', fontsize=9, ha='right', va='center')
        ax_xmin = self._add_control_axes([0.10, y, 0.09, h])
        self._add_text(0.20, y + h / 2, '-', fontsize=9, ha='center', va='center')
        ax_xmax = self._add_control_axes([0.21, y, 0.09, h])
        self._add_text(0.36, y + h / 2, 'Flux:', fontsize=9, ha='right', va='center')
        ax_ymin = self._add_control_axes([0.38, y, 0.10, h])
        self._add_text(0.49, y + h / 2, '-', fontsize=9, ha='center', va='center')
        ax_ymax = self._add_control_axes([0.50, y, 0.10, h])
        ax_reset = self._add_control_axes([0.65, y, 0.13, h])

        ax0 = axes_getter()[0]
        xl, yl = ax0.get_xlim(), ax0.get_ylim()
        tb_xmin = self._tabbable(self._track(TextBox(ax_xmin, '', initial=f'{xl[0]:.2f}')))
        tb_xmax = self._tabbable(self._track(TextBox(ax_xmax, '', initial=f'{xl[1]:.2f}')))
        tb_ymin = self._tabbable(self._track(TextBox(ax_ymin, '', initial=f'{yl[0]:.3f}')))
        tb_ymax = self._tabbable(self._track(TextBox(ax_ymax, '', initial=f'{yl[1]:.3f}')))
        btn_reset = self._track(Button(ax_reset, 'Reset Zoom'))
        for tb in (tb_xmin, tb_xmax, tb_ymin, tb_ymax):
            tb.text_disp.set_fontsize(9)
        btn_reset.label.set_fontsize(9)

        suppress = {'on': False}

        def _sync_textboxes():
            xl, yl = axes_getter()[0].get_xlim(), axes_getter()[0].get_ylim()
            suppress['on'] = True
            tb_xmin.set_val(f'{xl[0]:.2f}')
            tb_xmax.set_val(f'{xl[1]:.2f}')
            tb_ymin.set_val(f'{yl[0]:.3f}')
            tb_ymax.set_val(f'{yl[1]:.3f}')
            suppress['on'] = False

        def _submit(_text):
            if suppress['on']:
                return
            try:
                xmin, xmax = float(tb_xmin.text), float(tb_xmax.text)
                ymin, ymax = float(tb_ymin.text), float(tb_ymax.text)
            except ValueError:
                return
            if xmin >= xmax or ymin >= ymax:
                return
            self._zoom_locked = True
            self._redrawing = True
            for ax in axes_getter():
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
            self._redrawing = False
            self._fig.canvas.draw_idle()

        def _reset(_e):
            self._zoom_locked = False
            # Without this guard, autoscale_view() firing xlim_changed
            # would hit _on_axes_view_changed() and immediately set
            # _zoom_locked back to True, undoing the line above.
            self._redrawing = True
            for ax in axes_getter():
                ax.set_autoscale_on(True)
                ax.relim()
                ax.autoscale_view()
            self._redrawing = False
            _sync_textboxes()
            self._fig.canvas.draw_idle()

        for tb in (tb_xmin, tb_xmax, tb_ymin, tb_ymax):
            tb.on_submit(_submit)
        btn_reset.on_clicked(_reset)

        # Stored so every redraw function can keep the displayed numbers in
        # sync with reality (a fresh autoscale, or a restored zoom lock)
        # without the caller needing to know the zoom strip's internals.
        self._zoom_sync = _sync_textboxes

    # ------------------------------------------------------------------ #
    # Stage 1: Load                                                        #
    # ------------------------------------------------------------------ #

    def _build_stage_load(self):
        n_orders = len(self.spec.wavelength)
        wmin = min(w.min() for w in self.spec.wavelength)
        wmax = max(w.max() for w in self.spec.wavelength)

        msg = f"{self._n_exposures} exposure(s), {n_orders} orders, {wmin:.1f}-{wmax:.1f} Å"
        if self.spec.lines is not None:
            msg += f"  |  {len(self.spec.lines)} lines loaded"
            if self._linelist_path:
                msg += f" from {os.path.basename(self._linelist_path)}"
        else:
            msg += "  |  no linelist loaded (pass linelist= to enable Step 4)"
        print(msg)
        self._load_msg = msg

        self._ax_load = self._add_main_axes()
        # Same y as every other stage's order pager (0.135) -- keeps a
        # consistent, already-verified gap to the plot's own x-axis
        # label/ticks above it. The view-mode row goes BELOW this one
        # instead (see _build_load_multi_controls()), using the vertical
        # room freed up now that Combine Spectra lives in the header.
        self._order_pager(0.135, self._redraw_load)

        if self._n_exposures > 1:
            self._build_load_multi_controls()

        self._build_zoom_strip(lambda: [self._ax_load])
        self._redraw_load()

    def _build_load_multi_controls(self):
        """View-mode selector (one button per exposure + Overlay) and the
        Combine action -- only shown when more than one exposure was
        provided. Lets the user visually confirm the raw exposures line up
        (Overlay) before co-adding them, and inspect any single exposure on
        its own, rather than only ever seeing the already-combined result.

        Each exposure's button is colored to match that exposure's own
        line color in the Overlay plot (self._exposure_colors), and the
        currently-displayed view is marked with a bold border rather than
        a fill-color swap -- so the color-coding stays visible no matter
        which view is active. The row spans the full control width and
        sizes itself to however many exposures are actually present (2,
        3, or more), rather than assuming exactly 2.
        """
        labels = [f'Exposure {i + 1}' for i in range(self._n_exposures)]
        button_colors = list(self._exposure_colors)
        modes = list(range(self._n_exposures))
        labels.append('Overlay')
        button_colors.append('0.85')
        modes.append('overlay')
        if self._combined and getattr(self.spec, 'combine_debug', None):
            # Only offered once Combine Spectra has actually run --
            # combine_debug is populated by combine_spectra() itself (see
            # spectrum_data.py) and holds the aligned-but-not-yet-summed
            # inputs behind each combined order.
            # Same color as the combined trace itself in _redraw_load()'s
            # 'combined' branch (color='k') -- matches the Exposure-N
            # button <-> line color convention already used for Overlay.
            labels.append('Combined')
            button_colors.append('black')
            modes.append('combined')

        left, right, gap = 0.08, 0.94, 0.015
        width = (right - left - gap * (len(labels) - 1)) / len(labels)
        self._view_buttons = []
        for k, (lab, col, mode) in enumerate(zip(labels, button_colors, modes)):
            ax_b = self._add_control_axes([left + k * (width + gap), 0.07, width, 0.045])
            btn = self._track(Button(ax_b, lab, color=col, hovercolor=col))
            if col != '0.85':
                btn.label.set_color('white')
            btn.on_clicked(lambda _e, m=mode: self._set_view(m))
            self._view_buttons.append(btn)
        self._refresh_view_button_colors()

    def _view_label(self):
        if self._view_idx == 'overlay':
            return 'Overlay'
        if self._view_idx == 'combined':
            return 'Combined'
        return f'Exposure {self._view_idx + 1}'

    def _refresh_view_button_colors(self):
        active = self._view_label()
        for btn in self._view_buttons:
            is_active = (btn.label.get_text() == active)
            for spine in btn.ax.spines.values():
                spine.set_linewidth(1.75 if is_active else 0.8)
                # Not black -- indistinguishable from the Combined
                # button's own black fill when IT is the active one.
                # Same red used elsewhere in this codebase for emphasis
                # (e.g. plot_ew_fit()'s line-bound markers).
                spine.set_edgecolor('#e41a1c' if is_active else '0.4')
        self._fig.canvas.draw_idle()

    def _set_view(self, mode):
        self._view_idx = mode
        self._refresh_view_button_colors()
        self._redraw_load()

    def _do_combine(self):
        print(f"Normalizing all {self._n_exposures} exposures before combining...")
        for e in self._exposures:
            e.normalize_all()
        primary = self._exposures[0]
        for other in self._exposures[1:]:
            n = primary.combine_spectra(other, verbose=True)
            primary.update_combined()
            print(f"Combined {n} orders.")
        self._combined = True
        self._view_idx = 'combined'  # jump straight to the result just produced
        self._set_btn_visible(self._btn_combine, False)
        # Full rebuild, not just _redraw_load() -- the view-mode row needs
        # to add the new "Combined" button, not just redraw the plot.
        self._goto_stage(0)

    def _redraw_load(self):
        ax = self._ax_load
        saved_zoom = self._capture_zoom([ax])
        self._redrawing = True
        ax.clear()
        self._connect_view_watch(ax)
        o = self._order

        if self._view_idx == 'overlay':
            for i, e in enumerate(self._exposures):
                if o >= len(e.wavelength):
                    continue
                wave = e.wavelength[o]
                flux = self._raw_flux_snapshots[i][o]
                # Same color as that exposure's own view-mode button
                # (self._exposure_colors, set in _build_load_multi_controls())
                # so the two are visually tied together.
                ax.plot(wave, flux, lw=0.6, alpha=0.8, color=self._exposure_colors[i],
                        label=f'Exposure {i + 1} ({self._exposure_names[i]})')
            ax.legend(loc='upper right', fontsize=8)
            title_extra = f'overlay of {self._n_exposures} raw exposures'
        elif self._view_idx == 'combined':
            debug = self.spec.combine_debug
            if o >= len(debug):
                wave = self.spec.shifted_wavelength[o]
                ax.text(0.5, 0.5, 'No alignment data for this order.',
                        transform=ax.transAxes, ha='center', va='center')
                title_extra = 'combined result -- no alignment data for this order'
            else:
                d = debug[o]
                wave = d['wave']
                ax.plot(wave, d['flux_A'], lw=0.7, alpha=0.8, color=self._exposure_colors[0],
                        label='Exposure 1 (aligned)')
                if d['flux_B_aligned'] is not None:
                    b_color = self._exposure_colors[1] if len(self._exposure_colors) > 1 else '#ff7f0e'
                    ax.plot(wave, d['flux_B_aligned'], lw=0.7, alpha=0.8, color=b_color,
                            label='Exposure 2 (aligned, resampled)')
                ax.plot(wave, d['combined'], lw=1.0, color='k', label='Combined (post spike-correction)')
                ax.legend(loc='upper right', fontsize=8)
                title_extra = 'aligned inputs vs. combined result, pre-commit'
        else:
            i = self._view_idx
            wave = self._exposures[i].wavelength[o]
            flux = self._raw_flux_snapshots[i][o]
            ax.plot(wave, flux, '-', color=self._exposure_colors[i], lw=0.6)
            title_extra = f'Exposure {i + 1} ({self._exposure_names[i]}), raw, unfit'

        ax.set_xlabel('Wavelength (Å)')
        ax.set_ylabel('Raw flux')
        # fontsize=8, not 9 -- this is the only 2-line title of the 4
        # stages (the others are single-line), and at fontsize 9 its
        # rendered height collided with the Back/Next/Combine row above
        # (confirmed via get_window_extent()); 8 fits within the same
        # clearance the single-line titles already have.
        ax.set_title(f'{self._load_msg}\n'
                      f'Order {o}/{len(self.spec.wavelength) - 1}  '
                      f'({wave.min():.1f}-{wave.max():.1f} Å)  -- {title_extra}',
                      fontsize=8, loc='left')
        self._apply_zoom([ax], saved_zoom)
        self._redrawing = False
        self._zoom_sync()
        self._fig.canvas.draw_idle()

    # ------------------------------------------------------------------ #
    # Stage 2: Normalize                                                   #
    # ------------------------------------------------------------------ #

    def _build_stage_normalize(self):
        if not self._normalized:
            self.spec.normalize_all(lam=self._lam, p=self._p)
            self._normalized = True
            # Every order was just fit with these same lam/p -- see
            # _order_norm_params's docstring note in __init__.
            self._order_norm_params = {o: (self._lam, self._p)
                                        for o in range(len(self.spec.wavelength))}

        ax = self._add_main_axes()
        self._ax_norm = ax
        self._order_pager(0.135, self._on_normalize_order_change)

        ax_lam = self._add_control_axes([0.58, 0.135, 0.20, 0.03])
        # Same x/width as lam (0.58, inline in the same column) rather
        # than its own previous x=0.30 -- that used to sit close enough
        # to the order-pager row that its relabeled-above label read as
        # attached to "Next Order" instead of to the slider itself.
        # y=0.055, not 0.08 -- directly beneath lam now, p's own
        # relabeled-above label needs clearance from LAM'S axes above it
        # (confirmed via audit: 0.08 wasn't enough gap to lam at 0.135).
        ax_p = self._add_control_axes([0.58, 0.055, 0.20, 0.03])
        self._sl_lam = self._track(self._relabel_above(Slider(ax_lam, 'log10(lam)', 1.0, 6.0,
                                           valinit=np.log10(self._lam), valfmt='%1.3f')))
        self._sl_p = self._track(self._relabel_above(Slider(ax_p, 'p', 0.001, 0.10, valinit=self._p,
                                                              valfmt='%1.3f')))
        self._sl_lam.on_changed(lambda _v: self._refit_current_order())
        self._sl_p.on_changed(lambda _v: self._refit_current_order())

        # To the right of both sliders (which span 0.58-0.78), narrow and
        # tall enough for its label to wrap to two lines, spanning the
        # same overall height as lam+p stacked (0.055 to 0.165).
        ax_apply_all = self._add_control_axes([0.83, 0.055, 0.12, 0.11])
        self._btn_apply_all = self._track(Button(ax_apply_all, 'Apply to\nALL orders'))
        self._btn_apply_all.on_clicked(lambda _e: self._on_apply_all_clicked())

        # Raw+fit vs. the actual normalized result -- same toggle-button
        # pattern as the RV stage's "Showing: shifted/raw", so the user
        # can review what normalize_all()/Apply to ALL orders actually
        # produced before moving on to RV Shift (mirrors Stage 1's
        # Combine-then-inspect-the-result flow; there's no single
        # discrete "commit" action to gate Next on here the way Combine
        # had, so this adds the review capability rather than blocking
        # navigation until it's used). Directly under Previous Order.
        # y=0.075, not 0.055 -- too close to the zoom strip's own boxes
        # right below it at that position.
        ax_norm_toggle = self._add_control_axes([0.08, 0.075, 0.16, 0.045])
        color, hovercolor = _NORM_VIEW_COLORS[self._norm_view]
        self._btn_norm_toggle = self._track(Button(
            ax_norm_toggle, 'Raw+Fit' if self._norm_view == 'raw' else 'Normalized',
            color=color, hovercolor=hovercolor))
        self._btn_norm_toggle.on_clicked(lambda _e: self._toggle_norm_view())

        self._build_zoom_strip(lambda: [self._ax_norm])
        # Not just self._redraw_normalize() -- also covers revisiting this
        # stage while self._order was last changed somewhere ELSE (e.g.
        # paged to a different order on the RV stage, then clicked Back),
        # which would otherwise leave the sliders showing stale values
        # from whatever order was tuned before leaving Normalize.
        self._sync_norm_sliders_to_order()
        self._redraw_normalize()

    def _toggle_norm_view(self):
        self._norm_view = 'normalized' if self._norm_view == 'raw' else 'raw'
        self._btn_norm_toggle.label.set_text(
            'Raw+Fit' if self._norm_view == 'raw' else 'Normalized')
        color, hovercolor = _NORM_VIEW_COLORS[self._norm_view]
        self._btn_norm_toggle.color = color
        self._btn_norm_toggle.hovercolor = hovercolor
        self._btn_norm_toggle.ax.set_facecolor(color)
        self._redraw_normalize()

    def _refit_current_order(self):
        self._lam = 10 ** self._sl_lam.val
        self._p = self._sl_p.val
        self.spec.normalize(self._order, lam=self._lam, p=self._p)
        self._order_norm_params[self._order] = (self._lam, self._p)
        self._redraw_normalize()

    def _on_apply_all_clicked(self):
        # normalize_all() over every order can take a few seconds --
        # same busy-flash pattern as Back/Next/Combine, but this button
        # doesn't go through _goto_stage() (no stage change happens), so
        # it resets its own color here instead of relying on that.
        if self._busy:
            return
        self._busy = True
        self._flash_busy(self._btn_apply_all)
        self._apply_normalize_all()
        self._btn_apply_all.color = _DEFAULT_BTN_COLOR
        self._btn_apply_all.ax.set_facecolor(_DEFAULT_BTN_COLOR)
        self._fig.canvas.draw_idle()
        self._busy = False

    def _apply_normalize_all(self):
        self.spec.normalize_all(lam=self._lam, p=self._p)
        self._order_norm_params = {o: (self._lam, self._p)
                                    for o in range(len(self.spec.wavelength))}
        print(f"normalize_all(lam={self._lam:.1f}, p={self._p:.4f}) applied to all orders")
        self._redraw_normalize()

    def _sync_norm_sliders_to_order(self):
        """Set the lam/p sliders to whatever actually produced the CURRENT
        order's own fit, rather than leaving them showing values left
        over from whichever order was tuned last -- see
        self._order_norm_params's docstring in __init__."""
        lam, p = self._order_norm_params.get(self._order, (self._lam, self._p))
        self._lam, self._p = lam, p
        self._sl_lam.eventson = False
        self._sl_lam.set_val(np.log10(lam))
        self._sl_lam.eventson = True
        self._sl_p.eventson = False
        self._sl_p.set_val(p)
        self._sl_p.eventson = True

    def _on_normalize_order_change(self):
        """The order-pager's on_change callback for this stage -- syncs
        the sliders (see _sync_norm_sliders_to_order()) before redrawing,
        so they never silently show a different order's values."""
        self._sync_norm_sliders_to_order()
        self._redraw_normalize()

    def _redraw_normalize(self):
        ax = self._ax_norm
        saved_zoom = self._capture_zoom([ax])
        self._redrawing = True
        ax.clear()
        self._connect_view_watch(ax)
        o = self._order
        wave = self.spec.wavelength[o]

        if self._norm_view == 'normalized':
            ax.plot(wave, self.spec.normalized_flux[o], 'k-', lw=0.6, label='normalized flux')
            ax.axhline(1.0, color='#377eb8', lw=1.0, ls='--', label='continuum (=1)')
            ax.set_ylabel('Normalized flux')
            view_extra = 'normalized result'
        else:
            flux = self.spec.flux[o]
            pred = self.spec.pred_all[o]
            ax.plot(wave, flux, 'k-', lw=0.6, label='raw flux')
            ax.plot(wave, pred, '-', color='#377eb8', lw=1.2, label='continuum fit')
            ax.set_ylabel('Flux')
            view_extra = 'raw + continuum fit'

        ax.set_xlabel('Wavelength (Å)')
        ax.set_title(f'Order {o}/{len(self.spec.wavelength) - 1} '
                      f'({wave.min():.1f}-{wave.max():.1f} Å)  '
                      f'(lam={self._lam:.1f}, p={self._p:.4f})  -- {view_extra}',
                      fontsize=10, loc='left')
        ax.legend(loc='lower right', fontsize=8)
        self._apply_zoom([ax], saved_zoom)
        self._redrawing = False
        self._zoom_sync()
        self._fig.canvas.draw_idle()

    # ------------------------------------------------------------------ #
    # Stage 3: RV shift                                                    #
    # ------------------------------------------------------------------ #

    def _build_stage_rv(self):
        if not self._rv_applied:
            self.spec.apply_rv_shift(verbose=True)
            self._rv_applied = True

        ax = self._add_main_axes()
        self._ax_rv = ax
        self._order_pager(0.135, self._redraw_rv)

        # Directly under Previous Order -- same spot as the Normalize
        # stage's own Raw+Fit/Normalized toggle, for a consistent pattern.
        ax_toggle = self._add_control_axes([0.08, 0.075, 0.16, 0.045])
        color, hovercolor = _RV_VIEW_COLORS[self._show_shifted]
        self._btn_toggle = self._track(Button(
            ax_toggle, 'Shifted' if self._show_shifted else 'Observed',
            color=color, hovercolor=hovercolor))
        self._btn_toggle.on_clicked(lambda _e: self._toggle_shifted())

        # Widened now that "Remeasure RV" (removed -- unreliable RVs
        # anyway, per the user) no longer needs the space next to it.
        ax_rv_tb = self._add_control_axes([0.46, 0.08, 0.20, 0.03])
        rv0 = self.spec.rv[0] if self.spec.rv is not None else 0.0
        # NOT _relabel_above() here -- this row (y=0.08) sits directly
        # below the order-pager row with no vertical clearance, so
        # relabeling above collided with it (confirmed via audit). The
        # default left-bleeding position is safe: nothing else occupies
        # this row to its left now that the toggle moved away.
        self._tb_rv = self._track(
            TextBox(ax_rv_tb, 'RV override (km/s)  ', initial=f'{rv0:.3f}'))
        self._tb_rv.on_submit(self._on_rv_submit)

        self._build_zoom_strip(lambda: [self._ax_rv])
        self._redraw_rv()

    def _toggle_shifted(self):
        self._show_shifted = not self._show_shifted
        self._btn_toggle.label.set_text('Shifted' if self._show_shifted else 'Observed')
        color, hovercolor = _RV_VIEW_COLORS[self._show_shifted]
        self._btn_toggle.color = color
        self._btn_toggle.hovercolor = hovercolor
        self._btn_toggle.ax.set_facecolor(color)
        self._redraw_rv()

    def _on_rv_submit(self, text):
        try:
            rv = float(text)
        except ValueError:
            print(f"Could not parse RV override '{text}' as a number.")
            return
        self.spec.apply_rv_shift(rv=rv, verbose=True)
        self._redraw_rv()

    def _redraw_rv(self):
        ax = self._ax_rv
        saved_zoom = self._capture_zoom([ax])
        self._redrawing = True
        ax.clear()
        self._connect_view_watch(ax)
        o = self._order
        wave = self.spec.shifted_wavelength[o] if self._show_shifted else self.spec.wavelength[o]
        ax.plot(wave, self.spec.normalized_flux[o], 'k-', lw=0.6)
        if self.spec.lines is not None:
            lo, hi = wave.min(), wave.max()
            for line in self.spec.lines:
                if lo <= line <= hi:
                    ax.axvline(line, color='#e41a1c', lw=0.7, ls='--', alpha=0.6)
        rv_txt = f'{self.spec.rv[0]:.3f} ± {self.spec.rv[1]:.3f} km/s' if self.spec.rv else 'not measured'
        ax.set_xlabel('Wavelength (Å)')
        ax.set_ylabel('Normalized flux')
        ax.set_title(f'Order {o}  ({"rest-frame" if self._show_shifted else "observed"} wavelength)  '
                      f'|  RV = {rv_txt}', fontsize=10, loc='left')
        self._apply_zoom([ax], saved_zoom)
        self._redrawing = False
        self._zoom_sync()
        self._fig.canvas.draw_idle()

    # ------------------------------------------------------------------ #
    # Stage 4: Measure EW                                                  #
    # ------------------------------------------------------------------ #

    def _find_order_for_line(self, rest_wave, search_radius=0.15):
        for o, w in enumerate(self.spec.shifted_wavelength):
            if w[0] - search_radius <= rest_wave <= w[-1] + search_radius:
                return o
        return None

    def _build_stage_ew(self):
        if self.spec.lines is None or len(self.spec.lines) == 0:
            ax = self._add_main_axes()
            ax.text(0.5, 0.5, 'No linelist loaded -- pass linelist= to the\n'
                              'widget (or call spec.load_lines() yourself) '
                              'to enable EW measurement.',
                    transform=ax.transAxes, ha='center', va='center', fontsize=12)
            ax.set_xticks([])
            ax.set_yticks([])
            return

        if not self._rv_applied:
            print("Note: RV shift has not been applied yet (visit Step 3 first) -- "
                  "measuring against the unshifted wavelength solution.")

        self._ax_fit, self._ax_data = self._add_main_axes(2)
        n = len(self.spec.lines)

        ax_prev = self._add_control_axes([0.08, 0.135, 0.16, 0.04])
        ax_next = self._add_control_axes([0.26, 0.135, 0.16, 0.04])
        self._btn_line_prev = self._track(Button(ax_prev, '◀ Previous Line'))
        self._btn_line_next = self._track(Button(ax_next, 'Next Line ▶'))
        self._btn_line_prev.on_clicked(lambda _e: self._step_line(-1))
        self._btn_line_next.on_clicked(lambda _e: self._step_line(1))
        self._line_label = self._add_text(0.46, 0.155, '', fontsize=9, va='center')

        labels = ['cont. shift', 'left bound (Å)', 'right bound (Å)', 'center (Å)']
        ranges = [(-0.05, 0.05), (-0.3, 0.3), (-0.3, 0.3), (-0.2, 0.2)]
        self._sl_ex = []
        for k, (lab, (lo, hi)) in enumerate(zip(labels, ranges)):
            # y=0.055, not 0.08 -- _relabel_above() puts each slider's own
            # label above its axes, and at 0.08 that landed inside the
            # line-pager row (y=0.135) for the first slider (confirmed
            # via audit); 0.055 gives enough clearance underneath it.
            ax_s = self._add_control_axes([0.32 + k * 0.145, 0.055, 0.13, 0.03])
            sl = self._track(self._relabel_above(Slider(ax_s, lab, lo, hi, valinit=0.0)))
            # Slider's default numeric readout sits just to the RIGHT of
            # its axes (mirror image of the label's left-bleed problem);
            # packed this tightly (0.015 gap between axes), each one
            # landed on the NEXT slider's axes (confirmed via audit --
            # pre-existing, not introduced by today's other layout
            # changes). No room to relocate it without recreating the
            # same crowding, so it's hidden here; the slider's own
            # position still shows the current value interactively.
            sl.valtext.set_visible(False)
            sl.on_changed(lambda _v: self._remeasure_ew())
            self._sl_ex.append(sl)

        # (was previously placed at x > 1.0, off the right edge of the
        # figure -- fixed alongside the pager-button relabeling above)
        ax_reset = self._add_control_axes([0.90, 0.055, 0.06, 0.05])
        self._btn_ex_reset = self._track(Button(ax_reset, 'Reset'))
        self._btn_ex_reset.on_clicked(lambda _e: self._reset_ex_params())

        self._build_zoom_strip(lambda: [self._ax_fit, self._ax_data])
        self._remeasure_ew()

    def _step_line(self, delta):
        self._line_idx = (self._line_idx + delta) % len(self.spec.lines)
        self._zoom_locked = False  # a different line has a different natural range
        self._reset_ex_params(redraw=False)
        self._remeasure_ew()

    def _reset_ex_params(self, redraw=True):
        for sl in self._sl_ex:
            sl.eventson = False
            sl.set_val(0.0)
            sl.eventson = True
        if redraw:
            self._remeasure_ew()

    def _remeasure_ew(self):
        rest_wave = self.spec.lines[self._line_idx]
        order = self._find_order_for_line(rest_wave)
        self._line_label.set_text(
            f'Line {self._line_idx}/{len(self.spec.lines) - 1}: {rest_wave:.3f} Å')
        if order is None:
            self._redrawing = True
            self._ax_fit.clear()
            self._ax_data.clear()
            self._connect_view_watch(self._ax_fit)
            self._connect_view_watch(self._ax_data)
            self._redrawing = False
            self._ax_fit.text(0.5, 0.5, 'No order covers this line\'s wavelength.',
                               transform=self._ax_fit.transAxes, ha='center', va='center')
            self._fig.canvas.draw_idle()
            return

        axes = [self._ax_fit, self._ax_data]
        saved_zoom = self._capture_zoom(axes)
        self._redrawing = True
        ex_params = [sl.val for sl in self._sl_ex]
        try:
            self.spec.measure_ew(self._line_idx, order, ex_params=ex_params,
                                  axes=(self._ax_fit, self._ax_data))
            # measure_ew() -> plot_ew_fit() does its own fit_view.clear()/
            # data_view.clear() internally (plotting.py), which just as
            # much wipes the callback registry as our own ax.clear() calls
            # do -- must reconnect here too.
            self._connect_view_watch(self._ax_fit)
            self._connect_view_watch(self._ax_data)
            self._apply_zoom(axes, saved_zoom)
            self._zoom_sync()
        except Exception as exc:
            # measure_ew()'s fit can be degenerate for some ex_params
            # combinations (e.g. a center/bound shift narrow enough to
            # leave too few points in the fit window) -- a slider drag
            # must not crash the whole session, so report it in-place
            # and let the user back off the adjustment instead.
            self._ax_fit.clear()
            self._ax_data.clear()
            self._connect_view_watch(self._ax_fit)
            self._connect_view_watch(self._ax_data)
            self._ax_fit.text(0.5, 0.5, f'Fit failed with these ex_params:\n{exc}',
                               transform=self._ax_fit.transAxes, ha='center',
                               va='center', wrap=True, color='#e41a1c')
        finally:
            self._redrawing = False
        self._fig.canvas.draw_idle()


def ew_interactive(filename=None, spectrum=None, spectra=None, linelist=None,
                    KECK_file=True, **reader_kwargs):
    """
    Launch the interactive xspect-ew widget. Returns the EWWidget instance --
    keep a reference to it (e.g. `w = ew_interactive(...)`), since letting it
    be garbage-collected disconnects the figure's callbacks.

    Pass spectra=[exposure1, exposure2, ...] (paths and/or already-
    constructed Spectrum_Data) when more than one exposure of the same
    star is available -- Stage 1 (Load) then lets you inspect them
    individually or overlaid and combine them interactively. Use
    filename=/spectrum= for the single-exposure case.
    """
    return EWWidget(filename=filename, spectrum=spectrum, spectra=spectra,
                     linelist=linelist, KECK_file=KECK_file, **reader_kwargs)
