"""
Interactive matplotlib widget for MOOG synth mode.

Usage (terminal):
    from SPAE.moog.interactive import synth_interactive
    synth_interactive('batch.par')

Usage (Jupyter notebook):
    %matplotlib widget
    from SPAE.moog.interactive import synth_interactive
    synth_interactive('batch.par')

The widget supports:
  - Multi-pass synthesis with per-pass abundance controls
  - Live re-smoothing (Gaussian / rotation / macroturbulence / Lorentzian)
  - Wavelength shift (Δλ) applied live without re-synthesis
  - Auto-loaded observed spectrum from batch.par (if present)
  - Saving updated abundances back to batch.par
"""
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons, TextBox

from .state      import State
from .params     import params as _params
from .inmodel    import inmodel as _inmodel
from .inlines    import inlines as _inlines
from .eqlib      import eqlib as _eqlib
from .nearly     import nearly as _nearly
from .synspec    import synspec as _synspec
from .smooth     import smooth as _smooth
from .plot       import read_obs
from .atomic_data import ELEMENT_NAMES
from .state import NELEM

# Smoothing check-button labels (order must match _resmooth / _derive_smtype)
_SMOOTH_LABELS = ['Gaussian', 'Rotation', 'Macro', 'Lorentz']

_PASS_COLORS   = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
_MAX_SLOTS     = 6        # element slider rows visible at once
_STALE_COLOR   = '#cc3300'  # synthesize button color when re-synthesis needed
_NO_FILE       = 'no_filename_given'


def _iatom_from_name(name: str) -> int:
    """Return 1-indexed atomic number for an element symbol; -1 if not found."""
    name = name.strip().capitalize()
    for i, sym in enumerate(ELEMENT_NAMES):
        if sym.strip() == name:
            return i + 1
    return -1


class SynthWidget:
    """
    Interactive matplotlib GUI for exploring MOOG synth results.

    Keep a reference to the returned object — garbage-collecting it will
    destroy the figure callbacks.
    """

    def __init__(self, fparam: str = 'batch.par'):
        self._fparam  = os.path.abspath(fparam)
        self._rundir  = os.path.dirname(self._fparam)

        # ---- Load parameters and model ----
        self._state = State()
        self._state.fparam = os.path.basename(self._fparam)

        old_dir = os.getcwd()
        os.chdir(self._rundir)
        try:
            _params(self._state, self._state.fparam)
            _inmodel(self._state)

            # Save model log ε BEFORE synthesis may alter xabund
            self._model_logeps = {}
            for z in range(1, NELEM + 1):
                xa = self._state.xabund[z - 1]
                self._model_logeps[z] = (math.log10(xa) + 12.0) if xa > 0 else 0.0

            # Number of synthesis passes
            self._npass = (max(1, self._state.numatomsyn)
                           if self._state.numpecatom > 0 else 1)

            # Initial synthesis
            self._run_synthesis()

            # Optional observed spectrum
            self._obs_wave = self._obs_flux = None
            fobs = getattr(self._state, 'fobs', '') or ''
            if fobs and fobs != _NO_FILE:
                try:
                    self._obs_wave, self._obs_flux = read_obs(fobs)
                except Exception:
                    pass
        finally:
            os.chdir(old_dir)

        # ---- Build element → per-pass delta list ----
        self._pec_atoms = [z for z in range(1, NELEM + 1)
                           if self._state.pec[z - 1] == 1]

        # pecabund[z-1, k] is a log offset from xabu (the model abundance).
        # The slider delta IS that log offset directly — no conversion needed.
        self._delta_abund = {}
        for z in self._pec_atoms:
            self._delta_abund[z] = [
                self._state.pecabund[z - 1, k]
                for k in range(self._npass)
            ]

        # Widget state
        self._active_pass  = 0
        self._current_page = 0
        self._dlambda      = 0.0
        self._continuum    = 1.0     # multiplicative continuum scale (display only)
        self._needs_synthesis = False
        self._add_mode     = False
        self._zoom_locked  = False   # True once user sets explicit axis limits
        self._suppress_zoom = False  # Guard against set_val triggering callbacks
        self._suppress_abund_tb  = False  # Guard against circular abund TB updates
        self._suppress_smooth_tb = False  # Guard against circular smooth TB updates

        # Snapshot of initial smoothing values (from batch.par) for Reset
        st = self._state
        self._init_smooth = {
            'fwhmgauss': st.fwhmgauss,
            'vmac':       st.vmac,
            'vsini':      st.vsini,
            'limbdark':   st.limbdark,
            'fwhmloren':  st.fwhmloren,
            '_dlambda':   0.0,
            '_continuum': 1.0,
        }

        self._build_figure()

    # ------------------------------------------------------------------ #
    # Synthesis                                                            #
    # ------------------------------------------------------------------ #

    def _run_synthesis(self):
        """Run all synthesis passes and store raw depth arrays."""
        st = self._state
        waves, depths = [], []

        if self._npass <= 1 or st.numpecatom == 0:
            st.isynth  = 1
            st.isorun  = 1
            st.nlines  = 0
            st.waveold = 0.0
            _inlines(st, 1)
            _eqlib(st)
            _nearly(st, 1)
            w, d = _synspec(st)
            waves.append(w)
            depths.append(d)
        else:
            for n in range(self._npass):
                st.isynth  = n + 1
                st.isorun  = n + 1
                st.start   = st.oldstart
                st.sstop   = st.oldstop
                st.mode    = 3
                st.waveold = 0.0
                _inlines(st, 1)
                _eqlib(st)
                _nearly(st, 1)
                w, d = _synspec(st)
                waves.append(w)
                depths.append(d)

        self._wave       = waves[0]
        self._raw_depths = depths
        self._resmooth()

    @staticmethod
    def _derive_smtype(gauss: bool, rot: bool, mac: bool, lor: bool) -> str:
        """Map checkbox states to the closest MOOG smtype code (for batch.par saving)."""
        if lor and not (gauss or rot or mac):
            return 'l'
        if rot and mac and gauss:
            return 'r'
        if rot and gauss:
            return 'c'
        if mac and gauss:
            return 'd'
        if rot:
            return 'v'
        if mac:
            return 'm'
        if gauss:
            return 'g'
        return 'n'

    def _smooth_flags(self):
        """Return (gauss, rot, mac, lor) booleans from the CheckButtons state."""
        if hasattr(self, '_check_smooth'):
            return tuple(self._check_smooth.get_status())
        # Fallback during __init__ before the widget is built
        sm = self._state.smtype
        return (sm in ('g', 'c', 'd', 'r'),
                sm in ('v', 'c', 'r'),
                sm in ('m', 'd', 'r'),
                sm in ('l',))

    def _resmooth(self):
        """Apply current smoothing parameters to stored raw depths."""
        gauss, rot, mac, lor = self._smooth_flags()
        st = self._state
        self._flux_smooth = [
            _smooth(
                self._wave, 1.0 - d, st.step,
                vsini     = st.vsini     if rot   else 0.0,
                limbdark  = st.limbdark  if rot   else 0.0,
                vmac      = st.vmac      if mac   else 0.0,
                fwhmgauss = st.fwhmgauss if gauss else 0.0,
                fwhmloren = st.fwhmloren if lor   else 0.0,
                addflux   = st.addflux,
            )
            for d in self._raw_depths
        ]

    # ------------------------------------------------------------------ #
    # Figure layout                                                        #
    # ------------------------------------------------------------------ #

    def _build_figure(self):
        has_obs = self._obs_wave is not None

        self._fig = plt.figure('pymoog synth', figsize=(14, 9))
        self._fig.subplots_adjust(left=0, bottom=0, right=1, top=1)

        # Spectrum (and optional residuals).
        # Bottom of the lowest axes must sit high enough to leave room for x-axis
        # tick labels above the pass-button / zoom strip (top ≈ 0.490).
        if has_obs:
            self._ax_spec  = self._fig.add_axes([0.07, 0.62, 0.90, 0.35])
            self._ax_resid = self._fig.add_axes([0.07, 0.55, 0.90, 0.055],
                                                sharex=self._ax_spec)
            self._ax_resid.set_ylabel('Obs − Syn', fontsize=8)
            self._ax_resid.axhline(0, color='k', lw=0.5, ls='--')
            self._ax_resid.set_xlabel('Wavelength (Å)', fontsize=9)
            self._ax_resid.tick_params(labelsize=8)
            self._ax_spec.tick_params(labelbottom=False)   # hide duplicate x labels
        else:
            self._ax_spec  = self._fig.add_axes([0.07, 0.56, 0.90, 0.41])
            self._ax_resid = None
            self._ax_spec.set_xlabel('Wavelength (Å)', fontsize=9)
            self._ax_spec.tick_params(labelsize=8)

        self._ax_spec.set_ylabel('Relative flux', fontsize=9)
        self._ax_spec.tick_params(axis='y', labelsize=8)

        # Cursor readout — updated by motion_notify_event
        self._cursor_txt = self._ax_spec.text(
            0.01, 0.97, '', transform=self._ax_spec.transAxes,
            fontsize=8, va='top', ha='left',
            bbox=dict(facecolor='white', alpha=0.75, edgecolor='none', pad=2),
            zorder=10,
        )
        self._fig.canvas.mpl_connect('motion_notify_event', self._on_mouse_move)

        # Pass selector buttons (only if npass > 1)
        self._pass_btns = []
        if self._npass > 1:
            bw = min(0.07, 0.75 / self._npass)
            for k in range(self._npass):
                ax_pb = self._fig.add_axes(
                    [0.07 + k * (bw + 0.005), 0.462, bw, 0.028]
                )
                btn = Button(ax_pb, f'Pass {k + 1}',
                             color=_PASS_COLORS[k % len(_PASS_COLORS)],
                             hovercolor='0.80')
                btn.label.set_fontsize(8)
                btn.on_clicked(lambda _e, k=k: self._on_pass_select(k))
                self._pass_btns.append(btn)
            self._fig.text(
                0.07 + self._npass * (bw + 0.005) + 0.01, 0.474,
                'editing →', fontsize=8, va='center', color='0.40'
            )

        # ---- Zoom controls (λ range and flux range) ----
        # Placed in the right portion of the pass-button row.
        # Gaps chosen so text labels never touch adjacent box borders.
        _zy = 0.462
        _zh = 0.027
        _zc = _zy + _zh / 2   # vertical centre for text labels
        self._fig.text(0.509, _zc, 'λ:', fontsize=8, ha='right', va='center')
        ax_xmin = self._fig.add_axes([0.513, _zy, 0.062, _zh])   # ends 0.575
        self._fig.text(0.582, _zc, '–',  fontsize=8, ha='center', va='center')
        ax_xmax = self._fig.add_axes([0.590, _zy, 0.062, _zh])   # ends 0.652
        self._fig.text(0.670, _zc, 'F:', fontsize=8, ha='right', va='center')
        ax_ymin = self._fig.add_axes([0.674, _zy, 0.050, _zh])   # ends 0.724
        self._fig.text(0.731, _zc, '–',  fontsize=8, ha='center', va='center')
        ax_ymax = self._fig.add_axes([0.738, _zy, 0.050, _zh])   # ends 0.788
        ax_zrst = self._fig.add_axes([0.798, _zy, 0.085, _zh])   # ends 0.883

        st0 = self._state
        self._tb_xmin = TextBox(ax_xmin, '', initial=f'{st0.start:.2f}')
        self._tb_xmax = TextBox(ax_xmax, '', initial=f'{st0.sstop:.2f}')
        self._tb_ymin = TextBox(ax_ymin, '', initial='0.000')
        self._tb_ymax = TextBox(ax_ymax, '', initial='1.050')
        self._btn_zoom_reset = Button(ax_zrst, 'Reset zoom',
                                      color='0.85', hovercolor='0.70')

        for _tb in (self._tb_xmin, self._tb_xmax, self._tb_ymin, self._tb_ymax):
            _tb.text_disp.set_fontsize(8)
        self._btn_zoom_reset.label.set_fontsize(8)

        self._tb_xmin.on_submit(self._on_zoom_submit)
        self._tb_xmax.on_submit(self._on_zoom_submit)
        self._tb_ymin.on_submit(self._on_zoom_submit)
        self._tb_ymax.on_submit(self._on_zoom_submit)
        self._btn_zoom_reset.on_clicked(lambda _e: self._on_reset_zoom())

        # ---- Smoothing section ----
        self._fig.text(0.07, 0.448, 'Smoothing', fontsize=9, fontweight='bold')

        sm = self._state.smtype
        init_checks = [
            sm in ('g', 'c', 'd', 'r'),   # Gaussian
            sm in ('v', 'c', 'r'),          # Rotation
            sm in ('m', 'd', 'r'),          # Macro
            sm in ('l',),                   # Lorentz
        ]
        ax_smtype = self._fig.add_axes([0.07, 0.305, 0.18, 0.125])
        self._check_smooth = CheckButtons(ax_smtype, _SMOOTH_LABELS, init_checks)
        for lbl in self._check_smooth.labels:
            lbl.set_fontsize(8)
        self._check_smooth.on_clicked(self._on_smooth_check)

        # Smoothing parameter sliders.
        # Labels: fig.text() right-aligned just left of each slider axes.
        # Value text: moved INSIDE each slider axes so it never spills into
        # adjacent columns.  White bbox keeps it readable over the slider bar.
        # Left col: slider 0.40–0.53, TB 0.535–0.595
        # Right col: slider 0.72–0.85, TB 0.855–0.915  (label at xl=0.715)
        _sliders_cfg = [
            # (rect,                         fig-label,       xl,    yl,     attr,       vmin,  vmax,  fmt)
            ([0.40, 0.415, 0.13, 0.022], 'FWHM Gauss (Å)', 0.395, 0.4255, 'fwhmgauss',  0.0,   2.0, '%.3f'),
            ([0.72, 0.415, 0.13, 0.022], 'vmac (km/s)',     0.715, 0.4255, 'vmac',        0.0,  10.0, '%.2f'),
            ([0.40, 0.375, 0.13, 0.022], 'vsini (km/s)',    0.395, 0.3855, 'vsini',       0.0,  50.0, '%.1f'),
            ([0.72, 0.375, 0.13, 0.022], 'limb dark',       0.715, 0.3855, 'limbdark',    0.0,   1.0, '%.2f'),
            ([0.40, 0.335, 0.13, 0.022], 'Δλ shift (Å)',    0.395, 0.3455, '_dlambda',   -2.0,   2.0, '%+.3f'),
            ([0.72, 0.335, 0.13, 0.022], 'FWHM Loren (Å)', 0.715, 0.3455, 'fwhmloren',   0.0,   2.0, '%.3f'),
            ([0.40, 0.295, 0.13, 0.022], 'continuum',       0.395, 0.3055, '_continuum',  0.80,  1.20, '%.3f'),
        ]
        self._smooth_sliders = {}
        self._smooth_tb      = {}   # editable TextBoxes for each smoothing param
        self._smooth_fmt     = {}   # format strings keyed by attr
        for rect, fig_lbl, xl, yl, attr, vmin, vmax, fmt in _sliders_cfg:
            if attr == '_dlambda':
                valinit = self._dlambda
            elif attr == '_continuum':
                valinit = self._continuum
            else:
                valinit = getattr(self._state, attr, 0.0)
            ax_sl = self._fig.add_axes(rect)
            sl = Slider(ax_sl, '', vmin, vmax, valinit=valinit, valfmt=fmt)
            sl.valtext.set_visible(False)   # replaced by editable TextBox
            sl.on_changed(lambda val, a=attr: self._on_smooth_slider(a, val))
            self._smooth_sliders[attr] = sl
            self._smooth_fmt[attr]     = fmt
            self._fig.text(xl, yl, fig_lbl, fontsize=8, ha='right', va='center')
            # Editable TextBox immediately to the right of the slider
            tb_left = rect[0] + rect[2] + 0.005
            ax_tb = self._fig.add_axes([tb_left, rect[1], 0.065, rect[3]])
            tb = TextBox(ax_tb, '', initial=fmt % valinit, textalignment='center')
            tb.label.set_fontsize(8)
            tb.text_disp.set_fontsize(8)
            tb.on_submit(lambda text, a=attr, lo=vmin, hi=vmax: self._on_smooth_tb(a, text, lo, hi))
            self._smooth_tb[attr] = tb

        # ---- Abundance section ----
        self._fig.text(0.07, 0.284, 'Abundances', fontsize=9, fontweight='bold')

        ax_prev = self._fig.add_axes([0.82, 0.280, 0.040, 0.022])
        ax_next = self._fig.add_axes([0.87, 0.280, 0.040, 0.022])
        self._btn_prev = Button(ax_prev, '◀', color='0.85', hovercolor='0.70')
        self._btn_next = Button(ax_next, '▶', color='0.85', hovercolor='0.70')
        self._btn_prev.label.set_fontsize(8)
        self._btn_next.label.set_fontsize(8)
        self._btn_prev.on_clicked(lambda _e: self._on_page(-1))
        self._btn_next.on_clicked(lambda _e: self._on_page(+1))
        self._page_text = self._fig.text(0.76, 0.287, '', fontsize=8, va='center',
                                         color='0.40')
        # Hidden until there is more than one page of elements
        ax_prev.set_visible(False)
        ax_next.set_visible(False)

        # Element slider slots
        self._slot_axes      = []
        self._slot_sl        = []
        self._slot_texts     = []   # kept for compatibility (unused after refactor)
        self._slot_tb        = []   # editable abundance TextBoxes
        self._slot_tb_axes   = []
        self._slot_name_texts = []  # left-side element name labels
        self._slot_rm_axes   = []   # × remove-element buttons
        self._slot_rm_btn    = []
        slot_h   = 0.030
        slot_top = 0.268

        for i in range(_MAX_SLOTS):
            bot  = slot_top - (i + 1) * slot_h
            # Element name label (to the left of slider)
            name_txt = self._fig.text(0.08, bot + slot_h * 0.45, '',
                                      fontsize=9, va='center', ha='left',
                                      fontweight='bold')
            # Delta-abundance slider (shortened to leave room for entry box)
            ax_sl = self._fig.add_axes([0.18, bot + 0.005, 0.50, 0.018])
            sl = Slider(ax_sl, '', -1.5, 1.5, valinit=0.0, valfmt='%+.3f')
            sl.label.set_fontsize(8)
            sl.valtext.set_transform(ax_sl.transAxes)
            sl.valtext.set_position((0.98, 0.5))
            sl.valtext.set_ha('right')
            sl.valtext.set_va('center')
            sl.valtext.set_fontsize(8)
            sl.valtext.set_bbox(dict(facecolor='white', alpha=0.75, edgecolor='none', pad=1))
            sl.on_changed(lambda val, idx=i: self._on_abund_slider(idx, val))
            # Editable absolute log ε TextBox (to the right of slider)
            ax_tb = self._fig.add_axes([0.695, bot + 0.003, 0.105, 0.022])
            tb = TextBox(ax_tb, 'ε=', initial='', textalignment='left')
            tb.label.set_fontsize(8)
            tb.text_disp.set_fontsize(8)
            tb.on_submit(lambda text, idx=i: self._on_abund_tb(idx, text))
            # × remove-element button (right of TextBox)
            ax_rm = self._fig.add_axes([0.812, bot + 0.003, 0.030, 0.022])
            btn_rm = Button(ax_rm, '×', color='0.85', hovercolor='0.70')
            btn_rm.label.set_fontsize(8)
            btn_rm.on_clicked(lambda _e, idx=i: self._on_remove_element(idx))
            ax_rm.set_visible(False)
            self._slot_axes.append(ax_sl)
            self._slot_sl.append(sl)
            self._slot_texts.append(None)     # placeholder (unused)
            self._slot_tb.append(tb)
            self._slot_tb_axes.append(ax_tb)
            self._slot_name_texts.append(name_txt)
            self._slot_rm_axes.append(ax_rm)
            self._slot_rm_btn.append(btn_rm)

        # ---- Button row (five buttons) ----
        ax_syn  = self._fig.add_axes([0.07, 0.010, 0.20, 0.038])
        ax_rst  = self._fig.add_axes([0.28, 0.010, 0.10, 0.038])
        ax_sav  = self._fig.add_axes([0.39, 0.010, 0.13, 0.038])
        ax_ssp  = self._fig.add_axes([0.53, 0.010, 0.13, 0.038])
        ax_add  = self._fig.add_axes([0.67, 0.010, 0.14, 0.038])
        # Add-element input row — positioned ABOVE the button row, no overlap
        self._ax_tb   = self._fig.add_axes([0.10, 0.057, 0.40, 0.032])
        self._ax_conf = self._fig.add_axes([0.51, 0.057, 0.10, 0.032])

        self._btn_synth   = Button(ax_syn,  'Synthesize',    color='0.85', hovercolor='0.70')
        self._btn_reset   = Button(ax_rst,  'Reset',         color='0.85', hovercolor='0.70')
        self._btn_save    = Button(ax_sav,  'Save batch.par', color='0.85', hovercolor='0.70')
        self._btn_savesp  = Button(ax_ssp,  'Save spectra',   color='0.85', hovercolor='0.70')
        self._btn_add     = Button(ax_add,  '+ Add element',  color='0.85', hovercolor='0.70')
        self._tb_element  = TextBox(self._ax_tb,  'Element symbol: ')
        self._btn_confirm = Button(self._ax_conf, 'Add',        color='0.85', hovercolor='0.70')

        self._ax_tb.set_visible(False)
        self._ax_conf.set_visible(False)

        for btn in (self._btn_synth, self._btn_reset, self._btn_save,
                    self._btn_savesp, self._btn_add, self._btn_confirm):
            btn.label.set_fontsize(9)

        self._btn_synth.on_clicked(lambda _e: self._on_synthesize())
        self._btn_reset.on_clicked(lambda _e: self._on_reset())
        self._btn_save.on_clicked(lambda _e: self._on_save())
        self._btn_savesp.on_clicked(lambda _e: self._on_save_spectra())
        self._btn_add.on_clicked(lambda _e: self._on_add_element())
        self._btn_confirm.on_clicked(lambda _e: self._on_confirm_add())
        self._tb_element.on_submit(lambda _t: self._on_confirm_add())

        # ---- Initial plot ----
        self._lines_syn   = []
        self._line_obs    = None
        self._lines_resid = []
        self._init_plot()
        self._refresh_slots()
        self._update_plot()

    # ------------------------------------------------------------------ #
    # Plot management                                                      #
    # ------------------------------------------------------------------ #

    def _init_plot(self):
        for k in range(self._npass):
            lbl = f'Pass {k + 1}' if self._npass > 1 else 'Synth'
            ln, = self._ax_spec.plot(
                [], [], color=_PASS_COLORS[k % len(_PASS_COLORS)],
                lw=1.2, label=lbl
            )
            self._lines_syn.append(ln)

        if self._obs_wave is not None:
            self._line_obs, = self._ax_spec.plot(
                self._obs_wave, self._obs_flux,
                color='k', lw=1.0, alpha=0.65, label='Observed', zorder=0
            )
            for k in range(self._npass):
                lr, = self._ax_resid.plot(
                    [], [], color=_PASS_COLORS[k % len(_PASS_COLORS)], lw=1.0
                )
                self._lines_resid.append(lr)

        if self._npass > 1 or self._obs_wave is not None:
            self._ax_spec.legend(fontsize=8, loc='lower right')

    def _update_plot(self):
        """Redraw all synthetic spectrum lines with current smoothing and Δλ."""
        wave_shifted = self._wave + self._dlambda

        for k, ln in enumerate(self._lines_syn):
            flux = self._flux_smooth[k] if k < len(self._flux_smooth) else self._flux_smooth[0]
            ln.set_xdata(wave_shifted)
            ln.set_ydata(np.asarray(flux) * self._continuum)

        if self._obs_wave is not None and self._ax_resid is not None:
            for k, lr in enumerate(self._lines_resid):
                flux_k = self._flux_smooth[k] if k < len(self._flux_smooth) else self._flux_smooth[0]
                syn_i  = np.interp(self._obs_wave, wave_shifted,
                                   np.asarray(flux_k) * self._continuum)
                lr.set_xdata(self._obs_wave)
                lr.set_ydata(self._obs_flux - syn_i)

        if not self._zoom_locked:
            self._ax_spec.relim()
            self._ax_spec.autoscale_view()
            if self._ax_resid:
                self._ax_resid.relim()
                self._ax_resid.autoscale_view()
            self._update_zoom_textboxes()

        self._fig.canvas.draw_idle()

    # ------------------------------------------------------------------ #
    # Element slider slots                                                 #
    # ------------------------------------------------------------------ #

    def _refresh_slots(self):
        """Populate element slider slots for current page and active pass."""
        n_total = len(self._pec_atoms)
        n_pages = max(1, math.ceil(n_total / _MAX_SLOTS))
        self._current_page = min(self._current_page, max(0, n_pages - 1))

        multipage = n_pages > 1
        self._btn_prev.ax.set_visible(multipage)
        self._btn_next.ax.set_visible(multipage)
        self._page_text.set_text(
            f'pg {self._current_page + 1}/{n_pages}' if multipage else ''
        )

        page_start = self._current_page * _MAX_SLOTS
        page_atoms = self._pec_atoms[page_start: page_start + _MAX_SLOTS]

        self._suppress_abund_tb = True
        for i in range(_MAX_SLOTS):
            ax_sl    = self._slot_axes[i]
            sl       = self._slot_sl[i]
            tb       = self._slot_tb[i]
            ax_tb    = self._slot_tb_axes[i]
            ax_rm    = self._slot_rm_axes[i]
            name_txt = self._slot_name_texts[i]

            if i < len(page_atoms):
                z     = page_atoms[i]
                sym   = ELEMENT_NAMES[z - 1].strip()
                delta = self._delta_abund[z][self._active_pass]
                base  = self._model_logeps.get(z, 0.0)

                ax_sl.set_visible(True)
                ax_tb.set_visible(True)
                ax_rm.set_visible(True)
                name_txt.set_text(sym)

                sl.eventson = False
                sl.set_val(delta)
                sl.eventson = True

                tb.set_val(f'{base + delta:.3f}')
            else:
                ax_sl.set_visible(False)
                ax_tb.set_visible(False)
                ax_rm.set_visible(False)
                name_txt.set_text('')
                tb.set_val('')
        self._suppress_abund_tb = False

        self._fig.canvas.draw_idle()

    # ------------------------------------------------------------------ #
    # Callbacks                                                            #
    # ------------------------------------------------------------------ #

    # ------------------------------------------------------------------ #
    # Zoom helpers                                                         #
    # ------------------------------------------------------------------ #

    def _update_zoom_textboxes(self):
        """Sync TextBox text to current axis limits without triggering zoom."""
        if not hasattr(self, '_tb_xmin'):
            return
        xl = self._ax_spec.get_xlim()
        yl = self._ax_spec.get_ylim()
        self._suppress_zoom = True
        self._tb_xmin.set_val(f'{xl[0]:.2f}')
        self._tb_xmax.set_val(f'{xl[1]:.2f}')
        self._tb_ymin.set_val(f'{yl[0]:.3f}')
        self._tb_ymax.set_val(f'{yl[1]:.3f}')
        self._suppress_zoom = False

    def _on_zoom_submit(self, _text: str):
        """Apply wavelength and flux limits typed in the TextBoxes (Enter key)."""
        if self._suppress_zoom:
            return
        try:
            xmin = float(self._tb_xmin.text)
            xmax = float(self._tb_xmax.text)
            ymin = float(self._tb_ymin.text)
            ymax = float(self._tb_ymax.text)
        except ValueError:
            return
        if xmin >= xmax or ymin >= ymax:
            return
        self._zoom_locked = True
        self._ax_spec.set_xlim(xmin, xmax)
        self._ax_spec.set_ylim(ymin, ymax)
        # ax_resid shares x via sharex, so xlim propagates automatically
        self._fig.canvas.draw_idle()

    def _on_reset_zoom(self):
        """Restore autoscaled limits and sync TextBoxes."""
        self._zoom_locked = False
        # set_xlim/set_ylim disable the internal autoscale flag — re-enable first.
        self._ax_spec.set_autoscale_on(True)
        self._ax_spec.relim()
        self._ax_spec.autoscale_view()
        if self._ax_resid:
            self._ax_resid.set_autoscale_on(True)
            self._ax_resid.relim()
            self._ax_resid.autoscale_view()
        self._update_zoom_textboxes()
        self._fig.canvas.draw_idle()

    def _on_mouse_move(self, event):
        if event.inaxes is self._ax_spec and event.xdata is not None:
            self._cursor_txt.set_text(f'λ = {event.xdata:.3f} Å    F = {event.ydata:.4f}')
        else:
            self._cursor_txt.set_text('')
        self._fig.canvas.draw_idle()

    def _on_smooth_check(self, _label: str):
        gauss, rot, mac, lor = self._smooth_flags()
        self._state.smtype = self._derive_smtype(gauss, rot, mac, lor)
        self._resmooth()
        self._update_plot()

    def _on_smooth_slider(self, attr: str, val: float):
        if attr == '_dlambda':
            self._dlambda = val
        elif attr == '_continuum':
            self._continuum = val
        else:
            setattr(self._state, attr, val)
            self._resmooth()
        # Keep TextBox in sync with slider
        if attr in self._smooth_tb:
            self._suppress_smooth_tb = True
            self._smooth_tb[attr].set_val(self._smooth_fmt[attr] % val)
            self._suppress_smooth_tb = False
        self._update_plot()

    def _on_smooth_tb(self, attr: str, text: str, vmin: float, vmax: float):
        """User typed a value into a smoothing TextBox — update the slider."""
        if self._suppress_smooth_tb:
            return
        try:
            val = float(text)
        except ValueError:
            return
        val = max(vmin, min(vmax, val))
        self._smooth_sliders[attr].set_val(val)

    def _on_pass_select(self, k: int):
        self._active_pass = k
        for i, btn in enumerate(self._pass_btns):
            btn.ax.set_facecolor(
                _PASS_COLORS[i % len(_PASS_COLORS)] if i == k else '0.85'
            )
        self._refresh_slots()

    def _on_abund_slider(self, slot_idx: int, val: float):
        page_start = self._current_page * _MAX_SLOTS
        atom_idx   = page_start + slot_idx
        if atom_idx >= len(self._pec_atoms):
            return
        z = self._pec_atoms[atom_idx]
        self._delta_abund[z][self._active_pass] = val

        base = self._model_logeps.get(z, 0.0)
        self._suppress_abund_tb = True
        self._slot_tb[slot_idx].set_val(f'{base + val:.3f}')
        self._suppress_abund_tb = False

        if not self._needs_synthesis:
            self._needs_synthesis = True
            self._btn_synth.ax.set_facecolor(_STALE_COLOR)
            self._btn_synth.label.set_color('white')

        self._fig.canvas.draw_idle()

    def _on_abund_tb(self, slot_idx: int, text: str):
        """User typed an absolute log ε value into the TextBox — update the slider."""
        if self._suppress_abund_tb:
            return
        try:
            val = float(text)
        except ValueError:
            return
        page_start = self._current_page * _MAX_SLOTS
        atom_idx   = page_start + slot_idx
        if atom_idx >= len(self._pec_atoms):
            return
        z    = self._pec_atoms[atom_idx]
        base = self._model_logeps.get(z, 0.0)
        delta = val - base
        sl = self._slot_sl[slot_idx]
        delta = max(sl.valmin, min(sl.valmax, delta))
        # Updating the slider triggers _on_abund_slider, which updates the TextBox
        sl.set_val(delta)

    def _on_synthesize(self):
        """Push abundance deltas into state and re-run all synthesis passes."""
        for z in self._pec_atoms:
            for k in range(self._npass):
                self._state.pecabund[z - 1, k] = self._delta_abund[z][k]
            self._state.pec[z - 1] = 1

        old_dir = os.getcwd()
        os.chdir(self._rundir)
        try:
            self._run_synthesis()
        finally:
            os.chdir(old_dir)

        self._update_plot()

        self._needs_synthesis = False
        self._btn_synth.ax.set_facecolor('0.85')
        self._btn_synth.label.set_color('black')
        self._fig.canvas.draw_idle()

    def _on_reset(self):
        """Zero all abundance deltas and restore smoothing to batch.par values."""
        for z in self._pec_atoms:
            self._delta_abund[z] = [0.0] * self._npass
        self._refresh_slots()

        # Restore each smoothing slider; set_val triggers _on_smooth_slider which
        # updates state, resmooths, syncs the TextBox, and redraws.
        for attr, val in self._init_smooth.items():
            if attr in self._smooth_sliders:
                self._smooth_sliders[attr].set_val(val)

        self._on_synthesize()

    def _on_save(self):
        """Write current abundances back to batch.par."""
        try:
            self._save_batch_par()
            self._btn_save.label.set_text('Saved ✓')
        except Exception as exc:
            msg = str(exc)[:18]
            self._btn_save.label.set_text(f'Err: {msg}')
        self._fig.canvas.draw_idle()

    def _on_save_spectra(self):
        """Write abundances summary and smoothed spectra files."""
        try:
            self._save_results()
            self._btn_savesp.label.set_text('Saved ✓')
        except Exception as exc:
            self._btn_savesp.label.set_text(f'Err: {str(exc)[:14]}')
        self._fig.canvas.draw_idle()

    def _save_results(self):
        """
        Write files to the run directory:
          <base>_abund.txt  — per-pass absolute log ε for every peculiar element
          <base>_spec.dat   — wavelength + one flux column per pass
          <base>_obs.dat    — observed wavelength + flux (if observed spectrum loaded)
          <base>_plot.py    — stand-alone script that reproduces the publication figure
        <base> is derived from smoothed_out (f3out) or the batch.par stem.
        """
        f3 = getattr(self._state, 'f3out', '') or ''
        if f3 and f3 != _NO_FILE:
            base = os.path.splitext(f3)[0]
        else:
            base = os.path.splitext(os.path.basename(self._fparam))[0]

        abund_path = os.path.join(self._rundir, base + '_abund.txt')
        spec_path  = os.path.join(self._rundir, base + '_spec.dat')
        obs_path   = os.path.join(self._rundir, base + '_obs.dat')
        plot_path  = os.path.join(self._rundir, base + '_plot.py')

        # ---- Abundances summary ----
        with open(abund_path, 'w') as fh:
            fh.write(f'# Synthesis abundances: {os.path.basename(self._fparam)}\n')
            fh.write(f'# {"Pass":>4}  {"Element":>7}  {"delta(logeps)":>13}  {"logeps":>8}\n')
            for k in range(self._npass):
                for z in self._pec_atoms:
                    sym   = ELEMENT_NAMES[z - 1].strip()
                    delta = self._delta_abund[z][k]
                    abso  = self._model_logeps.get(z, 0.0) + delta
                    fh.write(f'  {k + 1:4d}  {sym:>7}  {delta:+13.3f}  {abso:8.3f}\n')

        # ---- Smoothed synthetic spectra ----
        wave_shifted = self._wave + self._dlambda
        with open(spec_path, 'w') as fh:
            pass_hdrs = ''.join(f'{"flux_p" + str(k + 1):>11}' for k in range(self._npass))
            fh.write(f'# Synthesized spectra: {os.path.basename(self._fparam)}\n')
            fh.write(f'# {"wavelength":>10}{pass_hdrs}\n')
            for i, w in enumerate(wave_shifted):
                fluxes = ''.join(
                    f'{(np.asarray(self._flux_smooth[k] if k < len(self._flux_smooth) else self._flux_smooth[0])[i]) * self._continuum:11.5f}'
                    for k in range(self._npass)
                )
                fh.write(f'  {w:10.3f}{fluxes}\n')

        # ---- Observed spectrum (if present) ----
        has_obs = self._obs_wave is not None and self._obs_flux is not None
        if has_obs:
            with open(obs_path, 'w') as fh:
                fh.write(f'# Observed spectrum: {os.path.basename(self._fparam)}\n')
                fh.write(f'# {"wavelength":>10}  {"flux":>10}\n')
                for w, f in zip(self._obs_wave, self._obs_flux):
                    fh.write(f'  {w:10.3f}  {f:10.5f}\n')

        # ---- Publication plotting script ----
        # Build per-pass label strings (e.g. "log ε(Eu) = 0.52")
        pass_labels = []
        for k in range(self._npass):
            parts = []
            for z in self._pec_atoms:
                sym  = ELEMENT_NAMES[z - 1].strip()
                abso = self._model_logeps.get(z, 0.0) + self._delta_abund[z][k]
                parts.append(f'log ε({sym}) = {abso:.2f}')
            pass_labels.append(', '.join(parts) if parts else f'Pass {k + 1}')

        # Current axis limits as defaults
        xl = self._ax_spec.get_xlim()
        yl = self._ax_spec.get_ylim()
        rl = self._ax_resid.get_ylim() if self._ax_resid else (-0.05, 0.05)

        colors_repr = "['#377eb8', '#e41a1c', '#4daf4a', '#984ea3', '#ff7f00']"

        obs_block = (
            f"obs      = np.loadtxt('{os.path.basename(obs_path)}', comments='#')\n"
            f"obs_wave = obs[:, 0]\n"
            f"obs_flux = obs[:, 1]\n"
        ) if has_obs else (
            "# No observed spectrum was loaded; comment in the lines below if you add one\n"
            "# obs      = np.loadtxt('<obs_file>', comments='#')\n"
            "# obs_wave = obs[:, 0]\n"
            "# obs_flux = obs[:, 1]\n"
        )

        obs_plot_block = (
            "ax.plot(obs_wave, obs_flux, color='k', lw=0.8, alpha=0.8,\n"
            "        label='Observed', zorder=0)\n"
        ) if has_obs else (
            "# ax.plot(obs_wave, obs_flux, color='k', lw=0.8, alpha=0.8,\n"
            "#         label='Observed', zorder=0)\n"
        )

        resid_block = (
            "for k, col in enumerate(COLORS[:n_pass]):\n"
            "    syn_i = np.interp(obs_wave, syn_wave, syn_flux[:, k])\n"
            "    axr.plot(obs_wave, obs_flux - syn_i, color=col, lw=0.8)\n"
            "axr.axhline(0, color='k', lw=0.5, ls='--')\n"
        ) if has_obs else (
            "# Residual panel requires an observed spectrum\n"
            "# axr.set_visible(False)\n"
        )

        script = f"""\
#!/usr/bin/env python3
\"\"\"
Publication figure for synthesis fit: {os.path.basename(self._fparam)}
Generated by SPAE interactive synthesis widget.
Run: python {os.path.basename(plot_path)}
\"\"\"
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ── Style ──────────────────────────────────────────────────────────────────
mpl.rcParams.update({{
    'font.family':       'serif',
    'font.size':         10,
    'axes.linewidth':    0.8,
    'xtick.direction':   'in',
    'ytick.direction':   'in',
    'xtick.top':         True,
    'ytick.right':       True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'legend.framealpha': 0.9,
    'legend.fontsize':   8,
}})

# ── Configuration (edit these) ─────────────────────────────────────────────
XLIM   = ({xl[0]:.3f}, {xl[1]:.3f})   # wavelength range (Å)
YLIM   = ({yl[0]:.4f}, {yl[1]:.4f})   # flux range
RLIM   = ({rl[0]:.4f}, {rl[1]:.4f})   # residual range
COLORS = {colors_repr}
OUTFILE = '{base}_figure.pdf'          # set to None to display interactively

# ── Pass labels (log ε for each synthesized pass) ──────────────────────────
PASS_LABELS = [
{''.join(f"    {repr(lbl)},{'  # pass ' + str(k+1)}{chr(10)}" for k, lbl in enumerate(pass_labels))}\
]

# ── Load data ──────────────────────────────────────────────────────────────
syn      = np.loadtxt('{os.path.basename(spec_path)}', comments='#')
syn_wave = syn[:, 0]
syn_flux = syn[:, 1:]          # shape (n_pixels, n_pass)
n_pass   = syn_flux.shape[1]

{obs_block}
# ── Build figure ───────────────────────────────────────────────────────────
fig, (ax, axr) = plt.subplots(
    2, 1, figsize=(7, 5), sharex=True,
    gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05),
)

{obs_plot_block}
for k, col in enumerate(COLORS[:n_pass]):
    ax.plot(syn_wave, syn_flux[:, k], color=col, lw=1.2,
            label=PASS_LABELS[k] if k < len(PASS_LABELS) else f'Pass {{k+1}}')

ax.set_xlim(*XLIM)
ax.set_ylim(*YLIM)
ax.set_ylabel('Relative flux')
ax.tick_params(labelbottom=False)
ax.legend(loc='lower right')

{resid_block}
axr.set_xlim(*XLIM)
axr.set_ylim(*RLIM)
axr.set_xlabel('Wavelength (Å)')
axr.set_ylabel('Obs − Syn')

if OUTFILE:
    fig.savefig(OUTFILE, bbox_inches='tight', dpi=300)
    print(f'Saved {{OUTFILE}}')
else:
    plt.show()
"""
        with open(plot_path, 'w') as fh:
            fh.write(script)

    def _on_add_element(self):
        """Toggle the add-element input row above the button bar."""
        self._add_mode = not self._add_mode
        self._ax_tb.set_visible(self._add_mode)
        self._ax_conf.set_visible(self._add_mode)
        if self._add_mode:
            self._btn_add.label.set_text('✕ Cancel')
            self._tb_element.set_val('')
        else:
            self._btn_add.label.set_text('+ Add element')
        self._fig.canvas.draw_idle()

    def _on_confirm_add(self):
        """Add the typed element to the peculiar atom list."""
        name  = self._tb_element.text.strip()
        iatom = _iatom_from_name(name)
        if iatom < 1:
            self._tb_element.set_val('? not found')
            return

        if iatom not in self._pec_atoms:
            self._pec_atoms.append(iatom)
            self._delta_abund[iatom] = [0.0] * self._npass
            self._state.pec[iatom - 1] = 1
            self._state.numpecatom += 1
            for k in range(self._npass):
                self._state.pecabund[iatom - 1, k] = 0.0  # zero log offset = xabu unchanged
            idx = self._pec_atoms.index(iatom)
            self._current_page = idx // _MAX_SLOTS

        # Close add-element mode
        self._add_mode = False
        self._ax_tb.set_visible(False)
        self._ax_conf.set_visible(False)
        self._btn_add.label.set_text('+ Add element')
        self._refresh_slots()

        if not self._needs_synthesis:
            self._needs_synthesis = True
            self._btn_synth.ax.set_facecolor(_STALE_COLOR)
            self._btn_synth.label.set_color('white')
        self._fig.canvas.draw_idle()

    def _on_remove_element(self, slot_idx: int):
        """Remove the element in the given slot from the peculiar atom list."""
        page_start = self._current_page * _MAX_SLOTS
        atom_idx   = page_start + slot_idx
        if atom_idx >= len(self._pec_atoms):
            return
        z = self._pec_atoms.pop(atom_idx)
        del self._delta_abund[z]
        self._state.pec[z - 1] = 0
        self._state.numpecatom = max(0, self._state.numpecatom - 1)
        self._state.pecabund[z - 1, :] = 0.0
        # Keep current page in bounds
        n_pages = max(1, math.ceil(len(self._pec_atoms) / _MAX_SLOTS))
        self._current_page = min(self._current_page, n_pages - 1)
        self._refresh_slots()
        if not self._needs_synthesis:
            self._needs_synthesis = True
            self._btn_synth.ax.set_facecolor(_STALE_COLOR)
            self._btn_synth.label.set_color('white')
        self._fig.canvas.draw_idle()

    def _on_page(self, direction: int):
        n_pages = max(1, math.ceil(len(self._pec_atoms) / _MAX_SLOTS))
        self._current_page = (self._current_page + direction) % n_pages
        self._refresh_slots()

    # ------------------------------------------------------------------ #
    # batch.par serialization                                              #
    # ------------------------------------------------------------------ #

    def _save_batch_par(self):
        with open(self._fparam) as fh:
            orig = fh.readlines()

        out = []
        i   = 0
        wrote_abundances = False

        while i < len(orig):
            line  = orig[i]
            parts = line.strip().split()
            kw    = parts[0] if parts else ''

            if kw == 'abundances':
                orig_n = int(parts[1]) if len(parts) > 1 else 0
                out.append(f"abundances {len(self._pec_atoms)} {self._npass}\n")
                for z in self._pec_atoms:
                    vals = '  '.join(
                        f'{self._delta_abund[z][k]:.3f}'
                        for k in range(self._npass)
                    )
                    out.append(f"  {z}  {vals}\n")
                i += orig_n + 1
                wrote_abundances = True
            else:
                out.append(line)
                i += 1

        if not wrote_abundances and self._pec_atoms:
            out.append(f"abundances {len(self._pec_atoms)} {self._npass}\n")
            for z in self._pec_atoms:
                vals = '  '.join(
                    f'{self._delta_abund[z][k]:.3f}'
                    for k in range(self._npass)
                )
                out.append(f"  {z}  {vals}\n")

        with open(self._fparam, 'w') as fh:
            fh.writelines(out)


def synth_interactive(fparam: str = 'batch.par') -> SynthWidget:
    """
    Launch the interactive synth widget.

    Parameters
    ----------
    fparam : str
        Path to the MOOG batch parameter file (default: 'batch.par').

    Returns
    -------
    SynthWidget
        Keep a reference to prevent garbage collection of callbacks.

    Notes
    -----
    In Jupyter notebooks use ``%matplotlib widget`` before calling this
    function (requires the ipympl package: ``pip install ipympl``).
    """
    widget = SynthWidget(fparam)
    plt.show()
    return widget
