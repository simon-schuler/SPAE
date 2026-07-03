"""
State dataclass replacing all Fortran COMMON blocks.

Sources:
    Atmos.com  — atmosphere arrays, file paths, control flags
    Linex.com  — line data, line-profile work arrays
    Mol.com    — molecular equilibrium
    Quants.com — internally stored atomic data
    Factor.com — synthesis abundance factors, isotopes
    Kappa.com  — continuous opacity work arrays
    Pstuff.com — output / smoothing profile parameters (non-X11 subset)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import numpy as np

# Fixed array dimensions matching MOOG Fortran declarations
NTAU_MAX  = 100    # max atmosphere depth points
NLINES_MAX = 2500  # max spectral lines
NELEM     = 95     # number of elements
NMOL_MAX  = 110    # max molecular species
NEQ_MAX   = 30     # max equilibrium equations
NSYN_MAX  = 5      # max synthesis abundance variations
NISO_MAX  = 20     # max isotopes
NCOG_MAX  = 3000   # curve-of-growth table size
NDEPTH_MAX = 5000  # synthesis depth array


def _zeros(*shape):
    return field(default_factory=lambda: np.zeros(shape))

def _izeros(*shape):
    return field(default_factory=lambda: np.zeros(shape, dtype=int))

def _empty_str():
    return field(default_factory=str)

def _str(default=''):
    return field(default=default)


@dataclass
class State:
    """All shared MOOG state, replacing Fortran COMMON blocks."""

    # ------------------------------------------------------------------ #
    # ATMOSPHERE  (Atmos.com)                                              #
    # ------------------------------------------------------------------ #

    # Per-depth atmosphere arrays, shape (NTAU_MAX,)
    t:          np.ndarray = _zeros(NTAU_MAX)   # temperature [K]
    theta:      np.ndarray = _zeros(NTAU_MAX)   # 5040/T
    tkev:       np.ndarray = _zeros(NTAU_MAX)   # T [keV] = k*T in eV
    tlog:       np.ndarray = _zeros(NTAU_MAX)   # log10(T)
    pgas:       np.ndarray = _zeros(NTAU_MAX)   # gas pressure [dyn/cm²]
    ne:         np.ndarray = _zeros(NTAU_MAX)   # electron number density
    nhtot:      np.ndarray = _zeros(NTAU_MAX)   # total H number density
    numdens:    np.ndarray = _zeros(8, 2, NTAU_MAX)  # number densities
    molweight:  np.ndarray = _zeros(NTAU_MAX)   # mean molecular weight
    vturb:      np.ndarray = _zeros(NTAU_MAX)   # microturbulence [cm/s]
    scont:      np.ndarray = _zeros(NTAU_MAX)   # continuum source function
    kapref:     np.ndarray = _zeros(NTAU_MAX)   # reference opacity
    kaplam:     np.ndarray = _zeros(NTAU_MAX)   # continuum opacity at λ
    tauref:     np.ndarray = _zeros(NTAU_MAX)   # reference optical depth
    taulam:     np.ndarray = _zeros(NTAU_MAX)   # optical depth at λ
    kaplamabs:  np.ndarray = _zeros(NTAU_MAX)   # absorptive continuum opacity
    kaplamsca:  np.ndarray = _zeros(NTAU_MAX)   # scattering continuum opacity
    rho:        np.ndarray = _zeros(NTAU_MAX)   # mass density
    rhox:       np.ndarray = _zeros(NTAU_MAX)   # column mass density
    xref:       np.ndarray = _zeros(NTAU_MAX)   # integration variable (rhox or tauref)
    xdepth:     np.ndarray = _zeros(NTAU_MAX)   # geometric depth

    # Per-element arrays, shape (NELEM,) or (NELEM, 4, NTAU_MAX)
    elem:       np.ndarray = _zeros(NELEM)           # element atomic numbers
    xabund:     np.ndarray = _zeros(NELEM)           # number fractions (N_el/N_H)
    xabu:       np.ndarray = _zeros(NELEM)           # saved copy of xabund
    u:          np.ndarray = _zeros(NELEM, 4, NTAU_MAX)  # partition functions

    # Atmosphere scalars
    ntau:       int   = 0       # actual number of depth points
    jtau5:      int   = 0       # depth index where taulam ≈ 0.5
    flux:       float = 0.0     # continuum flux
    fudge:      float = 0.0     # extra continuous opacity fudge factor
    wavref:     float = 0.0     # reference wavelength for kapref
    abscale:    float = 0.0     # abundance scale factor
    deltaabund: float = 0.0     # abundance offset (abandy mode)
    iunits:     int   = 0       # wavelength units (0=Å, 1=μm)
    itru:       int   = 0       # model atmosphere format flag
    iraf:       int   = 0       # IRAF output flag
    modelnum:   int   = 0       # model number counter

    # Print / option flags
    modprintopt:  int = 0   # atmosphere print verbosity
    linprintopt:  int = 0   # line print verbosity
    linprintalt:  int = 0   # alternate line print flag
    fluxintopt:   int = 0   # flux (0) vs intensity (1) integration
    plotopt:      int = 0   # plot option
    dampingopt:   int = 0   # damping option (0=Unsold, 1=Barklem)
    specfileopt:  int = 0   # observed spectrum file format
    linfileopt:   int = 0   # line list file format
    printstrong:  int = 0   # print strong lines flag
    linecount:    int = 0   # line counter
    oldcount:     int = 0   # previous line counter
    scatopt:      int = 0   # scattering option

    # File paths (replacing Fortran unit numbers; physics layer never opens these)
    fmodel:    str = 'star.mod'
    flines:    str = 'no_filename_given'
    fslines:   str = 'no_filename_given'
    fobs:      str = 'no_filename_given'
    ftable:    str = 'no_filename_given'
    fparam:    str = 'batch.par'
    fbarklem:  str = ''
    fbarklemUV: str = ''
    f1out:     str = 'moog_out.1'
    f2out:     str = 'moog_out.2'
    f3out:     str = 'no_filename_given'
    f4out:     str = 'no_filename_given'
    f5out:     str = 'no_filename_given'

    # String control variables
    modtype:   str = ''        # model atmosphere type (e.g. 'KURUCZ')
    control:   str = ''        # driver mode ('abfind ', 'synth  ', etc.)
    moditle:   str = ''        # model title line
    linitle:   str = ''        # linelist title line

    # Element name symbols, length-2 strings
    names: list = field(default_factory=lambda: [''] * NELEM)

    # ------------------------------------------------------------------ #
    # LINES  (Linex.com)                                                   #
    # ------------------------------------------------------------------ #

    # Per-line × per-depth arrays, shape (NLINES_MAX, NTAU_MAX)
    a:      np.ndarray = _zeros(NLINES_MAX, NTAU_MAX)    # damping parameter
    dopp:   np.ndarray = _zeros(NLINES_MAX, NTAU_MAX)    # Doppler width [cm/s]
    kapnu0: np.ndarray = _zeros(NLINES_MAX, NTAU_MAX)    # line-center opacity

    # Per-line arrays, shape (NLINES_MAX,)
    gf:       np.ndarray = _zeros(NLINES_MAX)  # gf value (input)
    wave1:    np.ndarray = _zeros(NLINES_MAX)  # wavelength [Å]
    atom1:    np.ndarray = _zeros(NLINES_MAX)  # species code (e.g. 26.0=Fe I)
    e:        np.ndarray = _zeros(NLINES_MAX, 2)   # lower/upper excitation pot [eV]
    chi:      np.ndarray = _zeros(NLINES_MAX, 3)   # ionization potentials [eV]
    amass:    np.ndarray = _zeros(NLINES_MAX)  # atomic mass [amu]
    charge:   np.ndarray = _zeros(NLINES_MAX)  # ionization state (1=neutral)
    d0:       np.ndarray = _zeros(NLINES_MAX)  # dissociation energy (molecular)
    dampnum:  np.ndarray = _zeros(NLINES_MAX)  # input damping constant
    gf1:      np.ndarray = _zeros(NLINES_MAX)  # working gf (adjusted during fit)
    width:    np.ndarray = _zeros(NLINES_MAX)  # observed equivalent width [Å]
    abundout: np.ndarray = _zeros(NLINES_MAX)  # derived abundance per line
    widout:   np.ndarray = _zeros(NLINES_MAX)  # computed equivalent width
    strength: np.ndarray = _zeros(NLINES_MAX)  # line strength at jtau5
    rdmass:   np.ndarray = _zeros(NLINES_MAX)  # reduced mass (molecular)
    gambark:  np.ndarray = _zeros(NLINES_MAX)  # Barklem γ_6
    alpbark:  np.ndarray = _zeros(NLINES_MAX)  # Barklem α
    gamrad:   np.ndarray = _zeros(NLINES_MAX)  # radiative damping
    wid1comp: np.ndarray = _zeros(NLINES_MAX)  # EW from first iteration
    group:    np.ndarray = _izeros(NLINES_MAX) # line group index
    damptype: list = field(default_factory=lambda: ['        '] * NLINES_MAX)

    # Per-depth work arrays, shape (NTAU_MAX,)
    kapnu:  np.ndarray = _zeros(NTAU_MAX)  # total line opacity at current λ
    taunu:  np.ndarray = _zeros(NTAU_MAX)  # line optical depth
    cd:     np.ndarray = _zeros(NTAU_MAX)  # contribution function
    sline:  np.ndarray = _zeros(NTAU_MAX)  # line source function

    # Synthesis / profile arrays
    d:      np.ndarray = _zeros(NDEPTH_MAX)  # spectrum depths
    dellam: np.ndarray = _zeros(400)         # wavelength offsets from line center
    w:      np.ndarray = _zeros(100)         # computed EW per iteration

    # Curve-of-growth lookup table
    rwtab:  np.ndarray = _zeros(NCOG_MAX)   # log(RW) values
    gftab:  np.ndarray = _zeros(NCOG_MAX)   # log(gf) values
    gfhold: float = 0.0
    ntabtot: int  = 0   # number of entries in COG table

    # Synthesis wavelength range (set by params)
    start:    float = 0.0
    sstop:    float = 0.0
    step:     float = 0.0
    delta:    float = 0.0
    contnorm: float = 0.0
    oldstart: float = 0.0
    oldstop:  float = 0.0
    oldstep:  float = 0.0
    olddelta: float = 0.0

    # Wavelength stepping
    wave:     float = 0.0
    waveold:  float = 0.0
    wavestep: float = 0.0
    delwave:  float = 0.0
    st1:      float = 0.0

    # COG limits
    rwlow:  float = 0.0
    rwhigh: float = 0.0
    rwstep: float = 0.0
    cogatom: float = 0.0
    xratio:  float = 0.001   # weedout strength/continuum opacity ratio cutoff

    # Damping totals
    gammatot: float = 0.0
    gammav:   float = 0.0
    gammas:   float = 0.0
    gammar:   float = 0.0

    # Line index / mode integers
    dostrong:  int = 0   # number of strong lines
    gfstyle:   int = 0   # gf input style
    lineflag:  int = 0   # line present in interval flag
    molflag:   int = 0   # molecular equilibrium needed flag
    lim1:      int = 0   # current lower line index
    lim2:      int = 0   # current upper line index
    mode:      int = 0   # abfind mode (2=abfind, 3=synth)
    nlines:    int = 0   # number of regular lines
    nstrong:   int = 0   # number of strong lines
    ndepths:   int = 0   # number of depth points in line profile
    ncurve:    int = 0   # iteration counter in lineabund
    lim1line:  int = 0   # first line of current species
    lim2line:  int = 0   # last line of current species
    lim1obs:   int = 0   # first observed line index
    lim2obs:   int = 0   # last observed line index
    n1marker:  int = 0   # line list position marker
    iabatom:   int = 0   # element index of current species
    iaa:       int = 0   # first atom in molecule
    ibb:       int = 0   # second atom in molecule

    # ------------------------------------------------------------------ #
    # MOLECULAR EQUILIBRIUM  (Mol.com)                                     #
    # ------------------------------------------------------------------ #

    pmol:        np.ndarray = _zeros(NMOL_MAX)              # log partial pressures
    xmol:        np.ndarray = _zeros(NMOL_MAX, NTAU_MAX)   # molecule number densities
    xamol:       np.ndarray = _zeros(NEQ_MAX, NTAU_MAX)    # neutral atom number densities
    xatom:       np.ndarray = _zeros(NEQ_MAX)              # working atom densities
    patom:       np.ndarray = _zeros(NEQ_MAX)              # atom partial pressures
    amol:        np.ndarray = _zeros(NMOL_MAX)             # species codes
    smallmollist: np.ndarray = _zeros(NMOL_MAX)            # small default molecule list
    largemollist: np.ndarray = _zeros(NMOL_MAX)            # large default molecule list
    datmol:      np.ndarray = _zeros(7, NMOL_MAX)          # molecular constants
    const:       np.ndarray = _zeros(6, NMOL_MAX)          # eq. constants per layer
    h2ocoeff:    np.ndarray = _zeros(5)                    # H₂O partition coefficients
    co2coeff:    np.ndarray = _zeros(5)                    # CO₂ partition coefficients
    xnh2o:       np.ndarray = _zeros(NTAU_MAX)             # H₂O number density
    xnco2:       np.ndarray = _zeros(NTAU_MAX)             # CO₂ number density
    uh2o:        np.ndarray = _zeros(NTAU_MAX)             # H₂O partition function
    uco2:        np.ndarray = _zeros(NTAU_MAX)             # CO₂ partition function
    iorder:      np.ndarray = _izeros(NEQ_MAX)             # element order in eq. system
    neq:         int = 0   # number of equilibrium equations
    lev:         int = 0   # current depth level in eqlib loop
    nmol:        int = 0   # number of molecular species
    natoms:      int = 0   # number of atomic species
    molopt:      int = 0   # molecular equilibrium verbosity
    molset:      int = 0   # molecule set (0=none, 1=small, 2=large)

    # ------------------------------------------------------------------ #
    # ATOMIC DATA  (Quants.com)                                            #
    # ------------------------------------------------------------------ #

    xsolar:      np.ndarray = _zeros(NELEM)        # solar abundances log ε
    xam:         np.ndarray = _zeros(NELEM)         # atomic masses
    newpartdata: np.ndarray = _zeros(50, 6)         # partition function data
    xchi1:       np.ndarray = _zeros(NELEM)         # 1st ionization potential [eV]
    xchi2:       np.ndarray = _zeros(NELEM)         # 2nd ionization potential [eV]
    xchi3:       np.ndarray = _zeros(NELEM)         # 3rd ionization potential [eV]
    nudata:      np.ndarray = _izeros(6, 380)       # partition function integer data
    partflag:    np.ndarray = _izeros(NELEM, 4)     # which partition fn to use
    nu:          int = 0

    # ------------------------------------------------------------------ #
    # SYNTHESIS FACTORS  (Factor.com)                                      #
    # ------------------------------------------------------------------ #

    pecabund:     np.ndarray = _zeros(NELEM, NSYN_MAX)   # peculiar abundances
    newpecabund:  np.ndarray = _zeros(NELEM, NSYN_MAX)
    abfactor:     np.ndarray = _zeros(NSYN_MAX)           # abundance offsets per synth
    pec:          np.ndarray = _izeros(NELEM)             # peculiar element flags
    newpec:       np.ndarray = _izeros(NELEM)
    numpecatom:   int = 0
    newnumpecatom: int = 0
    numatomsyn:   int = 0
    newnumatomsyn: int = 0
    isynth:       int = 0    # current synthesis index
    ninetynineflag: int = 0  # 99-marker flag in linelist

    isotope:      np.ndarray = _zeros(NISO_MAX)
    newisotope:   np.ndarray = _zeros(NISO_MAX)
    isoabund:     np.ndarray = _zeros(NISO_MAX, NSYN_MAX)
    newisoabund:  np.ndarray = _zeros(NISO_MAX, NSYN_MAX)
    numiso:       int = 0
    newnumiso:    int = 0
    numisosyn:    int = 0
    newnumisosyn: int = 0
    isorun:       int = 0

    # ------------------------------------------------------------------ #
    # CONTINUOUS OPACITY WORK ARRAYS  (Kappa.com)                          #
    # ------------------------------------------------------------------ #

    aH1:     np.ndarray = _zeros(NTAU_MAX)
    aHminus: np.ndarray = _zeros(NTAU_MAX)
    aHeminus: np.ndarray = _zeros(NTAU_MAX)
    sigel:   np.ndarray = _zeros(NTAU_MAX)
    sigH:    np.ndarray = _zeros(NTAU_MAX)
    sigH2:   np.ndarray = _zeros(NTAU_MAX)
    sigHe:   np.ndarray = _zeros(NTAU_MAX)
    aC1:     np.ndarray = _zeros(NTAU_MAX)
    aMg1:    np.ndarray = _zeros(NTAU_MAX)
    aMg2:    np.ndarray = _zeros(NTAU_MAX)
    aAl1:    np.ndarray = _zeros(NTAU_MAX)
    aSi1:    np.ndarray = _zeros(NTAU_MAX)
    aSi2:    np.ndarray = _zeros(NTAU_MAX)
    aFe1:    np.ndarray = _zeros(NTAU_MAX)
    evhkt:   np.ndarray = _zeros(NTAU_MAX)
    freq:    float = 0.0
    freqlg:  float = 0.0

    # ------------------------------------------------------------------ #
    # OUTPUT / SMOOTHING PARAMETERS  (Pstuff.com, non-X11 subset)         #
    # ------------------------------------------------------------------ #

    # Smoothing profile parameters (used by smooth.py / synth driver)
    smtype:      str   = 'n'   # smoothing type: n/g/l/v/c/m/d/r/p
    addflux:     float = 0.0   # veiling (extra continuum fraction)
    vsini:       float = 0.0   # projected rotation velocity [km/s]
    limbdark:    float = 0.0   # limb darkening coefficient
    vmac:        float = 0.0   # macroturbulence [km/s]
    fwhmgauss:   float = 0.0   # FWHM of Gaussian instrumental profile [Å]
    fwhmloren:   float = 0.0   # FWHM of Lorentzian profile [Å]

    # Smoothed profile work arrays (length 1000)
    p:    np.ndarray = _zeros(1000)   # raw spectrum depths
    prot: np.ndarray = _zeros(1000)   # rotation-broadened
    pmac: np.ndarray = _zeros(1000)   # macroturbulence-broadened

    # abfind statistics output (Pstuff.com)
    average:  float = 0.0
    deviate:  float = 0.0
    xxm1: float = 0.0;  xxb1: float = 0.0;  xxr1: float = 0.0  # EP trend
    xxm2: float = 0.0;  xxb2: float = 0.0;  xxr2: float = 0.0  # RW trend
    xxm3: float = 0.0;  xxb3: float = 0.0;  xxr3: float = 0.0  # wavelength trend
    deltaep:  float = 0.0
    deltarw:  float = 0.0
    deltawv:  float = 0.0

    # Observed spectrum (for synth overlay)
    xobs: np.ndarray = field(default_factory=lambda: np.zeros(500000, dtype=np.float32))
    yobs: np.ndarray = field(default_factory=lambda: np.zeros(500000, dtype=np.float32))
    obsitle: str = ''

    # Line counter for output
    lount: int = 0
    kount: int = 0
    nsyn:  int = 0
