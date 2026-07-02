"""
He⁻ free-free and bound-free continuous opacity.
Translated from OpacHelium.f.

Note: the original Fortran uses `c1` (an uninitialized implicit real*8 variable)
where it should use `cc` (the correctly computed third coefficient). That is a
bug in the source code. Here we use `cc`, which is the physically correct value
from the Doughty & Fraser (1966) formula for He⁻ opacity.
"""
import numpy as np


def opac_heminus(state) -> None:
    """
    He⁻ free-free opacity into state.aHeminus.

    Frequency-dependent coefficients from Doughty & Fraser (1966), as
    implemented in ATLAS (Kurucz 1970).
    """
    ntau = state.ntau
    freq = state.freq

    a1 =  3.397e-46 + (-5.216e-31 + 7.039e-15 / freq) / freq
    b1 = -4.116e-42 + ( 1.067e-26 + 8.135e-11 / freq) / freq
    cc =  5.081e-37 + (-8.724e-23 - 5.659e-08 / freq) / freq

    t   = state.t[:ntau]
    ne  = state.ne[:ntau]
    nhe = state.numdens[1, 0, :ntau]
    uhe = state.u[1, 0, :ntau]

    state.aHeminus[:ntau] = (a1 * t + b1 + cc / t) * ne * nhe / uhe
