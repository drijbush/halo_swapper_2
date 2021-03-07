module constants

  !! Module containing constants to be used by Ewald_minres

  Use, Intrinsic :: iso_fortran_env, Only:  real64, int64
  Implicit None

  Private

  Integer, Parameter, Public :: wp = real64 !! Real working precision
  Integer, Parameter, Public :: li = int64  !! Long integer

  Real( wp ), Parameter, Public :: zero = 0.0_wp, one = 1.0_wp !! Basic real constants
  Real( wp ), Parameter, Public :: pi = 3.141592653589793238462643383279502884197_wp !! Pi
  Real( wp ), Parameter, Public :: r4pie0 = 138935.48350000000_wp !! Electric prefactor



end module constants
