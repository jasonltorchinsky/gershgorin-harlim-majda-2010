! Defines the standard Fortran types, for compatibility across platforms

! Integer precision types: single byte, double bytes, quad bytes, oct bytes.
! E.g., sb can hold up to 10^2, db can hold up to 10^4.
INTEGER, PARAMETER :: sb = SELECTED_INT_KIND(2)
INTEGER, PARAMETER :: db = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: qb = SELECTED_INT_KIND(8)
INTEGER, PARAMETER :: ob = SELECTED_INT_KIND(16)

! Real precision types, a la Metcal et. al (2004, p. 71).
INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(33, 4931)
