# `time_stepper`

This is the directory containing the time-stepper module for the simulation.

The file `time_stepper.f90` contains the `TIME_STEPPER` module and all associated subroutines aside from the particular time-stepping schemes.

All other `.f90` files contain particular time-stepping schemes that may be used for the simulation. In particular, here are the `timeStepScheme` values that correspond to each time-stepping scheme:
  - `timeStepScheme` = 0: Forward Euler (issues with stability for large time-step sizes).
  - `timeStepScheme` = 1: Backward Euler (issues with stability for large time-step sizes).
  - `timeStepScheme` = 2: Forward Euler with Adaptive Time-Stepping (ensures absolute stability).
  - `timeStepScheme` = 3: Backward Euler with Adaptive Time-Stepping (ensures absolute stability).
  - `timeStepScheme` = 4: Analytic Statistics with Trapezoidal Numerical Integration (issues with stability for large time-step sizes).
