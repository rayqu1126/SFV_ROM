# SFV_ROM
Reduced order models (ROMs) of stochastic finite volume (SFV) method for 1D Burgers' equation and 1D Compressible Euler equations with 1D/2D stochastic variables of uncertainty.

## Main files

- **`SFVM2D_burgers.m`**: FOM and ROM of 1D Burgers' equation with one stochastic variable
- **`SFVM2D_euler.m`**: FOM and ROM of 1D compressible Euler equation with one stochastic variable
- **`SFVM3D_burgers.m`**: FOM and ROM of 1D Burgers' equation with two stochastic variables


## Publication tables and figures
This sections gives details for creating the tables and figures associaite with our paper titled ''Model order reduction techniques for the stochastic finite volume method''. The MATLAB version used is R2024a.
### Table 1
- Go to **`SFVM3D_burgers.m`** where the initial condition is set as default
- Vary `Ny` from 4 to 64 and run the code before ROM. The relative difference computation is displayed from the code.
- To measure runtime, uncomment the following lines:
  ```matlab
  % f1 = @() ode45(@(t,U) rhs_3D_burgers_state(t,U,params), tspan, ics, options);
  % timeit(f1)


