# SFV_ROM
Reduced order models (ROMs) of stochastic finite volume (SFV) method for 1D Burgers' equation and 1D Compressible Euler equations with 1D/2D stochastic variables of uncertainty.

## Main files

- **`SFVM2D_burgers.m`**: FOM and ROM of 1D Burgers' equation with one stochastic variable
- **`SFVM2D_euler.m`**: FOM and ROM of 1D compressible Euler equation with one stochastic variable
- **`SFVM3D_burgers.m`**: FOM and ROM of 1D Burgers' equation with two stochastic variables


## Paper tables and figures
This sections gives details for creating the tables and figures associaite with our paper titled ''Model order reduction techniques for the stochastic finite volume method''. The MATLAB version used is R2024a.

### Table 1
- Open **`SFVM3D_burgers.m`**, where the initial condition is defined by default.  
- Vary `Ny` from 4 to 64 and run the code **before** the ROM stage.  
  The relative difference between the two reconstruction methods will be displayed automatically.  
- To measure runtime, uncomment the following lines:
  ```matlab
  % f1 = @() ode45(@(t,U) rhs_3D_burgers_state(t,U,params), tspan, ics, options);
  % timeit(f1)
  % f2 = @() ode45(@(t,U) rhs_3D_burgers_flux(t,U,params), tspan, ics, options);
  % timeit(f2)
  
### Figure 1

- Open **`SFVM3D_burgers.m`**, where the initial condition is defined by default.  
- To generate the **left plot**, vary `Nmode` from 10 to 50 (in increments of 5) and set `NHR = 4 * Ny^2` (no hyper-reduction).  
  The relative ROM error will be displayed automatically by the code.  
- To generate the **right plot**, fix `Nmode = 50` and vary `NHR` as  
  `NHR = [5, 64, 128, 256, 512, 1024, 2048, 4096]`.  
  The corresponding relative ROM error will also be displayed automatically.

  
### Table 2

