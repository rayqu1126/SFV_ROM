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
- To measure runtime, uncomment and run the following lines:
  ```matlab
  % f1 = @() ode45(@(t,U) rhs_3D_burgers_state(t,U,params), tspan, ics, options);
  % timeit(f1)
  % f2 = @() ode45(@(t,U) rhs_3D_burgers_flux(t,U,params), tspan, ics, options);
  % timeit(f2)
  
### Figure 1

- Open **`SFVM3D_burgers.m`**, where the initial condition is defined by default.  
- To generate the **left plot**, vary `Nmode` from 10 to 50 (in increments of 5) and set `NHR = 4 * Ny^2` (no hyper-reduction).  
  Run the code, and the relative ROM error will be displayed automatically by the code.  
- To generate the **right plot**, fix `Nmode = 50` and vary `NHR` as  
  `NHR = [5, 64, 128, 256, 512, 1024, 2048, 4096]`.  
  Run the code, and the corresponding relative ROM error will also be displayed automatically.

  
### Figures 2 & 3

- Open **`SFVM3D_burgers.m`**, where the default initial condition is defined.  
- Set `Nmode` to either `10` or `20` and run the code.  
  The code will automatically generate and display the solution plots and statistical comparisons for both the full SFV method (using WENO with reconstructed fluxes) and the ROM.

### Table 2
- Open **`SFVM2D_euler.m`**, where the initial condition is defined by default.
- Vary `Ny` from 4 to 64 and run the code **before** the ROM stage.  
  The relative difference between the two reconstruction methods will be displayed automatically.
- To measure runtime, uncomment and run the following lines:
  ```matlab
  % f1 = @() ode45(@(t,U) rhs_2D_euler_state(t,U,params), tspan, ics, options);
  % timeit(f1)
  % f2 = @() ode45(@(t,U) rhs_2D_euler_flux(t,U,params), tspan, ics, options);
  % timeit(f2)
  ```

  
### Figure 4

- Open **`SFVM2D_euler.m`**, where the initial condition is defined by default.  
- Vary `Nmode` from 15 to 55 (in increments of 5).  
  Run the code, and the relative ROM error will be displayed automatically by the code.


### Figures 5 & 6

- Open **`SFVM2D_euler.m`**, where the default initial condition is defined.  
- Set `Nmode` to either `15` or `30` and run the code.  
  The code will automatically generate and display the solution plots and statistical comparisons for both the full SFV method (using WENO with reconstructed fluxes) and the ROM.

### Figures 7

- Open **`SFVM2D_euler.m`**, where the default initial condition is defined.  
- Set `Nmode =20` to and `NHR = Nmode`. Run the code.  
  The code will automatically generate and display the solution plots and statistical comparisons for both the full SFV method (using WENO with reconstructed fluxes) and the ROM.

### Figures 8
- Open **`SFVM2D_euler.m`**.
- Change initial conditions to
 ```matlab
rho(X < 0.3 + 0.3*Y) = 1; 
p(X < 0.3 + 0.3*Y) = 1;
```
- 

