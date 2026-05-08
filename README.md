# SFV_ROM
Reduced order models (ROMs) of stochastic finite volume (SFV) method for 1D Burgers' equation and 1D Compressible Euler equations with 1D/2D stochastic variables of uncertainty.

## Main files
- **`SFVM2D_euler.m`**: Full SFV method and corresponding ROM of 1D compressible Euler equation with one stochastic variable
- **`SFVM3D_euler.m`**: Full SFV method and corresponding ROM of 1D compressible Euler equation with two stochastic variables


## Paper tables and figures
This sections gives details for creating the tables and figures associaite with our paper titled ''Model order reduction techniques for the stochastic finite volume method''. The MATLAB version used is R2024a. 

## Example 1: stochastic Sod shock tube with one stochastic variable
### Table 1
- Open **`SFVM2D_euler.m`**, where the initial condition and tspan are defined by default. (If not, uncomment the tspan and initial conditions to Sod shock tube.)
- Comment out the preset `Ny`
- Uncomment (or add) the following:
  ```matlab
  % Ny_list = [4;8;16;32;64]
  % state_time_list = [];
  % state_step_list = [];
  % flux_time_list = [];
  % flux_step_list = [];
  % difference_list = [];
  % 
  % for r = 1:length(Ny_list)
  % Ny = Ny_list(r);
  % ...
  % f1 = @() ode45(@(t,U) rhs_2D_euler_state(t,U,params), tspan, ics, options);
  % state_time_list = [state_time_list;timeit(f1)];
  % state_step_list = [state_step_list;size(U_state,1)];
  % ...
  % f2 = @() ode45(@(t,U) rhs_2D_euler_flux(t,U,params), tspan, ics, options);
  % flux_time_list = [flux_time_list;timeit(f2)];
  % flux_step_list = [flux_step_list;size(U_flux,1)];
  % ...
  % difference_list = [difference_list;er_flux];
  % ...
  % end
- Run the code **before** the ROM stage.
  The relative difference between the two reconstruction methods, total time-stepping runtime, and total number of time steps are stored in the corresponding vectors.

  
### Figure 2
- Open **`SFVM2D_euler.m`**, where the initial condition and tspan are defined by default. (If not, uncomment the tspan and initial conditions to Sod shock tube.)
- Run the entire code and the singular value plot should display. The singular values are stored in `s`.

### Figure 3
- Open **`SFVM2D_euler.m`**, where the initial condition and tspan are defined by default. (If not, uncomment the tspan and initial conditions to Sod shock tube.)
- Comment out the preset `Nmode`
- Uncomment the following:
  ```matlab
  % ROM_er_list = [];
  % ROM_time_list = [];
  % ROM_step_list = [];
  %
  % N_list = 6:11;
  % for m = 1:length(N_list)
  % Nmode = N_list(m);
  % ...
  % f3 = @() ode45(@(t,U) rhs_2D_euler_interp_HR(t,U,params), tspan, ics, options);
  % ROM_time_list = [ROM_time_list; timeit(f3)];
  % ROM_step_list = [ROM_step_list; size(U_ROM,1)];
  % ...
  % % ROM_er_list = [ROM_er_list; er_ROM];
  % end
- Run the code, and ROM error, total time-stepping runtime, total time steps are stored in the corresponding vectors.
- Plot them by uncommenting the provided figure command. 
  

### Figures 4 & 5

- Open **`SFVM2D_euler.m`**, where the initial condition and tspan are defined by default. (If not, uncomment the tspan and initial conditions to Sod shock tube.)  
- Set `Nmode` to either `5` or `13` and run the code.  
  The code will automatically generate and display the solution plots and statistical comparisons for both the full SFV method (using WENO with reconstructed fluxes) and the ROM.


