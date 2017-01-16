The Stratified Ocean Model with Adaptive Refinement (SOMAR) v1.0 DOI:10.5281/zenodo.247741
=====

Welcome to the SOMAR repository!

SOMAR is [free software](https://www.gnu.org/licenses/lgpl-2.1.html 'The GNU Lesser General Public License, version 2.1 applies.') provided jointly by the [Marine Sciences](http://marine.unc.edu/ 'UNC Marine Sciences website') and [Physics](http://physics.unc.edu/ 'UNC Physics website') departments of the [University of North Carolina at Chapel Hill](http://unc.edu/ 'UNC at Chapel Hill website').


Features
-----
- **Nonhydrostatic -** Our model solves the Boussinesq Navier-Stokes equations *without* the hydrostatic approximation in order to properly model the internal waves and tides that are ubiquitous in the ocean.


- **Complex topography -** We have maintained general covariance so that the domain can be described in curvilinear coordinates. This allows us to model irregular boundaries with stretched grids and a logically rectangular coordinate system.


- **Separation of background density and its deviation -** By splitting the density field into a vertical background stratification and a deviation, we relieve the Poisson solver of computing the associated hydrostatic component of the pressure. This treatment, already implemented in some regional models including [MITgcm](http://mitgcm.org/ 'The MITgcm website'), also prevents diffusion of oceanic features that are maintained by unmodeled phenomena.


- **Stable integration of stiff forcing terms -** The forcing terms that lead to buoyancy oscillations often impose instabilities in the form of fast, high-frequency waves. To quell these unphysical modes, we developed an integration method that updates all of the state variables semi-implicitly without requiring additional costly Poisson solves, as is the case with implicit Runge-Kutta schemes.


- **Anisotropic grid refinement.** A coarse underlying grid along with dynamic local refinement over transient features eliminates unnecessary computation in large portions of the domain. Since the background stratification requires some level of vertical resolution, it is often the case that we only need further resolution in the horizontal. Our anisotropic refinement methods are capable of providing additional cells only in those directions that are under-resolved. Furthermore, our coarse grids operate on larger timesteps than the finer grids. This refinement in both time and space provides a drastic speedup of computation, minimizes the number of required Poisson solves, and ensures all levels are evolving at a Courant number close to one, reducing numerical dissipation.


- **Anisotropic Poisson solvers -** In the absence of the hydrostatic approximation, we are faced with an ill-conditioned and exceedingly expensive Poisson problem for the pressure. To that end, we provide a leptic Poisson solver as well as a sophisticated semicoarsening multigrid solver to efficiently enforce the incompressibility condition. The leptic iterative method is a well-behaved, perturbative solution to highly anisotropic elliptic problems that has been demonstrated to work over several different geometries and degrees of anisotropy (See: [Scotti & Mitran 2008](http://www.sciencedirect.com/science/article/pii/S1463500308001005 'An approximated method for the solution of elliptic problems in thin domains: Application to nonlinear internal waves'), [Santilli & Scotti, 2011](http://dx.doi.org/10.1016/j.jcp.2011.06.022 'An efficient method for solving highly anisotropic elliptic equations')).


Software prerequisites
-----
SOMAR is built upon Chombo 3.1, an AMR framework that has been developed and is being distributed by the [Applied Numerical Algorithms Group](http://crd.lbl.gov/groups-depts/ANAG/ 'ANAG website') of [Lawrence Berkeley National Lab](http://www.lbl.gov/ 'LBNL website'). Before you can compile SOMAR, you must first download the Chombo software from [https://commons.lbl.gov/display/chombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations) and build its libraries. We strongly advise that you look at the [Chombo Design Document](https://seesar.lbl.gov/anag/chombo/ChomboDesign-3.1.pdf) for compilation instructions and prerequisites.


Compilation
-----
If you have successfully created the Chombo libraries, you should be ready to compile SOMAR. Open `exec/GNUmakefile` and edit the `CHOMBO_HOME` variable to point to your Chombo libraries. Then, open a terminal and run `make all` from the `exec` folder. An executable with the suffix `.ex` should be created. You have successfully compiled SOMAR.

If you plan to do any real work with SOMAR, it is likely that you'll want to tweak some compilation parameters. For example, you may want to do one or all of the following:

- change the dimensionality of the solver
- toggle debugging mode
- change the optimization level
- toggle code profiling

All of these tweaks can be performed in the `Chombo/lib/mk/Make.defs.local` file. You will find yourself altering this file quite often.

To get ready for your first demo run, try compiling SOMAR in 2D mode with debugging off. This is done by setting the following lines in `Make.defs.local` before running make.
```Makefile
DIM   = 2
DEBUG = FALSE
OPT   = HIGH
```
With these settings in place, make should produce an executable with a `somar2d` prefix and an `OPTHIGH.MPI.ex` suffix.


Taking a test drive
-----
Before creating your own simulation, you should successfully test one of the packaged simulations that come with the solver. In the `exec` folder, you will see several input files. Each of these files are targeted to run a specific problem on a specific machine and are identified accordingly (`inputs.problem.machine`). As a rule, an input file for machine A should never be altered to work on machine B. So if we want to run, say, the lock exchange demo problem on a machine called HAL, we should first copy one of the lock exchange input files to `inputs.LockExchange.HAL`, and then edit it as we please.

Open this new file and look for the following lines (it's okay if some of them are missing).
```
# amr.restart_file = chkpt_000010.2d.hdf5
# plot.plot_prefix = plot_
# plot.checkpoint_prefix = chkpt_
# plot.plot_period = 0.1
plot.plot_interval = 1
plot.checkpoint_interval = 100
```
In this listing, several of the lines are commented out with a `#`. This usually means we either do not want a particular feature or we are happy with the default values. In this case, we do not want to restart a simulation from a saved state and we are happy with the default file name prefixes. As you can see, all plot files will be prefixed with `plot_` and all checkpoint files (files used to restart a simulation) will be prefixed with `chkpt_`. If we want the plot files to go into a different directory, we can uncomment the plot\_prefix line to read `plot.plot_prefix = /home/user/myPlotFolder/plot_`. Other input parameters can similarly be altered as you see fit.

With your code compiled and input file prepared, you are now ready to run the lock exchange demo. To run the demo in serial, use `./somar2d.[config].OPTHIGH.MPI.ex inputs.LockExchange.HAL`. To run this simulation in parallel over 8 processors, use `mpirun -np 8 ./somar2d.[config].OPTHIGH.MPI.ex inputs.LockExchange.HAL`. The result should be two sets of [HDF5](http://www.hdfgroup.org/HDF5/ 'The HDF group website') files -- one set of checkpoint files that are used to restart a simulation, and one set of plot files that can be viewed in [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit 'The VisIt webpage').
