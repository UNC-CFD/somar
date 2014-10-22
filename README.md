The Stratified Ocean Model with Adaptive Refinement (SOMAR)
=====

Introduction
-----


Software prerequisites
-----
SOMAR is built upon Chombo 3.1, an AMR framework that has been developed and is being distributed by the Applied Numerical Algorithms Group of Lawrence Berkeley National Lab. Before you can compile SOMAR, you must first download the Chombo software from [https://commons.lbl.gov/display/chombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations) and build its libraries. We strongly advise that you look at the [Chombo Design Document](https://seesar.lbl.gov/anag/chombo/ChomboDesign-3.1.pdf) for compilation instructions and prerequisites.


Compilation
-----
If you have successfully created the Chombo libraries, you should be ready to compile SOMAR. Open `exec/GNUmakefile` and edit the `CHOMBO_HOME` variable to point to your chombo libraries. It is not likely, but you may also need to alter `cxxcppflags` to include system libraries. Once all variables are properly set, run `make all` from the `exec` folder. An executable with the prefix `somar` should be created.

If you edit the `Chombo/lib/mk/Make.defs.local` file, you will need to rerun make for the changes to take effect. You will be altering this file quite often for various reasons.

- To change the dimensionality of the solver.
- To toggle debugging mode.
- To change the optimization level.
- To toggle code profiling.

For example, to compile a 2D simulation with debugging off, we would set the following lines in `Make.defs.local` before running make.
```Makefile
DIM   = 2
DEBUG = FALSE
OPT   = HIGH
```
With these settings, make should produce an executable with an `OPTHIGH.MPI.ex` suffix. 


Taking a test drive
-----
Before creating your own simulation, you should successfully test one of the packaged simulations that come with the solver. In the `exec` folder, you will see several input files. Each of these files are targeted to run a specific problem on a specific machine and are identified accordingly (`inputs.problem.machine`). As a rule, an input file for machine A should never be altered to work on machine B. So if we want to run, say, the lock exchange demo problem on a machine called Topsail, we should first copy one of the lock exchange input files to `inputs.LockExchange.Topsail`, and then edit it as we please.

Open this new file and look for the following lines (it's okay if some of them are missing).
```
# amr.restart_file = chkpt_000010.2d.hdf5
# plot.plot_prefix = plot_
# plot.checkpoint_prefix = chkpt_
# plot.plot_period = 0.1
plot.plot_interval = 1
plot.checkpoint_interval = 100
```
In this listing, several of the lines are commented out with a `#`. This usually means we do not want a particular feature turned on or we are happy with the default values. In this case, we do not want to restart a simulation from a saved state and we are happy with the default file name prefixes. As you can see, all plot files will be prefixed with `plot_` and all checkpoint files (files used to restart a simulation) will be prefixed with `chkpt_`. If we want the plot files to go into a different directory, we can uncomment the plot\_prefix line to read `plot.plot_prefix = /home/user/myPlotFolder/plot_`. Other input parameters such as `amr.final` and `amr.maxsteps` can be altered as you see fit.

Now it is time to compile a 2D simulation with debugging off. Make sure that the Chombo makefile `Chombo/lib/mk/Make.defs.local` knows this before running make.
```Makefile
DIM   = 2
DEBUG = FALSE
OPT   = HIGH
```
Next, run `make all` from the `exec` directory. This will create an executable with an `OPTHIGH.MPI.ex` suffix. This executable must be given an input file as a parameter. For example, to run a simulation on 8 processors, use `mpirun -np 8 amrins2d.OPTHIGH.MPI.ex inputs.lockexchange.topsail`. The result should be two sets of HDF5 files -- one set of checkpoint files that are used to restart a simulation, and one set of plot files that can be viewed in [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit).






**TODO**
