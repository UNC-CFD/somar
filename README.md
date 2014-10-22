The Stratified Ocean Model with Adaptive Refinement (SOMAR)
=====

Introduction
-----
TODO


Software prerequisites
-----
SOMAR is built upon Chombo 3.1, an AMR framework that has been developed and is being distributed by the Applied Numerical Algorithms Group of Lawrence Berkeley National Lab. Before you can compile SOMAR, you must first download the Chombo software from [https://commons.lbl.gov/display/chombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations) and build its libraries. We strongly advise that you look at the [Chombo Design Document](https://seesar.lbl.gov/anag/chombo/ChomboDesign-3.1.pdf) for compilation instructions and prerequisites.


Compilation
-----
If you have successfully created the Chombo libraries, you should be ready to compile SOMAR. Open `exec/GNUmakefile` and edit the `CHOMBO_HOME` variable to point to your Chombo libraries. Then, open a terminal and run `make all` from the `exec` folder. An executable with the suffix `.ex` should be created. You have successfully compiled SOMAR.

If you plan to do any real work with SOMAR, it is likely that you'll want to tweak some compilation parameters. For example, you may want to do one or all of the following:

- Change the dimensionality of the solver
- Toggle debugging mode
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
In this listing, several of the lines are commented out with a `#`. This usually means we either do not want a particular feature or we are happy with the default values. In this case, we do not want to restart a simulation from a saved state and we are happy with the default file name prefixes. As you can see, all plot files will be prefixed with `plot_` and all checkpoint files (files used to restart a simulation) will be prefixed with `chkpt_`. If we want the plot files to go into a different directory, we can uncomment the plot\_prefix line to read `plot.plot_prefix = /home/user/myPlotFolder/plot_`. Other input parameters such as `amr.final` and `amr.maxsteps` can be altered as you see fit.

With your code compiled and input file prepared, you are now ready to run the lock exchange demo. To run the demo in serial, use `./somar2d.[config].OPTHIGH.MPI.ex inputs.LockExchange.HAL`. To run this simulation in parallel over 8 processors, use `mpirun -np 8 ./somar2d.[config].OPTHIGH.MPI.ex inputs.LockExchange.HAL`. The result should be two sets of HDF5 files -- one set of checkpoint files that are used to restart a simulation, and one set of plot files that can be viewed in [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit).

