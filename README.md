# README #

## Dependencies ##

Boost, GSL, OpenMPI

## Build ##

Build without treatment:

    $ cd <source directory>
    $ ./build.sh

Build *with* treatment:

    $ cd <source directory>
    $ ./build.sh --treatment

These scripts create a build/ directory and use CMake
to build executables into that folder.

The build log is stored to build/build.log, and should
be checked in the case that the build fails.

NOTE: If building on the PMCI cluster, you must start an interactive PBS job
with `qsub -I`. The build process won't work on bioinf-head.

## Run (local) ##

1.  Store parameter sets, one per line, in a config file. [Non-treatment only.]
2.  Choose a descriptive name for your group of tests. [Non-treatment only.]
3.  Decide how many times to run the simulation for each param set.
4.  Finally, run:

        $ ./run_local_simulations.sh <test group name> <runs per param set> <config file>

    for a non-treatment simulation, or

        $ ./run_local_treatment_sim.sh <num runs>

    for a simulation with treatment.

## Run (cluster) ##

Follow steps 1-3 as for a local simulation. Then

1.  Run:

        $ ./setup_batch_job.sh [options] <test group name> <runs per param set> <config file>

    for a non-treatment simulation, or

        $ ./setup_batch_job_treatment.sh [options] <num runs> <num processors>

    for a treatment simulation.

    The `options` are

    *   `-m <email address>` to have results of cluster job emailed to you
    *   `-w <walltime>` to set the job walltime (in HH:MM:SS format)
    *   `-q <queue name>` to send job to a particular cluster queue (e.g. cloud, bacg).

    In either case, a new PBS script will have been created.

2.  Ensure you are on the bioinf-head node (otherwise qsub won't work) and run

        $ qsub <PBS script name>

    Your job should now be running on the cluster. To check, run

        $ qstat -t

    to see all your active jobs.

## Clean up ##

Just run

    $ ./clean.sh

to remove the build folder.

## Notes ##

The garbage collector is coded as unique pointers.
These version do not require yet the HPC layer.
The models are described in the pdf

## Contact ##

Do not hesitate to contact Luis Lara or Yoshua Wakeham with questions.
