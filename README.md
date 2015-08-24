# README #

## Dependencies ##

Boost, GSL, OpenMP

## Build ##

    $ cd <source directory>
    $ ./build.sh

This script creates a build/ directory and uses CMake
to build executables into that folder.

The build log is stored to build/build.log, and should
be checked in the case that the build fails.

## Run (local) ##

1.  Store parameter sets, one per line, in a config file.
2.  Choose a descriptive name for your group of tests.
3.  Decide how many times to run the simulation for each param set.
4.  Finally, run:

        $ ./run_local_simulations.sh <test group name> <runs per param set> <config file>

## Run (cluster) ##

Follow steps 1-3 as for a local simulation. Then

1.  Run:

        $ ./setup_batch_job.sh [options] <test group name> <runs per param set> <config file>

    where `options` are

    *   `-m <email address>` to have results of cluster job emailed to you
    *   `-w <walltime>` to set the job walltime (in HH:MM:SS format)
    *   `-q <queue name>` to send job to a particular cluster queue (e.g. cloud, bacg).

    This script creates a new PBS script named `PBS-<test group name>.sh`.

2.  Ensure you are on the head node (otherwise qsub won't work) and run

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
