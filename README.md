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
to build executables into that folder. No need to worry too much about
the directory structure, the run scripts will take care
of running the executable.

The build log is stored to build/build.log, and should
be checked in the case that the build fails.

NOTE: If building on the PMCI cluster, you must first start an interactive PBS job
with `qsub -I`. The build process won't work on bioinf-head.

## Run (local) ##

1.  Store parameter sets, one per line, in a config file. (For treatment simulation,
    only one parameter set can be stored in the conf file. See also
    [Config file formats](#config-formats).)
2.  Choose a descriptive name for your group of tests. (Non-treatment only. Treatment
    sim will instead create time-stamped folder to store all sim results.)
3.  Decide how many times to run the simulation for each param set. (If running treatment
    simulation, decide how many processors (treatment config files) are needed.)
4.  Finally, run:

        $ ./run_local_simulations.sh <test group name> <num runs> <config file>

    for a non-treatment simulation, or

        $ ./run_local_treatment_sim.sh <num runs> <num processors> <config file>

    for a simulation with treatment.

## Run (cluster) ##

Follow steps 1-3 as for a local simulation. Then

1.  Run:

        $ ./setup_batch_job.sh [options] <test group name> <runs per param set> <config file>

    for a non-treatment simulation, or

        $ ./setup_batch_job_treatment.sh [options] <num runs> <num processors> <config file>

    for a treatment simulation.

    The `options` are

    *   `-m <email address>` to have results of cluster job emailed to you
    *   `-w <walltime>` to set the job walltime (in HH:MM:SS format)
    *   `-q <queue name>` to send job to a particular cluster queue (e.g. cloud, bacg).

    In either case, a new PBS script will have been created. Take note
    of the name of the script.

2.  Ensure you are on the bioinf-head node (otherwise qsub won't work) and run

        $ qsub <PBS script name>

    using the name of the script created in the previous step. Your job should
    now be running on the cluster. To check, run

        $ qstat -t

    to see all your active jobs. You should see one job for every param set
    in your conf file (non-treatment sim) or for every replicate run (treatment sim).

## Clean up ##

If you need to rebuild the simulation, run

    $ ./clean.sh

to remove the build folder before rebuilding.

## Analysis ##

To generate some simple summaries of a set of results, run

    $ analysis/summarise_results.py <results directory> [--compress]

where the directory points to results from a single group of sims.
If the --compress flag is supplied, the script will attempt to compress
the raw results in the sim directory after summarising them.

## <a name="config-formats"></a>Config file formats ##

General / initial parameters should be stored in a CSV-style config file,
with param names in the header row:

    param1,param2,param3, ...
    value11,value21,value31, ...
    value12,value22,value32, ...

The file should have a '.conf' extension. If running a treatment simulation,
only one row of values may be specified.

Drug parameters should be stored in an INI-like config file:

    param1 = value1
    param2 = value2
    param3 = value3
    ...

The file should have a '.drug' extension. (In fact, as currently implemented,
the file must be named 'CTX_Scheme_ID_<em>&lt;num&gt;</em>.drug' where <em>&lt;num&gt;</em> is the processor
number that will access the file.) The spaces around the '=' are optional.

## Contact ##

For questions, contact:

Luis Lara
llara@student.unimelb.edu.au

Yoshua Wakeham
ywakeham@student.unimelb.edu.au
