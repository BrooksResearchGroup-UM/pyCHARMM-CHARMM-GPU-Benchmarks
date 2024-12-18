This benchmark was adapted from the apoa1 benchmark from
namd. It was chaged to use the modern force fields for lipids and proteins
by CHarles L. Brooks III.

To run this benchmark simply cd to either the openmm or bladelib folders and
execute "python ../apoa1_bench.py > apoa1_bench.log". The python benchmark script runs pyCHARMM
recording the time for the MD steps in a loop for five times and reports the
average and sd of the runs.

Note that I am using pyCHARMM build following the perscription from the
pyCHARMM Worokshop.