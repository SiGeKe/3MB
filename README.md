<p align="center">
  <img src="logo.png" width="200px">
</p>

# TRP - Meta-Basin-Analysis for Polydisperse 2D Systems

This project involves performing MD simulations of 2D systems and subsequent analyses to determine the lifetime of meta-basins. The simulation is based largely on the work of Burkard Doliwa and Andreas Heuer. Example output files can be found in the ``00_Code/01_Python`` folder (see also below).

> The scripts need the following directory structure to work (otherwise you can always change all direct paths):
> ``00_Configs``, ``01_Logs``, ``02_Restarts``, ``03_PEL``, ``04_IBM`` and ``05_Results``

> ``${index}`` denotes one system's number and will be referenced in the following text. 

## Overview

* ``00_Code/01_Python``: files for analysis.
* ``00_Code/02_LAMMPS``: LAMMPS input files for the MD simulations.
* ``00_Code/03_Bash``: Slurm scripts containing the procedures.

## Molecular Dynamics Simulations

The work begins with an initial LAMMPS run. This is outlined in ``run.sh`` and ``run.lmp`` and requires an initial structure in ``./00_Config/conf_${index}.dat``. During this run, restart files are written to ``02_Restarts/${index}_restart.${nr_relax}`` every ``${nr_relax}`` steps.

To subsequently find the minimal energy of every structure a conjugate gradient minimization is run with LAMMPS. Here, ``min.sh`` and ``min.lmp`` are used. The energies are written to ``./03_PEL/${index}_IS.dat``. 

## Interval Bisection Method

The MD simulations yield a time series of minimized energies, with a resolution of ``${nr_relax}``. However, it is not possible to determine which minima are visited between the ones that are found. This means that the exact timeline of each minimum is not known.

To achieve a resolution of one MD step, we use the IBM. This is done with ``ibm.sh``, ``MB_Analysis.py`` and ``ibm.lmp``. The procedure works as follows: 

1. The MD simulations yield inherent structures as identified by their energy $\xi$. First, the time points $t_i$ where $\xi(t_i) \neq \xi(t_j = t_i + N_r)$ are identified. $N_r$ denotes ``${nr_relax}``.

2. Starting from these points, trajectories at $t_k=(t_i+t_j)/2$ are reconstructed.

3. If $\xi(t_k)~=\xi(t_i)$, then $t_i=t_k$ is set, else $t_j=t_k$.

4. Steps 2. and 3. are repeated until the condition $t_i-t_j=1$ is met. 

In practice this means that the exact jump is found by repeatedly halving the time between two structures with different energies, and finding the minimal energy at this time point.

As this procedure is performed with Python and LAMMPS scripts, the condition can not be determined automatically. However, it is possible to calculate the exact number of iterations required, $L$, via $N_r$.

$L = \log{N_r}/\log{2}$

> Note that this procedure creates a large number of restart files, that are only cleaned up during the subsequent analysis.

## Meta Basin Analysis

Finally, ``MB_Analysis.py`` and ``ana.sh`` run the analysis of the time series of minima $\xi$. The exact procedure is explained in detail in Burkhard Doliwa's PhD thesis. Here, only a brief summary is provided.

1. The intervals between the first occurance $t^\ast$ and the last occurance $t^\dagger$ of a specific $\xi$ is determined.

2. If two intervals overlap by less than 50% one of them is shortened at random.

3. Two intervals that overlap by more than 50% are summarized into a new interval.

4. Intervals within other intervals are deleted.

The final intervals with their $t^\ast$ and $t^\dagger$ and the lowest visited $\xi$ characterizes the MB and its lifetime. 