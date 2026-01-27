# PLIP-HPC
Code and documentation for running the Protein-Ligand Interaction Profiler (PLIP) on the MSU HPCC.

PLIP can be run via the command line using the script `plip_singluarlity.py` for a sinlge sample. Run `python3 plip_singluarity.py --help` for usage.

Example implementations for my personal workflow are provided in `batch_plip_runs` (for running multiple samples sequentially in a single SLURM job) and `slurm_array_plip_runs.sb` (for running multiple samples in parallel via a SLURM job array). These implementations are specific to my directory structure and use case and are intended to provide a loose guide for batch implementation on the MSU HPCC.
