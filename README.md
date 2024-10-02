# TCR sharing in hypersensitivity pneumonitis: code to reproduce analyses
For the repertoire metrics, the R file in the /immunarch folder will produce the charts in *Figure 1*. That folder also contains the individual repertoires (and some metadata to make sense of which file is which) to use for the individual repertoire analyses to compare neighbours.

The Python script to perform those individual analyses is in the /repertoires folder (individual_repertoire_analysis.py). That script should be run with a command line argument to specify which repertoire to analyse. Details in the comments in the file.

The other two Python scripts in the /repertoire folder take the concatenated HP repertoire file (all_hp_bal_tcrb.tsv) as input to produce a TCR publicity report and comparisons against the VDJdb dataset in the folder as suggested by their names. They need no command line arguments.