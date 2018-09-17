# Protein analysis using pytraj
Some miscellaneous scripts I made to do protein analysis using [pytraj](http://amber-md.github.io/pytraj/latest/index.html)
 with Python 2.7.5 (since it's unavailable for Python 3).

### cppy.py
This class loads the MD production files ready for processing. Before calculating any quantity it centers the molecule at the origin and minimizes RMSD.

### getpdbs.py
Simple script to get PDB file from RCSB Protein Data Bank.
