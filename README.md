# Protein analysis using pytraj
Some miscellaneous scripts I made to do protein analysis using [pytraj](http://amber-md.github.io/pytraj/latest/index.html)

_work in progress..._

#### cppy.py
Loads the MD production files ready for processing. Before calculating any quantity it centers the molecule at the origin and minimizes RMSD.

#### getpdbs.py
Simple script to get PDB file from RCSB Protein Data Bank.

#### dataPrep.sh
Runs afterProd.in using cpptraj, organizes data in folders as it's convenient for cppy.py

#### afterProd.in
Removes water molecules for later analysis. The most important stuff is done with pytraj when running cppy.py.
