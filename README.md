# PDB_stacker
takes a time directory of PDB files and shows stacking percentage as a function of time.  This is great for packages like AMBER.  This program was first used to study [RNA tetramers](http://pubs.acs.org/doi/abs/10.1021/ct501025q) and has been used also to study [compaction of duplex nucleic acids](http://pubs.acs.org/doi/full/10.1021/acscentsci.7b00084)

## examples
execution is as simple as it gets:
```shell
perl pdb_stacking.pl
```
and then the script scans all PDB files in your directory.  The PDB files should be named like 1701ns.pdb.  The Perl script will automatically create stacking plots for you, so you can see how all of the stacking variables, defined [here](http://pubs.acs.org/doi/abs/10.1021/ct501025q) vary over time.  Examples can be found in that link.
