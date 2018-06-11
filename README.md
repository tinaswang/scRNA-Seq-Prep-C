# scRNA-Seq-Prep-C: 10x scRNA-Seq v2 Preprocessing Tool


This is a C++ version of <https://github.com/pachterlab/scRNA-Seq-TCC-prep/tree/10xTCCprep-SC3Pv2>.

Does everything except split the reads and create the kallisto input.

HOW TO COMPILE AND RUN:

Run the Makefile included

Necessary dependencies:

boost, zlib

Command line example:

```
./preprocess --dmin (dmin value) --expected (expected number of cells) --filepath 
(path to directory with the read files and cell barcodes)
```




Allowed options:
  -  [ --help ]         Load in the samples required
  -  [ --dmin ] arg     set dmin for each sample like --dmin 4 (4 is the default)
  -  [ --expected ] arg set expected number of cells, ex. --expected 5000
  -  [ --filepath ] arg set input directories
