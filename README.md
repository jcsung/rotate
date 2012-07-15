rotate
======

Quick and dirty program for rotating atoms around origin in PDB format

Compile using g++ or compiler of choice:

  g++ rotate.cpp -o rotate.exe

Modify file 'input' with indexes (starting from 1) to rotate
and coordinates to rotate by.

Then run:
  ./rotate.exe INITIAL_PDB FINAL_PDB

Only works with PDB files.
Only rotates ATOM and HETATM type elements.
Not tested for extreme cases.
No periodic boundary stuff.

