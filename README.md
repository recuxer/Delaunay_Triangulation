# Delaunay_Triangulation
DT code, algorithms, research, etc. for KB

This folder contains all the files related to the C code implementation of the Delaunay Triangulation.

kristi.vtk: Do not modify this file. This is the output file that is generated. This is a file that you can open up with visit. It just shows the current triangulation and is a way to check your work.

visit_writer.c: Do not modify this file. This is necessary for the function which outputs the file above.

DTcode.C: This is the version of the DT code file that should be modified. It can be compiled with g++, but make sure the above file is present in the same folder or you will get a compilation error.

To compile:
g++ DTcode.C

To run:
./a.out      (or whatever executable name you give it)

To commit/push/pull:
Use git commands. If you'd like you can create a new branch to work from to limit any conflicts. I have a copy of this repo that is separate, so don't worry about overwritting code that shouldn't be - there's a backup.
