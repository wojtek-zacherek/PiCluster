# Parallel Computing Final Project
## Link to Project
https://github.com/wojtek-zacherek/PiCluster/tree/main/FractalProject
## Files
- FractalMain.c: Source code for my custom fractal generator. Takes command line input from 0-4, which determines which type of fractal will be generated.
- GnomeSort.c: Source code for a Gnome sort. The array size can be changed by modify the `N` value in the source code. Changes take effect after re-compiling.
- GnomeSortAnalysis.c: Source code for my Gnome Sort anlysis program. This various the number of processors assigned to my sorting algorithm, from 2 up to the number used to run this code (comm_sz).

## Compiling
I have included a `Makefile`. Running `make` should work. If that doesn't work, run the following lines:
```
mpicc -c FractalMain.c -o FractalMain.bin -lm
mpicc -c GnomeSort.c -o GnomeSort.bin -lm
mpicc -c GnomeSortAnalysis.c -o GnomeSortAnalysis.bin -lm
```

## Running
Running on the cluster should work like everyone else's runs, with `mpirun` or `mpiexec`. Which one works seems to depend on permission levels. For me, running `mpirun FractalMain.bin` works.

When using VSCode with WSL MPI, I needed to run the following command
```
mpirun -np 4 --mca btl_vader_single_copy_mechanism none FractalMain.bin
```
where the `np` flag sets the number of processors (up to the number of core my processor has?), and the ` --mca btl_vader_single_copy_mechanism none` phrase surpresses some useless and annoying errors.

## Issues
If there are any issues, please reach out to me at zacherw@rose-hulman.edu or 708-369-4270 (leaving a voicemessage - I don't call back to unknown numbers).