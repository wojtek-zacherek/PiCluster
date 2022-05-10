# Homework 6
## To run
You *should* be able to run `make` and it should compile the program.
If not ... do `mpicc GEPP.c -o GEPP.bin -lm`

## Fundemental Logic
For this one, I decided to extra the logic to avoid copying over memory lines durnig each pivoting step. This required a lot of thought (for me at least) and hours of debugging. In hindsight, I should've stuck with the normal memory-memory rotation, but alas ... only hindsight is 20/20.

## Warning
This code is messy (it was a lot messier before). I've reduced it down "safely" to the minimum it needs to run. It isn't necessarily the most efficient, but my limited trials seemed to work.

## Modify input
For now, there are 4 variables to modify at the top of the GEPP.c source code.
- response - set this to 'y' or 'Y' to use your user input
- N_user - the size of the square matrix. This is seperate from 'N' to simulate distributed loads
- matInputs - this is your row-major matrix written in an array/list format.
- bInputs - same as matInputs, but for your b vector.

## Other Stuff
Currently the program will print the matrix and b vector...
- at the beginning
- during each column step
- at the very end
The matrix is not actively reorganized into a descending traingular matrix, but the results are correct (I verified some results with online calculators).
