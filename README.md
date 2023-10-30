# Marching Squares

Marching Squares is a graphics algorithm introduced in the 1980s that can be used to delimit contours in an image. It is commonly used for drawing elevation lines on topographic maps, temperatures on thermal maps, pressure points on pressure field maps, and more.

This project is an optimization of the original algorithm, which leverages pthreads to parallelize the execution of the algorithm. The scalability of the algorithm is improved by allowing a dynamic number of threads via command line input. On the author's machine, the speedup for 2 threads is 1.89x, and for 4 threads, it's 3.71x.

## How to Run

To run the full test suite, use the following command in the `src` directory:

./local.sh checker


Additionally, there's a script that can be used to run a single test case:

./run_single.sh <test_number>


## Implementation Details

The parallelization in this project is achieved in several key areas:

1. **Rescaling of the Image**: The image rescaling is parallelized.
2. **Sample Grid (First Part of Marching Squares Algorithm)**: The first part of the Marching Squares algorithm is parallelized.
3. **March (Second Part of Marching Squares Algorithm)**: The second part of the Marching Squares algorithm is parallelized.

Input and output (I/O) operations are not parallelized to avoid extra overhead.

To synchronize the threads, two barriers are used:

- The first barrier synchronizes the threads after the rescaling of the image.
- The second barrier is used after the first part of the Marching Squares algorithm (sample grid).

For images larger than 2048x2048, they are downscaled to 2048x2048 before executing the algorithm. This is handled in the thread function by using a flag that is received as a parameter (`needs_rescaled`).

The formula used for parallelization is as follows:

```c
int start = ID * (double)N / P;
int end = min((ID + 1) * (double)N / P, N);
