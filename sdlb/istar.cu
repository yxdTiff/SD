#include <cuda.h>
#include <cuda_runtime.h>				// Stops underlining of __global__
#include <device_launch_parameters.h>	// Stops underlining of threadIdx etc.

#include <stdio.h>

__global__ void compute_istar_cuda(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}