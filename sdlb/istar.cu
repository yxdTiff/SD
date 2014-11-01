// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// CUDA runtime
#include <cuda_runtime.h>

// helper functions and utilities to work with CUDA
//#include <helper_cuda.h>
//#include <helper_functions.h>
#define DBL_MAX         1.7976931348623158e+308 /* max value */
extern "C" {
#include "cuda_wrapper.h"
}



__global__ void compute_istar_cuda(sigma_type *sigma, delta_type *delta, i_type *istar, double *argmax, vector pi_Tbar_x, vector Xvect, int rv_cols, int ictr, int count)
{

	int sig_pi=0, del_pi=0, c=0;
	double largest_temp;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < count)
	{
		double largestSoFar = -DBL_MAX;
		for (sig_pi = 0; sig_pi < sigma->cnt; sig_pi++)
		{
			if (sigma->ck[sig_pi] <= ictr)
			{
				/* Find the row in delta corresponding to this row in sigma */
				del_pi = sigma->lamb[sig_pi];
				/* Start with (Pi x Rbar) + (Pi x Romega) + (Pi x Tbar) x X */
				largest_temp = sigma->val[sig_pi].R + delta->val[del_pi][idx].R - pi_Tbar_x[sig_pi];

				/* Subtract (Pi x Tomega) x X. Multiply only non-zero VxT values */
				for (c = 1; c <= rv_cols; c++)
					largest_temp -= delta->val[del_pi][idx].T[c] * Xvect[delta->col[c]];

				if (largest_temp > largestSoFar)
				{
					largestSoFar = largest_temp;
					istar[idx].sigma = sig_pi;
					istar[idx].delta = del_pi;
					argmax[idx] = largestSoFar;
				}

			}
		}
	}
}

extern "C"
void launch_kernel(sigma_type *sigma, delta_type *delta, i_type *istar, double *argmax, vector pi_Tbar_x, vector Xvect, int rv_cols, int ictr, int count)
{
	compute_istar_cuda<<<(count / 640) + 1, 640>>>(sigma, delta, istar, argmax, pi_Tbar_x, Xvect, rv_cols, ictr, count);
	cudaDeviceSynchronize();
}

__global__ void foo()
{
}

extern "C" 
void CudaMain(void)
{
	foo<<<1, 1>>>();
}
