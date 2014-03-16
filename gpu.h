#ifndef GPU_H
#define GPU_H

#include "defines.h"
#include <cuda.h>
#include <stdio.h>

#define CUPRINTF(fmt, ...) printf("[%d, %d , %d]:\t" fmt, \
                                  blockIdx.x*gridDim.x+threadIdx.x,\
                                  blockIdx.y*gridDim.y+threadIdx.y,\
                                  blockIdx.z*gridDim.z+threadIdx.z,\
                                  __VA_ARGS__)

// Point is active
#define GPU_ACTIVE_POINT ((i < (gpu_def->locNx)) && (j < (gpu_def->locNy)) && (k < (gpu_def->locNz)) \
		&& (!( ((((gpu_def->rankx) != 0 && i == 0) || ((gpu_def->rankx) != (gpu_def->sizex) - 1 && i == (gpu_def->locNx) - 1)) && (gpu_def->Nx) >= 2) \
			|| ((((gpu_def->ranky) != 0 && j == 0) || ((gpu_def->ranky) != (gpu_def->sizey) - 1 && j == (gpu_def->locNy) - 1)) && (gpu_def->Ny) >= 2) \
			|| ((((gpu_def->rankz) != 0 && k == 0) || ((gpu_def->rankz) != (gpu_def->sizez) - 1 && k == (gpu_def->locNz) - 1)) && (gpu_def->Nz) >= 2))))

// Point is internal condition
#define GPU_INTERNAL_POINT (((((i != 0) && (i != (gpu_def->locNx) - 1)) || ((gpu_def->locNx) < 2)) && (j != 0) && (j != (gpu_def->locNy) - 1) \
		&& (((k != 0) && (k != (gpu_def->locNz) - 1)) || ((gpu_def->locNz) < 2))) && GPU_ACTIVE_POINT)

// Point is on boundary
#define GPU_BOUNDARY_POINT (((((i == 0) || (i == (gpu_def->locNx) - 1)) && ((gpu_def->locNx) >= 2)) || (j == 0) || (j == (gpu_def->locNy) - 1) \
		|| (((k == 0) || (k == (gpu_def->locNz) - 1)) && ((gpu_def->locNz) >= 2))) && GPU_ACTIVE_POINT)

__device__ int device_local_to_global(int local_index, char axis);
__device__ double device_ro_eff_gdy(int local);
__device__ int device_set_boundary_basic_coordinate(int i, int j, int k);
__device__ int device_is_injection_well(int i, int j, int k);
__device__ int device_is_output_well(int i, int j, int k);
__device__ void device_assing_k(double* k_w, double* k_n, double S_w);
__device__ void device_assign_S(int local);
__device__ void device_assign_ro(int local);
__device__ double device_left_difference (double* ptr, char axis);
__device__ double device_right_difference (double* ptr, char axis);
#ifdef ENERGY
__device__ void device_assign_H(int local);
__device__ void device_assign_E_current(int local);
#endif

#endif
