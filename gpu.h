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

__device__ int device_is_active_point(int i, int j, int k);
__device__ int device_local_to_global(int local_index, char axis);
__device__ double device_ro_eff_gdy(ptr_Arrays DevArraysPtr, int local);
__device__ int device_set_boundary_basic_coordinate(int i, int j, int k);
__device__ int device_is_injection_well(int i, int j, int k);
__device__ int device_is_output_well(int i, int j, int k);
__device__ void device_assing_k(double* k_w, double* k_n, double S_w);

#endif
