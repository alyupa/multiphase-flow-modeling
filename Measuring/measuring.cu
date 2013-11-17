#include "measuring.h"
#include <cuda.h>

#define CUPRINTF(fmt, ...) printf("[%d, %d , %d]:\t" fmt, \
	blockIdx.x*gridDim.x+threadIdx.x,\
	blockIdx.y*gridDim.y+threadIdx.y,\
	blockIdx.z*gridDim.z+threadIdx.z,\
	__VA_ARGS__)

#include "../../shared_test.cu"


void device_memory_allocation(double** DevBuffer, int buffer_size)
{
	cudaMalloc((void**) DevBuffer,  buffer_size * sizeof(double));
}

void device_memory_free(double* DevBuffer)
{
	cudaFree(DevBuffer);
}

// Инициализация ускорителя
// Расчет происходит на ускорителе, номер которого равен
// номеру запускающего процессора
void device_initialization(int rank)
{
	// Считаем, что ядер на узле не меньше, чем ускорителей
	int device = rank % GPU_PER_NODE;
	cudaSetDevice(device);

	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, device);

	if (devProp.major < 2)
	{
		printf("\nError! Compute capability < 2, rank=%d\n", rank);
	}
	/*
	if (!rank)
	{
		printf("Name : %s\n", devProp.name);
		printf("Compute capability : %d.%d\n", devProp.major, devProp.minor);
		printf("Total Global Memory : %ld\n", devProp.totalGlobalMem);
		printf("Shared memory per block: %d\n", devProp.sharedMemPerBlock);
		printf("Registers per block : %d\n", devProp.regsPerBlock);
		printf("Warp size : %d\n", devProp.warpSize);
		printf("Max threads per block : %d\n", devProp.maxThreadsPerBlock);
		printf("Total constant memory : %d\n", devProp.totalConstMem);
		printf("Number of multiprocessors: %d\n",  devProp.multiProcessorCount);
		//printf("Kernel execution timeout: %s\n\n", (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
		for (int i = 0; i < 3; ++i)
		{
			printf("Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
		}
		for (int i = 0; i < 3; ++i)
		{
			printf("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
		}
	}*/

	fflush(stdout);

	// Инициализируем библиотеку cuPrintf для вывода текста на консоль
	// прямо из kernel
	cudaPrintfInit();
}

// Финализация ускорителя
void device_finalization(void)
{
	// Останавливаем библиотеку cuPrintf для вывода текста на консоль
	// прямо из kernel
	cudaPrintfEnd();
}

__global__ void load_exchange_data_part_kernel(double* DevBuffer, int buffer_size)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i < buffer_size)
	{
		DevBuffer[i] = DevBuffer[buffer_size - i - 1];
	}
}


void load_exchange_data_part(double* HostBuffer, double* DevBuffer, int buffer_size)
{
	load_exchange_data_part_kernel <<< dim3(buffer_size / BlockNX, 1), dim3(BlockNX, 1)>>>(DevBuffer, buffer_size);
	checkErrors("load_exchange_data_part", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy(HostBuffer, DevBuffer, buffer_size * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

__global__ void save_exchange_data_part_kernel(double* DevBuffer, int buffer_size)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i < buffer_size)
	{
		DevBuffer[i] = DevBuffer[buffer_size - i - 1];
	}
}

void save_exchange_data_part(double* HostBuffer, double* DevBuffer, int buffer_size)
{
	cudaMemcpy(DevBuffer, HostBuffer, buffer_size * sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_part_kernel <<< dim3(buffer_size / BlockNX, 1), dim3(BlockNX, 1)>>>(DevBuffer, buffer_size);
	checkErrors("save_exchange_data_part", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

__global__ void set_devbuffer_values_kernel(double *DevBuffer, int buffer_size)
{
	for (int i = 0; i < buffer_size; i++)
	{
		DevBuffer[i] = buffer_size / i;
	}
}

void set_devbuffer_values(double *DevBuffer, int buffer_size)
{
	set_devbuffer_values_kernel <<< dim3(buffer_size / BlockNX, 1), dim3(BlockNX, 1)>>>(DevBuffer, buffer_size);
	checkErrors("set_devbuffer_values", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}
