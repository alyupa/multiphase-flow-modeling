#include "gpu.h"
#include "shared_test.cu"

__constant__ consts gpu_def [1];
__device__ ptr_Arrays DevArraysPtr[1];
__device__ double *DevBuffer;
ptr_Arrays DevArraysPtrLoc[1];
double *DevBufferLoc;
extern double *HostBuffer;

#include "three-phase.cu"

#ifdef ENERGY
#include "gauss.cu"
#include "energy.cu"
#endif

__device__ void device_assign_ro(int local)
{
#ifdef ENERGY
	// !!! Вынести коэффициенты теплового расширения в const consts &def и использовать T_0 оттуда же
	double alfa_w = 1.32E-7; // 1/K !!! E-4
	double alfa_n = 9.2E-7;
	double T_0 = 273;

	DevArraysPtr->ro_w[local] = gpu_def->ro0_w * (1. + (gpu_def->beta_w) * (DevArraysPtr->P_w[local] - gpu_def->P_atm) - alfa_w * (DevArraysPtr->T[local] - T_0));
	DevArraysPtr->ro_n[local] = gpu_def->ro0_n * (1. + (gpu_def->beta_n) * (DevArraysPtr->P_n[local] - gpu_def->P_atm) - alfa_n * (DevArraysPtr->T[local] - T_0));
	DevArraysPtr->ro_g[local] = gpu_def->ro0_g * (DevArraysPtr->P_g[local] / gpu_def->P_atm) * (T_0 / DevArraysPtr->T[local]);
#else
	DevArraysPtr->ro_w[local] = def.ro0_w * (1. + (def.beta_w) * (DevArraysPtr->P_w[local] - def.P_atm));
	DevArraysPtr->ro_n[local] = def.ro0_n * (1. + (def.beta_n) * (DevArraysPtr->P_n[local] - def.P_atm));
	DevArraysPtr->ro_g[local] = def.ro0_g * DevArraysPtr->P_g[local] / def.P_atm;
#endif
	device_test_positive(DevArraysPtr->ro_g[local], __FILE__, __LINE__);
	device_test_positive(DevArraysPtr->ro_w[local], __FILE__, __LINE__);
	device_test_positive(DevArraysPtr->ro_n[local], __FILE__, __LINE__);
}

// Вычисление координаты точки, через которую будет вычисляться значение на границе (i1, j1, k1)
__device__ int device_set_boundary_basic_coordinate(int i, int j, int k)
{
	int i1, j1, k1;

	i1 = i;
	j1 = j;
	k1 = k;

	if (i == 0)
	{
		i1 ++;
	}
	if (i == (gpu_def->locNx) - 1)
	{
		i1 --;
	}
	if (j == 0)
	{
		j1 ++;
	}
	if (j == (gpu_def->locNy) - 1)
	{
		j1 --;
	}
	if ((k == 0) && ((gpu_def->locNz) > 2))
	{
		k1 ++;
	}
	if ((k == (gpu_def->locNz) - 1) && ((gpu_def->locNz) > 2))
	{
		k1 --;
	}

	return (i1 + j1 * (gpu_def->locNx) + k1 * (gpu_def->locNx) * (gpu_def->locNy));
}

// Расчет центральной разности
__device__ double central_difference (double* ptr, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (*(ptr+1) - *(ptr-1) )/ (2. * (gpu_def->hx));	
		}
	case 'y':
		{
			return (*(ptr + gpu_def->locNx) - *(ptr - gpu_def->locNx) )/ (2. * (gpu_def->hy));
		}
	case 'z':
		{
			return (*(ptr + gpu_def->locNx * (gpu_def->locNy)) - *(ptr - gpu_def->locNx * (gpu_def->locNy)) )/ (2. * (gpu_def->hz));
		}
	default:
		{
			device_print_error("Axis of [central_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет центральной разности для произведения двух элементов структуры
__device__ double multi_central_difference (double* ptr1, double* ptr2, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return ((*(ptr1+1)) * (*(ptr2+1)) - (*(ptr1-1)) * (*(ptr2-1)) )/ (2. * (gpu_def->hx));
		}
	case 'y':
		{
			return ((*(ptr1+gpu_def->locNx)) * (*(ptr2+gpu_def->locNx)) - (*(ptr1-gpu_def->locNx)) * (*(ptr2-gpu_def->locNx)) )/ (2. * (gpu_def->hy));
		}
	case 'z':
		{
			return ((*(ptr1+gpu_def->locNx * (gpu_def->locNy))) * (*(ptr2+gpu_def->locNx * (gpu_def->locNy)))
				- (*(ptr1-gpu_def->locNx * (gpu_def->locNy))) * (*(ptr2-gpu_def->locNx * (gpu_def->locNy))) )/ (2. * (gpu_def->hz));
		}
	default:
		{
			device_print_error("Axis of [central_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет направленной разности
__device__ double directed_difference (double x1, double x2, double* Xi, double* ro, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
				(x1 + fabs(x1)) / 2. * (-1.) * (*(Xi-1)) * (*(ro-1)) +
				(x2 - fabs(x2)) / 2. * (-1.) * (*(Xi+1)) * (*(ro+1))) / gpu_def->hx;
		}
	case 'y':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
				(x1 + fabs(x1)) / 2. * (-1.) * (*(Xi - gpu_def->locNx)) * (*(ro - gpu_def->locNx)) +
				(x2 - fabs(x2)) / 2. * (-1.) * (*(Xi + gpu_def->locNx)) * (*(ro + gpu_def->locNx))) / gpu_def->hy;
		}
	case 'z':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
				(x1 + fabs(x1)) / 2. * (-1.) * (*(Xi - gpu_def->locNx * (gpu_def->locNy))) * (*(ro - gpu_def->locNx * (gpu_def->locNy))) +
				(x2 - fabs(x2)) / 2. * (-1.) * (*(Xi + gpu_def->locNx * (gpu_def->locNy))) * (*(ro + gpu_def->locNx * (gpu_def->locNy)))) / gpu_def->hz;
		}
	default:
		{
			device_print_error("Axis of [directed_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет левой разности
__device__ double left_difference (double* ptr, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (*ptr - *(ptr-1) )/ gpu_def->hx;	
		}
	case 'y':
		{
			return (*ptr - *(ptr-gpu_def->locNx) )/ gpu_def->hy;
		}
	case 'z':
		{
			return (*ptr - *(ptr - gpu_def->locNx * (gpu_def->locNy)) )/ gpu_def->hz;
		}
	default:
		{
			device_print_error("Axis of [left_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет правой разности
__device__ double right_difference (double* ptr, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (*(ptr+1) - *ptr )/ gpu_def->hx;	
		}
	case 'y':
		{
			return (*(ptr + gpu_def->locNx) - *ptr )/ gpu_def->hy;
		}
	case 'z':
		{
			return (*(ptr + gpu_def->locNx * (gpu_def->locNy)) - *ptr )/ gpu_def->hz;
		}
	default:
		{
			device_print_error("Axis of [right_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет divgrad для произведения двух элементов структуры
__device__ double multi_divgrad (double* ptr1, double* ptr2, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return ((*(ptr1+1)) * (*(ptr2+1)) - 2 * (*ptr1) * (*ptr2)
				+ (*(ptr1-1)) * (*(ptr2-1))) / ((gpu_def->hx) * (gpu_def->hx));
		}
	case 'y':
		{
			return ((*(ptr1+gpu_def->locNx)) * (*(ptr2+gpu_def->locNx)) - 2 * (*ptr1) * (*ptr2)
				+ (*(ptr1-gpu_def->locNx)) * (*(ptr2-gpu_def->locNx))) / ((gpu_def->hy) * (gpu_def->hy));
		}
	case 'z':
		{
			return ((*(ptr1+gpu_def->locNx * (gpu_def->locNy))) * (*(ptr2+gpu_def->locNx * (gpu_def->locNy))) - 2 * (*ptr1) * (*ptr2)
				+ (*(ptr1-gpu_def->locNx * (gpu_def->locNy))) * (*(ptr2-gpu_def->locNx * (gpu_def->locNy)))) / ((gpu_def->hz) * (gpu_def->hz));
		}
	default:
		{
			device_print_error("Axis of [right_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Преобразование локальных координат процессора к глобальным
// Каждый процессор содержит дополнительную точку в массиве для
// обмена данными, если имеет соседа
// (если 2 соседа с обеих сторон,то +2 точки).
// Глобальные границы хранятся как обычные точки (отсюда и условие на rank==0)
__device__ int device_local_to_global(int local_index, char axis)
{
	int global_index = local_index;
	switch (axis)
	{
		case 'x':
		{
			global_index += gpu_def->rankx * (gpu_def->Nx) / gpu_def->sizex + min(gpu_def->rankx, gpu_def->Nx % gpu_def->sizex);
			break;
		}
		case 'y':
		{
			global_index += gpu_def->ranky * (gpu_def->Ny) / gpu_def->sizey + min(gpu_def->ranky, gpu_def->Ny % gpu_def->sizey);
			break;
		}
		case 'z':
		{
			global_index += gpu_def->rankz * (gpu_def->Nz) / gpu_def->sizez + min(gpu_def->rankz, gpu_def->Nz % gpu_def->sizez);
			break;
		}
		default:
		{
			//CUPRINTF("Error!");
		}
	}
	//some_test(global_index);
	return global_index;
}

// Функция вычисления "эффективной" плотности
__device__ double device_ro_eff_gdy(int local)
{
	double ro_g_dy = (DevArraysPtr->ro_g[local] * (1. - DevArraysPtr->S_w[local] - DevArraysPtr->S_n[local]) + DevArraysPtr->ro_w[local] * DevArraysPtr->S_w[local]
	                  + DevArraysPtr->ro_n[local] * DevArraysPtr->S_n[local]) * (DevArraysPtr->m[local]) * (gpu_def->g_const) * (gpu_def->hy);
	return ro_g_dy;
}

void prepare_all_vars()
{
	prepare_local_vars_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>();
	checkErrors("assign P, Xi", __FILE__, __LINE__);
}

void solve_nonlinear_system()
{
	Newton_method_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>();
	checkErrors("assign Pw and Sn", __FILE__, __LINE__);
}

// Расчет скорости в каждой точке сетки
__global__ void assign_u_kernel()
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if (GPU_ACTIVE_POINT)
	{
		int local=i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
		if ((gpu_def->Nx) > 2)
		{
			if (i == 0)
			{
				DevArraysPtr->ux_w[local] = DevArraysPtr->Xi_w[local] * right_difference(DevArraysPtr->P_w+local, 'x');
				DevArraysPtr->ux_n[local] = DevArraysPtr->Xi_n[local] * right_difference(DevArraysPtr->P_n+local, 'x');
				DevArraysPtr->ux_g[local] = DevArraysPtr->Xi_g[local] * right_difference(DevArraysPtr->P_g+local, 'x');
			}
			else
			{
				if (i == (gpu_def->locNx) - 1)
				{
					DevArraysPtr->ux_w[local] = DevArraysPtr->Xi_w[local] * left_difference(DevArraysPtr->P_w+local, 'x');
					DevArraysPtr->ux_n[local] = DevArraysPtr->Xi_n[local] * left_difference(DevArraysPtr->P_n+local, 'x');
					DevArraysPtr->ux_g[local] = DevArraysPtr->Xi_g[local] * left_difference(DevArraysPtr->P_g+local, 'x');
				}
				else
				{
					DevArraysPtr->ux_w[local] = DevArraysPtr->Xi_w[local] * central_difference (DevArraysPtr->P_w+local, 'x');
					DevArraysPtr->ux_n[local] = DevArraysPtr->Xi_n[local] * central_difference (DevArraysPtr->P_n+local, 'x');
					DevArraysPtr->ux_g[local] = DevArraysPtr->Xi_g[local] * central_difference (DevArraysPtr->P_g+local, 'x');
				}
			}
		}
		else
		{
			DevArraysPtr->ux_w[local] = 0.;
			DevArraysPtr->ux_n[local] = 0.;
			DevArraysPtr->ux_g[local] = 0.;
		}

		if ((gpu_def->Ny) > 2)
		{
			if (j == 0)
			{
				DevArraysPtr->uy_w[local] = DevArraysPtr->Xi_w[local] * (right_difference (DevArraysPtr->P_w+local, 'y') - DevArraysPtr->ro_w[local] * (gpu_def->g_const));
				DevArraysPtr->uy_n[local] = DevArraysPtr->Xi_n[local] * (right_difference (DevArraysPtr->P_n+local, 'y') - DevArraysPtr->ro_n[local] * (gpu_def->g_const));
				DevArraysPtr->uy_g[local] = DevArraysPtr->Xi_g[local] * (right_difference (DevArraysPtr->P_g+local, 'y') - DevArraysPtr->ro_g[local] * (gpu_def->g_const));
			}
			else
			{
				if (j == (gpu_def->locNy) - 1)
				{
					DevArraysPtr->uy_w[local] = DevArraysPtr->Xi_w[local] * (left_difference (DevArraysPtr->P_w+local, 'y') - DevArraysPtr->ro_w[local] * (gpu_def->g_const));
					DevArraysPtr->uy_n[local] = DevArraysPtr->Xi_n[local] * (left_difference (DevArraysPtr->P_n+local, 'y') - DevArraysPtr->ro_n[local] * (gpu_def->g_const));
					DevArraysPtr->uy_g[local] = DevArraysPtr->Xi_g[local] * (left_difference (DevArraysPtr->P_g+local, 'y') - DevArraysPtr->ro_g[local] * (gpu_def->g_const));
				}
				else
				{
					DevArraysPtr->uy_w[local] = DevArraysPtr->Xi_w[local] * (central_difference (DevArraysPtr->P_w+local, 'y')	- DevArraysPtr->ro_w[local] * (gpu_def->g_const));
					DevArraysPtr->uy_n[local] = DevArraysPtr->Xi_n[local] * (central_difference (DevArraysPtr->P_n+local, 'y')	- DevArraysPtr->ro_n[local] * (gpu_def->g_const));
					DevArraysPtr->uy_g[local] = DevArraysPtr->Xi_g[local] * (central_difference (DevArraysPtr->P_g+local, 'y')	- DevArraysPtr->ro_g[local] * (gpu_def->g_const));
				}
			}
		}
		else
		{
			DevArraysPtr->uy_w[local] = 0.;
			DevArraysPtr->uy_n[local] = 0.;
			DevArraysPtr->uy_g[local] = 0.;
		}

		if ((gpu_def->Nz) > 2)
		{
			if (k == 0)
			{
				DevArraysPtr->uz_w[local] = DevArraysPtr->Xi_w[local] * right_difference (DevArraysPtr->P_w+local, 'z');
				DevArraysPtr->uz_n[local] = DevArraysPtr->Xi_n[local] * right_difference (DevArraysPtr->P_n+local, 'z');
				DevArraysPtr->uz_g[local] = DevArraysPtr->Xi_g[local] * right_difference (DevArraysPtr->P_g+local, 'z');
			}
			else
			{
				if (k == (gpu_def->locNz) - 1)
				{
					DevArraysPtr->uz_w[local] = DevArraysPtr->Xi_w[local] * left_difference (DevArraysPtr->P_w+local, 'z');
					DevArraysPtr->uz_n[local] = DevArraysPtr->Xi_n[local] * left_difference (DevArraysPtr->P_n+local, 'z');
					DevArraysPtr->uz_g[local] = DevArraysPtr->Xi_g[local] * left_difference (DevArraysPtr->P_g+local, 'z');
				}
				else
				{
					DevArraysPtr->uz_w[local] = DevArraysPtr->Xi_w[local] * central_difference (DevArraysPtr->P_w+local, 'z');
					DevArraysPtr->uz_n[local] = DevArraysPtr->Xi_n[local] * central_difference (DevArraysPtr->P_n+local, 'z');
					DevArraysPtr->uz_g[local] = DevArraysPtr->Xi_g[local] * central_difference (DevArraysPtr->P_g+local, 'z');
				}
			}
		}
		else
		{
			DevArraysPtr->uz_w[local] = 0.;
			DevArraysPtr->uz_n[local] = 0.;
			DevArraysPtr->uz_g[local] = 0.;
		}

		device_test_u(DevArraysPtr->ux_w[local], __FILE__, __LINE__);
		device_test_u(DevArraysPtr->ux_n[local], __FILE__, __LINE__);
		device_test_u(DevArraysPtr->uy_w[local], __FILE__, __LINE__);
		device_test_u(DevArraysPtr->uy_n[local], __FILE__, __LINE__);
		device_test_u(DevArraysPtr->uz_w[local], __FILE__, __LINE__);
		device_test_u(DevArraysPtr->uz_n[local], __FILE__, __LINE__);
		device_test_u(DevArraysPtr->ux_g[local], __FILE__, __LINE__);
		device_test_u(DevArraysPtr->uy_g[local], __FILE__, __LINE__);
		device_test_u(DevArraysPtr->uz_g[local], __FILE__, __LINE__);
	}
}

// Расчет скоростей во всех точках сетки
void u_calculation()
{
	assign_u_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>();
	checkErrors("assign u", __FILE__, __LINE__);
}

// Расчет вспомогательной насыщенности в каждой точке сетки
__device__ void device_assign_S(int local)
{
	DevArraysPtr->S_g[local] = 1. - DevArraysPtr->S_w[local] - DevArraysPtr->S_n[local];
}

// Расчет ro*S в каждой точке сетки методом направленных разностей
__global__ void assign_roS_kernel_nr(double t)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if (GPU_INTERNAL_POINT)
	{
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		if(! DevArraysPtr->m[local])
			return;

		double q_w = 0.;
		double q_n = 0.;
		double q_g = 0.;

		// Значения q на скважинах
		device_wells_q(i, j, k, &q_w, &q_n, &q_g);

		DevArraysPtr->roS_w[local] = DevArraysPtr->ro_w[local] * DevArraysPtr->S_w[local];
		DevArraysPtr->roS_g[local] = DevArraysPtr->ro_g[local]
		        * (1. - DevArraysPtr->S_w[local] - DevArraysPtr->S_n[local]);
		DevArraysPtr->roS_n[local] = DevArraysPtr->ro_n[local] * DevArraysPtr->S_n[local];

		double fx_g, fy_g, fz_g, A3 = 0.;

		double x1, x2, y1, y2, z1, z2, fx_w, fy_w, fz_w, fx_n, fy_n, fz_n, A1 = 0., A2 = 0.;

		if ((gpu_def->Nz) < 2)
		{
			fz_w = 0.;
			fz_n = 0.;
			fz_g = 0.;
		}
		else
		{
			z2 = -1. * right_difference (DevArraysPtr->P_w+local, 'z');
			z1 = -1. * left_difference (DevArraysPtr->P_w+local, 'z');
			fz_w = directed_difference (z1, z2, DevArraysPtr->Xi_w+local, DevArraysPtr->ro_w+local, 'z');

			z2 = -1. * right_difference (DevArraysPtr->P_n+local, 'z');
			z1 = -1. * left_difference (DevArraysPtr->P_n+local, 'z');
			fz_n = directed_difference (z1, z2, DevArraysPtr->Xi_n+local, DevArraysPtr->ro_n+local, 'z');

			z2 = -1. * right_difference (DevArraysPtr->P_g+local, 'z');
			z1 = -1. * left_difference (DevArraysPtr->P_g+local, 'z');
			fz_g = directed_difference (z1, z2, DevArraysPtr->Xi_g+local, DevArraysPtr->ro_g+local, 'z');

		}

		if ((gpu_def->Nx) < 2)
		{
			fx_w = 0.;
			fx_n = 0.;
			fx_g = 0.;
		}
		else
		{
			x2 = -1. * right_difference (DevArraysPtr->P_w+local, 'x'); //-(DevArraysPtr->P_w[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - Pw) / gpu_def->hx;
			x1 = -1. * left_difference (DevArraysPtr->P_w+local, 'x'); //-(Pw - DevArraysPtr->P_w[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) / gpu_def->hx;
			fx_w = directed_difference (x1, x2, DevArraysPtr->Xi_w+local, DevArraysPtr->ro_w+local, 'x');

			x2 = -1. * right_difference (DevArraysPtr->P_n+local, 'x');
			x1 = -1. * left_difference (DevArraysPtr->P_n+local, 'x');
			fx_n = directed_difference (x1, x2, DevArraysPtr->Xi_n+local, DevArraysPtr->ro_n+local, 'x');

			x2 = -1. * right_difference (DevArraysPtr->P_g+local, 'x');
			x1 = -1. * left_difference (DevArraysPtr->P_g+local, 'x');
			fx_g = directed_difference (x1, x2, DevArraysPtr->Xi_g+local, DevArraysPtr->ro_g+local, 'x');
		}

		y2 = -1. * right_difference (DevArraysPtr->P_w+local, 'y') + gpu_def->g_const * (DevArraysPtr->ro_w[local]);
		y1 = -1. * left_difference (DevArraysPtr->P_w+local, 'y') + gpu_def->g_const * (DevArraysPtr->ro_w[local]);
		fy_w = directed_difference (y1, y2, DevArraysPtr->Xi_w+local, DevArraysPtr->ro_w+local, 'y');

		y2 = -1. * right_difference (DevArraysPtr->P_n+local, 'y') + gpu_def->g_const * (DevArraysPtr->ro_n[local]);
		y1 = -1. * left_difference (DevArraysPtr->P_n+local, 'y') + gpu_def->g_const * (DevArraysPtr->ro_n[local]);
		fy_n = directed_difference (y1, y2, DevArraysPtr->Xi_n+local, DevArraysPtr->ro_n+local, 'y');

		A1 = DevArraysPtr->roS_w[local] - (gpu_def->dt / DevArraysPtr->m[local]) * (-q_w + fx_w + fy_w + fz_w);
		A2 = DevArraysPtr->roS_n[local] - (gpu_def->dt / DevArraysPtr->m[local]) * (-q_n + fx_n + fy_n + fz_n);

		DevArraysPtr->roS_w_old[local] = DevArraysPtr->roS_w[local];
		DevArraysPtr->roS_n_old[local] = DevArraysPtr->roS_n[local];
		DevArraysPtr->roS_w[local] = A1;
		DevArraysPtr->roS_n[local] = A2;

		device_test_positive(DevArraysPtr->roS_w[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr->roS_n[local], __FILE__, __LINE__);

		y2 = -1. * right_difference (DevArraysPtr->P_g+local, 'y') + gpu_def->g_const * (DevArraysPtr->ro_g[local]);
		y1 = -1. * left_difference (DevArraysPtr->P_g+local, 'y') + gpu_def->g_const * (DevArraysPtr->ro_g[local]);
		fy_g = directed_difference (y1, y2, DevArraysPtr->Xi_g+local, DevArraysPtr->ro_g+local, 'y');

		A3 = DevArraysPtr->roS_g[local] - (gpu_def->dt / DevArraysPtr->m[local]) * (q_g + fx_g + fy_g + fz_g);

		DevArraysPtr->roS_g_old[local] = DevArraysPtr->roS_g[local];
		DevArraysPtr->roS_g[local] = A3;

		device_test_positive(DevArraysPtr->roS_g[local], __FILE__, __LINE__);
	}
}

// Расчет ro*S в каждой точке сетки
__global__ void assign_roS_kernel(double t)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if (GPU_INTERNAL_POINT)
	{
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
		double divgrad1, divgrad2, Tx1, Ty1, Tz1, Tx2, Ty2, Tz2, A1 = 0, A2 = 0;

		double divgrad3, Tx3, Ty3, Tz3, A3 = 0;

		DevArraysPtr->roS_w[local] = DevArraysPtr->ro_w[local] * DevArraysPtr->S_w[local];
		DevArraysPtr->roS_n[local] = DevArraysPtr->ro_n[local] * DevArraysPtr->S_n[local];
		DevArraysPtr->roS_g[local] = DevArraysPtr->ro_g[local] * DevArraysPtr->S_g[local];

		if ((gpu_def->Nz) < 2)
		{
			divgrad1 = 0.;
			divgrad2 = 0.;
			divgrad3 = 0.;
			Tz1 = 0.;
			Tz2 = 0.;
			Tz3 = 0.;
		}
		else
		{
			divgrad1 = multi_divgrad (DevArraysPtr->ro_w + local, DevArraysPtr->S_w + local, 'z');
			divgrad2 = multi_divgrad (DevArraysPtr->ro_n + local, DevArraysPtr->S_n + local, 'z');
			divgrad3 = multi_divgrad (DevArraysPtr->ro_g + local, DevArraysPtr->S_g + local, 'z');

			Tz1 = multi_central_difference (DevArraysPtr->ro_w + local, DevArraysPtr->uz_w + local, 'z');
			Tz2 = multi_central_difference (DevArraysPtr->ro_n + local, DevArraysPtr->uz_n + local, 'z');
			Tz3 = multi_central_difference (DevArraysPtr->ro_g + local, DevArraysPtr->uz_g + local, 'z');
		}

		if ((gpu_def->Nx) < 2)
		{
			Tx1 = 0.;
			Tx2 = 0.;
			Tx3 = 0.;
			divgrad3 = 0.;
		}
		else
		{
			divgrad1 += multi_divgrad (DevArraysPtr->ro_w + local, DevArraysPtr->S_w + local, 'x');
			divgrad2 += multi_divgrad (DevArraysPtr->ro_n + local, DevArraysPtr->S_n + local, 'x');
			divgrad3 += multi_divgrad (DevArraysPtr->ro_g + local, DevArraysPtr->S_g + local, 'x');

			Tx1 = multi_central_difference (DevArraysPtr->ro_w + local, DevArraysPtr->uz_w + local, 'x');
			Tx2 = multi_central_difference (DevArraysPtr->ro_n + local, DevArraysPtr->uz_n + local, 'x');
			Tx3 = multi_central_difference (DevArraysPtr->ro_g + local, DevArraysPtr->uz_g + local, 'x');
		}

		divgrad1 += multi_divgrad (DevArraysPtr->ro_w + local, DevArraysPtr->S_w + local, 'y');
		divgrad1 *= DevArraysPtr->m[local] * (gpu_def->l) * (gpu_def->c_w);

		divgrad2 += multi_divgrad (DevArraysPtr->ro_n + local, DevArraysPtr->S_n + local, 'y');
		divgrad2 *= DevArraysPtr->m[local] * (gpu_def->l) * (gpu_def->c_n);

		divgrad3 += multi_divgrad (DevArraysPtr->ro_g + local, DevArraysPtr->S_g + local, 'y');
		divgrad3 *= DevArraysPtr->m[local] * (gpu_def->l) * (gpu_def->c_g);

		Ty1 = multi_central_difference (DevArraysPtr->ro_w + local, DevArraysPtr->uz_w + local, 'y');
		Ty2 = multi_central_difference (DevArraysPtr->ro_n + local, DevArraysPtr->uz_n + local, 'y');
		Ty3 = multi_central_difference (DevArraysPtr->ro_g + local, DevArraysPtr->uz_g + local, 'y');

		device_test_arrowhead(Tx1 + Ty1 + Tz1, divgrad1, __FILE__, __LINE__);
		device_test_arrowhead(Tx2 + Ty2 + Tz2, divgrad2, __FILE__, __LINE__);
		device_test_arrowhead(Tx3 + Ty3 + Tz3, divgrad3, __FILE__, __LINE__);

		double q_w = 0.;
		double q_n = 0.;
		double q_g = 0.;

		// Значения q на скважинах
		device_wells_q(i, j, k, &q_w, &q_n, &q_g);

		if ((t < 2 * (gpu_def->dt)) || TWO_LAYERS)
		{
			A1 = DevArraysPtr->roS_w[local] + ((gpu_def->dt) / DevArraysPtr->m[local]) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = DevArraysPtr->roS_n[local] + ((gpu_def->dt) / DevArraysPtr->m[local]) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2);
			A3 = DevArraysPtr->roS_g[local] + ((gpu_def->dt) / DevArraysPtr->m[local]) * (q_g + divgrad3 - Tx3 - Ty3 - Tz3);
		}
		else
		{
			A1 = (1. / ((DevArraysPtr->m[local]) * (gpu_def->dt) + 2. * (gpu_def->tau))) * (2. * (gpu_def->dt) * (gpu_def->dt) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1)
			        + ((DevArraysPtr->m[local]) * (gpu_def->dt) - 2. * (gpu_def->tau)) * DevArraysPtr->roS_w_old[local]
			        + 4. * (gpu_def->tau) * DevArraysPtr->roS_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]);
			A2 = (1. / ((DevArraysPtr->m[local]) * (gpu_def->dt) + 2. * (gpu_def->tau))) * (2. * (gpu_def->dt) * (gpu_def->dt) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2)
			        + ((DevArraysPtr->m[local]) * (gpu_def->dt) - 2. * (gpu_def->tau)) * DevArraysPtr->roS_n_old[local]
			        + 4. * (gpu_def->tau) * DevArraysPtr->roS_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]);

			A3 = (1. / ((DevArraysPtr->m[local]) * (gpu_def->dt) + 2. * (gpu_def->tau))) * (2. * (gpu_def->dt) * (gpu_def->dt) * (q_g + divgrad3 - Tx3 - Ty3 - Tz3)
			        + ((DevArraysPtr->m[local]) * (gpu_def->dt) - 2. * (gpu_def->tau)) * DevArraysPtr->roS_g_old[local]
			        + 4. * (gpu_def->tau) * DevArraysPtr->roS_g[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]);
		}

		DevArraysPtr->roS_w_old[local] = DevArraysPtr->roS_w[local];
		DevArraysPtr->roS_n_old[local] = DevArraysPtr->roS_n[local];
		DevArraysPtr->roS_g_old[local] = DevArraysPtr->roS_g[local];
		DevArraysPtr->roS_w[local] = A1;
		DevArraysPtr->roS_n[local] = A2;
		DevArraysPtr->roS_g[local] = A3;

		device_test_positive(DevArraysPtr->roS_w[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr->roS_n[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr->roS_g[local], __FILE__, __LINE__);
	}
}

void find_values_from_partial_equations(double t)
{
#ifdef NR
	assign_roS_kernel_nr <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(t);
#else
	assign_roS_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(t);
#endif
	checkErrors("assign roS", __FILE__, __LINE__);
#ifdef ENERGY
	assign_E_new_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>();
	checkErrors("assign E new", __FILE__, __LINE__);
#endif
}

// Применение граничных условий
void boundary_conditions()
{
	Border_S_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>();
	checkErrors("assign S", __FILE__, __LINE__);

	Border_P_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>();
	checkErrors("assign Pw", __FILE__, __LINE__);

#ifdef ENERGY
	Border_T_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>();
	checkErrors("assign T", __FILE__, __LINE__);
#endif
}

// Функция загрузки данных в память хоста
void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr)
{
	cudaMemcpy(HostArrayPtr, DevArrayPtr, (def.locNx) * (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
}

// Функция загрузки данных типа double в память ускорителя
void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr)
{
	cudaMemcpy(DevArrayPtr, HostArrayPtr, (def.locNx) * (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy double data to device", __FILE__, __LINE__);
}

// Функция загрузки данных типа int в память ускорителя
void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr)
{
	cudaMemcpy(DevArrayPtr, HostArrayPtr, (def.locNx) * (def.locNy) * (def.locNz)*sizeof(int), cudaMemcpyHostToDevice);
	checkErrors("copy int data to device", __FILE__, __LINE__);
}

// Выделение памяти ускорителя под массив точек расчетной области
void device_memory_allocation()
{
	int buffer_size = 0;

	if(def.sizex > 1)
		buffer_size = (def.locNy) * (def.locNz);
	if(def.sizey > 1 && (def.locNx) * (def.locNz) > buffer_size)
		buffer_size = (def.locNx) * (def.locNz);
	if(def.sizez > 1 && (def.locNx) * (def.locNy) > buffer_size)
		buffer_size = (def.locNx) * (def.locNy);

	if(buffer_size) {
	//	cudaMalloc((void**) &DevBufferLoc, buffer_size * sizeof(double));
	//	cudaMemcpyToSymbol(DevBuffer, &DevBufferLoc,  buffer_size * sizeof(double), 0, cudaMemcpyHostToDevice);
	}

	int sz = (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double);

	cudaMalloc((void**) & (DevArraysPtrLoc->P_w), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->P_n), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->S_n), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->ro_w), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->ro_n), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->ux_w), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->uy_w), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->uz_w), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->ux_n), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->uy_n), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->uz_n), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->Xi_w), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->Xi_n), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->roS_w), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->roS_n), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->roS_w_old), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->roS_n_old), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->m), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->K), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->S_w), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->P_g), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->S_g), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->ro_g), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->ux_g), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->uy_g), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->uz_g), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->Xi_g), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->roS_g), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->roS_g_old), sz);
#ifdef ENERGY
	cudaMalloc((void**) & (DevArraysPtrLoc->T), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->H_w), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->H_n), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->H_g), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->H_r), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->E), sz);
	cudaMalloc((void**) & (DevArraysPtrLoc->E_new), sz);
#endif
	ptr_Arrays *DevArraysTmp = new ptr_Arrays[1];
	DevArraysTmp[0] = *DevArraysPtrLoc;
	cudaMemcpyToSymbol(DevArraysPtr, DevArraysTmp, sizeof(ptr_Arrays));

	cudaMemcpy((double*)(DevArraysPtrLoc->P_w), HostArraysPtr.P_w, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->P_n), HostArraysPtr.P_n, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->S_n), HostArraysPtr.S_n, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->ro_w), HostArraysPtr.ro_w, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->ro_n), HostArraysPtr.ro_n, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->ux_w), HostArraysPtr.ux_w, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->uy_w), HostArraysPtr.uy_w, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->uz_w), HostArraysPtr.uz_w, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->ux_n), HostArraysPtr.ux_n, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->uy_n), HostArraysPtr.uy_n, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->uz_n), HostArraysPtr.uz_n, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->Xi_w), HostArraysPtr.Xi_w, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->Xi_n), HostArraysPtr.Xi_n, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->roS_w), HostArraysPtr.roS_w, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->roS_n), HostArraysPtr.roS_n, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->roS_w_old), HostArraysPtr.roS_w_old, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->roS_n_old), HostArraysPtr.roS_n_old, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->m), HostArraysPtr.m, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->K), HostArraysPtr.K, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->S_w), HostArraysPtr.S_w, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->P_g), HostArraysPtr.P_g, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->S_g), HostArraysPtr.S_g, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->ro_g), HostArraysPtr.ro_g, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->ux_g), HostArraysPtr.ux_g, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->uy_g), HostArraysPtr.uy_g, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->uz_g), HostArraysPtr.uz_g, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->Xi_g), HostArraysPtr.Xi_g, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->roS_g), HostArraysPtr.roS_g, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->roS_g_old), HostArraysPtr.roS_g_old, sz, cudaMemcpyHostToDevice);
#ifdef ENERGY
	cudaMemcpy((double*)(DevArraysPtrLoc->T), HostArraysPtr.T, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->H_w), HostArraysPtr.H_w, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->H_n), HostArraysPtr.H_n, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->H_g), HostArraysPtr.H_g, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->H_r), HostArraysPtr.H_r, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->E), HostArraysPtr.E, sz, cudaMemcpyHostToDevice);
	cudaMemcpy((double*)(DevArraysPtrLoc->E_new), HostArraysPtr.E_new, sz, cudaMemcpyHostToDevice);
#endif

	checkErrors("memory allocation", __FILE__, __LINE__);
}

// Освобожение памяти ускорителя из под массива точек расчетной области
void device_memory_free()
{
//	cudaFree(DevBufferLoc);
	cudaFree((double*)(DevArraysPtrLoc->P_w));
	cudaFree((double*)(DevArraysPtrLoc->P_n));
	cudaFree((double*)(DevArraysPtrLoc->S_n));
	cudaFree((double*)(DevArraysPtrLoc->ro_w));
	cudaFree((double*)(DevArraysPtrLoc->ro_n));
	cudaFree((double*)(DevArraysPtrLoc->ux_w));
	cudaFree((double*)(DevArraysPtrLoc->uy_w));
	cudaFree((double*)(DevArraysPtrLoc->uz_w));
	cudaFree((double*)(DevArraysPtrLoc->ux_n));
	cudaFree((double*)(DevArraysPtrLoc->uy_n));
	cudaFree((double*)(DevArraysPtrLoc->uz_n));
	cudaFree((double*)(DevArraysPtrLoc->Xi_w));
	cudaFree((double*)(DevArraysPtrLoc->Xi_n));
	cudaFree((double*)(DevArraysPtrLoc->roS_w));
	cudaFree((double*)(DevArraysPtrLoc->roS_n));
	cudaFree((double*)(DevArraysPtrLoc->roS_w_old));
	cudaFree((double*)(DevArraysPtrLoc->roS_n_old));
	cudaFree((double*)(DevArraysPtrLoc->m));
	cudaFree((double*)(DevArraysPtrLoc->K));
	cudaFree((double*)(DevArraysPtrLoc->S_w));
	cudaFree((double*)(DevArraysPtrLoc->P_g));
	cudaFree((double*)(DevArraysPtrLoc->S_g));
	cudaFree((double*)(DevArraysPtrLoc->ro_g));
	cudaFree((double*)(DevArraysPtrLoc->ux_g));
	cudaFree((double*)(DevArraysPtrLoc->uy_g));
	cudaFree((double*)(DevArraysPtrLoc->uz_g));
	cudaFree((double*)(DevArraysPtrLoc->Xi_g));
	cudaFree((double*)(DevArraysPtrLoc->roS_g));
	cudaFree((double*)(DevArraysPtrLoc->roS_g_old));
#ifdef ENERGY
	cudaFree((double*)(DevArraysPtrLoc->T));
	cudaFree((double*)(DevArraysPtrLoc->H_w));
	cudaFree((double*)(DevArraysPtrLoc->H_n));
	cudaFree((double*)(DevArraysPtrLoc->H_g));
	cudaFree((double*)(DevArraysPtrLoc->H_r));
	cudaFree((double*)(DevArraysPtrLoc->E));
	cudaFree((double*)(DevArraysPtrLoc->E_new));
#endif
	checkErrors("memory release", __FILE__, __LINE__);
}

// Инициализация ускорителя
// Расчет происходит на ускорителе, номер которого равен
// номеру запускающего процессора
void device_initialization()
{
	// Было бы очень неплохо вместо GPU_PER_NODE использовать cudaGetDeviceCount
	//int deviceCount;
	//cudaGetDeviceCount ( &deviceCount );

	// Считаем, что ядер на узле не меньше, чем ускорителей
	int device = def.rank % GPU_PER_NODE;
	cudaSetDevice(device);

	// Количество запускаемых блоков
	// Если число точек сетки не кратно размеру блока,
	// то количество блоков будет на 1 больше.
	def.blocksX = def.locNx / BlockNX;
	if ((def.locNx % BlockNX) != 0)
	{
		(def.blocksX)++;
	}
	def.blocksY = (def.locNy) / BlockNY;
	if ((def.locNy) % BlockNY != 0)
	{
		(def.blocksY)++;
	}
	def.blocksZ = (def.locNz) / BlockNZ;
	if ((def.locNz) % BlockNZ != 0)
	{
		(def.blocksZ)++;
	}

	consts* deff = new consts[1];
	deff[0] = def;
	cudaMemcpyToSymbol(gpu_def, deff, sizeof(consts));
	checkErrors("constant memory copy", __FILE__, __LINE__);

	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, device);

	if (devProp.major < 2)
	{
		printf("\nError! Compute capability < 2, rank=%d\n", def.rank);
	}

	if (!def.rank)
	{
		//printf ( "Device %d\n", device );
		printf("Name : %s\n", devProp.name);
		printf("Compute capability : %d.%d\n", devProp.major, devProp.minor);
		printf("Total Global Memory : %ld\n", devProp.totalGlobalMem);
		printf("Shared memory per block: %ld\n", devProp.sharedMemPerBlock);
		printf("Registers per block : %d\n", devProp.regsPerBlock);
		printf("Warp size : %d\n", devProp.warpSize);
		printf("Max threads per block : %d\n", devProp.maxThreadsPerBlock);
		printf("Total constant memory : %ld\n", devProp.totalConstMem);
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


		// Максимальный размер расчетной сетки для ускорителя
		// sizeof(ptr_Arrays)/4 - количество параметров в точке, т.к. 4 -размер одного указателя
		printf("\nTotal NAPL_Filtration grid size : %ld\n\n", devProp.totalGlobalMem / (sizeof(ptr_Arrays)*sizeof(double) / 4));
	}

	// (def.locNx)+2 потому что 2NyNz на буфер обмена выделяется
	// Нужно переписать!!! Учесть размер буфера правильно!!!
	if ((def.locNx + 2) * (def.locNy) * (def.locNz) > (devProp.totalGlobalMem / (sizeof(ptr_Arrays)*sizeof(double) / 4)))
	{
		printf("\nError! Not enough memory at GPU, rank=%d\n", def.rank);
	}
	fflush(stdout);
}

// Финализация ускорителя
void device_finalization(void)
{
}

__global__ void load_exchange_data_part_xl_kernel(double* DevArray)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < gpu_def->locNy && k < (gpu_def->locNz))
	{
		DevBuffer[j + (gpu_def->locNy)*k] = DevArray[1 + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k];
		device_test_nan(DevBuffer[j + (gpu_def->locNy)*k], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_xr_kernel(double* DevArray)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < gpu_def->locNy && k < (gpu_def->locNz))
	{
		DevBuffer[j + (gpu_def->locNy)*k] = DevArray[(gpu_def->locNx) - 2 + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k];
		device_test_nan(DevBuffer[j + (gpu_def->locNy)*k], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_yl_kernel(double* DevArray)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && k < (gpu_def->locNz))
	{
		DevBuffer[i + (gpu_def->locNx)*k] = DevArray[i + (gpu_def->locNx) + (gpu_def->locNx) * (gpu_def->locNy) * k];
		device_test_nan(DevBuffer[i + (gpu_def->locNx)*k], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_yr_kernel(double* DevArray)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && k < (gpu_def->locNz))
	{
		DevBuffer[i + (gpu_def->locNx)*k] = DevArray[i + (gpu_def->locNx) * (gpu_def->locNy - 2) + (gpu_def->locNx) * (gpu_def->locNy) * k];
		device_test_nan(DevBuffer[i + (gpu_def->locNx)*k], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_zl_kernel(double* DevArray)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && j < (gpu_def->locNy))
	{
		DevBuffer[i + (gpu_def->locNx)*j] = DevArray[i + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy)];
		device_test_nan(DevBuffer[i + (gpu_def->locNx)*j], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_zr_kernel(double* DevArray)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && j < (gpu_def->locNy))
	{
		DevBuffer[i + (gpu_def->locNx)*j] = DevArray[i + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * (gpu_def->locNz - 2)];
		device_test_nan(DevBuffer[i + (gpu_def->locNx)*j], __FILE__, __LINE__);
	}
}

void load_exchange_data_part_xl(double* DevArray)
{
/*	load_exchange_data_part_xl_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArray);
	checkErrors("load_exchange_data_part_xl", __FILE__, __LINE__);

	cudaMemcpy(HostBuffer, DevBufferLoc, (def.locNy) * (def.locNz) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
	*/
}

void load_exchange_data_part_xr(double* DevArray)
{
/*	load_exchange_data_part_xr_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArray);
	checkErrors("load_exchange_data_part_xr", __FILE__, __LINE__);

	cudaMemcpy(HostBuffer, DevBufferLoc, (def.locNy) * (def.locNz) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
*/}

void load_exchange_data_part_yl(double* DevArray)
{
/*	load_exchange_data_part_yl_kernel <<< dim3(def.blocksX, def.blocksZ), dim3(BlockNX, BlockNZ)>>>(DevArray);
	checkErrors("load_exchange_data_part_yl", __FILE__, __LINE__);

	cudaMemcpy(HostBuffer, DevBufferLoc, (def.locNx) * (def.locNz) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
*/}

void load_exchange_data_part_yr(double* DevArray)
{
/*	load_exchange_data_part_yr_kernel <<< dim3(def.blocksX, def.blocksZ), dim3(BlockNX, BlockNZ)>>>(DevArray);
	checkErrors("load_exchange_data_part_yr", __FILE__, __LINE__);

	cudaMemcpy(HostBuffer, DevBufferLoc, (def.locNx) * (def.locNz) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
*/}

void load_exchange_data_part_zl(double* DevArray)
{
/*	load_exchange_data_part_zl_kernel <<< dim3(def.blocksX, def.blocksY), dim3(BlockNX, BlockNY)>>>(DevArray);
	checkErrors("load_exchange_data_part_zl", __FILE__, __LINE__);
	
	cudaMemcpy(HostBuffer, DevBufferLoc, (def.locNx) * (def.locNy) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
*/}

void load_exchange_data_part_zr(double* DevArray)
{
/*	load_exchange_data_part_zr_kernel <<< dim3(def.blocksX, def.blocksY), dim3(BlockNX, BlockNY)>>>(DevArray);
	checkErrors("load_exchange_data_part_zr", __FILE__, __LINE__);

	cudaMemcpy(HostBuffer, DevBufferLoc, (def.locNx) * (def.locNy) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
*/}

__global__ void save_exchange_data_part_xl_kernel(double* DevArray)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < gpu_def->locNy && k < (gpu_def->locNz))
	{
		DevArray[(gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k] = DevBuffer[j + (gpu_def->locNy)*k];
		device_test_nan(DevArray[(gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_xr_kernel(double* DevArray)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < gpu_def->locNy && k < (gpu_def->locNz))
	{
		DevArray[(gpu_def->locNx) - 1 + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k] = DevBuffer[j + (gpu_def->locNy)*k];
		device_test_nan(DevArray[(gpu_def->locNx) - 1 + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_yl_kernel(double* DevArray)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && k < (gpu_def->locNz))
	{
		DevArray[i + (gpu_def->locNx) * (gpu_def->locNy) * k] = DevBuffer[i + (gpu_def->locNx)*k];
		device_test_nan(DevArray[i + (gpu_def->locNx) * (gpu_def->locNy) * k], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_yr_kernel(double* DevArray)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && k < (gpu_def->locNz))
	{
		DevArray[i + (gpu_def->locNx) * (gpu_def->locNy - 1) + (gpu_def->locNx) * (gpu_def->locNy) * k] = DevBuffer[i + (gpu_def->locNx)*k];
		device_test_nan(DevArray[i + (gpu_def->locNx) * (gpu_def->locNy - 1) + (gpu_def->locNx) * (gpu_def->locNy) * k], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_zl_kernel(double* DevArray)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && j < (gpu_def->locNy))
	{
		DevArray[i + (gpu_def->locNx) * j] = DevBuffer[i + (gpu_def->locNx)*j];
		device_test_nan(DevArray[i + (gpu_def->locNx) * j], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_zr_kernel(double* DevArray)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && j < (gpu_def->locNy))
	{
		DevArray[i + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * (gpu_def->locNz - 1)] = DevBuffer[i + (gpu_def->locNx)*j];
		device_test_nan(DevArray[i + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * (gpu_def->locNz - 1)], __FILE__, __LINE__);
	}
}

void save_exchange_data_part_xl(double* DevArray)
{
/*	cudaMemcpy(DevBufferLoc, HostBuffer, (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);

	save_exchange_data_part_xl_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArray);
	checkErrors("save_exchange_data_part_xl", __FILE__, __LINE__);
*/}

void save_exchange_data_part_xr(double* DevArray)
{
/*	cudaMemcpy(DevBufferLoc, HostBuffer, (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);

	save_exchange_data_part_xr_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArray);
	checkErrors("save_exchange_data_part_xr", __FILE__, __LINE__);
*/}

void save_exchange_data_part_yl(double* DevArray)
{
/*	cudaMemcpy(DevBufferLoc, HostBuffer, (def.locNx) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);

	save_exchange_data_part_yl_kernel <<< dim3(def.blocksX, def.blocksZ), dim3(BlockNX, BlockNZ)>>>(DevArray);
	checkErrors("save_exchange_data_part_yl", __FILE__, __LINE__);
*/}

void save_exchange_data_part_yr(double* DevArray)
{
/*	cudaMemcpy(DevBufferLoc, HostBuffer, (def.locNx) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);

	save_exchange_data_part_yr_kernel <<< dim3(def.blocksX, def.blocksZ), dim3(BlockNX, BlockNZ)>>>(DevArray);
	checkErrors("save_exchange_data_part_yr", __FILE__, __LINE__);
*/}

void save_exchange_data_part_zl(double* DevArray)
{
/*	cudaMemcpy(DevBufferLoc, HostBuffer, (def.locNx) * (def.locNy)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);

	save_exchange_data_part_zl_kernel <<< dim3(def.blocksX, def.blocksY), dim3(BlockNX, BlockNY)>>>(DevArray);
	checkErrors("save_exchange_data_part_zl", __FILE__, __LINE__);
*/}

void save_exchange_data_part_zr(double* DevArray)
{
/*	cudaMemcpy(DevBufferLoc, HostBuffer, (def.locNx) * (def.locNy)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);

	save_exchange_data_part_zr_kernel <<< dim3(def.blocksX, def.blocksY), dim3(BlockNX, BlockNY)>>>(DevArray);
	checkErrors("save_exchange_data_part_zr", __FILE__, __LINE__);
*/}
