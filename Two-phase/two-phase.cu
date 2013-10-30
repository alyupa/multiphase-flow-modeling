//#include "../defines.h"
#include "../gpu.h"
//#include "two-phase.h"

// Заглушка! Убрать как функция будет перенесена
void data_initialization(const ptr_Arrays &HostArraysPtr, long int* t, const consts &def)
{
	*t = 0;
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				if (is_active_point(i, j, k, def))
				{
					// Преобразование локальных координат процессора к глобальным
					int I = local_to_global(i, 'x', def);
					int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

					HostArraysPtr.m[local]=def.porosity[0];

					// Если точка на верхней границе, не далее (def.source) точек от центра,
					// то в ней начальная насыщенность. Иначе, нулевая
					if ((j == 0) && (I >= (def.Nx) / 2 - (def.source)) && (I <= (def.Nx) / 2 + (def.source)) && (k >= (def.Nz) / 2 - (def.source)) && (k <= (def.Nz) / 2 + (def.source)))
					{
						HostArraysPtr.S_n[local] = def.S_n_gr;
					}
					else
					{
						HostArraysPtr.S_n[local] = 0;
					}

					if (j == 0)
					{
						HostArraysPtr.P_w[local] = def.P_atm;
					}
					else
					{
						HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local - (def.locNx)] + ro_eff_gdy(HostArraysPtr, local - (def.locNx), def);
					}

					HostArraysPtr.ro_w[local] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm));

					///!!!! Не учитываются капиллярные силы! Или надо считать перед этим шагом P_n
					HostArraysPtr.ro_n[local] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[local] - def.P_atm));

					/*
					if ((HostArraysPtr.x[local]>=(def.NX)/2.*(def.h1)) && (HostArraysPtr.x[local]<=4.*(def.NX)/5.*(def.h1)))
						if ((HostArraysPtr.y[local]<=2./5.*def.locNy*(def.h2)) && (HostArraysPtr.y[local]>=(-1.)*HostArraysPtr.x[local]/4.+2./5.*def.locNy*(def.h2)))
							HostArraysPtr.media[local]=1;

					if ((HostArraysPtr.x[local]>=(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[local]<=2.*(def.NX)/5.*(def.h1)))
						if ((HostArraysPtr.y[local]<=4./5.*def.locNy*(def.h2)) && (HostArraysPtr.y[local]>=3./5.*def.locNy*(def.h2)))
							HostArraysPtr.media[local]=1;
							*/

					/*
					if ((HostArraysPtr.x[local]>=2.*(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[local]<=3.*(def.NX)/5.*(def.h1)))
						if ((HostArraysPtr.y[local]>=1./10.*def.locNy*(def.h2)) && (HostArraysPtr.y[local]<=3./10.*def.locNy*(def.h2)))
							HostArraysPtr.media[local]=1;
					*/

					test_nan(HostArraysPtr.S_n[local], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[local], __FILE__, __LINE__);
					test_nan(HostArraysPtr.m[local], __FILE__, __LINE__);
				}
}

// Расчет плотностей, давления NAPL P2 и Xi в каждой точке сетки (независимо от остальных точек)
__global__ void assign_P_Xi_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx)) && (j < (gpu_def->locNy)) && (k < (gpu_def->locNz)) && (device_is_active_point(i, j, k) == 1))
	{
		int media = 0;
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
		double S_n = DevArraysPtr.S_n[local];
		double P_w = DevArraysPtr.P_w[local];

		double S_e = (1. - S_n - gpu_def->S_wr[media]) / (1. - gpu_def->S_wr[media]);
		double k_w = pow(S_e, (2. + 3. * gpu_def->lambda[media]) / gpu_def->lambda[media]);
		double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + gpu_def->lambda[media]) / gpu_def->lambda[media]));
		double P_k = gpu_def->P_d[media] * pow((1. - S_n - gpu_def->S_wr[media]) / (1. - gpu_def->S_wr[media]), -1. / gpu_def->lambda[media]);

		DevArraysPtr.P_n[local] = P_w + P_k;
		DevArraysPtr.Xi_w[local] = -1 * gpu_def->K[media] * k_w / gpu_def->mu_w;
		DevArraysPtr.Xi_n[local] = -1 * gpu_def->K[media] * k_n / gpu_def->mu_n;

		device_test_positive(DevArraysPtr.P_n[local], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_w[local], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_n[local], __FILE__, __LINE__);
	}
}


// Метод Ньютона для каждой точки сетки (независимо от остальных точек)
__global__ void Newton_method_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx) - 1) && (j < gpu_def->locNy - 1) && (k < (gpu_def->locNz)) && (i != 0) && (j != 0) && (((k != 0) && (k != (gpu_def->locNz) - 1)) || ((gpu_def->locNz) < 2)))
	{
		int media = 0;
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
		double S_e, S_n, P_w, P_k, AAA, F1, F2, PkS, F1P, F2P, F1S, F2S, det;

		for (int w = 1; w <= gpu_def->newton_iterations; w++)
		{
			S_n = DevArraysPtr.S_n[local];
			P_w = DevArraysPtr.P_w[local];

			S_e = (1 - S_n - gpu_def->S_wr[media]) / (1 - gpu_def->S_wr[media]);
			P_k = gpu_def->P_d[media] * pow(S_e, -1 / gpu_def->lambda[media]);
			AAA = pow(S_e, ((-1 / gpu_def->lambda[media]) - 1));
			F1 = gpu_def->ro0_w * (1 + (gpu_def->beta_w) * (P_w - gpu_def->P_atm)) * (1 - S_n) - DevArraysPtr.roS_w[local];
			F2 = gpu_def->ro0_n * (1 + (gpu_def->beta_n) * (P_w + P_k - gpu_def->P_atm)) * S_n - DevArraysPtr.roS_n[local];

			PkS = AAA * gpu_def->P_d[media] / (gpu_def->lambda[media] * (1 - gpu_def->S_wr[media]));
			F1P = gpu_def->ro0_w * (gpu_def->beta_w) * (1 - S_n);
			F2P = gpu_def->ro0_n * (gpu_def->beta_n) * S_n;
			F1S = (-1) * gpu_def->ro0_w * (1 + (gpu_def->beta_w) * (P_w - gpu_def->P_atm));
			F2S = gpu_def->ro0_n * (1 + (gpu_def->beta_n) * (P_w + P_k - gpu_def->P_atm + (S_n * PkS)));

			det = F1P * F2S - F1S * F2P;

			DevArraysPtr.P_w[local] = P_w - (1 / det) * (F2S * F1 - F1S * F2);
			DevArraysPtr.S_n[local] = S_n - (1 / det) * (F1P * F2 - F2P * F1);
		}

		device_test_positive(DevArraysPtr.P_w[local], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_n[local], __FILE__, __LINE__);
	}
}

// Задание граничных условий
__global__ void Border_S_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < gpu_def->locNx) && (j < gpu_def->locNy) && (k < gpu_def->locNz))
		if (((i == 0) || (i == (gpu_def->locNx) - 1) || (j == 0) || (j == (gpu_def->locNy) - 1) ||
			(((k == 0) || (k == (gpu_def->locNz) - 1)) && ((gpu_def->locNz) >= 2))) && (device_is_active_point(i, j, k) == 1))
	{
		int local1 = device_set_boundary_basic_coordinate(i, j, k);
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		if ((j != 0) || ((gpu_def->source) <= 0))
		{
			DevArraysPtr.S_n[local] = DevArraysPtr.S_n[local1];
		}

		if ((j == 0) && ((gpu_def->source) > 0))
		{
			int I = device_local_to_global(i, 'x');
			if ((I >= (gpu_def->Nx) / 2 - (gpu_def->source)) && (I <= (gpu_def->Nx) / 2 + (gpu_def->source)) && (k >= (gpu_def->Nz) / 2 - (gpu_def->source)) && (k <= (gpu_def->Nz) / 2 + (gpu_def->source)))
			{
				DevArraysPtr.S_n[local] = gpu_def->S_n_gr;
			}
			else
				//DevArraysPtr.S_n[local] = 0;
			{
				DevArraysPtr.S_n[local] = DevArraysPtr.S_n[local1];
			}
		}

		device_test_S(DevArraysPtr.S_n[local], __FILE__, __LINE__);
	}
}

__global__ void Border_P_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < gpu_def->locNx) && (j < gpu_def->locNy) && (k < gpu_def->locNz))
		if (((i == 0) || (i == (gpu_def->locNx) - 1) || (j == 0) || (j == (gpu_def->locNy) - 1) ||
			(((k == 0) || (k == (gpu_def->locNz) - 1)) && ((gpu_def->locNz) >= 2))) && (device_is_active_point(i, j, k) == 1))
	{
		int local1 = device_set_boundary_basic_coordinate(i, j, k);
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		if ((j != 0) && (j != (gpu_def->locNy) - 1))
		{
			DevArraysPtr.P_w[local] = DevArraysPtr.P_w[local1];
		}
		else if (j == 0)
		{
			DevArraysPtr.P_w[local] = gpu_def->P_atm;
		}
		else
		{
			DevArraysPtr.P_w[local] = DevArraysPtr.P_w[local1] + device_ro_eff_gdy(DevArraysPtr, local1);
		}

		device_test_positive(DevArraysPtr.P_w[local], __FILE__, __LINE__);
	}
}

// Является ли точка нагнетательной скважиной
__device__ int device_is_injection_well(int i, int j, int k)
{
	return 0;
}

// Является ли точка добывающей скважиной
__device__ int device_is_output_well(int i, int j, int k)
{
	return 0;
}

// Устанавливает значения втекаемых/вытекаемых жидкостей q_i на скважинах
__device__ void device_wells_q(ptr_Arrays DevArraysPtr, int i, int j, int k, double* q_w, double* q_n, double* q_g)
{
	*q_w = 0.0;
	*q_g = 0.0;
	*q_n = 0.0;
}