#include "gpu.h"

// Коэффициенты удельных теплоемкостей при постоянном давлении  для water, napl, gas and rock в Вт/(м*К)
__device__ double c_w (double T)
{
	return gpu_def->c0_w - gpu_def->C_w * (T - gpu_def->T_0) + gpu_def->C_w2 * (T - gpu_def->T_0) * (T - gpu_def->T_0);
}

__device__ double c_n (double T)
{
	return gpu_def->c0_n + gpu_def->C_n * (T - gpu_def->T_0);
}

__device__ double c_g (double T)
{
	return gpu_def->c0_g + gpu_def->C_g * (T - gpu_def->T_0);
}

__device__ double c_r (double T)
{
	return gpu_def->c0_r + gpu_def->C_r * (T - gpu_def->T_0);
}

// Коэффициенты теплопроводности для water, napl, gas and rock
__device__ double lambda_w (double T)
{
	return gpu_def->lambda0_w * (1 - gpu_def->lambdaA_w * (T - gpu_def->T_0));
}

__device__ double lambda_n (double T)
{
	return gpu_def->lambda0_n * (1 - gpu_def->lambdaA_n * (T - gpu_def->T_0));
}

__device__ double lambda_g (double T)
{
	return gpu_def->lambda0_g * pow((T / gpu_def->T_0), gpu_def->lambdaA_g);
}

__device__ double lambda_r (double T)
{
	return gpu_def->lambda0_r;
}

// Эффективный коэффициент теплопроводности в точке (будет использоваться при расчете теплового потока)
__device__ double assign_lambda_eff (int local)
{
	return DevArraysPtr->m[local] * (DevArraysPtr->S_w[local] * lambda_w (DevArraysPtr->T[local])
		+ DevArraysPtr->S_n[local] * lambda_n (DevArraysPtr->T[local])
		+ DevArraysPtr->S_g[local] * lambda_g (DevArraysPtr->T[local]))
		+ (1. - DevArraysPtr->m[local]) * lambda_r (DevArraysPtr->T[local]);
}

// Расчет энтальпии по температуре и теплоемкости
// !!! Переписать, задав точность для метода Симпсона и передавая указатель на функцию, чтобы не копировать одно и то же
__device__ double assign_H_w (double T)
{
	/* Возможно, что Симпсон понадобиться позже, а пока все равно нужно знать явную зависимость энергии от температуры
	double integral = 0, sum = 0, h_temp;
	int N_temp = 50;

	h_temp = (T - gpu_def->T_0) / N_temp;
	
	integral += (gpu_def->P_atm / gpu_def->ro0_w);
	integral += с_w(gpu_def->T_0);
	integral += с_w(T);

	for(int i = 2; i < N_temp; i+=2)
		sum += с_w(gpu_def->T_0 + i * h_temp);

	sum *= 2;
	integral += sum;
	sum = 0;

	for(int i = 1; i < N_temp; i+=2)
		sum += с_w(gpu_def->T_0 + i * h_temp);

	sum *= 4;
	integral += sum;

	h_temp /= 3;
	integral *= h_temp;

	return integral;
	*/
	return (gpu_def->P_atm / gpu_def->ro0_w) + (T - gpu_def->T_0) * (gpu_def->c0_w - (T - gpu_def->T_0) \
			* (gpu_def->C_w / 2 + gpu_def->C_w2 * (T - gpu_def->T_0) / 3));
}

__device__ double assign_H_n (double T)
{
	return (gpu_def->P_atm / gpu_def->ro0_n) + (T - gpu_def->T_0) * (gpu_def->c0_n + gpu_def->C_n * (T - gpu_def->T_0) / 2);
}

__device__ double assign_H_g (double T)
{
	return (gpu_def->P_atm / gpu_def->ro0_g) + (T - gpu_def->T_0) * (gpu_def->c0_g + gpu_def->C_g * (T - gpu_def->T_0) / 2);
}

__device__ double assign_H_r (double T)
{
	return (gpu_def->P_atm / gpu_def->ro_r) + (T - gpu_def->T_0) * (gpu_def->c0_r + gpu_def->C_r * (T - gpu_def->T_0) / 2);
}

__device__ void device_assign_H (int local)
{
	DevArraysPtr->H_w[local] = assign_H_w (DevArraysPtr->T[local]);
	DevArraysPtr->H_n[local] = assign_H_n (DevArraysPtr->T[local]);
	DevArraysPtr->H_g[local] = assign_H_g (DevArraysPtr->T[local]);
	DevArraysPtr->H_r[local] = assign_H_r (DevArraysPtr->T[local]);

	device_test_nan(DevArraysPtr->H_w[local], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr->H_n[local], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr->H_g[local], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr->H_r[local], __FILE__, __LINE__);
}

// Возвращает значение плотности в точке фазы phase как функции от P, T
__device__ double ro(double P, double T, char phase)
{
	double result_ro;
	switch (phase)
	{
	case 'w':
		result_ro = gpu_def->ro0_w * (1. + (gpu_def->beta_w) * (P - gpu_def->P_atm) - gpu_def->alfa_w * (T - gpu_def->T_0));
		break;
	case 'n':
		result_ro = gpu_def->ro0_n * (1. + (gpu_def->beta_n) * (P - gpu_def->P_atm) - gpu_def->alfa_n * (T - gpu_def->T_0));
		break;
	case 'g':
		result_ro = gpu_def->ro0_g * (P / gpu_def->P_atm) * (gpu_def->T_0 / T);
		break;
	default:
		printf ("Wrong phase in function ro!\n");
		break;
	}

//	device_test_ro(result_ro, __FILE__, __LINE__);
//	device_test_ro(result_ro, __FILE__, __LINE__);
//	device_test_ro(result_ro, __FILE__, __LINE__);

	return result_ro;
}

// Возвращает значение частной производной плотности в точке фазы phase по переменной var
__device__ double d_ro(double P, double T, char phase, char var)
{
	double result_d_ro = 0;
	switch (phase)
	{
	case 'w':
		if (var == 'P')
		{
			result_d_ro = gpu_def->ro0_w * (gpu_def->beta_w);
		} 
		else if (var == 'T')
		{
			result_d_ro = (-1) * gpu_def->ro0_w * gpu_def->alfa_w;
		}
		else 
		{
			printf("Wrong var in function d_ro!\n");
		}
		break;
	case 'n':
		if (var == 'P')
		{
			result_d_ro = gpu_def->ro0_n * (gpu_def->beta_n);
		} 
		else if (var == 'T')
		{
			result_d_ro = (-1) * gpu_def->ro0_n * gpu_def->alfa_n;
		}
		else 
		{
			printf("Wrong var in function d_ro!\n");
		}
		break;
	case 'g':
		if (var == 'P')
		{
			result_d_ro = (gpu_def->ro0_g / gpu_def->P_atm) * (gpu_def->T_0 / T);
		} 
		else if (var == 'T')
		{
			result_d_ro = gpu_def->ro0_g * (P / gpu_def->P_atm) * (-1) * (gpu_def->T_0 / T) / T;
		}
		else 
		{
			printf("Wrong var in function d_ro!\n");
		}
		break;
	default:
		printf ("Wrong phase in function d_ro!\n");
		break;
	}

	//	device_test_ro(result_d_ro, __FILE__, __LINE__);
	//	device_test_ro(result_d_ro, __FILE__, __LINE__);
	//	device_test_ro(result_d_ro, __FILE__, __LINE__);

	return result_d_ro;
}

// Коэффициенты вязкости для water, napl, gas and rock

// Расчет теплового потока в точке
__device__ double assign_T_flow (int i, int j, int k)
{
	if (GPU_INTERNAL_POINT)
	{
		double T_flow = 0;
		int local=i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		if ((gpu_def->locNx) > 2)
		{
			T_flow += (assign_lambda_eff(local + 1) * DevArraysPtr->T[local + 1]
			- 2 * assign_lambda_eff(local) * DevArraysPtr->T[local]
			+ assign_lambda_eff(local - 1) * DevArraysPtr->T[local - 1]) / ((gpu_def->hx) * (gpu_def->hx));
		}
		if ((gpu_def->locNy) > 2)
		{
			T_flow += (assign_lambda_eff(local + gpu_def->locNx) * DevArraysPtr->T[local + gpu_def->locNx]
			- 2 * assign_lambda_eff(local) * DevArraysPtr->T[local]
			+ assign_lambda_eff(local - gpu_def->locNx) * DevArraysPtr->T[local - gpu_def->locNx]) / ((gpu_def->hy) * (gpu_def->hy));
		}

		if ((gpu_def->locNz) > 2)
		{
			T_flow = (assign_lambda_eff(local + (gpu_def->locNx) * (gpu_def->locNy)) * DevArraysPtr->T[local + (gpu_def->locNx) * (gpu_def->locNy)]
			- 2 * assign_lambda_eff(local) * DevArraysPtr->T[local]
			+ assign_lambda_eff(local - (gpu_def->locNx) * (gpu_def->locNy)) * DevArraysPtr->T[local - (gpu_def->locNx) * (gpu_def->locNy)]) / ((gpu_def->hz) * (gpu_def->hz));
		}

		device_test_u(T_flow, __FILE__, __LINE__);
		return T_flow;
	}
	else
		return 0;
}

// Расчет направленной разности
__device__ double directed_difference_E (double* P, double* Xi, double* ro, double* H, char axis)
{
	double x1 = 0, x2 = 0;
	switch (axis)
	{
	case 'x':
		{
			x2 = -device_right_difference (P, 'x');
			x1 = -device_left_difference (P, 'x');
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (*Xi) * (*ro) * (*H) -
		      (x1 + fabs(x1)) / 2. * (*(Xi-1)) * (*(ro-1)) * (*(H-1)) +
		      (x2 - fabs(x2)) / 2. * (*(Xi+1)) * (*(ro+1)) * (*(H+1))) / gpu_def->hx * (-1.0);
		}
	case 'y':
		{
			x2 = -device_right_difference (P, 'y') + gpu_def->g_const * (*ro);
			x1 = -device_left_difference (P, 'y') + gpu_def->g_const * (*ro);
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (*Xi) * (*ro) * (*H) -
		      (x1 + fabs(x1)) / 2. * (*(Xi-gpu_def->locNx)) * (*(ro-gpu_def->locNx)) * (*(H-gpu_def->locNx)) +
		      (x2 - fabs(x2)) / 2. * (*(Xi+gpu_def->locNx)) * (*(ro+gpu_def->locNx)) * (*(H+gpu_def->locNx))) / gpu_def->hy * (-1.0);
		}
	case 'z':
		{
			x2 = -device_right_difference (P, 'z');
			x1 = -device_left_difference (P, 'z');
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (*Xi) * (*ro) * (*H) -
		      (x1 + fabs(x1)) / 2. * (*(Xi-gpu_def->locNx * (gpu_def->locNy))) * (*(ro-gpu_def->locNx * (gpu_def->locNy))) * (*(H-gpu_def->locNx * (gpu_def->locNy))) +
		      (x2 - fabs(x2)) / 2. * (*(Xi+gpu_def->locNx * (gpu_def->locNy))) * (*(ro+gpu_def->locNx * (gpu_def->locNy))) * (*(H-gpu_def->locNx * (gpu_def->locNy)))) / gpu_def->hz * (-1.0);
		}
	default:
		{
			device_print_error("wrong axis", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет потока энергии в точке
__device__ double assign_E_flow (int i, int j, int k)
{
	if (GPU_INTERNAL_POINT)
	{
		double E_flow = 0;
		int local=i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
#ifdef NR
		if ((gpu_def->locNx) > 2)
		{
			E_flow += (directed_difference_E (DevArraysPtr->P_w+local, DevArraysPtr->Xi_w+local, DevArraysPtr->ro_w+local, DevArraysPtr->H_w+local, 'x')
					 + directed_difference_E (DevArraysPtr->P_n+local, DevArraysPtr->Xi_n+local, DevArraysPtr->ro_n+local, DevArraysPtr->H_n+local, 'x')
					 + directed_difference_E (DevArraysPtr->P_g+local, DevArraysPtr->Xi_g+local, DevArraysPtr->ro_g+local, DevArraysPtr->H_g+local, 'x'));
		}
		if ((gpu_def->locNy) > 2)
		{
			E_flow += (directed_difference_E (DevArraysPtr->P_w+local, DevArraysPtr->Xi_w+local, DevArraysPtr->ro_w+local, DevArraysPtr->H_w+local, 'y')
					 + directed_difference_E (DevArraysPtr->P_n+local, DevArraysPtr->Xi_n+local, DevArraysPtr->ro_n+local, DevArraysPtr->H_n+local, 'y')
					 + directed_difference_E (DevArraysPtr->P_g+local, DevArraysPtr->Xi_g+local, DevArraysPtr->ro_g+local, DevArraysPtr->H_g+local, 'y'));
		}
		if ((gpu_def->locNz) > 2)
		{
			E_flow += (directed_difference_E (DevArraysPtr->P_w+local, DevArraysPtr->Xi_w+local, DevArraysPtr->ro_w+local, DevArraysPtr->H_w+local, 'z')
					 + directed_difference_E (DevArraysPtr->P_n+local, DevArraysPtr->Xi_n+local, DevArraysPtr->ro_n+local, DevArraysPtr->H_n+local, 'z')
					 + directed_difference_E (DevArraysPtr->P_g+local, DevArraysPtr->Xi_g+local, DevArraysPtr->ro_g+local, DevArraysPtr->H_g+local, 'z'));
		}
#else
		if ((gpu_def->locNx) > 2)
		{
			E_flow += (DevArraysPtr->ro_w[local + 1] * DevArraysPtr->H_w[local + 1] * DevArraysPtr->ux_w[local + 1]
				- DevArraysPtr->ro_w[local - 1] * DevArraysPtr->H_w[local - 1] * DevArraysPtr->ux_w[local - 1]
				+ DevArraysPtr->ro_n[local + 1] * DevArraysPtr->H_n[local + 1] * DevArraysPtr->ux_n[local + 1]
				- DevArraysPtr->ro_n[local - 1] * DevArraysPtr->H_n[local - 1] * DevArraysPtr->ux_n[local - 1]
				+ DevArraysPtr->ro_g[local + 1] * DevArraysPtr->H_g[local + 1] * DevArraysPtr->ux_g[local + 1]
				- DevArraysPtr->ro_g[local - 1] * DevArraysPtr->H_g[local - 1] * DevArraysPtr->ux_g[local - 1]
				) / (2. * (gpu_def->hx));
		}
		if ((gpu_def->locNy) > 2)
		{
			E_flow += (DevArraysPtr->ro_w[local + gpu_def->locNx] * DevArraysPtr->H_w[local + gpu_def->locNx] * DevArraysPtr->uy_w[local + gpu_def->locNx]
			- DevArraysPtr->ro_w[local - gpu_def->locNx] * DevArraysPtr->H_w[local - gpu_def->locNx] * DevArraysPtr->uy_w[local - gpu_def->locNx]
			+ DevArraysPtr->ro_n[local + gpu_def->locNx] * DevArraysPtr->H_n[local + gpu_def->locNx] * DevArraysPtr->uy_n[local + gpu_def->locNx]
			- DevArraysPtr->ro_n[local - gpu_def->locNx] * DevArraysPtr->H_n[local - gpu_def->locNx] * DevArraysPtr->uy_n[local - gpu_def->locNx]
			+ DevArraysPtr->ro_g[local + gpu_def->locNx] * DevArraysPtr->H_g[local + gpu_def->locNx] * DevArraysPtr->uy_g[local + gpu_def->locNx]
			- DevArraysPtr->ro_g[local - gpu_def->locNx] * DevArraysPtr->H_g[local - gpu_def->locNx] * DevArraysPtr->uy_g[local - gpu_def->locNx]
			)/ (2. * (gpu_def->hy));	
		}

		if ((gpu_def->locNz) > 2)
		{
			E_flow += (DevArraysPtr->ro_w[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->H_w[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->uy_w[local + (gpu_def->locNx) * (gpu_def->locNy)]
			- DevArraysPtr->ro_w[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->H_w[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->uy_w[local - (gpu_def->locNx) * (gpu_def->locNy)]
			+ DevArraysPtr->ro_n[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->H_n[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->uy_n[local + (gpu_def->locNx) * (gpu_def->locNy)]
			- DevArraysPtr->ro_n[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->H_n[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->uy_n[local - (gpu_def->locNx) * (gpu_def->locNy)]
			+ DevArraysPtr->ro_g[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->H_g[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->uy_g[local + (gpu_def->locNx) * (gpu_def->locNy)]
			- DevArraysPtr->ro_g[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->H_g[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr->uy_g[local - (gpu_def->locNx) * (gpu_def->locNy)]
			)/ (2. * (gpu_def->hz));	
		}
#endif
		device_test_u(E_flow, __FILE__, __LINE__);
		return E_flow;
	}
	else
		return 0;
}

// Расчет внутренней энергии всей системы в точке
__device__ void device_assign_E_current (int local)
{
	DevArraysPtr->E[local] = (DevArraysPtr->m[local] * (DevArraysPtr->S_w[local] * (DevArraysPtr->ro_w[local] * DevArraysPtr->H_w[local] - DevArraysPtr->P_w[local])
		+ DevArraysPtr->S_n[local] * (DevArraysPtr->ro_n[local] * DevArraysPtr->H_n[local] - DevArraysPtr->P_n[local])
		+ DevArraysPtr->S_g[local] * (DevArraysPtr->ro_g[local] * DevArraysPtr->H_g[local] - DevArraysPtr->P_g[local]))
		+ (1. - DevArraysPtr->m[local]) * (gpu_def->ro_r * DevArraysPtr->H_r[local] - DevArraysPtr->P_w[local]));

	device_test_nan(DevArraysPtr->E[local], __FILE__, __LINE__);
}

// Расчет внутренней энергии всей системы в точке на следующем шаге по времени
__global__ void assign_E_new_kernel ()
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (GPU_INTERNAL_POINT)
	{
		double Q_hw = 0, Q_hr = 0; // источниковые члены

		int local=i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		DevArraysPtr->E_new[local] = DevArraysPtr->E[local] + (gpu_def->dt) * (assign_T_flow(i, j, k) + Q_hw + Q_hr - assign_E_flow(i, j, k));

		device_test_nan(DevArraysPtr->E_new[local], __FILE__, __LINE__);
	}
}

// Задание граничных условий на температуру
__global__ void Border_T_kernel()
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (GPU_BOUNDARY_POINT)
	{
		int local1 = device_set_boundary_basic_coordinate(i, j, k);
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		if (j == 0)
		{
			DevArraysPtr->T[local] = 400;
		}
		else if(j == (gpu_def->locNy) - 1)
		{
			DevArraysPtr->T[local] = 273;
		}
		else
		{
			// Будем считать границы области не теплопроводящими
			DevArraysPtr->T[local] = DevArraysPtr->T[local1];
		}

		device_test_positive(DevArraysPtr->T[local], __FILE__, __LINE__);
	}
}

// Расчет методом Ньютона значений переменных на новом шаге по времени, когда учитываем изменение энергии (случай 4х переменных)
// !!! Пока "выбросим" капиллярные давления
__global__ void Newton_method_kernel()
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (GPU_INTERNAL_POINT)
	{
		enum { n = 5 }; // Размерность системы
		double F[n]; // Вектор значений функций (из системы уравнений)
		double correction[n]; // Вектор поправок к функциям
		double dF[n*n]; // Матрица Якоби (в виде одномерного массива)

		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		for (int w = 1; w <= gpu_def->newton_iterations; w++)
		{
			F[0] = DevArraysPtr->S_g[local] + DevArraysPtr->S_w[local] + DevArraysPtr->S_n[local] - 1.;
			F[1] = ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'w') * DevArraysPtr->S_w[local] - DevArraysPtr->roS_w[local];
			F[2] = ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'n') * DevArraysPtr->S_n[local] - DevArraysPtr->roS_n[local];
			F[3] = ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'g') * DevArraysPtr->S_g[local] - DevArraysPtr->roS_g[local];
			F[4] = DevArraysPtr->m[local] * (DevArraysPtr->S_w[local] * (ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'w') * DevArraysPtr->H_w[local] - DevArraysPtr->P_w[local])
				+ DevArraysPtr->S_n[local] * (ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'n') * DevArraysPtr->H_n[local] - DevArraysPtr->P_w[local])
				+ DevArraysPtr->S_g[local] * (ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'g') * DevArraysPtr->H_g[local] - DevArraysPtr->P_w[local]))
				+ (1. - DevArraysPtr->m[local]) * (gpu_def->ro_r * DevArraysPtr->H_r[local] - DevArraysPtr->P_w[local])
				- DevArraysPtr->E_new[local];

			// Матрица частных производных. Строки: dF/dSw, dF/dSn, dF/dSg, dF/dP, dF/dT

			dF[0] = 1.;
			dF[1] = 1.;
			dF[2] = 1.;
			dF[3] = 0.;
			dF[4] = 0.;

			dF[0 + n] = ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'w');
			dF[1 + n] = 0.;
			dF[2 + n] = 0.;
			dF[3 + n] = d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'w', 'P') * DevArraysPtr->S_w[local];
			dF[4 + n] = d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'w', 'T') * DevArraysPtr->S_w[local];

			dF[0 + n * 2] = 0.;
			dF[1 + n * 2] = ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'n');
			dF[2 + n * 2] = 0.;
			dF[3 + n * 2] = d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'n', 'P') * DevArraysPtr->S_n[local];
			dF[4 + n * 2] = d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'n', 'T') * DevArraysPtr->S_n[local];

			dF[0 + n * 3] = 0.;
			dF[1 + n * 3] = 0.;
			dF[2 + n * 3] = ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'g');
			dF[3 + n * 3] = d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'g', 'P') * DevArraysPtr->S_g[local];
			dF[4 + n * 3] = d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'g', 'T') * DevArraysPtr->S_g[local];

			dF[0 + n * 4] = DevArraysPtr->m[local] * (ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'w') * DevArraysPtr->H_w[local] - DevArraysPtr->P_w[local]);
			dF[1 + n * 4] = DevArraysPtr->m[local] * (ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'n') * DevArraysPtr->H_n[local] - DevArraysPtr->P_w[local]);
			dF[2 + n * 4] = DevArraysPtr->m[local] * (ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'g') * DevArraysPtr->H_g[local] - DevArraysPtr->P_w[local]);
			dF[3 + n * 4] = DevArraysPtr->m[local] * (DevArraysPtr->S_w[local] * (d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'w', 'P') * DevArraysPtr->H_w[local] - 1.)
				+ DevArraysPtr->S_n[local] * (d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'n', 'P') * DevArraysPtr->H_n[local] - 1.)
				+ DevArraysPtr->S_g[local] * (d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'g', 'P') * DevArraysPtr->H_g[local] - 1.))
				+ (1. - DevArraysPtr->m[local]) * (-1);
			dF[4 + n * 4] = DevArraysPtr->m[local] * (DevArraysPtr->S_w[local] * (d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'w', 'T') * DevArraysPtr->H_w[local]
				+ ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'w') * c_w(DevArraysPtr->T[local]))
				+ DevArraysPtr->S_n[local] * (d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'n', 'T') * DevArraysPtr->H_n[local]
				+ ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'n') * c_n(DevArraysPtr->T[local]))
				+ DevArraysPtr->S_g[local] * (d_ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'g', 'T') * DevArraysPtr->H_g[local]
				+ ro(DevArraysPtr->P_w[local], DevArraysPtr->T[local], 'g') * c_g(DevArraysPtr->T[local])))
				+ (1. - DevArraysPtr->m[local]) * gpu_def->ro_r * c_r(DevArraysPtr->T[local]);

			device_reverse_matrix(dF, n);
			device_mult_matrix_vector(correction, dF, F, n);

			DevArraysPtr->S_w[local] = DevArraysPtr->S_w[local] - correction[0];
			DevArraysPtr->S_n[local] = DevArraysPtr->S_n[local] - correction[1];
			DevArraysPtr->S_g[local] = DevArraysPtr->S_g[local] - correction[2];
			DevArraysPtr->P_w[local] = DevArraysPtr->P_w[local] - correction[3];
			DevArraysPtr->T[local] = DevArraysPtr->T[local] - correction[4];
			device_assign_H(local);
		}

		// Обновление значения суммарной энергии, т.к. оно больше не понадобится
		// !!! Лучше вынести в отдельную функцию (просто обменять указатели).
		DevArraysPtr->E[local] = DevArraysPtr->E_new[local];

		device_test_S(DevArraysPtr->S_w[local], __FILE__, __LINE__);
		device_test_S(DevArraysPtr->S_n[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr->P_w[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr->T[local], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr->E[local], __FILE__, __LINE__);
	}
}
