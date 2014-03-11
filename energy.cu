#include "gpu.h"

// !!! Нужно потом будет вынести в структуру констант
// Базовая температура
__constant__ double T_0 = 273.; // К
// Плотность породы
__constant__ double ro_r = 2000.; // кг/м^3
// Теплопроводность
__constant__ double lambda0_w = 0.553; // Вт/(м*К)
__constant__ double lambda0_n = 0.14;
__constant__ double lambda0_g = 0.0237;
__constant__ double lambda0_r = 1.;
__constant__ double lambdaA_w = 3E-3; // 1/K
__constant__ double lambdaA_n = 1E-3;
__constant__ double lambdaA_g = 0.82;
// Теплоемкость
__constant__ double c0_w = 4.194E3; // Дж/(кг*К)
__constant__ double c0_n = 1.7E3;
__constant__ double c0_g = 1E3;
__constant__ double c0_r = 0.8E3;
__constant__ double C_w = 1.15;
__constant__ double C_w2 = 0.015;
__constant__ double C_n = 3.4;
__constant__ double C_g = 0.119;
__constant__ double C_r = 0.75;
// 1/K !!! E-4 Коэффициент теплового расширения (для плотности)
__constant__ double alfa_w = 1.32E-7; 
__constant__ double alfa_n = 9.2E-7;


// Коэффициенты удельных теплоемкостей при постоянном давлении  для water, napl, gas and rock в Вт/(м*К)
__device__ double c_w (double T)
{
	return c0_w - C_w * (T - T_0) + C_w2 * (T - T_0) * (T - T_0);
}

__device__ double c_n (double T)
{
	return c0_n + C_n * (T - T_0);
}

__device__ double c_g (double T)
{
	return c0_g + C_g * (T - T_0);
}

__device__ double c_r (double T)
{
	return c0_r + C_r * (T - T_0);
}

// Коэффициенты теплопроводности для water, napl, gas and rock
__device__ double lambda_w (double T)
{
	return lambda0_w * (1 - lambdaA_w * (T - T_0));
}

__device__ double lambda_n (double T)
{
	return lambda0_n * (1 - lambdaA_n * (T - T_0));
}

__device__ double lambda_g (double T)
{
	return lambda0_g * pow((T / T_0), lambdaA_g);
}

__device__ double lambda_r (double T)
{
	return lambda0_r;
}

// Эффективный коэффициент теплопроводности в точке (будет использоваться при расчете теплового потока)
__device__ double assign_lambda_eff (ptr_Arrays DevArraysPtr, int local)
{
	return DevArraysPtr.m[local] * (DevArraysPtr.S_w[local] * lambda_w (DevArraysPtr.T[local])
		+ DevArraysPtr.S_n[local] * lambda_n (DevArraysPtr.T[local])
		+ DevArraysPtr.S_g[local] * lambda_g (DevArraysPtr.T[local])) 
		+ (1. - DevArraysPtr.m[local]) * lambda_r (DevArraysPtr.T[local]);
}

// Расчет энтальпии по температуре и теплоемкости
// !!! Переписать, задав точность для метода Симпсона и передавая указатель на функцию, чтобы не копировать одно и то же
__device__ double assign_H_w (double T)
{
	/* Возможно, что Симпсон понадобиться позже, а пока все равно нужно знать явную зависимость энергии от температуры
	double integral = 0, sum = 0, h_temp;
	int N_temp = 50;

	h_temp = (T - T_0) / N_temp;
	
	integral += (gpu_def->P_atm / gpu_def->ro0_w);
	integral += с_w(T_0);
	integral += с_w(T);

	for(int i = 2; i < N_temp; i+=2)
		sum += с_w(T_0 + i * h_temp);

	sum *= 2;
	integral += sum;
	sum = 0;

	for(int i = 1; i < N_temp; i+=2)
		sum += с_w(T_0 + i * h_temp);

	sum *= 4;
	integral += sum;

	h_temp /= 3;
	integral *= h_temp;

	return integral;
	*/
	return (gpu_def->P_atm / gpu_def->ro0_w) + (T - T_0) * (c0_w - (T - T_0) * (C_w / 2 + C_w2 * (T - T_0) / 3));
}

__device__ double assign_H_n (double T)
{
	return (gpu_def->P_atm / gpu_def->ro0_n) + (T - T_0) * (c0_n + C_n * (T - T_0) / 2);
}

__device__ double assign_H_g (double T)
{
	return (gpu_def->P_atm / gpu_def->ro0_g) + (T - T_0) * (c0_g + C_g * (T - T_0) / 2);
}

__device__ double assign_H_r (double T)
{
	return (gpu_def->P_atm / ro_r) + (T - T_0) * (c0_r + C_r * (T - T_0) / 2);
}

__device__ void device_assign_H (ptr_Arrays DevArraysPtr, int local)
{
	DevArraysPtr.H_w[local] = assign_H_w (DevArraysPtr.T[local]);
	DevArraysPtr.H_n[local] = assign_H_n (DevArraysPtr.T[local]);
	DevArraysPtr.H_g[local] = assign_H_g (DevArraysPtr.T[local]);
	DevArraysPtr.H_r[local] = assign_H_r (DevArraysPtr.T[local]);

	device_test_nan(DevArraysPtr.H_w[local], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.H_n[local], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.H_g[local], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.H_r[local], __FILE__, __LINE__);
}

// Возвращает значение плотности в точке фазы phase как функции от P, T
__device__ double ro(double P, double T, char phase)
{
	double result_ro;
	switch (phase)
	{
	case 'w':
		result_ro = gpu_def->ro0_w * (1. + (gpu_def->beta_w) * (P - gpu_def->P_atm) - alfa_w * (T - T_0));
		break;
	case 'n':
		result_ro = gpu_def->ro0_n * (1. + (gpu_def->beta_n) * (P - gpu_def->P_atm) - alfa_n * (T - T_0));
		break;
	case 'g':
		result_ro = gpu_def->ro0_g * (P / gpu_def->P_atm) * (T_0 / T);
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
			result_d_ro = (-1) * gpu_def->ro0_w * alfa_w;
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
			result_d_ro = (-1) * gpu_def->ro0_n * alfa_n;
		}
		else 
		{
			printf("Wrong var in function d_ro!\n");
		}
		break;
	case 'g':
		if (var == 'P')
		{
			result_d_ro = (gpu_def->ro0_g / gpu_def->P_atm) * (T_0 / T);
		} 
		else if (var == 'T')
		{
			result_d_ro = gpu_def->ro0_g * (P / gpu_def->P_atm) * (-1) * (T_0 / T) / T;
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
__device__ double assign_T_flow (ptr_Arrays DevArraysPtr, int i, int j, int k)
{
	if (GPU_INTERNAL_POINT)
	{
		double T_flow = 0;
		int local=i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		if ((gpu_def->locNx) > 2)
		{
			T_flow += (assign_lambda_eff(DevArraysPtr, local + 1) * DevArraysPtr.T[local + 1]
			- 2 * assign_lambda_eff(DevArraysPtr, local) * DevArraysPtr.T[local]
			+ assign_lambda_eff(DevArraysPtr, local - 1) * DevArraysPtr.T[local - 1]) / ((gpu_def->hx) * (gpu_def->hx));
		}
		if ((gpu_def->locNy) > 2)
		{
			T_flow += (assign_lambda_eff(DevArraysPtr, local + gpu_def->locNx) * DevArraysPtr.T[local + gpu_def->locNx]
			- 2 * assign_lambda_eff(DevArraysPtr, local) * DevArraysPtr.T[local]
			+ assign_lambda_eff(DevArraysPtr, local - gpu_def->locNx) * DevArraysPtr.T[local - gpu_def->locNx]) / ((gpu_def->hy) * (gpu_def->hy));
		}

		if ((gpu_def->locNz) > 2)
		{
			T_flow = (assign_lambda_eff(DevArraysPtr, local + (gpu_def->locNx) * (gpu_def->locNy)) * DevArraysPtr.T[local + (gpu_def->locNx) * (gpu_def->locNy)]
			- 2 * assign_lambda_eff(DevArraysPtr, local) * DevArraysPtr.T[local]
			+ assign_lambda_eff(DevArraysPtr, local - (gpu_def->locNx) * (gpu_def->locNy)) * DevArraysPtr.T[local - (gpu_def->locNx) * (gpu_def->locNy)]) / ((gpu_def->hz) * (gpu_def->hz));
		}

		device_test_u(T_flow, __FILE__, __LINE__);
		return T_flow;
	}
	else
		return 0;
}

// Расчет потока энергии в точке
__device__ double assign_E_flow (ptr_Arrays DevArraysPtr, int i, int j, int k)
{
	
	if (GPU_INTERNAL_POINT)
	{
		double E_flow = 0;
		int local=i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		if ((gpu_def->locNx) > 2)
		{
			E_flow += (DevArraysPtr.ro_w[local + 1] * DevArraysPtr.H_w[local + 1] * DevArraysPtr.ux_w[local + 1]
				- DevArraysPtr.ro_w[local - 1] * DevArraysPtr.H_w[local - 1] * DevArraysPtr.ux_w[local - 1]
				+ DevArraysPtr.ro_n[local + 1] * DevArraysPtr.H_n[local + 1] * DevArraysPtr.ux_n[local + 1]
				- DevArraysPtr.ro_n[local - 1] * DevArraysPtr.H_n[local - 1] * DevArraysPtr.ux_n[local - 1]
				+ DevArraysPtr.ro_g[local + 1] * DevArraysPtr.H_g[local + 1] * DevArraysPtr.ux_g[local + 1]
				- DevArraysPtr.ro_g[local - 1] * DevArraysPtr.H_g[local - 1] * DevArraysPtr.ux_g[local - 1]
				) / (2. * (gpu_def->hx));
		}
		if ((gpu_def->locNy) > 2)
		{
			E_flow += (DevArraysPtr.ro_w[local + gpu_def->locNx] * DevArraysPtr.H_w[local + gpu_def->locNx] * DevArraysPtr.uy_w[local + gpu_def->locNx]
			- DevArraysPtr.ro_w[local - gpu_def->locNx] * DevArraysPtr.H_w[local - gpu_def->locNx] * DevArraysPtr.uy_w[local - gpu_def->locNx]
			+ DevArraysPtr.ro_n[local + gpu_def->locNx] * DevArraysPtr.H_n[local + gpu_def->locNx] * DevArraysPtr.uy_n[local + gpu_def->locNx]
			- DevArraysPtr.ro_n[local - gpu_def->locNx] * DevArraysPtr.H_n[local - gpu_def->locNx] * DevArraysPtr.uy_n[local - gpu_def->locNx]
			+ DevArraysPtr.ro_g[local + gpu_def->locNx] * DevArraysPtr.H_g[local + gpu_def->locNx] * DevArraysPtr.uy_g[local + gpu_def->locNx]
			- DevArraysPtr.ro_g[local - gpu_def->locNx] * DevArraysPtr.H_g[local - gpu_def->locNx] * DevArraysPtr.uy_g[local - gpu_def->locNx]
			)/ (2. * (gpu_def->hy));	
		}

		if ((gpu_def->locNz) > 2)
		{
			E_flow += (DevArraysPtr.ro_w[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.H_w[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_w[local + (gpu_def->locNx) * (gpu_def->locNy)]
			- DevArraysPtr.ro_w[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.H_w[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_w[local - (gpu_def->locNx) * (gpu_def->locNy)]
			+ DevArraysPtr.ro_n[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.H_n[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_n[local + (gpu_def->locNx) * (gpu_def->locNy)]
			- DevArraysPtr.ro_n[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.H_n[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_n[local - (gpu_def->locNx) * (gpu_def->locNy)]
			+ DevArraysPtr.ro_g[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.H_g[local + (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_g[local + (gpu_def->locNx) * (gpu_def->locNy)]
			- DevArraysPtr.ro_g[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.H_g[local - (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_g[local - (gpu_def->locNx) * (gpu_def->locNy)]
			)/ (2. * (gpu_def->hz));	
		}

		device_test_u(E_flow, __FILE__, __LINE__);
		return E_flow;
	}
	else
		return 0;
}

// Расчет внутренней энергии всей системы в точке
__device__ void device_assign_E_current (ptr_Arrays DevArraysPtr, int local)
{
	DevArraysPtr.E[local] = (DevArraysPtr.m[local] * (DevArraysPtr.S_w[local] * (DevArraysPtr.ro_w[local] * DevArraysPtr.H_w[local] - DevArraysPtr.P_w[local])
		+ DevArraysPtr.S_n[local] * (DevArraysPtr.ro_n[local] * DevArraysPtr.H_n[local] - DevArraysPtr.P_n[local])
		+ DevArraysPtr.S_g[local] * (DevArraysPtr.ro_g[local] * DevArraysPtr.H_g[local] - DevArraysPtr.P_g[local])) 
		+ (1. - DevArraysPtr.m[local]) * (ro_r * DevArraysPtr.H_r[local] - DevArraysPtr.P_w[local]));

	device_test_nan(DevArraysPtr.E[local], __FILE__, __LINE__);
}

// Расчет внутренней энергии всей системы в точке на следующем шаге по времени
__global__ void assign_E_new_kernel (ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (GPU_INTERNAL_POINT)
	{
		double Q_hw = 0, Q_hr = 0; // источниковые члены

		int local=i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		DevArraysPtr.E_new[local] = DevArraysPtr.E[local] + (gpu_def->dt) * (assign_T_flow(DevArraysPtr, i, j, k) + Q_hw + Q_hr - assign_E_flow(DevArraysPtr, i, j, k));

		device_test_nan(DevArraysPtr.E_new[local], __FILE__, __LINE__);
	}
}

// Задание граничных условий на температуру
__global__ void Border_T_kernel(ptr_Arrays DevArraysPtr)
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
			DevArraysPtr.T[local] = 400;
		}
		else if(j == (gpu_def->locNy) - 1)
		{
			DevArraysPtr.T[local] = 273;
		}
		else
		{
			// Будем считать границы области не теплопроводящими
			DevArraysPtr.T[local] = DevArraysPtr.T[local1];
		}

		device_test_positive(DevArraysPtr.T[local], __FILE__, __LINE__);
	}
}

// Расчет методом Ньютона значений переменных на новом шаге по времени, когда учитываем изменение энергии (случай 4х переменных)
// !!! Пока "выбросим" капиллярные давления
__global__ void Newton_method_kernel(ptr_Arrays DevArraysPtr)
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
			F[0] = DevArraysPtr.S_g[local] + DevArraysPtr.S_w[local] + DevArraysPtr.S_n[local] - 1.;
			F[1] = ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'w') * DevArraysPtr.S_w[local] - DevArraysPtr.roS_w[local];
			F[2] = ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'n') * DevArraysPtr.S_n[local] - DevArraysPtr.roS_n[local];
			F[3] = ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'g') * DevArraysPtr.S_g[local] - DevArraysPtr.roS_g[local];
			F[4] = DevArraysPtr.m[local] * (DevArraysPtr.S_w[local] * (ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'w') * DevArraysPtr.H_w[local] - DevArraysPtr.P_w[local])
				+ DevArraysPtr.S_n[local] * (ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'n') * DevArraysPtr.H_n[local] - DevArraysPtr.P_w[local])
				+ DevArraysPtr.S_g[local] * (ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'g') * DevArraysPtr.H_g[local] - DevArraysPtr.P_w[local])) 
				+ (1. - DevArraysPtr.m[local]) * (ro_r * DevArraysPtr.H_r[local] - DevArraysPtr.P_w[local]) 
				- DevArraysPtr.E_new[local];

			// Матрица частных производных. Строки: dF/dSw, dF/dSn, dF/dSg, dF/dP, dF/dT

			dF[0] = 1.;
			dF[1] = 1.;
			dF[2] = 1.;
			dF[3] = 0.;
			dF[4] = 0.;

			dF[0 + n] = ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'w');
			dF[1 + n] = 0.;
			dF[2 + n] = 0.;
			dF[3 + n] = d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'w', 'P') * DevArraysPtr.S_w[local];
			dF[4 + n] = d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'w', 'T') * DevArraysPtr.S_w[local];

			dF[0 + n * 2] = 0.;
			dF[1 + n * 2] = ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'n');
			dF[2 + n * 2] = 0.;
			dF[3 + n * 2] = d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'n', 'P') * DevArraysPtr.S_n[local];
			dF[4 + n * 2] = d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'n', 'T') * DevArraysPtr.S_n[local];

			dF[0 + n * 3] = 0.;
			dF[1 + n * 3] = 0.;
			dF[2 + n * 3] = ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'g');
			dF[3 + n * 3] = d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'g', 'P') * DevArraysPtr.S_g[local];
			dF[4 + n * 3] = d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'g', 'T') * DevArraysPtr.S_g[local];

			dF[0 + n * 4] = DevArraysPtr.m[local] * (ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'w') * DevArraysPtr.H_w[local] - DevArraysPtr.P_w[local]);
			dF[1 + n * 4] = DevArraysPtr.m[local] * (ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'n') * DevArraysPtr.H_n[local] - DevArraysPtr.P_w[local]);
			dF[2 + n * 4] = DevArraysPtr.m[local] * (ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'g') * DevArraysPtr.H_g[local] - DevArraysPtr.P_w[local]);
			dF[3 + n * 4] = DevArraysPtr.m[local] * (DevArraysPtr.S_w[local] * (d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'w', 'P') * DevArraysPtr.H_w[local] - 1.) 
				+ DevArraysPtr.S_n[local] * (d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'n', 'P') * DevArraysPtr.H_n[local] - 1.)  
				+ DevArraysPtr.S_g[local] * (d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'g', 'P') * DevArraysPtr.H_g[local] - 1.))
				+ (1. - DevArraysPtr.m[local]) * (-1);
			dF[4 + n * 4] = DevArraysPtr.m[local] * (DevArraysPtr.S_w[local] * (d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'w', 'T') * DevArraysPtr.H_w[local] 
				+ ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'w') * c_w(DevArraysPtr.T[local])) 
				+ DevArraysPtr.S_n[local] * (d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'n', 'T') * DevArraysPtr.H_n[local]
				+ ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'n') * c_n(DevArraysPtr.T[local]))
				+ DevArraysPtr.S_g[local] * (d_ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'g', 'T') * DevArraysPtr.H_g[local]
				+ ro(DevArraysPtr.P_w[local], DevArraysPtr.T[local], 'g') * c_g(DevArraysPtr.T[local])))
				+ (1. - DevArraysPtr.m[local]) * ro_r * c_r(DevArraysPtr.T[local]);

			device_reverse_matrix(dF, n);
			device_mult_matrix_vector(correction, dF, F, n);

			DevArraysPtr.S_w[local] = DevArraysPtr.S_w[local] - correction[0];
			DevArraysPtr.S_n[local] = DevArraysPtr.S_n[local] - correction[1];
			DevArraysPtr.S_g[local] = DevArraysPtr.S_g[local] - correction[2];
			DevArraysPtr.P_w[local] = DevArraysPtr.P_w[local] - correction[3];
			DevArraysPtr.T[local] = DevArraysPtr.T[local] - correction[4];
			device_assign_H(DevArraysPtr, local);
		}

		// Обновление значения суммарной энергии, т.к. оно больше не понадобится
		// !!! Лучше вынести в отдельную функцию (просто обменять указатели).
		DevArraysPtr.E[local] = DevArraysPtr.E_new[local];

		device_test_S(DevArraysPtr.S_w[local], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_n[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.P_w[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.T[local], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.E[local], __FILE__, __LINE__);
	}
}
