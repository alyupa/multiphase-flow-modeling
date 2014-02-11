#include "defines.h"

// !!! Нужно потом будет вынести в структуру констант
// Базовая температура
double T_0 = 273.; // К
// Плотность породы
double ro_r = 2000.; // кг/м^3
// Теплопроводность
double lambda0_w = 0.553; // Вт/(м*К)
double lambda0_n = 0.14;
double lambda0_g = 0.0237;
double lambda0_r = 1.;
double lambdaA_w = 3E-3; // 1/K
double lambdaA_n = 1E-3;
double lambdaA_g = 0.82;
// Теплоемкость
double c0_w = 4.194E3; // Дж/(кг*К)
double c0_n = 1.7E3;
double c0_g = 1E3;
double c0_r = 0.8E3;
double C_w = 1.15;
double C_w2 = 0.015;
double C_n = 3.4;
double C_g = 0.119;
double C_r = 0.75;
// 1/K !!! E-4 Коэффициент теплового расширения (для плотности)
double alfa_w = 1.32E-7; 
double alfa_n = 9.2E-7;


// Коэффициенты удельных теплоемкостей при постоянном давлении  для water, napl, gas and rock в Вт/(м*К)
double c_w (double T, const consts &def)
{
	return c0_w - C_w * (T - T_0) + C_w2 * (T - T_0) * (T - T_0);
}

double c_n (double T, const consts &def)
{
	return c0_n + C_n * (T - T_0);
}

double c_g (double T, const consts &def)
{
	return c0_g + C_g * (T - T_0);
}

double c_r (double T, const consts &def)
{
	return c0_r + C_r * (T - T_0);
}

// Коэффициенты теплопроводности для water, napl, gas and rock
double lambda_w (double T, const consts &def)
{
	return lambda0_w * (1 - lambdaA_w * (T - T_0));
}

double lambda_n (double T, const consts &def)
{
	return lambda0_n * (1 - lambdaA_n * (T - T_0));
}

double lambda_g (double T, const consts &def)
{
	return lambda0_g * pow((T / T_0), lambdaA_g);
}

double lambda_r (double T, const consts &def)
{
	return lambda0_r;
}

// Эффективный коэффициент теплопроводности в точке (будет использоваться при расчете теплового потока)
double assign_lambda_eff (const ptr_Arrays &HostArraysPtr, int local, const consts &def)
{
	return HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * lambda_w (HostArraysPtr.T[local], def)
		+ HostArraysPtr.S_n[local] * lambda_n (HostArraysPtr.T[local], def)
		+ HostArraysPtr.S_g[local] * lambda_g (HostArraysPtr.T[local], def)) 
		+ (1. - HostArraysPtr.m[local]) * lambda_r (HostArraysPtr.T[local], def);
}

// Расчет энтальпии по температуре и теплоемкости
// !!! Переписать, задав точность для метода Симпсона и передавая указатель на функцию, чтобы не копировать одно и то же
double assign_H_w (double T, const consts &def)
{
	/* Возможно, что Симпсон понадобиться позже, а пока все равно нужно знать явную зависимость энергии от температуры
	double integral = 0, sum = 0, h_temp;
	int N_temp = 50;

	h_temp = (T - T_0) / N_temp;
	
	integral += (def.P_atm / def.ro0_w);
	integral += с_w(T_0, def);
	integral += с_w(T, def);

	for(int i = 2; i < N_temp; i+=2)
		sum += с_w(T_0 + i * h_temp, def);

	sum *= 2;
	integral += sum;
	sum = 0;

	for(int i = 1; i < N_temp; i+=2)
		sum += с_w(T_0 + i * h_temp, def);

	sum *= 4;
	integral += sum;

	h_temp /= 3;
	integral *= h_temp;

	return integral;
	*/
	return (def.P_atm / def.ro0_w) + (T - T_0) * (c0_w - (T - T_0) * (C_w / 2 + C_w2 * (T - T_0) / 3));
}

double assign_H_n (double T, const consts &def)
{
	return (def.P_atm / def.ro0_n) + (T - T_0) * (c0_n + C_n * (T - T_0) / 2);
}

double assign_H_g (double T, const consts &def)
{
	return (def.P_atm / def.ro0_g) + (T - T_0) * (c0_g + C_g * (T - T_0) / 2);
}

double assign_H_r (double T, const consts &def)
{
	return (def.P_atm / ro_r) + (T - T_0) * (c0_r + C_r * (T - T_0) / 2);
}

void assign_H (const ptr_Arrays &HostArraysPtr, int local, const consts &def)
{
	HostArraysPtr.H_w[local] = assign_H_w (HostArraysPtr.T[local], def);
	HostArraysPtr.H_n[local] = assign_H_n (HostArraysPtr.T[local], def);
	HostArraysPtr.H_g[local] = assign_H_g (HostArraysPtr.T[local], def);
	HostArraysPtr.H_r[local] = assign_H_r (HostArraysPtr.T[local], def);

	test_nan(HostArraysPtr.H_w[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.H_n[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.H_g[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.H_r[local], __FILE__, __LINE__);
}

// Возвращает значение плотности в точке фазы phase как функции от P, T
double ro(double P, double T, char phase, const consts &def)
{
	double result_ro;
	switch (phase)
	{
	case 'w':
		result_ro = def.ro0_w * (1. + (def.beta_w) * (P - def.P_atm) - alfa_w * (T - T_0));
		break;
	case 'n':
		result_ro = def.ro0_n * (1. + (def.beta_n) * (P - def.P_atm) - alfa_n * (T - T_0));
		break;
	case 'g':
		result_ro = def.ro0_g * (P / def.P_atm) * (T_0 / T);
		break;
	default:
		printf ("Wrong phase in function ro!\n");
		break;
	}

//	test_ro(result_ro, __FILE__, __LINE__);
//	test_ro(result_ro, __FILE__, __LINE__);
//	test_ro(result_ro, __FILE__, __LINE__);

	return result_ro;
}

// Возвращает значение частной производной плотности в точке фазы phase по переменной var
double d_ro(double P, double T, char phase, char var, const consts &def)
{
	double result_d_ro = 0;
	switch (phase)
	{
	case 'w':
		if (var == 'P')
		{
			result_d_ro = def.ro0_w * (def.beta_w);
		} 
		else if (var == 'T')
		{
			result_d_ro = (-1) * def.ro0_w * alfa_w;
		}
		else 
		{
			printf("Wrong var in function d_ro!\n");
		}
		break;
	case 'n':
		if (var == 'P')
		{
			result_d_ro = def.ro0_n * (def.beta_n);
		} 
		else if (var == 'T')
		{
			result_d_ro = (-1) * def.ro0_n * alfa_n;
		}
		else 
		{
			printf("Wrong var in function d_ro!\n");
		}
		break;
	case 'g':
		if (var == 'P')
		{
			result_d_ro = (def.ro0_g / def.P_atm) * (T_0 / T);
		} 
		else if (var == 'T')
		{
			result_d_ro = def.ro0_g * (P / def.P_atm) * (-1) * (T_0 / T) / T;
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

	//	test_ro(result_d_ro, __FILE__, __LINE__);
	//	test_ro(result_d_ro, __FILE__, __LINE__);
	//	test_ro(result_d_ro, __FILE__, __LINE__);

	return result_d_ro;
}

// Коэффициенты вязкости для water, napl, gas and rock

// Расчет теплового потока в точке
double assign_T_flow (const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{	
	if (INTERNAL_POINT)
	{
		double T_flow = 0;
		int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if ((def.locNx) > 2)
		{
			T_flow += (assign_lambda_eff(HostArraysPtr, local + 1, def) * HostArraysPtr.T[local + 1]
			- 2 * assign_lambda_eff(HostArraysPtr, local, def) * HostArraysPtr.T[local]
			+ assign_lambda_eff(HostArraysPtr, local - 1, def) * HostArraysPtr.T[local - 1]) / ((def.hx) * (def.hx));
		}
		if ((def.locNy) > 2)
		{
			T_flow += (assign_lambda_eff(HostArraysPtr, local + def.locNx, def) * HostArraysPtr.T[local + def.locNx]
			- 2 * assign_lambda_eff(HostArraysPtr, local, def) * HostArraysPtr.T[local]
			+ assign_lambda_eff(HostArraysPtr, local - def.locNx, def) * HostArraysPtr.T[local - def.locNx]) / ((def.hy) * (def.hy));
		}

		if ((def.locNz) > 2)
		{
			T_flow = (assign_lambda_eff(HostArraysPtr, local + (def.locNx) * (def.locNy), def) * HostArraysPtr.T[local + (def.locNx) * (def.locNy)]
			- 2 * assign_lambda_eff(HostArraysPtr, local, def) * HostArraysPtr.T[local]
			+ assign_lambda_eff(HostArraysPtr, local - (def.locNx) * (def.locNy), def) * HostArraysPtr.T[local - (def.locNx) * (def.locNy)]) / ((def.hz) * (def.hz));
		}

		test_u(T_flow, __FILE__, __LINE__);
		return T_flow;
	}
	else
		return 0;
}

// Расчет потока энергии в точке
double assign_E_flow (const ptr_Arrays &HostArraysPtr, int i, int j, int k,  const consts &def)
{
	if (INTERNAL_POINT)
	{
		double E_flow = 0;
		int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if ((def.locNx) > 2)
		{
			E_flow += (HostArraysPtr.ro_w[local + 1] * HostArraysPtr.H_w[local + 1] * HostArraysPtr.ux_w[local + 1]
				- HostArraysPtr.ro_w[local - 1] * HostArraysPtr.H_w[local - 1] * HostArraysPtr.ux_w[local - 1]
				+ HostArraysPtr.ro_n[local + 1] * HostArraysPtr.H_n[local + 1] * HostArraysPtr.ux_n[local + 1]
				- HostArraysPtr.ro_n[local - 1] * HostArraysPtr.H_n[local - 1] * HostArraysPtr.ux_n[local - 1]
				+ HostArraysPtr.ro_g[local + 1] * HostArraysPtr.H_g[local + 1] * HostArraysPtr.ux_g[local + 1]
				- HostArraysPtr.ro_g[local - 1] * HostArraysPtr.H_g[local - 1] * HostArraysPtr.ux_g[local - 1]
				) / (2. * (def.hx));
		}
		if ((def.locNy) > 2)
		{
			E_flow += (HostArraysPtr.ro_w[local + def.locNx] * HostArraysPtr.H_w[local + def.locNx] * HostArraysPtr.uy_w[local + def.locNx]
			- HostArraysPtr.ro_w[local - def.locNx] * HostArraysPtr.H_w[local - def.locNx] * HostArraysPtr.uy_w[local - def.locNx]
			+ HostArraysPtr.ro_n[local + def.locNx] * HostArraysPtr.H_n[local + def.locNx] * HostArraysPtr.uy_n[local + def.locNx]
			- HostArraysPtr.ro_n[local - def.locNx] * HostArraysPtr.H_n[local - def.locNx] * HostArraysPtr.uy_n[local - def.locNx]
			+ HostArraysPtr.ro_g[local + def.locNx] * HostArraysPtr.H_g[local + def.locNx] * HostArraysPtr.uy_g[local + def.locNx]
			- HostArraysPtr.ro_g[local - def.locNx] * HostArraysPtr.H_g[local - def.locNx] * HostArraysPtr.uy_g[local - def.locNx]
			)/ (2. * (def.hy));	
		}

		if ((def.locNz) > 2)
		{
			E_flow += (HostArraysPtr.ro_w[local + (def.locNx) * (def.locNy)] * HostArraysPtr.H_w[local + (def.locNx) * (def.locNy)] * HostArraysPtr.uy_w[local + (def.locNx) * (def.locNy)]
			- HostArraysPtr.ro_w[local - (def.locNx) * (def.locNy)] * HostArraysPtr.H_w[local - (def.locNx) * (def.locNy)] * HostArraysPtr.uy_w[local - (def.locNx) * (def.locNy)]
			+ HostArraysPtr.ro_n[local + (def.locNx) * (def.locNy)] * HostArraysPtr.H_n[local + (def.locNx) * (def.locNy)] * HostArraysPtr.uy_n[local + (def.locNx) * (def.locNy)]
			- HostArraysPtr.ro_n[local - (def.locNx) * (def.locNy)] * HostArraysPtr.H_n[local - (def.locNx) * (def.locNy)] * HostArraysPtr.uy_n[local - (def.locNx) * (def.locNy)]
			+ HostArraysPtr.ro_g[local + (def.locNx) * (def.locNy)] * HostArraysPtr.H_g[local + (def.locNx) * (def.locNy)] * HostArraysPtr.uy_g[local + (def.locNx) * (def.locNy)]
			- HostArraysPtr.ro_g[local - (def.locNx) * (def.locNy)] * HostArraysPtr.H_g[local - (def.locNx) * (def.locNy)] * HostArraysPtr.uy_g[local - (def.locNx) * (def.locNy)]
			)/ (2. * (def.hz));	
		}

		test_u(E_flow, __FILE__, __LINE__);
		return E_flow;
	}
	else
		return 0;
}

// Расчет внутренней энергии всей системы в точке
void assign_E_current (const ptr_Arrays &HostArraysPtr, int local, const consts &def)
{
	HostArraysPtr.E[local] = (HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (HostArraysPtr.ro_w[local] * HostArraysPtr.H_w[local] - HostArraysPtr.P_w[local])
		+ HostArraysPtr.S_n[local] * (HostArraysPtr.ro_n[local] * HostArraysPtr.H_n[local] - HostArraysPtr.P_n[local])
		+ HostArraysPtr.S_g[local] * (HostArraysPtr.ro_g[local] * HostArraysPtr.H_g[local] - HostArraysPtr.P_g[local])) 
		+ (1. - HostArraysPtr.m[local]) * (ro_r * HostArraysPtr.H_r[local] - HostArraysPtr.P_w[local]));

	test_nan(HostArraysPtr.E[local], __FILE__, __LINE__);
}

// Расчет внутренней энергии всей системы в точке на следующем шаге по времени
void assign_E_new (const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	if (INTERNAL_POINT)
	{
		double Q_hw = 0, Q_hr = 0; // источниковые члены

		int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		HostArraysPtr.E_new[local] = HostArraysPtr.E[local] + (def.dt) * (assign_T_flow(HostArraysPtr, i, j, k, def) + Q_hw + Q_hr - assign_E_flow(HostArraysPtr, i, j, k, def));

		test_nan(HostArraysPtr.E_new[local], __FILE__, __LINE__);
	}
}

// Задание граничных условий на температуру
void Border_T(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	if (BOUNDARY_POINT)
	{
		int local1 = set_boundary_basic_coordinate(i, j, k, def);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if (j == 0)
		{
			HostArraysPtr.T[local] = 400;
		}
		else if(j == (def.locNy) - 1)
		{
			HostArraysPtr.T[local] = 273;
		}
		else
		{
			// Будем считать границы области не теплопроводящими
			HostArraysPtr.T[local] = HostArraysPtr.T[local1];
		}

		test_positive(HostArraysPtr.T[local], __FILE__, __LINE__);
	}
}

// Расчет методом Ньютона значений переменных на новом шаге по времени, когда учитываем изменение энергии (случай 4х переменных)
// !!! Пока "выбросим" капиллярные давления
void Newton(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	if (INTERNAL_POINT)
	{
		int n = 5; // Размерность системы
		double *F; // Вектор значений функций (из системы уравнений)
		double *correction; // Вектор поправок к функциям
		double *dF; // Матрица Якоби (в виде одномерного массива)

		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
		
		F = new double [n];
		correction = new double [n];
		dF = new double [n * n];

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			F[0] = HostArraysPtr.S_g[local] + HostArraysPtr.S_w[local] + HostArraysPtr.S_n[local] - 1.;
			F[1] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', def) * HostArraysPtr.S_w[local] - HostArraysPtr.roS_w[local];
			F[2] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', def) * HostArraysPtr.S_n[local] - HostArraysPtr.roS_n[local];
			F[3] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', def) * HostArraysPtr.S_g[local] - HostArraysPtr.roS_g[local];
			F[4] = HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', def) * HostArraysPtr.H_w[local] - HostArraysPtr.P_w[local])
				+ HostArraysPtr.S_n[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', def) * HostArraysPtr.H_n[local] - HostArraysPtr.P_w[local])
				+ HostArraysPtr.S_g[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', def) * HostArraysPtr.H_g[local] - HostArraysPtr.P_w[local])) 
				+ (1. - HostArraysPtr.m[local]) * (ro_r * HostArraysPtr.H_r[local] - HostArraysPtr.P_w[local]) 
				- HostArraysPtr.E_new[local];

			// Матрица частных производных. Строки: dF/dSw, dF/dSn, dF/dSg, dF/dP, dF/dT

			dF[0] = 1.;
			dF[1] = 1.;
			dF[2] = 1.;
			dF[3] = 0.;
			dF[4] = 0.;

			dF[0 + n] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', def);
			dF[1 + n] = 0.;
			dF[2 + n] = 0.;
			dF[3 + n] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', 'P', def) * HostArraysPtr.S_w[local];
			dF[4 + n] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', 'T', def) * HostArraysPtr.S_w[local];

			dF[0 + n * 2] = 0.;
			dF[1 + n * 2] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', def);
			dF[2 + n * 2] = 0.;
			dF[3 + n * 2] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', 'P', def) * HostArraysPtr.S_n[local];
			dF[4 + n * 2] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', 'T', def) * HostArraysPtr.S_n[local];

			dF[0 + n * 3] = 0.;
			dF[1 + n * 3] = 0.;
			dF[2 + n * 3] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', def);
			dF[3 + n * 3] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', 'P', def) * HostArraysPtr.S_g[local];
			dF[4 + n * 3] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', 'T', def) * HostArraysPtr.S_g[local];

			dF[0 + n * 4] = HostArraysPtr.m[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', def) * HostArraysPtr.H_w[local] - HostArraysPtr.P_w[local]);
			dF[1 + n * 4] = HostArraysPtr.m[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', def) * HostArraysPtr.H_n[local] - HostArraysPtr.P_w[local]);
			dF[2 + n * 4] = HostArraysPtr.m[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', def) * HostArraysPtr.H_g[local] - HostArraysPtr.P_w[local]);
			dF[3 + n * 4] = HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', 'P', def) * HostArraysPtr.H_w[local] - 1.) 
				+ HostArraysPtr.S_n[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', 'P', def) * HostArraysPtr.H_n[local] - 1.)  
				+ HostArraysPtr.S_g[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', 'P', def) * HostArraysPtr.H_g[local] - 1.))
				+ (1. - HostArraysPtr.m[local]) * (-1);
			dF[4 + n * 4] = HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', 'T', def) * HostArraysPtr.H_w[local] 
				+ ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', def) * c_w(HostArraysPtr.T[local], def)) 
				+ HostArraysPtr.S_n[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', 'T', def) * HostArraysPtr.H_n[local]
				+ ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', def) * c_n(HostArraysPtr.T[local], def))
				+ HostArraysPtr.S_g[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', 'T', def) * HostArraysPtr.H_g[local]
				+ ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', def) * c_g(HostArraysPtr.T[local], def)))
				+ (1. - HostArraysPtr.m[local]) * ro_r * c_r(HostArraysPtr.T[local], def);

			reverse_matrix(dF, n);
			mult_matrix_vector(correction, dF, F, n);

			HostArraysPtr.S_w[local] = HostArraysPtr.S_w[local] - correction[0];
			HostArraysPtr.S_n[local] = HostArraysPtr.S_n[local] - correction[1];
			HostArraysPtr.S_g[local] = HostArraysPtr.S_g[local] - correction[2];
			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local] - correction[3];
			HostArraysPtr.T[local] = HostArraysPtr.T[local] - correction[4];
			assign_H(HostArraysPtr, local, def);
		}

		// Обновление значения суммарной энергии, т.к. оно больше не понадобится
		// !!! Лучше вынести в отдельную функцию (просто обменять указатели).
		HostArraysPtr.E[local] = HostArraysPtr.E_new[local];

		delete[] F;
		delete[] correction;
		delete[] dF;

		test_S(HostArraysPtr.S_w[local], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.T[local], __FILE__, __LINE__);
		test_nan(HostArraysPtr.E[local], __FILE__, __LINE__);
	}
}