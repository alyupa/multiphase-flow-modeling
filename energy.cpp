#include "defines.h"

extern consts def;

// Коэффициенты удельных теплоемкостей при постоянном давлении  для water, napl, gas and rock в Вт/(м*К)
static inline double c_w (double T)
{
	return def.c0_w - def.C_w * (T - def.T_0) + def.C_w2 * (T - def.T_0) * (T - def.T_0);
}

static inline double c_n (double T)
{
	return def.c0_n + def.C_n * (T - def.T_0);
}

static inline double c_g (double T)
{
	return def.c0_g + def.C_g * (T - def.T_0);
}

static inline double c_r (double T)
{
	return def.c0_r + def.C_r * (T - def.T_0);
}

// Коэффициенты теплопроводности для water, napl, gas and rock
static inline double lambda_w (double T)
{
	return def.lambda0_w * (1 - def.lambdaA_w * (T - def.T_0));
}

static inline double lambda_n (double T)
{
	return def.lambda0_n * (1 - def.lambdaA_n * (T - def.T_0));
}

static inline double lambda_g (double T)
{
	return def.lambda0_g * pow((T / def.T_0), def.lambdaA_g);
}

static inline double lambda_r (double T)
{
	return def.lambda0_r;
}

// Эффективный коэффициент теплопроводности в точке (будет использоваться при расчете теплового потока)
static inline double assign_lambda_eff (int local)
{
	return HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * lambda_w (HostArraysPtr.T[local])
		+ HostArraysPtr.S_n[local] * lambda_n (HostArraysPtr.T[local])
		+ HostArraysPtr.S_g[local] * lambda_g (HostArraysPtr.T[local]))
		+ (1. - HostArraysPtr.m[local]) * lambda_r (HostArraysPtr.T[local]);
}

// Расчет энтальпии по температуре и теплоемкости
// !!! Переписать, задав точность для метода Симпсона и передавая указатель на функцию, чтобы не копировать одно и то же
static inline double assign_H_w (double T)
{
	/* Возможно, что Симпсон понадобиться позже, а пока все равно нужно знать явную зависимость энергии от температуры
	double integral = 0, sum = 0, h_temp;
	int N_temp = 50;

	h_temp = (T - def.T_0) / N_temp;
	
	integral += (def.P_atm / def.ro0_w);
	integral += с_w(def.T_0);
	integral += с_w(T);

	for(int i = 2; i < N_temp; i+=2)
		sum += с_w(def.T_0 + i * h_temp);

	sum *= 2;
	integral += sum;
	sum = 0;

	for(int i = 1; i < N_temp; i+=2)
		sum += с_w(def.T_0 + i * h_temp);

	sum *= 4;
	integral += sum;

	h_temp /= 3;
	integral *= h_temp;

	return integral;
	*/
	return (def.P_atm / def.ro0_w) + (T - def.T_0) * (def.c0_w - (T - def.T_0) * (def.C_w / 2 + def.C_w2 * (T - def.T_0) / 3));
}

static inline double assign_H_n (double T)
{
	return (def.P_atm / def.ro0_n) + (T - def.T_0) * (def.c0_n + def.C_n * (T - def.T_0) / 2);
}

static inline double assign_H_g (double T)
{
	return (def.P_atm / def.ro0_g) + (T - def.T_0) * (def.c0_g + def.C_g * (T - def.T_0) / 2);
}

static inline double assign_H_r (double T)
{
	return (def.P_atm / def.ro_r) + (T - def.T_0) * (def.c0_r + def.C_r * (T - def.T_0) / 2);
}

void assign_H (int local)
{
	HostArraysPtr.H_w[local] = assign_H_w (HostArraysPtr.T[local]);
	HostArraysPtr.H_n[local] = assign_H_n (HostArraysPtr.T[local]);
	HostArraysPtr.H_g[local] = assign_H_g (HostArraysPtr.T[local]);
	HostArraysPtr.H_r[local] = assign_H_r (HostArraysPtr.T[local]);

	test_nan(HostArraysPtr.H_w[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.H_n[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.H_g[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.H_r[local], __FILE__, __LINE__);
}

// Возвращает значение плотности в точке фазы phase как функции от P, T
static inline double ro(double P, double T, char phase)
{
	double result_ro;
	switch (phase)
	{
	case 'w':
		result_ro = def.ro0_w * (1. + (def.beta_w) * (P - def.P_atm) - def.alfa_w * (T - def.T_0));
		break;
	case 'n':
		result_ro = def.ro0_n * (1. + (def.beta_n) * (P - def.P_atm) - def.alfa_n * (T - def.T_0));
		break;
	case 'g':
		result_ro = def.ro0_g * (P / def.P_atm) * (def.T_0 / T);
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
static inline double d_ro(double P, double T, char phase, char var)
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
			result_d_ro = (-1) * def.ro0_w * def.alfa_w;
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
			result_d_ro = (-1) * def.ro0_n * def.alfa_n;
		}
		else 
		{
			printf("Wrong var in function d_ro!\n");
		}
		break;
	case 'g':
		if (var == 'P')
		{
			result_d_ro = (def.ro0_g / def.P_atm) * (def.T_0 / T);
		} 
		else if (var == 'T')
		{
			result_d_ro = def.ro0_g * (P / def.P_atm) * (-1) * (def.T_0 / T) / T;
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
static double assign_T_flow (int i, int j, int k)
{	
	if (INTERNAL_POINT)
	{
		double T_flow = 0;
		int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if ((def.locNx) > 2)
		{
			T_flow += (assign_lambda_eff(local + 1) * HostArraysPtr.T[local + 1]
			- 2 * assign_lambda_eff(local) * HostArraysPtr.T[local]
			+ assign_lambda_eff(local - 1) * HostArraysPtr.T[local - 1]) / ((def.hx) * (def.hx));
		}
		if ((def.locNy) > 2)
		{
			T_flow += (assign_lambda_eff(local + def.locNx) * HostArraysPtr.T[local + def.locNx]
			- 2 * assign_lambda_eff(local) * HostArraysPtr.T[local]
			+ assign_lambda_eff(local - def.locNx) * HostArraysPtr.T[local - def.locNx]) / ((def.hy) * (def.hy));
		}

		if ((def.locNz) > 2)
		{
			T_flow = (assign_lambda_eff(local + (def.locNx) * (def.locNy)) * HostArraysPtr.T[local + (def.locNx) * (def.locNy)]
			- 2 * assign_lambda_eff(local) * HostArraysPtr.T[local]
			+ assign_lambda_eff(local - (def.locNx) * (def.locNy)) * HostArraysPtr.T[local - (def.locNx) * (def.locNy)]) / ((def.hz) * (def.hz));
		}

		test_u(T_flow, __FILE__, __LINE__);
		return T_flow;
	}
	else
		return 0;
}

// Расчет направленной разности
double directed_difference_E (double* P, double* Xi, double* ro, double* H, char axis)
{
	double x1 = 0, x2 = 0;
	switch (axis)
	{
	case 'x':
		{
			x2 = -right_difference (P, 'x');
			x1 = -left_difference (P, 'x');
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (*Xi) * (*ro) * (*H) -
		      (x1 + fabs(x1)) / 2. * (*(Xi-1)) * (*(ro-1)) * (*(H-1)) +
		      (x2 - fabs(x2)) / 2. * (*(Xi+1)) * (*(ro+1)) * (*(H+1))) / def.hx * (-1.0);
		}
	case 'y':
		{
			x2 = -right_difference (P, 'y') + def.g_const * (*ro);
			x1 = -left_difference (P, 'y') + def.g_const * (*ro);
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (*Xi) * (*ro) * (*H) -
		      (x1 + fabs(x1)) / 2. * (*(Xi-def.locNx)) * (*(ro-def.locNx)) * (*(H-def.locNx)) +
		      (x2 - fabs(x2)) / 2. * (*(Xi+def.locNx)) * (*(ro+def.locNx)) * (*(H+def.locNx))) / def.hy * (-1.0);
		}
	case 'z':
		{
			x2 = -right_difference (P, 'z');
			x1 = -left_difference (P, 'z');
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (*Xi) * (*ro) * (*H) -
		      (x1 + fabs(x1)) / 2. * (*(Xi-def.locNx * (def.locNy))) * (*(ro-def.locNx * (def.locNy))) * (*(H-def.locNx * (def.locNy))) +
		      (x2 - fabs(x2)) / 2. * (*(Xi+def.locNx * (def.locNy))) * (*(ro+def.locNx * (def.locNy))) * (*(H-def.locNx * (def.locNy)))) / def.hz * (-1.0);
		}
	default:
		{
			print_error("Axis of [directed_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет потока энергии в точке
static double assign_E_flow (int i, int j, int k)
{
	if (INTERNAL_POINT)
	{
		double E_flow = 0;
		int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
#ifdef NR
		if ((def.locNx) > 2)
		{
			E_flow += (directed_difference_E (HostArraysPtr.P_w+local, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, HostArraysPtr.H_w+local, 'x')
					 + directed_difference_E (HostArraysPtr.P_n+local, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, HostArraysPtr.H_n+local, 'x')
					 + directed_difference_E (HostArraysPtr.P_g+local, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, HostArraysPtr.H_g+local, 'x'));
		}
		if ((def.locNy) > 2)
		{
			E_flow += (directed_difference_E (HostArraysPtr.P_w+local, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, HostArraysPtr.H_w+local, 'y')
					 + directed_difference_E (HostArraysPtr.P_n+local, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, HostArraysPtr.H_n+local, 'y')
					 + directed_difference_E (HostArraysPtr.P_g+local, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, HostArraysPtr.H_g+local, 'y'));
		}
		if ((def.locNz) > 2)
		{
			E_flow += (directed_difference_E (HostArraysPtr.P_w+local, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, HostArraysPtr.H_w+local, 'z')
					 + directed_difference_E (HostArraysPtr.P_n+local, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, HostArraysPtr.H_n+local, 'z')
					 + directed_difference_E (HostArraysPtr.P_g+local, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, HostArraysPtr.H_g+local, 'z'));
		}
#else
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
#endif
		test_u(E_flow, __FILE__, __LINE__);
		return E_flow;
	}
	else
		return 0;
}

// Расчет внутренней энергии всей системы в точке
void assign_E_current (int local)
{
	HostArraysPtr.E[local] = (HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (HostArraysPtr.ro_w[local] * HostArraysPtr.H_w[local] - HostArraysPtr.P_w[local])
		+ HostArraysPtr.S_n[local] * (HostArraysPtr.ro_n[local] * HostArraysPtr.H_n[local] - HostArraysPtr.P_n[local])
		+ HostArraysPtr.S_g[local] * (HostArraysPtr.ro_g[local] * HostArraysPtr.H_g[local] - HostArraysPtr.P_g[local])) 
		+ (1. - HostArraysPtr.m[local]) * (def.ro_r * HostArraysPtr.H_r[local] - HostArraysPtr.P_w[local]));

	test_nan(HostArraysPtr.E[local], __FILE__, __LINE__);
}

// Расчет внутренней энергии всей системы в точке на следующем шаге по времени
void assign_E_new (int i, int j, int k)
{
	if (INTERNAL_POINT)
	{
		double Q_hw = 0, Q_hr = 0; // источниковые члены

		int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		HostArraysPtr.E_new[local] = HostArraysPtr.E[local] + (def.dt) * (assign_T_flow(i, j, k) + Q_hw + Q_hr - assign_E_flow(i, j, k));

		test_nan(HostArraysPtr.E_new[local], __FILE__, __LINE__);
	}
}

// Задание граничных условий на температуру
void Border_T(int i, int j, int k)
{
	if (BOUNDARY_POINT)
	{
		int local1 = set_boundary_basic_coordinate(i, j, k);
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
void Newton(int i, int j, int k)
{
	if (INTERNAL_POINT)
	{
		enum { n = 5 }; // Размерность системы
		double F[n]; // Вектор значений функций (из системы уравнений)
		double correction[n]; // Вектор поправок к функциям
		double dF[n*n]; // Матрица Якоби (в виде одномерного массива)

		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			F[0] = HostArraysPtr.S_g[local] + HostArraysPtr.S_w[local] + HostArraysPtr.S_n[local] - 1.;
			F[1] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w') * HostArraysPtr.S_w[local] - HostArraysPtr.roS_w[local];
			F[2] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n') * HostArraysPtr.S_n[local] - HostArraysPtr.roS_n[local];
			F[3] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g') * HostArraysPtr.S_g[local] - HostArraysPtr.roS_g[local];
			F[4] = HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w') * HostArraysPtr.H_w[local] - HostArraysPtr.P_w[local])
				+ HostArraysPtr.S_n[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n') * HostArraysPtr.H_n[local] - HostArraysPtr.P_w[local])
				+ HostArraysPtr.S_g[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g') * HostArraysPtr.H_g[local] - HostArraysPtr.P_w[local]))
				+ (1. - HostArraysPtr.m[local]) * (def.ro_r * HostArraysPtr.H_r[local] - HostArraysPtr.P_w[local])
				- HostArraysPtr.E_new[local];

			// Матрица частных производных. Строки: dF/dSw, dF/dSn, dF/dSg, dF/dP, dF/dT

			dF[0] = 1.;
			dF[1] = 1.;
			dF[2] = 1.;
			dF[3] = 0.;
			dF[4] = 0.;

			dF[0 + n] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w');
			dF[1 + n] = 0.;
			dF[2 + n] = 0.;
			dF[3 + n] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', 'P') * HostArraysPtr.S_w[local];
			dF[4 + n] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', 'T') * HostArraysPtr.S_w[local];

			dF[0 + n * 2] = 0.;
			dF[1 + n * 2] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n');
			dF[2 + n * 2] = 0.;
			dF[3 + n * 2] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', 'P') * HostArraysPtr.S_n[local];
			dF[4 + n * 2] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', 'T') * HostArraysPtr.S_n[local];

			dF[0 + n * 3] = 0.;
			dF[1 + n * 3] = 0.;
			dF[2 + n * 3] = ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g');
			dF[3 + n * 3] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', 'P') * HostArraysPtr.S_g[local];
			dF[4 + n * 3] = d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', 'T') * HostArraysPtr.S_g[local];

			dF[0 + n * 4] = HostArraysPtr.m[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w') * HostArraysPtr.H_w[local] - HostArraysPtr.P_w[local]);
			dF[1 + n * 4] = HostArraysPtr.m[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n') * HostArraysPtr.H_n[local] - HostArraysPtr.P_w[local]);
			dF[2 + n * 4] = HostArraysPtr.m[local] * (ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g') * HostArraysPtr.H_g[local] - HostArraysPtr.P_w[local]);
			dF[3 + n * 4] = HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', 'P') * HostArraysPtr.H_w[local] - 1.)
				+ HostArraysPtr.S_n[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', 'P') * HostArraysPtr.H_n[local] - 1.)
				+ HostArraysPtr.S_g[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', 'P') * HostArraysPtr.H_g[local] - 1.))
				+ (1. - HostArraysPtr.m[local]) * (-1);
			dF[4 + n * 4] = HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w', 'T') * HostArraysPtr.H_w[local]
				+ ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'w') * c_w(HostArraysPtr.T[local]))
				+ HostArraysPtr.S_n[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n', 'T') * HostArraysPtr.H_n[local]
				+ ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'n') * c_n(HostArraysPtr.T[local]))
				+ HostArraysPtr.S_g[local] * (d_ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g', 'T') * HostArraysPtr.H_g[local]
				+ ro(HostArraysPtr.P_w[local], HostArraysPtr.T[local], 'g') * c_g(HostArraysPtr.T[local])))
				+ (1. - HostArraysPtr.m[local]) * def.ro_r * c_r(HostArraysPtr.T[local]);

			reverse_matrix(dF, n);
			mult_matrix_vector(correction, dF, F, n);

			HostArraysPtr.S_w[local] = HostArraysPtr.S_w[local] - correction[0];
			HostArraysPtr.S_n[local] = HostArraysPtr.S_n[local] - correction[1];
			HostArraysPtr.S_g[local] = HostArraysPtr.S_g[local] - correction[2];
			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local] - correction[3];
			HostArraysPtr.T[local] = HostArraysPtr.T[local] - correction[4];
			assign_H(local);
		}

		// Обновление значения суммарной энергии, т.к. оно больше не понадобится
		// !!! Лучше вынести в отдельную функцию (просто обменять указатели).
		HostArraysPtr.E[local] = HostArraysPtr.E_new[local];

		test_S(HostArraysPtr.S_w[local], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.T[local], __FILE__, __LINE__);
		test_nan(HostArraysPtr.E[local], __FILE__, __LINE__);
	}
}
