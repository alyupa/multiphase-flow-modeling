#include "defines.h"
#include "three-phase.h"

// Номер среды
int media = 0;
// Переломные точки насыщенностей при вычислении капиллярных давлений
double S_w_range[2] = {0.1, 0.99};
double S_g_range[2] = {0.005, 0.95};

// Функции вычисления эффективных значений насыщенностей
double assign_S_w_e(const ptr_Arrays &HostArraysPtr, int local, const consts &def)
{
	return (HostArraysPtr.S_w[local] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
}

double assign_S_n_e(const ptr_Arrays &HostArraysPtr, int local, const consts &def)
{
	return (HostArraysPtr.S_n[local] - def.S_nr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
}

// Вычисление капиллярных давлений
// Функции кап. давлений и их производных для центральной части интервала
double P_k_nw(double S, const consts &def)
{
	return 0;
	double A = def.lambda[media];
	return def.P_d_nw[media] * pow((pow(S, A / (1. - A)) - 1.), 1. / A);
}

double P_k_gn(double S, const consts &def)
{
	return 0;
	double A = def.lambda[media];
	return def.P_d_gn[media] * pow(pow((1. - S), A / (1. - A)) - 1., 1. / A);
}

double P_k_nw_S(double S, const consts &def)
{
	return 0;
	double A = def.lambda[media];
	return def.P_d_nw[media] * pow(pow(S, A / (1. - A)) - 1., 1. / A - 1.) * pow(S, (A / (1. - A) - 1.)) / (1. - A)
		/ (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
}

double P_k_gn_S(double S, const consts &def)
{
	return 0;
	double A = def.lambda[media];
	return def.P_d_gn[media] * pow(pow(1. - S, A / (1. - A)) - 1., 1. / A - 1.) * pow(1. - S, A / (1. - A) - 1.) / (1. - A)
		/ (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
}

// Функции вычисления капиллярных давлений и производных на всем интервале
// По краям интервала [0, 1] функции капиллярных давлений гладко заменяем линейными, производные меняются соответственно.
// Описание можно посмотреть в файле mathcad.
double assign_P_k_nw(double S_w_e, const consts &def)
{
	return 0;

	double Pk_nw = 0;

	if (S_w_e <= S_w_range[0])
	{
		Pk_nw = P_k_nw_S(S_w_range[0], def) * (S_w_e - S_w_range[0]) + P_k_nw(S_w_range[0], def);
	}
	else if (S_w_e >= S_w_range[1])
	{
		Pk_nw = P_k_nw_S(S_w_range[1], def) * (S_w_e - S_w_range[1]) + P_k_nw(S_w_range[1], def);;
	}
	else
	{
		Pk_nw = P_k_nw(S_w_e, def);
	}

	return Pk_nw;
}

double assign_P_k_gn(double S_g_e, const consts &def)
{
	return 0;

	double Pk_gn = 0;

	if (S_g_e <= S_g_range[0])
	{
		Pk_gn = P_k_gn_S(S_g_range[0], def) * (S_g_e - S_g_range[0]) + P_k_gn(S_g_range[0], def);
	}
	else if (S_g_e >= S_g_range[1])
	{
		Pk_gn = P_k_gn_S(S_g_range[1], def) * (S_g_e - S_g_range[1]) + P_k_gn(S_g_range[1], def);
	}
	else
	{
		Pk_gn = P_k_gn(S_g_e, def);
	}

	return Pk_gn;
}

// Функции вычисления производных капиллярных давлений по насыщенностям
double assign_P_k_nw_S(double S_w_e, const consts &def)
{
	return 0;

	double PkSw = 0;

	if (S_w_e <= S_w_range[0])
	{
		PkSw = P_k_nw_S(S_w_range[0], def);
	}
	else if (S_w_e >= S_w_range[1])
	{
		PkSw = P_k_nw_S(S_w_range[1], def);
	}
	else
	{
		PkSw = P_k_nw_S(S_w_e, def);
	}

	return PkSw;
}

double assign_P_k_gn_S(double S_g_e, const consts &def)
{
	return 0;

	double PkSn = 0;

	if (S_g_e <= S_g_range[0])
	{
		PkSn = (-1) * P_k_gn_S(S_g_range[0], def);
	}
	else if (S_g_e >= S_g_range[1])
	{
		PkSn = (-1) * P_k_gn_S(S_g_range[1], def);
	}
	else
	{
		PkSn = P_k_gn_S(S_g_e, def);
	}

	return PkSn;
}

// Функции вычисления относительных проницаемостей
double assign_k_w(double S_w_e, const consts &def)
{
	double A = def.lambda[media];
	double k_w = 0;

	if (S_w_e >= 1e-3)
	{
		k_w = pow(S_w_e, 0.5) * pow(1. - pow(1. - pow(S_w_e, A / (A - 1.)), (A - 1.) / A), 2.);
	}

	return k_w;
}

double assign_k_g(double S_g_e, const consts &def)
{
	double A = def.lambda[media];
	double k_g = 0;

	if (S_g_e >= 1e-3)
	{
		k_g = pow(S_g_e, 0.5) * pow(1. - pow(1. - S_g_e, A / (A - 1.)), 2. * (A - 1.) / A);
	}

	return k_g;
}

double assign_k_n(double S_w_e, double S_n_e, const consts &def)
{
	double A = def.lambda[media];
	double k_n = 0;
	double S_g_e = 1. - S_w_e - S_n_e;

	if (S_n_e >= 1e-3)
	{
		double k_n_w = pow(1. - S_w_e, 0.5) * pow(1. - pow(S_w_e, A / (A - 1.)), 2. * (A - 1.) / A);
		double k_n_g = pow(S_n_e, 0.5) * pow(1. - pow(1. - pow(S_n_e, A / (A - 1.)), (A - 1.) / A), 2.);
		k_n = S_n_e * k_n_w * k_n_g / (1 - S_w_e) / (1 - S_g_e);
	}

	return k_n;
}

//Функция вычисления значений давлений, плотностей и коэффициентов в законе Дарси в точке (i,j,k) среды media,
//исходя из известных значений основных параметров (Pw,Sw,Sn)
//1. Запоминаем, с какой именно из сред работаем
//2. Вычисление значения насыщенности фазы n из условия равенства насыщенностей в сумме единице
//3. Вычисление эффективных насыщенностей по формулам модели трехфазной фильтрации
//4. Вычисление относительных фазовых проницаемостей в соответствии с приближением Стоуна в модификации Азиза и Сеттари
//5. Вычисление капиллярных давлений в соответствии с приближенной моделью Паркера
//6. Вычисление фазовых давлений c помощью капиллярных
//7. Вычисление коэффициентов закона Дарси

void assign_P_Xi(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
	double k_w, k_g, k_n, Pk_nw, Pk_gn;
	double S_w_e = assign_S_w_e(HostArraysPtr, local, def);
	double S_n_e = assign_S_n_e(HostArraysPtr, local, def);
	double S_g_e = 1. - S_w_e - S_n_e;

	k_w = assign_k_w(S_w_e, def);
	k_g = assign_k_g(S_g_e, def);
	k_n = assign_k_n(S_w_e, S_n_e, def);

#ifdef ENERGY
	// Вынести в константы!!!
	double mu_w = 1. / (29.21 * HostArraysPtr.T[local] - 7506.64);
	double mu_n = 7.256E-10 * exp(4141.9 / HostArraysPtr.T[local]);
	double mu_g = 1.717E-5 * pow((HostArraysPtr.T[local] / 273.), 0.683);

	HostArraysPtr.Xi_w[local] = (-1.) * (def.K[media]) * k_w / mu_w;
	HostArraysPtr.Xi_n[local] = (-1.) * (def.K[media]) * k_n / mu_n;
	HostArraysPtr.Xi_g[local] = (-1.) * (def.K[media]) * k_g / mu_g;
#else
	HostArraysPtr.Xi_w[local] = (-1.) * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[local] = (-1.) * (def.K[media]) * k_n / def.mu_n;
	HostArraysPtr.Xi_g[local] = (-1.) * (def.K[media]) * k_g / def.mu_g;
#endif

	if ((((i != 0) && (i != (def.locNx) - 1)) || ((def.locNx) < 2)) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		Pk_nw = assign_P_k_nw(S_w_e, def);
		Pk_gn = assign_P_k_gn(S_g_e, def);

		HostArraysPtr.P_n[local] = HostArraysPtr.P_w[local] + Pk_nw;
		HostArraysPtr.P_g[local] = HostArraysPtr.P_w[local] + Pk_nw + Pk_gn;

		test_positive(HostArraysPtr.P_n[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_g[local], __FILE__, __LINE__);
	}

	test_nan(HostArraysPtr.Xi_w[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_g[local], __FILE__, __LINE__);
}

// Вспомогательная функции для метода Ньютона:
// Нахождение обратной матрицы 3*3;
// Теперь будет использоваться функция reverse_matrix(double* aб int n) из gauss.cpp
/*void reverse_matrix(double* a)
{
	int n = 3;
	double b[9], det = 0;

	// Вычисление дополнительных миноров матрицы
	for(int j = 0; j < n; j++)
		for(int i = 0; i < n; i++)
		{
			b[i + n * j] = a[(i + 1) % n + n * ((j + 1) % n)] * a[(i + 2) % n + n * ((j + 2) % n)]
			- a[(i + 2) % n + n * ((j + 1) % n)] * a[(i + 1) % n + n * ((j + 2) % n)];
		}

	// Нахождение детерминанта матрицы 3*3;
	for(int i = 0; i < n; i++)
	{
		det += a[i] * b[i];
	}
	test_nan(det, __FILE__, __LINE__);

	// Транспонирование и деление на детерминант
	for(int j = 0; j < n; j++)
		for(int i = 0; i < n; i++)
		{
			a[i + n * j] = b[j + n * i] / det;
			test_nan(a[i + n * j], __FILE__, __LINE__);
		}
}
*/

//Функция решения системы 3*3 на основные параметры (Pn,Sw,Sg) методом Ньютона в точке (i,j,k) среды media
//1. Вычисление эффективных насыщенностей
//2. Переобозначение степенного коэффициента
//3. Вычисление капиллярных давлений
//4. Вычисление насыщенности фазы n
//5. Нахождение значения трех функций системы
//6. Вычисление частных производных производных капиллярных давлений по насыщенностям
//7. Вычисление матрицы частных производных
//8. Вычисление детерминанта матрицы частных производных
//9. Получение решения системы методом Крамера в явном виде
#ifndef ENERGY
void Newton(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	if ((((i != 0) && (i != (def.locNx) - 1)) || ((def.locNx) < 2)) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		double S_w_e, S_g_e, S_n_e, Pk_nw, Pk_gn, PkSw, PkSn, Sg, F1, F2, F3;
		double dF[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			S_w_e = assign_S_w_e(HostArraysPtr, local, def);
			S_n_e = assign_S_n_e(HostArraysPtr, local, def);
			S_g_e = 1. - S_w_e - S_n_e;

			Pk_nw = assign_P_k_nw(S_w_e, def);
			Pk_gn = assign_P_k_gn(S_g_e, def);
			PkSw = assign_P_k_nw_S(S_w_e, def);
			PkSn = assign_P_k_gn_S(S_g_e, def);

			Sg = 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local];

			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm))
			     * HostArraysPtr.S_w[local] - HostArraysPtr.roS_w[local];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[local] + Pk_nw - def.P_atm))
			     * HostArraysPtr.S_n[local] - HostArraysPtr.roS_n[local];
			F3 = def.ro0_g * (HostArraysPtr.P_w[local] + Pk_nw + Pk_gn) / def.P_atm
			     * Sg - HostArraysPtr.roS_g[local];

			// Матрица частных производных. Индексу от 0 до 8 соответствуют F1P, F1Sw, F1Sn, F2P, F2Sw, F2Sn, F3P, F3Sw, F3Sn
			dF[0] = def.ro0_w * def.beta_w * HostArraysPtr.S_w[local];
			dF[3] = def.ro0_n * def.beta_n * HostArraysPtr.S_n[local];
			dF[6] = def.ro0_g * Sg / def.P_atm;
			dF[1] = def.ro0_w * (1 + def.beta_w * (HostArraysPtr.P_w[local] - def.P_atm));
			dF[4] = def.ro0_n * (1. + (def.beta_n) * PkSw) * HostArraysPtr.S_n[local];
			dF[7] = (-1) * def.ro0_g * (HostArraysPtr.P_w[local] + Pk_nw + Pk_gn - Sg * (PkSn + PkSw)) / def.P_atm;
			dF[2] = 0;
			dF[5] = def.ro0_n * (1. + def.beta_n * (HostArraysPtr.P_w[local] + Pk_nw - def.P_atm));
			dF[8] = (-1) * def.ro0_g * (HostArraysPtr.P_w[local] + Pk_nw + Pk_gn - Sg * PkSn) / def.P_atm;

			reverse_matrix(dF, 3);

			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local]
			        - (dF[0] * F1 + dF[1] * F2 + dF[2] * F3);
			HostArraysPtr.S_w[local] = HostArraysPtr.S_w[local]
			        - (dF[3] * F1 + dF[4] * F2 + dF[5] * F3);
			HostArraysPtr.S_n[local] = HostArraysPtr.S_n[local]
			        - (dF[6] * F1 + dF[7] * F2 + dF[8] * F3);
		}

		test_S(HostArraysPtr.S_w[local], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
	}
}
#endif

//Задание граничных условий отдельно для (Sw,Sg),Pn

// Задание граничных условий с меньшим числом проверок, но с введением дополнительных переменных
void Border_S(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	if ((((i == 0) || (i == (def.locNx) - 1)) && ((def.locNx) >= 2)) || (j == 0) || (j == (def.locNy) - 1) || (((k == 0) || (k == (def.locNz) - 1)) && ((def.locNz) >= 2)))
	{
		int local1 = set_boundary_basic_coordinate(i, j, k, def);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if ((j != 0) || ((def.source) <= 0))
		{
			HostArraysPtr.S_w[local] = HostArraysPtr.S_w[local1];
			HostArraysPtr.S_n[local] = HostArraysPtr.S_n[local1];
		}

		if ((j == 0) && ((def.source) > 0))
		{
			HostArraysPtr.S_w[local] = def.S_w_gr;
			HostArraysPtr.S_n[local] = def.S_n_gr;
		}
		test_S(HostArraysPtr.S_w[local], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);
	}
}

void Border_P(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	if ((((i == 0) || (i == (def.locNx) - 1)) && ((def.locNx) >= 2)) || (j == 0) || (j == (def.locNy) - 1) || (((k == 0) || (k == (def.locNz) - 1)) && ((def.locNz) >= 2)))
	{
		int local1 = set_boundary_basic_coordinate(i, j, k, def);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		double S_w_e = assign_S_w_e(HostArraysPtr, local1, def);
		double S_n_e = assign_S_n_e(HostArraysPtr, local1, def);
		double S_g_e = 1. - S_w_e - S_n_e;

		double Pk_nw = assign_P_k_nw(S_w_e, def);
		double Pk_gn = assign_P_k_gn(S_g_e, def);

		// Если отдельно задаем значения на границах через градиент (условия непротекания)
		if ((j != 0) && (j != (def.locNy) - 1))
		{
			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local1];
			HostArraysPtr.P_n[local] = HostArraysPtr.P_w[local1] + Pk_nw;
			HostArraysPtr.P_g[local] = HostArraysPtr.P_w[local1] + Pk_nw + Pk_gn;
		
		}
		else if (j == 0)
		{
			/*if((i > (def.locNx) / 3) && (i < 2 * (def.locNx) / 3) && (((def.locNz) < 2) || (k > (def.locNz) / 3) && (k < 2 * (def.locNz) / 3)))
			{
				//Открытая верхняя граница
				HostArraysPtr.P_w[local] = def.P_atm;
				HostArraysPtr.P_n[local] = def.P_atm;
				HostArraysPtr.P_g[local] = def.P_atm;
			}
			else*/
			{
				// Условия непротекания
				HostArraysPtr.P_w[local] = (HostArraysPtr.P_w[local1]
				- (def.ro0_w) * (def.g_const) * (def.hy) * (1. - (def.beta_w) * (def.P_atm))) 
					/ (1. + (def.beta_w) * (def.ro0_w) * (def.g_const) * (def.hy));
				HostArraysPtr.P_n[local] = (HostArraysPtr.P_w[local1]
				+ Pk_nw - (def.ro0_n) * (def.g_const) * (def.hy) * (1. - (def.beta_n) * (def.P_atm))) 
					/ (1. + (def.beta_n) * (def.ro0_n) * (def.g_const) * (def.hy));
				HostArraysPtr.P_g[local] = (HostArraysPtr.P_w[local1]
				+ Pk_nw + Pk_gn) / (1. + (def.ro0_g) * (def.g_const) * (def.hy) / (def.P_atm));
			}
		}
		else
		{
			HostArraysPtr.P_w[local] = (HostArraysPtr.P_w[local1]
			+ (def.ro0_w) * (def.g_const) * (def.hy) * (1. - (def.beta_w) * (def.P_atm))) 
				/ (1. - (def.beta_w) * (def.ro0_w) * (def.g_const) * (def.hy));
			HostArraysPtr.P_n[local] = (HostArraysPtr.P_w[local1]
			+ Pk_nw + (def.ro0_n) * (def.g_const) * (def.hy) * (1. - (def.beta_n) * (def.P_atm))) 
				/ (1. - (def.beta_n) * (def.ro0_n) * (def.g_const) * (def.hy));
			HostArraysPtr.P_g[local] = (HostArraysPtr.P_w[local1]
			+ Pk_nw + Pk_gn) / (1. - (def.ro0_g) * (def.g_const) * (def.hy) / (def.P_atm));
		}
		test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_n[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_g[local], __FILE__, __LINE__);
	}
}

// Является ли точка нагнетательной скважиной
int is_injection_well(int i, int j, int k, const consts &def)
{

		return 0;
}

// Является ли точка добывающей скважиной
int is_output_well(int i, int j, int k, const consts &def)
{
		return 0;
}

// Устанавливает значения втекаемых/вытекаемых жидкостей q_i на скважинах
void wells_q(const ptr_Arrays &HostArraysPtr, int i, int j, int k, double* q_w, double* q_n, double* q_g, const consts &def)
{
	*q_w = 0.0;
	*q_g = 0.0;
	*q_n = 0.0;
}



