#include "defines.h"
#include "three-phase.h"

extern consts def;
extern ptr_Arrays HostArraysPtr;

// Номер среды
const int media = 0;
// Переломные точки насыщенностей при вычислении капиллярных давлений
const double S_w_range[2] = {0.1, 0.9};
const double S_g_range[2] = {0.1, 0.9};

// Функции вычисления эффективных значений насыщенностей
double assign_S_w_e(int local)
{
	return (HostArraysPtr.S_w[local] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
}

double assign_S_n_e(int local)
{
	return (HostArraysPtr.S_n[local] - def.S_nr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
}

double assign_S_g_e(int local)
{
	return (HostArraysPtr.S_g[local] - def.S_gr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
}

// Вычисление капиллярных давлений
// Функции кап. давлений и их производных для центральной части интервала
static inline double P_k_nw(double S)
{
	double A = def.lambda[media];
	return def.P_d_nw[media] * pow((pow(S, A / (1. - A)) - 1.), 1. / A);
}

static inline double P_k_gn(double S)
{
	double A = def.lambda[media];
	return def.P_d_gn[media] * pow(pow((1. - S), A / (1. - A)) - 1., 1. / A);
}

static inline double P_k_nw_S(double S)
{
	double A = def.lambda[media];
	return def.P_d_nw[media] * pow(pow(S, A / (1. - A)) - 1., 1. / A - 1.) * pow(S, (A / (1. - A) - 1.)) / (1. - A)
		/ (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
}

static inline double P_k_gn_S(double S)
{
	double A = def.lambda[media];
	return def.P_d_gn[media] * pow(pow(1. - S, A / (1. - A)) - 1., 1. / A - 1.) * pow(1. - S, A / (1. - A) - 1.) / (1. - A)
		/ (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]) * (-1.0);
}

// Функции вычисления капиллярных давлений и производных на всем интервале
// По краям интервала [0, 1] функции капиллярных давлений гладко заменяем линейными, производные меняются соответственно.
// Описание можно посмотреть в файле mathcad.
double assign_P_k_nw(double S_w_e)
{
	double Pk_nw = 0;

	if (S_w_e <= S_w_range[0])
	{
		Pk_nw = P_k_nw_S(S_w_range[0]) * (S_w_e - S_w_range[0]) + P_k_nw(S_w_range[0]);
	}
	else if (S_w_e >= S_w_range[1])
	{
		Pk_nw = P_k_nw_S(S_w_range[1]) * (S_w_e - S_w_range[1]) + P_k_nw(S_w_range[1]);;
	}
	else
	{
		Pk_nw = P_k_nw(S_w_e);
	}

	return Pk_nw;
}

double assign_P_k_gn(double S_g_e)
{
	double Pk_gn = 0;

	if (S_g_e <= S_g_range[0])
	{
		Pk_gn = P_k_gn_S(S_g_range[0]) * (S_g_e - S_g_range[0]) + P_k_gn(S_g_range[0]);
	}
	else if (S_g_e >= S_g_range[1])
	{
		Pk_gn = P_k_gn_S(S_g_range[1]) * (S_g_e - S_g_range[1]) + P_k_gn(S_g_range[1]);
	}
	else
	{
		Pk_gn = P_k_gn(S_g_e);
	}

	return Pk_gn;
}

// Функции вычисления производных капиллярных давлений по насыщенностям
static inline double assign_P_k_nw_S(double S_w_e)
{
	double PkSw = 0;

	if (S_w_e <= S_w_range[0])
	{
		PkSw = P_k_nw_S(S_w_range[0]);
	}
	else if (S_w_e >= S_w_range[1])
	{
		PkSw = P_k_nw_S(S_w_range[1]);
	}
	else
	{
		PkSw = P_k_nw_S(S_w_e);
	}

	return PkSw;
}

static inline double assign_P_k_gn_S(double S_g_e)
{
	double PkSn = 0;

	if (S_g_e <= S_g_range[0])
	{
		PkSn = P_k_gn_S(S_g_range[0]);
	}
	else if (S_g_e >= S_g_range[1])
	{
		PkSn = P_k_gn_S(S_g_range[1]);
	}
	else
	{
		PkSn = P_k_gn_S(S_g_e);
	}

	return PkSn;
}

// Функции вычисления относительных проницаемостей
static inline double assign_k_w(double S_w_e)
{
	double A = def.lambda[media];
	double k_w = 0;

	if (S_w_e >= 1e-3)
	{
		k_w = pow(S_w_e, 0.5) * pow(1. - pow(1. - pow(S_w_e, A / (A - 1.)), (A - 1.) / A), 2.);
	}

	return k_w;
}

static inline double assign_k_g(double S_g_e)
{
	double A = def.lambda[media];
	double k_g = 0;

	if (S_g_e >= 1e-3)
	{
		k_g = pow(S_g_e, 0.5) * pow(1. - pow(1. - S_g_e, A / (A - 1.)), 2. * (A - 1.) / A);
	}

	return k_g;
}

static inline double assign_k_n(double S_w_e, double S_n_e)
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

void prepare_local_vars(int i, int j, int k)
{
	int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
	double k_w, k_g, k_n, Pk_nw, Pk_gn;

	assign_S(local);

	double S_w_e = assign_S_w_e(local);
	double S_n_e = assign_S_n_e(local);
	double S_g_e = 1. - S_w_e - S_n_e;

	k_w = assign_k_w(S_w_e);
	k_g = assign_k_g(S_g_e);
	k_n = assign_k_n(S_w_e, S_n_e);

	Pk_nw = assign_P_k_nw(S_w_e);
	Pk_gn = assign_P_k_gn(S_g_e);

	HostArraysPtr.P_n[local] = HostArraysPtr.P_w[local] + Pk_nw;
	HostArraysPtr.P_g[local] = HostArraysPtr.P_w[local] + Pk_nw + Pk_gn;

	assign_ro(local);

#ifdef ENERGY
	assign_H(local);
	assign_E_current(local);

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

	test_positive(HostArraysPtr.P_n[local], __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_g[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_g[local], __FILE__, __LINE__);
}

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
void Newton(int i, int j, int k)
{
	if (INTERNAL_POINT)
	{
		double S_w_e, S_g_e, S_n_e, Pk_nw, Pk_gn, PkSw, PkSn, Sg, F1, F2, F3;
		double dF[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			S_w_e = assign_S_w_e(local);
			S_n_e = assign_S_n_e(local);
			S_g_e = 1. - S_w_e - S_n_e;

			Pk_nw = assign_P_k_nw(S_w_e);
			Pk_gn = assign_P_k_gn(S_g_e);
			PkSw = assign_P_k_nw_S(S_w_e);
			PkSn = assign_P_k_gn_S(S_g_e);

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
void Border_S(int i, int j, int k)
{
	if (BOUNDARY_POINT)
	{
		int local1 = set_boundary_basic_coordinate(i, j, k);
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

void Border_P(int i, int j, int k)
{
	if (BOUNDARY_POINT)
	{
		int local1 = set_boundary_basic_coordinate(i, j, k);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		double S_w_e = assign_S_w_e(local1);
		double S_n_e = assign_S_n_e(local1);
		double S_g_e = 1. - S_w_e - S_n_e;

		double Pk_nw = assign_P_k_nw(S_w_e);
		double Pk_gn = assign_P_k_gn(S_g_e);

		// Если отдельно задаем значения на границах через градиент (условия непротекания)
		if ((j != 0) && (j != (def.locNy) - 1))
		{
			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local1];
			HostArraysPtr.P_n[local] = HostArraysPtr.P_w[local1] + Pk_nw;
			HostArraysPtr.P_g[local] = HostArraysPtr.P_w[local1] + Pk_nw + Pk_gn;
		
		}
		else if (j == 0)
		{
			HostArraysPtr.P_w[local] = 1.1 * def.P_atm;
			HostArraysPtr.P_n[local] = 1.1 * def.P_atm;
			HostArraysPtr.P_g[local] = 1.1 * def.P_atm;

			/*if((i > (def.locNx) / 3) && (i < 2 * (def.locNx) / 3) && (((def.locNz) < 2) || (k > (def.locNz) / 3) && (k < 2 * (def.locNz) / 3)))
			{*/
				//Открытая верхняя граница
			/*	HostArraysPtr.P_w[local] = def.P_atm;
				HostArraysPtr.P_n[local] = def.P_atm;
				HostArraysPtr.P_g[local] = def.P_atm;
			*/
			/*}
			else*/
			/*{
				// Условия непротекания
				HostArraysPtr.P_w[local] = (HostArraysPtr.P_w[local1]
				- (def.ro0_w) * (def.g_const) * (def.hy) * (1. - (def.beta_w) * (def.P_atm))) 
					/ (1. + (def.beta_w) * (def.ro0_w) * (def.g_const) * (def.hy));
				HostArraysPtr.P_n[local] = (HostArraysPtr.P_w[local1]
				+ Pk_nw - (def.ro0_n) * (def.g_const) * (def.hy) * (1. - (def.beta_n) * (def.P_atm))) 
					/ (1. + (def.beta_n) * (def.ro0_n) * (def.g_const) * (def.hy));
				HostArraysPtr.P_g[local] = (HostArraysPtr.P_w[local1]
				+ Pk_nw + Pk_gn) / (1. + (def.ro0_g) * (def.g_const) * (def.hy) / (def.P_atm));
			}*/
		}
		else
		{
			HostArraysPtr.P_w[local] = def.P_atm;
			HostArraysPtr.P_n[local] = def.P_atm;
			HostArraysPtr.P_g[local] = def.P_atm;

			/*HostArraysPtr.P_w[local] = 1.1 * def.P_atm;
			HostArraysPtr.P_n[local] = 1.1 * def.P_atm;
			HostArraysPtr.P_g[local] = 1.1 * def.P_atm;*/

			// Условия непротекания (normal u = 0)
/*#ifdef ENERGY
			HostArraysPtr.P_w[local] = (HostArraysPtr.P_w[local1]
			+ (def.ro0_w) * (def.g_const) * (def.hy) * (1. - (def.beta_w) * (def.P_atm))
			- (def.alfa_w) * (HostArraysPtr.T[local] - def.T_0))
				/ (1. - (def.beta_w) * (def.ro0_w) * (def.g_const) * (def.hy));
			HostArraysPtr.P_n[local] = (HostArraysPtr.P_w[local1]
			+ Pk_nw + (def.ro0_n) * (def.g_const) * (def.hy) * (1. - (def.beta_n) * (def.P_atm))
			- (def.alfa_n) * (HostArraysPtr.T[local] - def.T_0))
				/ (1. - (def.beta_n) * (def.ro0_n) * (def.g_const) * (def.hy));
			HostArraysPtr.P_g[local] = (HostArraysPtr.P_w[local1]
			+ Pk_nw + Pk_gn) / (1. - (def.ro0_g) * (def.g_const) * (def.hy) * (def.T_0)
				/ ((def.P_atm) * HostArraysPtr.T[local]));
#else
			HostArraysPtr.P_w[local] = (HostArraysPtr.P_w[local1]
			+ (def.ro0_w) * (def.g_const) * (def.hy) * (1. - (def.beta_w) * (def.P_atm)))
				/ (1. + (def.beta_w) * (def.ro0_w) * (def.g_const) * (def.hy));
			HostArraysPtr.P_n[local] = (HostArraysPtr.P_w[local1]
			+ Pk_nw + (def.ro0_n) * (def.g_const) * (def.hy) * (1. - (def.beta_n) * (def.P_atm)))
				/ (1. + (def.beta_n) * (def.ro0_n) * (def.g_const) * (def.hy));
			HostArraysPtr.P_g[local] = (HostArraysPtr.P_w[local1]
			+ Pk_nw + Pk_gn) / (1. - (def.ro0_g) * (def.g_const) * (def.hy) / (def.P_atm));
#endif*/
		}
		test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_n[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_g[local], __FILE__, __LINE__);
	}
}

#ifdef ENERGY
// Задание граничных условий на температуру
void Border_T(int i, int j, int k)
{
	if (BOUNDARY_POINT)
	{
		int local1 = set_boundary_basic_coordinate(i, j, k);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if (j == 0)
		{
			//HostArraysPtr.T[local] = 400;
			//HostArraysPtr.T[local] = 285;
			HostArraysPtr.T[local] = 320;
		}
		//else if(j == (def.locNy) - 1)
		//{
			//HostArraysPtr.T[local] = 273;
			//HostArraysPtr.T[local] = 285;
		//}
		else
		{
			// Будем считать границы области не теплопроводящими
			HostArraysPtr.T[local] = HostArraysPtr.T[local1];
		}

		test_positive(HostArraysPtr.T[local], __FILE__, __LINE__);
	}
}
#endif

// Является ли точка нагнетательной скважиной
int is_injection_well(int i, int j, int k)
{

		return 0;
}

// Является ли точка добывающей скважиной
int is_output_well(int i, int j, int k)
{
		return 0;
}

// Устанавливает значения втекаемых/вытекаемых жидкостей q_i на скважинах
void wells_q(int i, int j, int k, double* q_w, double* q_n, double* q_g)
{
	*q_w = 0.0;
	*q_g = 0.0;
	*q_n = 0.0;

/*	int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

	if (local_to_global(j, 'y') == 1 && (HostArraysPtr.S_w[local] <= 0.6))
		*q_w = 1.0 * (def.dt);
*/
}



