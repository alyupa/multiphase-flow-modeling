#include "../defines.h"
#include "two-phase.h"

void assign_P_Xi(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	int media = 0;
	int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
	double S_e = (1. - HostArraysPtr.S_n[local] - def.S_wr[media]) / (1. - def.S_wr[media]);
	//if (S_e<0)
	//	S_e=0;
	double k_w = pow(S_e, (2. + 3. * (def.lambda[media])) / def.lambda[media]);
	double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + def.lambda[media]) / def.lambda[media]));
	double P_k = def.P_d[media] * pow((1. - HostArraysPtr.S_n[local] - def.S_wr[media]) / (1. - def.S_wr[media]), -1. / def.lambda[media]);

	HostArraysPtr.P_n[local] = HostArraysPtr.P_w[local] + P_k;
	HostArraysPtr.Xi_w[local] = -1. * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[local] = -1. * (def.K[media]) * k_n / def.mu_n;

	test_S(S_e, __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_n[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[local], __FILE__, __LINE__);
}

void Newton(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		int media = 0;
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
		double S_e, P_k, AAA, F1, F2, PkS, F1P, F2P, F1S, F2S, det;

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			S_e = (1 - HostArraysPtr.S_n[local] - def.S_wr[media]) / (1 - def.S_wr[media]);
			P_k = def.P_d[media] * pow(S_e, (-1.) / def.lambda[media]);
			AAA = pow(S_e, (((-1.) / def.lambda[media]) - 1.));
			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm)) * (1. - HostArraysPtr.S_n[local]) - HostArraysPtr.roS_w[local];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[local] + P_k - def.P_atm)) * HostArraysPtr.S_n[local] - HostArraysPtr.roS_n[local];
			PkS = AAA * (def.P_d[media]) / (def.lambda[media] * (1 - def.S_wr[media]));
			F1P = def.ro0_w * (def.beta_w) * (1 - HostArraysPtr.S_n[local]);
			F2P = def.ro0_n * (def.beta_n) * HostArraysPtr.S_n[local];
			F1S = (-1) * (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm));
			F2S = def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_w[local] + P_k - def.P_atm + (HostArraysPtr.S_n[local] * PkS)));

			det = F1P * F2S - F1S * F2P;

			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local] - (1 / det) * (F2S * F1 - F1S * F2);
			HostArraysPtr.S_n[local] = HostArraysPtr.S_n[local] - (1 / det) * (F1P * F2 - F2P * F1);
		}

		test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);

	}
}

// Задание граничных условий
void Border_S(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	if ((i == 0) || (i == (def.locNx) - 1) || (j == 0) || (j == (def.locNy) - 1) || (((k == 0) || (k == (def.locNz) - 1)) && ((def.locNz) >= 2)))
	{
		int local1=set_boundary_basic_coordinate(i, j, k, def);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if ((j != 0) || ((def.source) <= 0))
		{
			HostArraysPtr.S_n[local] = HostArraysPtr.S_n[local1];
		}

		if ((j == 0) && ((def.source) > 0))
		{
			int I = local_to_global(i, 'x', def);
			if ((I >= (def.Nx) / 2 - (def.source)) && (I <= (def.Nx) / 2 + (def.source)) && (k >= (def.Nz) / 2 - (def.source)) && (k <= (def.Nz) / 2 + (def.source)))
			{
				HostArraysPtr.S_n[local] = def.S_n_gr;
			}
			else
				//HostArraysPtr.S_n[local] = 0;
			{
				HostArraysPtr.S_n[local] = HostArraysPtr.S_n[local1];
			}
		}

		test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);
	}
}

void Border_P(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	if ((i == 0) || (i == (def.locNx) - 1) || (j == 0) || (j == (def.locNy) - 1) || (((k == 0) || (k == (def.locNz) - 1)) && ((def.locNz) >= 2)))
	{
		int local1=set_boundary_basic_coordinate(i, j, k, def);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if ((j != 0) && (j != (def.locNy) - 1))
		{
			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local1];
		}
		else if (j == 0)
		{
			HostArraysPtr.P_w[local] = def.P_atm;
		}
		else
		{
			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local1] + ro_eff_gdy(HostArraysPtr, local1, def);
		}

		test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
	}
}

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


