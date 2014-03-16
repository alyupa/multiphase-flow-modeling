#include "defines.h"

extern double *HostBuffer;

void prepare_all_vars()
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
			{
				prepare_local_vars(i, j, k);
			}
}

void u_calculation()
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (ACTIVE_POINT)
				{
					assign_u(i, j, k);
				}
}

void find_values_from_partial_equations(double t)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (ACTIVE_POINT)
				{
#ifdef NR
					assign_roS_nr(t, i, j, k);
#else
					assign_roS(t, i, j, k);
#endif
#ifdef ENERGY
					assign_E_new(i, j, k);
#endif
				}
}

void solve_nonlinear_system()
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (ACTIVE_POINT)
				{
					Newton(i, j, k);
				}
}

void boundary_conditions()
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (ACTIVE_POINT)
				{
					Border_S(i, j, k);
					Border_P(i, j, k);
#ifdef ENERGY
					Border_T(i, j, k);
#endif
				}
}

// Вычисление координаты точки, через которую будет вычисляться значение на границе (i1, j1, k1)
int set_boundary_basic_coordinate(int i, int j, int k)
{
	int i1, j1, k1;

	i1 = i;
	j1 = j;
	k1 = k;

	if ((i == 0) && ((def.locNx) > 2))
	{
		i1 ++;
	}
	if ((i == (def.locNx) - 1) && ((def.locNx) > 2))
	{
		i1 --;
	}
	if (j == 0)
	{
		j1 ++;
	}
	if (j == (def.locNy) - 1)
	{
		j1 --;
	}
	if ((k == 0) && ((def.locNz) > 2))
	{
		k1 ++;
	}
	if ((k == (def.locNz) - 1) && ((def.locNz) > 2))
	{
		k1 --;
	}

	return (i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy));
}

void assign_ro(int local)
{
#ifdef ENERGY
	// !!! Вынести коэффициенты теплового расширения в const consts &def и использовать T_0 оттуда же
	double alfa_w = 1.32E-7; // 1/K !!! E-4
	double alfa_n = 9.2E-7;
	double T_0 = 273;

	HostArraysPtr.ro_w[local] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm) - alfa_w * (HostArraysPtr.T[local] - T_0));
	HostArraysPtr.ro_n[local] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[local] - def.P_atm) - alfa_n * (HostArraysPtr.T[local] - T_0));
	HostArraysPtr.ro_g[local] = def.ro0_g * (HostArraysPtr.P_g[local] / def.P_atm) * (T_0 / HostArraysPtr.T[local]);
#else
	HostArraysPtr.ro_w[local] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm));
	HostArraysPtr.ro_n[local] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[local] - def.P_atm));
	HostArraysPtr.ro_g[local] = def.ro0_g * HostArraysPtr.P_g[local] / def.P_atm;
#endif
	test_ro(HostArraysPtr.ro_g[local], __FILE__, __LINE__);
	test_ro(HostArraysPtr.ro_w[local], __FILE__, __LINE__);
	test_ro(HostArraysPtr.ro_n[local], __FILE__, __LINE__);
}

void assign_S(int local)
{
	HostArraysPtr.S_g[local] = 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local];
	test_S(HostArraysPtr.S_g[local], __FILE__, __LINE__);
}

// Расчет центральной разности
double central_difference (double* ptr, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (*(ptr+1) - *(ptr-1) )/ (2. * (def.hx));	
		}
	case 'y':
		{
			return (*(ptr+def.locNx) - *(ptr-def.locNx) )/ (2. * (def.hy));
		}
	case 'z':
		{
			return (*(ptr + def.locNx * (def.locNy)) - *(ptr - def.locNx * (def.locNy)) )/ (2. * (def.hz));
		}
	default:
		{
			print_error("Axis of [central_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет центральной разности для произведения двух элементов структуры
double multi_central_difference (double* ptr1, double* ptr2, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return ((*(ptr1+1)) * (*(ptr2+1)) - (*(ptr1-1)) * (*(ptr2-1)) )/ (2. * (def.hx));	
		}
	case 'y':
		{
			return ((*(ptr1+def.locNx)) * (*(ptr2+def.locNx)) - (*(ptr1-def.locNx)) * (*(ptr2-def.locNx)) )/ (2. * (def.hy));
		}
	case 'z':
		{
			return ((*(ptr1+def.locNx * (def.locNy))) * (*(ptr2+def.locNx * (def.locNy)))
				- (*(ptr1-def.locNx * (def.locNy))) * (*(ptr2-def.locNx * (def.locNy))) )/ (2. * (def.hz));
		}
	default:
		{
			print_error("Axis of [central_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет направленной разности
double directed_difference (double* P, double* Xi, double* ro, char axis)
{
	double x1 = 0, x2 = 0;
	switch (axis)
	{
	case 'x':
		{
			x2 = right_difference (P, 'x');
			x1 = left_difference (P, 'x');
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (*Xi) * (*ro) -
		      (x1 + fabs(x1)) / 2. * (*(Xi-1)) * (*(ro-1)) +
		      (x2 - fabs(x2)) / 2. * (*(Xi+1)) * (*(ro+1))) / def.hx;
		}
	case 'y':
		{
			x2 = right_difference (P, 'y') + def.g_const * (*ro);
			x1 = left_difference (P, 'y') + def.g_const * (*ro);
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (*Xi) * (*ro) -
		      (x1 + fabs(x1)) / 2. * (*(Xi-def.locNx)) * (*(ro-def.locNx)) +
		      (x2 - fabs(x2)) / 2. * (*(Xi+def.locNx)) * (*(ro+def.locNx))) / def.hy;
		}
	case 'z':
		{
			x2 = right_difference (P, 'z');
			x1 = left_difference (P, 'z');
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (*Xi) * (*ro) -
		      (x1 + fabs(x1)) / 2. * (*(Xi-def.locNx * (def.locNy))) * (*(ro-def.locNx * (def.locNy))) +
		      (x2 - fabs(x2)) / 2. * (*(Xi+def.locNx * (def.locNy))) * (*(ro+def.locNx * (def.locNy)))) / def.hz;
		}
	default:
		{
			print_error("Axis of [directed_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет левой разности
double left_difference (double* ptr, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (*ptr - *(ptr-1) )/ def.hx;
		}
	case 'y':
		{
			return (*ptr - *(ptr-def.locNx) )/ def.hy;
		}
	case 'z':
		{
			return (*ptr - *(ptr - def.locNx * (def.locNy)) )/ def.hz;
		}
	default:
		{
			print_error("Axis of [left_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет правой разности
double right_difference (double* ptr, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (*(ptr+1) - *ptr )/ def.hx;
		}
	case 'y':
		{
			return (*(ptr+def.locNx) - *ptr )/ def.hy;
		}
	case 'z':
		{
			return (*(ptr + def.locNx * (def.locNy)) - *ptr )/ def.hz;
		}
	default:
		{
			print_error("Axis of [right_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет divgrad для произведения двух элементов структуры
double multi_divgrad (double* ptr1, double* ptr2, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return ((*(ptr1+1)) * (*(ptr2+1)) - 2 * (*ptr1) * (*ptr2)
				+ (*(ptr1-1)) * (*(ptr2-1))) / ((def.hx) * (def.hx));
		}
	case 'y':
		{
			return ((*(ptr1+def.locNx)) * (*(ptr2+def.locNx)) - 2 * (*ptr1) * (*ptr2)
				+ (*(ptr1-def.locNx)) * (*(ptr2-def.locNx))) / ((def.hy) * (def.hy));
		}
	case 'z':
		{
			return ((*(ptr1+def.locNx * (def.locNy))) * (*(ptr2+def.locNx * (def.locNy))) - 2 * (*ptr1) * (*ptr2)
				+ (*(ptr1-def.locNx * (def.locNy))) * (*(ptr2-def.locNx * (def.locNy)))) / ((def.hz) * (def.hz));
		}
	default:
		{
			print_error("Axis of [right_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}


// Расчет скоростей в точке
void assign_u(int i, int j, int k)
{
	int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
	if ((def.Nx) > 2)
	{
		if (i == 0)
		{
			HostArraysPtr.ux_w[local] = HostArraysPtr.Xi_w[local] * right_difference(HostArraysPtr.P_w+local, 'x');
			HostArraysPtr.ux_n[local] = HostArraysPtr.Xi_n[local] * right_difference(HostArraysPtr.P_n+local, 'x');
			HostArraysPtr.ux_g[local] = HostArraysPtr.Xi_g[local] * right_difference(HostArraysPtr.P_g+local, 'x');
		}
		else
		{
			if (i == (def.locNx) - 1)
			{
				HostArraysPtr.ux_w[local] = HostArraysPtr.Xi_w[local] * left_difference(HostArraysPtr.P_w+local, 'x');
				HostArraysPtr.ux_n[local] = HostArraysPtr.Xi_n[local] * left_difference(HostArraysPtr.P_n+local, 'x');
				HostArraysPtr.ux_g[local] = HostArraysPtr.Xi_g[local] * left_difference(HostArraysPtr.P_g+local, 'x');
			}
			else
			{
				HostArraysPtr.ux_w[local] = HostArraysPtr.Xi_w[local] * central_difference (HostArraysPtr.P_w+local, 'x');
				HostArraysPtr.ux_n[local] = HostArraysPtr.Xi_n[local] * central_difference (HostArraysPtr.P_n+local, 'x');
				HostArraysPtr.ux_g[local] = HostArraysPtr.Xi_g[local] * central_difference (HostArraysPtr.P_g+local, 'x');
			}
		}
	}
	else
	{
		HostArraysPtr.ux_w[local] = 0.;
		HostArraysPtr.ux_n[local] = 0.;
		HostArraysPtr.ux_g[local] = 0.;
	}

	if ((def.Ny) > 2)
	{
		if (j == 0)
		{
			HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (right_difference (HostArraysPtr.P_w+local, 'y') - HostArraysPtr.ro_w[local] * (def.g_const));
			HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (right_difference (HostArraysPtr.P_n+local, 'y') - HostArraysPtr.ro_n[local] * (def.g_const));
			HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (right_difference (HostArraysPtr.P_g+local, 'y') - HostArraysPtr.ro_g[local] * (def.g_const));
		}
		else
		{
			if (j == (def.locNy) - 1)
			{
				HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (left_difference (HostArraysPtr.P_w+local, 'y') - HostArraysPtr.ro_w[local] * (def.g_const));
				HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (left_difference (HostArraysPtr.P_n+local, 'y') - HostArraysPtr.ro_n[local] * (def.g_const));
				HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (left_difference (HostArraysPtr.P_g+local, 'y') - HostArraysPtr.ro_g[local] * (def.g_const));
			}
			else
			{
				HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (central_difference (HostArraysPtr.P_w+local, 'y')	- HostArraysPtr.ro_w[local] * (def.g_const));
				HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (central_difference (HostArraysPtr.P_n+local, 'y')	- HostArraysPtr.ro_n[local] * (def.g_const));
				HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (central_difference (HostArraysPtr.P_g+local, 'y')	- HostArraysPtr.ro_g[local] * (def.g_const));
			}
		}
	}
	else
	{
		HostArraysPtr.uy_w[local] = 0.;
		HostArraysPtr.uy_n[local] = 0.;
		HostArraysPtr.uy_g[local] = 0.;
	}

	if ((def.Nz) > 2)
	{
		if (k == 0)
		{
			HostArraysPtr.uz_w[local] = HostArraysPtr.Xi_w[local] * right_difference (HostArraysPtr.P_w+local, 'z');
			HostArraysPtr.uz_n[local] = HostArraysPtr.Xi_n[local] * right_difference (HostArraysPtr.P_n+local, 'z');
			HostArraysPtr.uz_g[local] = HostArraysPtr.Xi_g[local] * right_difference (HostArraysPtr.P_g+local, 'z');
		}
		else
		{
			if (k == (def.locNz) - 1)
			{
				HostArraysPtr.uz_w[local] = HostArraysPtr.Xi_w[local] * left_difference (HostArraysPtr.P_w+local, 'z');
				HostArraysPtr.uz_n[local] = HostArraysPtr.Xi_n[local] * left_difference (HostArraysPtr.P_n+local, 'z');
				HostArraysPtr.uz_g[local] = HostArraysPtr.Xi_g[local] * left_difference (HostArraysPtr.P_g+local, 'z');
			}
			else
			{
				HostArraysPtr.uz_w[local] = HostArraysPtr.Xi_w[local] * central_difference (HostArraysPtr.P_w+local, 'z');
				HostArraysPtr.uz_n[local] = HostArraysPtr.Xi_n[local] * central_difference (HostArraysPtr.P_n+local, 'z');
				HostArraysPtr.uz_g[local] = HostArraysPtr.Xi_g[local] * central_difference (HostArraysPtr.P_g+local, 'z');
			}
		}
	}
	else
	{
		HostArraysPtr.uz_w[local] = 0.;
		HostArraysPtr.uz_n[local] = 0.;
		HostArraysPtr.uz_g[local] = 0.;
	}

	test_u(HostArraysPtr.ux_w[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.ux_n[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uy_w[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uy_n[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uz_w[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uz_n[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.ux_g[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uy_g[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uz_g[local], __FILE__, __LINE__);
}

void assign_roS(double t, int i, int j, int k)
{
	if (INTERNAL_POINT)
	{
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
		double divgrad1, divgrad2, Tx1, Ty1, Tz1, Tx2, Ty2, Tz2, A1 = 0, A2 = 0;

		double divgrad3, Tx3, Ty3, Tz3, A3 = 0;

		HostArraysPtr.roS_w[local] = HostArraysPtr.ro_w[local] * HostArraysPtr.S_w[local];
		HostArraysPtr.roS_n[local] = HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local];
		HostArraysPtr.roS_g[local] = HostArraysPtr.ro_g[local] * HostArraysPtr.S_g[local];

		if ((def.Nz) < 2)
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
			divgrad1 = multi_divgrad (HostArraysPtr.ro_w + local, HostArraysPtr.S_w + local, 'z');
			divgrad2 = multi_divgrad (HostArraysPtr.ro_n + local, HostArraysPtr.S_n + local, 'z');
			divgrad3 = multi_divgrad (HostArraysPtr.ro_g + local, HostArraysPtr.S_g + local, 'z');

			Tz1 = multi_central_difference (HostArraysPtr.ro_w + local, HostArraysPtr.uz_w + local, 'z');
			Tz2 = multi_central_difference (HostArraysPtr.ro_n + local, HostArraysPtr.uz_n + local, 'z');
			Tz3 = multi_central_difference (HostArraysPtr.ro_g + local, HostArraysPtr.uz_g + local, 'z');
		}

		if ((def.Nx) < 2)
		{
			Tx1 = 0.;
			Tx2 = 0.;
			Tx3 = 0.;
			divgrad3 = 0.;
		}
		else
		{
			divgrad1 += multi_divgrad (HostArraysPtr.ro_w + local, HostArraysPtr.S_w + local, 'x');
			divgrad2 += multi_divgrad (HostArraysPtr.ro_n + local, HostArraysPtr.S_n + local, 'x');
			divgrad3 += multi_divgrad (HostArraysPtr.ro_g + local, HostArraysPtr.S_g + local, 'x');

			Tx1 = multi_central_difference (HostArraysPtr.ro_w + local, HostArraysPtr.uz_w + local, 'x');
			Tx2 = multi_central_difference (HostArraysPtr.ro_n + local, HostArraysPtr.uz_n + local, 'x');
			Tx3 = multi_central_difference (HostArraysPtr.ro_g + local, HostArraysPtr.uz_g + local, 'x');
		}

		divgrad1 += multi_divgrad (HostArraysPtr.ro_w + local, HostArraysPtr.S_w + local, 'y');
		divgrad1 *= HostArraysPtr.m[local] * (def.l) * (def.c_w);

		divgrad2 += multi_divgrad (HostArraysPtr.ro_n + local, HostArraysPtr.S_n + local, 'y');
		divgrad2 *= HostArraysPtr.m[local] * (def.l) * (def.c_n);

		divgrad3 += multi_divgrad (HostArraysPtr.ro_g + local, HostArraysPtr.S_g + local, 'y');
		divgrad3 *= HostArraysPtr.m[local] * (def.l) * (def.c_g);

		Ty1 = multi_central_difference (HostArraysPtr.ro_w + local, HostArraysPtr.uz_w + local, 'y');
		Ty2 = multi_central_difference (HostArraysPtr.ro_n + local, HostArraysPtr.uz_n + local, 'y');
		Ty3 = multi_central_difference (HostArraysPtr.ro_g + local, HostArraysPtr.uz_g + local, 'y');

		test_arrowhead(Tx1 + Ty1 + Tz1, divgrad1, __FILE__, __LINE__);
		test_arrowhead(Tx2 + Ty2 + Tz2, divgrad2, __FILE__, __LINE__);
		test_arrowhead(Tx3 + Ty3 + Tz3, divgrad3, __FILE__, __LINE__);

		double q_w = 0.;
		double q_n = 0.;
		double q_g = 0.;

		// Значения q на скважинах
		wells_q(i, j, k, &q_w, &q_n, &q_g);

		if ((t < 2 * (def.dt)) || TWO_LAYERS)
		{
			A1 = HostArraysPtr.roS_w[local] + ((def.dt) / HostArraysPtr.m[local]) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = HostArraysPtr.roS_n[local] + ((def.dt) / HostArraysPtr.m[local]) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2);
			A3 = HostArraysPtr.roS_g[local] + ((def.dt) / HostArraysPtr.m[local]) * (q_g + divgrad3 - Tx3 - Ty3 - Tz3);
		}
		else
		{
			A1 = (1. / ((HostArraysPtr.m[local]) * (def.dt) + 2. * (def.tau))) * (2. * (def.dt) * (def.dt) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1)
			        + ((HostArraysPtr.m[local]) * (def.dt) - 2. * (def.tau)) * HostArraysPtr.roS_w_old[local]
			        + 4. * (def.tau) * HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
			A2 = (1. / ((HostArraysPtr.m[local]) * (def.dt) + 2. * (def.tau))) * (2. * (def.dt) * (def.dt) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2)
			        + ((HostArraysPtr.m[local]) * (def.dt) - 2. * (def.tau)) * HostArraysPtr.roS_n_old[local]
			        + 4. * (def.tau) * HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);

			A3 = (1. / ((HostArraysPtr.m[local]) * (def.dt) + 2. * (def.tau))) * (2. * (def.dt) * (def.dt) * (q_g + divgrad3 - Tx3 - Ty3 - Tz3)
			        + ((HostArraysPtr.m[local]) * (def.dt) - 2. * (def.tau)) * HostArraysPtr.roS_g_old[local]
			        + 4. * (def.tau) * HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
		}

		HostArraysPtr.roS_w_old[local] = HostArraysPtr.roS_w[local];
		HostArraysPtr.roS_n_old[local] = HostArraysPtr.roS_n[local];
		HostArraysPtr.roS_g_old[local] = HostArraysPtr.roS_g[local];
		HostArraysPtr.roS_w[local] = A1;
		HostArraysPtr.roS_n[local] = A2;
		HostArraysPtr.roS_g[local] = A3;

		test_positive(HostArraysPtr.roS_w[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_n[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_g[local], __FILE__, __LINE__);
	}
}


void assign_roS_nr(double t, int i, int j, int k)
{
	if (INTERNAL_POINT)
	{
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if(! HostArraysPtr.m[local])
			return;

		double q_w = 0., q_n = 0., q_g = 0;

		// Значения q на скважинах
		wells_q(i, j, k, &q_w, &q_n, &q_g);

		HostArraysPtr.roS_w[local] = HostArraysPtr.ro_w[local] * HostArraysPtr.S_w[local];
		HostArraysPtr.roS_g[local] = HostArraysPtr.ro_g[local]
		        * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local]);
		HostArraysPtr.roS_n[local] = HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local];

		double f_w = 0., f_n = 0., f_g = 0., A1 = 0., A2 = 0., A3 = 0.;

		if ((def.Nx) > 2)
		{
			f_w += directed_difference (HostArraysPtr.P_w+local, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, 'x');
			f_n += directed_difference (HostArraysPtr.P_n+local, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, 'x');
			f_g += directed_difference (HostArraysPtr.P_g+local, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, 'x');
		}
		if ((def.Ny) > 2)
		{
			f_w += directed_difference (HostArraysPtr.P_w+local, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, 'y');
			f_n += directed_difference (HostArraysPtr.P_n+local, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, 'y');
			f_g += directed_difference (HostArraysPtr.P_g+local, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, 'y');
		}
		if ((def.Nz) > 2)
		{
			f_w += directed_difference (HostArraysPtr.P_w+local, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, 'z');
			f_n += directed_difference (HostArraysPtr.P_n+local, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, 'z');
			f_g += directed_difference (HostArraysPtr.P_g+local, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, 'z');
		}

		A1 = HostArraysPtr.roS_w[local] + (def.dt / HostArraysPtr.m[local]) * (q_w - f_w);
		A2 = HostArraysPtr.roS_n[local] + (def.dt / HostArraysPtr.m[local]) * (q_n - f_n);
		A3 = HostArraysPtr.roS_g[local] + (def.dt / HostArraysPtr.m[local]) * (q_g - f_g);

		HostArraysPtr.roS_w_old[local] = HostArraysPtr.roS_w[local];
		HostArraysPtr.roS_n_old[local] = HostArraysPtr.roS_n[local];
		HostArraysPtr.roS_g_old[local] = HostArraysPtr.roS_g[local];
		HostArraysPtr.roS_w[local] = A1;
		HostArraysPtr.roS_n[local] = A2;
		HostArraysPtr.roS_g[local] = A3;

		test_positive(HostArraysPtr.roS_w[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_n[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_g[local], __FILE__, __LINE__);
	}
}

// Функция загрузки данных в память хоста
void load_data_to_host(double *HostArray, double *DevArray)
{
}

// Функция загрузки данных типа double в память ускорителя
void load_data_to_device(double *HostArray, double *DevArray)
{
}

// Функция загрузки данных типа int в память ускорителя
void load_data_to_device_int(int *HostArray, int *DevArray)
{
}

// Выделение памяти ускорителя под массив точек расчетной области
void device_memory_allocation()
{
}

// Освобожение памяти ускорителя из под массива точек расчетной области
void device_memory_free()
{
}

// Инициализация ускорителя
void device_initialization()
{
}

// Финализация ускорителя
void device_finalization(void)
{
}

// Загрузка в буфер данных для обмена на границе. Для каждого из направлений своя функция. Направление - это ось координат и лево/право.
void load_exchange_data_part_xl(double *HostArray)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[j + (def.locNy)*k] = HostArray[1 + (def.locNx) * j + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[j + (def.locNy)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_xr(double *HostArray)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[j + (def.locNy)*k] = HostArray[(def.locNx) - 2 + (def.locNx) * j + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[j + (def.locNy)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_yl(double *HostArray)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[i + (def.locNx)*k] = HostArray[i + (def.locNx) + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[i + (def.locNx)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_yr(double *HostArray)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[i + (def.locNx)*k] = HostArray[i + (def.locNx) * ((def.locNy) - 2) + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[i + (def.locNx)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_zl(double *HostArray)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostBuffer[i + (def.locNx)*j] = HostArray[i + (def.locNx) * j + (def.locNx) * (def.locNy)];
			test_nan(HostBuffer[i + (def.locNx)*j], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_zr(double *HostArray)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostBuffer[i + (def.locNx)*j] = HostArray[i + (def.locNx) * j + (def.locNx) * (def.locNy) * ((def.locNz) - 2)];
			test_nan(HostBuffer[i + (def.locNx)*j], __FILE__, __LINE__);
		}
}

// Загрузка из буфера данных обмена на границе. Для каждого из направлений своя функция. Направление - это ось координат и лево/право.
void save_exchange_data_part_xl(double *HostArray)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArray[(def.locNx)*j + (def.locNx) * (def.locNy)*k] = HostBuffer[j + (def.locNy) * k];
			test_nan(HostArray[(def.locNx)*j + (def.locNx) * (def.locNy)*k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_xr(double *HostArray)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArray[(def.locNx) - 1 + (def.locNx)*j + (def.locNx) * (def.locNy)*k] = HostBuffer[j + (def.locNy) * k];
			test_nan(HostArray[(def.locNx) - 1 + (def.locNx)*j + (def.locNx) * (def.locNy)*k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_yl(double *HostArray)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArray[i + (def.locNx) * (def.locNy) * k] = HostBuffer[i + (def.locNx)*k];
			test_nan(HostArray[i + (def.locNx) * (def.locNy) * k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_yr(double *HostArray)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArray[i + (def.locNx) * ((def.locNy) - 1) + (def.locNx) * (def.locNy) * k] = HostBuffer[i + (def.locNx)*k];
			test_nan(HostArray[i + (def.locNx) * ((def.locNy) - 1) + (def.locNx) * (def.locNy) * k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_zl(double *HostArray)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostArray[i + (def.locNx) * j] = HostBuffer[i + (def.locNx)*j];
			test_nan(HostArray[i + (def.locNx) * j], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_zr(double *HostArray)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostArray[i + (def.locNx) * j + (def.locNx) * (def.locNy) * ((def.locNz) - 1)] = HostBuffer[i + (def.locNx)*j];
			test_nan(HostArray[i + (def.locNx) * j + (def.locNx) * (def.locNy) * ((def.locNz) - 1)], __FILE__, __LINE__);
		}
}

