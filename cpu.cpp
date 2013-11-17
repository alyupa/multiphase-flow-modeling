#include "defines.h"

void ro_P_Xi_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
			{
				if (is_active_point(i, j, k, def))
				{
					assign_P_Xi(HostArraysPtr, i, j, k, def);
					assign_ro(HostArraysPtr, i, j, k, def);
				}
			}
}

void u_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (is_active_point(i, j, k, def))
				{
					assign_u(HostArraysPtr, i, j, k, def);
				}
}

void S_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				{
					assign_S(HostArraysPtr, i, j, k, def);
				}
}

#ifdef ENERGY
void H_E_current_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
			{
				int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
				assign_H(HostArraysPtr, local, def);
				assign_E_current(HostArraysPtr, local, def);
			}
}
#endif

void roS_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double t, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (is_active_point(i, j, k, def))
				{
#ifdef NR
					assign_roS_nr(HostArraysPtr, t, i, j, k, def);
#else
					assign_roS(HostArraysPtr, t, i, j, k, def);
#endif
				}
}

#ifdef ENERGY
void E_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (is_active_point(i, j, k, def))
				{
					assign_E_new(HostArraysPtr, i, j, k, def);
				}
}
#endif

void P_S_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (is_active_point(i, j, k, def))
				{
					Newton(HostArraysPtr, i, j, k, def);
				}
}

void boundary_conditions(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (is_active_point(i, j, k, def))
				{
					Border_S(HostArraysPtr, i, j, k, def);
					Border_P(HostArraysPtr, i, j, k, def);
#ifdef ENERGY
					Border_T(HostArraysPtr, i, j, k, def);
#endif
				}
}

// Вычисление координаты точки, через которую будет вычисляться значение на границе (i1, j1, k1)
int set_boundary_basic_coordinate(int i, int j, int k, const consts &def)
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

void assign_ro(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

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

void assign_S(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

	HostArraysPtr.S_g[local] = 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local];
	test_S(HostArraysPtr.S_g[local], __FILE__, __LINE__);
}

// Расчет центральной разности
double central_difference (double* ptr, char axis, const consts &def)
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
double multi_central_difference (double* ptr1, double* ptr2, char axis, const consts &def)
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
double directed_difference (double x1, double x2, double* Xi, double* ro, char axis, const consts &def)
{
	switch (axis)
	{
	case 'x':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
		      (x1 + fabs(x1)) / 2. * (-1.) * (*(Xi-1)) * (*(ro-1)) +
		      (x2 - fabs(x2)) / 2. * (-1.) * (*(Xi+1)) * (*(ro+1))) / def.hx;
		}
	case 'y':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
		      (x1 + fabs(x1)) / 2. * (-1.) * (*(Xi-def.locNx)) * (*(ro-def.locNx)) +
		      (x2 - fabs(x2)) / 2. * (-1.) * (*(Xi+def.locNx)) * (*(ro+def.locNx))) / def.hy;
		}
	case 'z':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
		      (x1 + fabs(x1)) / 2. * (-1.) * (*(Xi-def.locNx * (def.locNy))) * (*(ro-def.locNx * (def.locNy))) +
		      (x2 - fabs(x2)) / 2. * (-1.) * (*(Xi+def.locNx * (def.locNy))) * (*(ro+def.locNx * (def.locNy)))) / def.hz;
		}
	default:
		{
			print_error("Axis of [directed_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// Расчет левой разности
double left_difference (double* ptr, char axis, const consts &def)
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
double right_difference (double* ptr, char axis, const consts &def)
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
double multi_divgrad (double* ptr1, double* ptr2, char axis, const consts &def)
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
void assign_u(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def)
{
	int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
	if ((def.Nx) > 2)
	{
		if (i == 0)
		{
			HostArraysPtr.ux_w[local] = HostArraysPtr.Xi_w[local] * right_difference(HostArraysPtr.P_w+local, 'x', def);
			HostArraysPtr.ux_n[local] = HostArraysPtr.Xi_n[local] * right_difference(HostArraysPtr.P_n+local, 'x', def);
			HostArraysPtr.ux_g[local] = HostArraysPtr.Xi_g[local] * right_difference(HostArraysPtr.P_g+local, 'x', def);
		}
		else
		{
			if (i == (def.locNx) - 1)
			{
				HostArraysPtr.ux_w[local] = HostArraysPtr.Xi_w[local] * left_difference(HostArraysPtr.P_w+local, 'x', def);
				HostArraysPtr.ux_n[local] = HostArraysPtr.Xi_n[local] * left_difference(HostArraysPtr.P_n+local, 'x', def);
				HostArraysPtr.ux_g[local] = HostArraysPtr.Xi_g[local] * left_difference(HostArraysPtr.P_g+local, 'x', def);
			}
			else
			{
				HostArraysPtr.ux_w[local] = HostArraysPtr.Xi_w[local] * central_difference (HostArraysPtr.P_w+local, 'x', def);
				HostArraysPtr.ux_n[local] = HostArraysPtr.Xi_n[local] * central_difference (HostArraysPtr.P_n+local, 'x', def);
				HostArraysPtr.ux_g[local] = HostArraysPtr.Xi_g[local] * central_difference (HostArraysPtr.P_g+local, 'x', def);
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
			HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (right_difference (HostArraysPtr.P_w+local, 'y', def) - HostArraysPtr.ro_w[local] * (def.g_const));
			HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (right_difference (HostArraysPtr.P_n+local, 'y', def) - HostArraysPtr.ro_n[local] * (def.g_const));
			HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (right_difference (HostArraysPtr.P_g+local, 'y', def) - HostArraysPtr.ro_g[local] * (def.g_const));
		}
		else
		{
			if (j == (def.locNy) - 1)
			{
				HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (left_difference (HostArraysPtr.P_w+local, 'y', def) - HostArraysPtr.ro_w[local] * (def.g_const));
				HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (left_difference (HostArraysPtr.P_n+local, 'y', def) - HostArraysPtr.ro_n[local] * (def.g_const));
				HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (left_difference (HostArraysPtr.P_g+local, 'y', def) - HostArraysPtr.ro_g[local] * (def.g_const));
			}
			else
			{
				HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (central_difference (HostArraysPtr.P_w+local, 'y', def)	- HostArraysPtr.ro_w[local] * (def.g_const));
				HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (central_difference (HostArraysPtr.P_n+local, 'y', def)	- HostArraysPtr.ro_n[local] * (def.g_const));
				HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (central_difference (HostArraysPtr.P_g+local, 'y', def)	- HostArraysPtr.ro_g[local] * (def.g_const));
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
			HostArraysPtr.uz_w[local] = HostArraysPtr.Xi_w[local] * right_difference (HostArraysPtr.P_w+local, 'z', def);
			HostArraysPtr.uz_n[local] = HostArraysPtr.Xi_n[local] * right_difference (HostArraysPtr.P_n+local, 'z', def);
			HostArraysPtr.uz_g[local] = HostArraysPtr.Xi_g[local] * right_difference (HostArraysPtr.P_g+local, 'z', def);
		}
		else
		{
			if (k == (def.locNz) - 1)
			{
				HostArraysPtr.uz_w[local] = HostArraysPtr.Xi_w[local] * left_difference (HostArraysPtr.P_w+local, 'z', def);
				HostArraysPtr.uz_n[local] = HostArraysPtr.Xi_n[local] * left_difference (HostArraysPtr.P_n+local, 'z', def);
				HostArraysPtr.uz_g[local] = HostArraysPtr.Xi_g[local] * left_difference (HostArraysPtr.P_g+local, 'z', def);
			}
			else
			{
				HostArraysPtr.uz_w[local] = HostArraysPtr.Xi_w[local] * central_difference (HostArraysPtr.P_w+local, 'z', def);
				HostArraysPtr.uz_n[local] = HostArraysPtr.Xi_n[local] * central_difference (HostArraysPtr.P_n+local, 'z', def);
				HostArraysPtr.uz_g[local] = HostArraysPtr.Xi_g[local] * central_difference (HostArraysPtr.P_g+local, 'z', def);
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

void assign_roS(const ptr_Arrays &HostArraysPtr, double t, int i, int j, int k, const consts &def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.Nz) < 2)))
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
			divgrad1 = multi_divgrad (HostArraysPtr.ro_w + local, HostArraysPtr.S_w + local, 'z', def);
			divgrad2 = multi_divgrad (HostArraysPtr.ro_n + local, HostArraysPtr.S_n + local, 'z', def);
			divgrad3 = multi_divgrad (HostArraysPtr.ro_g + local, HostArraysPtr.S_g + local, 'z', def);

			Tz1 = multi_central_difference (HostArraysPtr.ro_w + local, HostArraysPtr.uz_w + local, 'z', def);
			Tz2 = multi_central_difference (HostArraysPtr.ro_n + local, HostArraysPtr.uz_n + local, 'z', def);
			Tz3 = multi_central_difference (HostArraysPtr.ro_g + local, HostArraysPtr.uz_g + local, 'z', def);
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
			divgrad1 += multi_divgrad (HostArraysPtr.ro_w + local, HostArraysPtr.S_w + local, 'x', def);
			divgrad2 += multi_divgrad (HostArraysPtr.ro_n + local, HostArraysPtr.S_n + local, 'x', def);
			divgrad3 += multi_divgrad (HostArraysPtr.ro_g + local, HostArraysPtr.S_g + local, 'x', def);

			Tx1 = multi_central_difference (HostArraysPtr.ro_w + local, HostArraysPtr.uz_w + local, 'x', def);
			Tx2 = multi_central_difference (HostArraysPtr.ro_n + local, HostArraysPtr.uz_n + local, 'x', def);
			Tx3 = multi_central_difference (HostArraysPtr.ro_g + local, HostArraysPtr.uz_g + local, 'x', def);
		}

		divgrad1 += multi_divgrad (HostArraysPtr.ro_w + local, HostArraysPtr.S_w + local, 'y', def);
		divgrad1 *= HostArraysPtr.m[local] * (def.l) * (def.c_w);

		divgrad2 += multi_divgrad (HostArraysPtr.ro_n + local, HostArraysPtr.S_n + local, 'y', def);
		divgrad2 *= HostArraysPtr.m[local] * (def.l) * (def.c_n);

		divgrad3 += multi_divgrad (HostArraysPtr.ro_g + local, HostArraysPtr.S_g + local, 'y', def);
		divgrad3 *= HostArraysPtr.m[local] * (def.l) * (def.c_g);

		Ty1 = multi_central_difference (HostArraysPtr.ro_w + local, HostArraysPtr.uz_w + local, 'y', def);		
		Ty2 = multi_central_difference (HostArraysPtr.ro_n + local, HostArraysPtr.uz_n + local, 'y', def);
		Ty3 = multi_central_difference (HostArraysPtr.ro_g + local, HostArraysPtr.uz_g + local, 'y', def);

		test_arrowhead(Tx1 + Ty1 + Tz1, divgrad1, __FILE__, __LINE__);
		test_arrowhead(Tx2 + Ty2 + Tz2, divgrad2, __FILE__, __LINE__);
		test_arrowhead(Tx3 + Ty3 + Tz3, divgrad3, __FILE__, __LINE__);

		double q_w = 0.;
		double q_n = 0.;
		double q_g = 0.;

		// Значения q на скважинах
		wells_q(HostArraysPtr, i, j, k, &q_w, &q_n, &q_g, def);

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

			test_tau(HostArraysPtr.roS_w_old[local], HostArraysPtr.roS_w[local], A1, local, def, __FILE__, __LINE__);
			test_tau(HostArraysPtr.roS_n_old[local], HostArraysPtr.roS_n[local], A2, local, def, __FILE__, __LINE__);
			test_tau(HostArraysPtr.roS_g_old[local], HostArraysPtr.roS_g[local], A3, local, def, __FILE__, __LINE__);
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


void assign_roS_nr(const ptr_Arrays &HostArraysPtr, double t, int i, int j, int k, const consts &def)
{
	if ((((i != 0) && (i != (def.locNx) - 1)) || ((def.locNx) < 2)) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if(! HostArraysPtr.m[local])
			return;

		double q_w = 0.;
		double q_n = 0.;
		double q_g = 0.;

		// Значения q на скважинах
		wells_q(HostArraysPtr, i, j, k, &q_w, &q_n, &q_g, def);

		HostArraysPtr.roS_w[local] = HostArraysPtr.ro_w[local] * HostArraysPtr.S_w[local];
		HostArraysPtr.roS_g[local] = HostArraysPtr.ro_g[local]
		        * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local]);

		double Pg = HostArraysPtr.P_g[local];
		double fx_g, fy_g, fz_g, A3 = 0.;

		HostArraysPtr.roS_n[local] = HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local];
		double Pw = HostArraysPtr.P_w[local];
		double Pn = HostArraysPtr.P_n[local];

		double x1, x2, y1, y2, z1, z2, fx_w, fy_w, fz_w, fx_n, fy_n, fz_n, A1 = 0., A2 = 0.;

		if ((def.Nz) < 2)
		{
			fz_w = 0.;
			fz_n = 0.;
			fz_g = 0.;
		}
		else
		{
			z2 = -1. * right_difference (HostArraysPtr.P_w+local, 'z', def);
			z1 = -1. * left_difference (HostArraysPtr.P_w+local, 'z', def);
			fz_w = directed_difference (z1, z2, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, 'z', def);

			z2 = -1. * right_difference (HostArraysPtr.P_n+local, 'z', def); 
			z1 = -1. * left_difference (HostArraysPtr.P_n+local, 'z', def); 
			fz_n = directed_difference (z1, z2, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, 'z', def);

			z2 = -1. * right_difference (HostArraysPtr.P_g+local, 'z', def); 
			z1 = -1. * left_difference (HostArraysPtr.P_g+local, 'z', def); 
			fz_g = directed_difference (z1, z2, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, 'z', def);

		}

		if ((def.Nx) < 2)
		{
			fx_w = 0.;
			fx_n = 0.;
			fx_g = 0.;
		}
		else
		{
			x2 = -1. * right_difference (HostArraysPtr.P_w+local, 'x', def); //-(HostArraysPtr.P_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - Pw) / def.hx;
			x1 = -1. * left_difference (HostArraysPtr.P_w+local, 'x', def); //-(Pw - HostArraysPtr.P_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hx;
			fx_w = directed_difference (x1, x2, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, 'x', def);

			x2 = -1. * right_difference (HostArraysPtr.P_n+local, 'x', def); 
			x1 = -1. * left_difference (HostArraysPtr.P_n+local, 'x', def);
			fx_n = directed_difference (x1, x2, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, 'x', def);

			x2 = -1. * right_difference (HostArraysPtr.P_g+local, 'x', def); 
			x1 = -1. * left_difference (HostArraysPtr.P_g+local, 'x', def); 
			fx_g = directed_difference (x1, x2, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, 'x', def);
		}

		y2 = -1. * right_difference (HostArraysPtr.P_w+local, 'y', def) + def.g_const * (HostArraysPtr.ro_w[local]);
		y1 = -1. * left_difference (HostArraysPtr.P_w+local, 'y', def) + def.g_const * (HostArraysPtr.ro_w[local]);
		fy_w = directed_difference (y1, y2, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, 'y', def);
 
		y2 = -1. * right_difference (HostArraysPtr.P_n+local, 'y', def) + def.g_const * (HostArraysPtr.ro_n[local]);
		y1 = -1. * left_difference (HostArraysPtr.P_n+local, 'y', def) + def.g_const * (HostArraysPtr.ro_n[local]);
		fy_n = directed_difference (y1, y2, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, 'y', def);

		A1 = HostArraysPtr.roS_w[local] - (def.dt / HostArraysPtr.m[local]) * (-q_w + fx_w + fy_w + fz_w);
		A2 = HostArraysPtr.roS_n[local] - (def.dt / HostArraysPtr.m[local]) * (-q_n + fx_n + fy_n + fz_n);

		HostArraysPtr.roS_w_old[local] = HostArraysPtr.roS_w[local];
		HostArraysPtr.roS_n_old[local] = HostArraysPtr.roS_n[local];
		HostArraysPtr.roS_w[local] = A1;
		HostArraysPtr.roS_n[local] = A2;

		test_positive(HostArraysPtr.roS_w[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_n[local], __FILE__, __LINE__);

		y2 = -1. * right_difference (HostArraysPtr.P_g+local, 'y', def) + def.g_const * (HostArraysPtr.ro_g[local]);
		y1 = -1. * left_difference (HostArraysPtr.P_g+local, 'y', def) + def.g_const * (HostArraysPtr.ro_g[local]);
		fy_g = directed_difference (y1, y2, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, 'y', def);

		A3 = HostArraysPtr.roS_g[local] - (def.dt / HostArraysPtr.m[local]) * (q_g + fx_g + fy_g + fz_g);

		HostArraysPtr.roS_g_old[local] = HostArraysPtr.roS_g[local];
		HostArraysPtr.roS_g[local] = A3;

		test_positive(HostArraysPtr.roS_g[local], __FILE__, __LINE__);
	}
}

// Функция загрузки данных в память хоста
void load_data_to_host(double *HostArrayPtr, double *DevArrayPtr, const consts &def)
{
}

// Функция загрузки данных типа double в память ускорителя
void load_data_to_device(double *HostArrayPtr, double *DevArrayPtr, const consts &def)
{
}

// Функция загрузки данных типа int в память ускорителя
void load_data_to_device_int(int *HostArrayPtr, int *DevArrayPtr, const consts &def)
{
}

// Выделение памяти ускорителя под массив точек расчетной области
void device_memory_allocation(ptr_Arrays *ArraysPtr, double **DevBuffer, const consts &def)
{
}

// Освобожение памяти ускорителя из под массива точек расчетной области
void device_memory_free(ptr_Arrays ptDev, double *DevBuffer)
{
}

// Инициализация ускорителя
void device_initialization(consts *def)
{
}

// Финализация ускорителя
void device_finalization(void)
{
}

// Загрузка в буфер данных для обмена на границе. Для каждого из направлений своя функция. Направление - это ось координат и лево/право.
void load_exchange_data_part_xl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[j + (def.locNy)*k] = HostArrayPtr[1 + (def.locNx) * j + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[j + (def.locNy)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_xr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[j + (def.locNy)*k] = HostArrayPtr[(def.locNx) - 2 + (def.locNx) * j + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[j + (def.locNy)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_yl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[i + (def.locNx)*k] = HostArrayPtr[i + (def.locNx) + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[i + (def.locNx)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_yr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[i + (def.locNx)*k] = HostArrayPtr[i + (def.locNx) * ((def.locNy) - 2) + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[i + (def.locNx)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_zl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostBuffer[i + (def.locNx)*j] = HostArrayPtr[i + (def.locNx) * j + (def.locNx) * (def.locNy)];
			test_nan(HostBuffer[i + (def.locNx)*j], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_zr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostBuffer[i + (def.locNx)*j] = HostArrayPtr[i + (def.locNx) * j + (def.locNx) * (def.locNy) * ((def.locNz) - 2)];
			test_nan(HostBuffer[i + (def.locNx)*j], __FILE__, __LINE__);
		}
}

// Загрузка из буфера данных обмена на границе. Для каждого из направлений своя функция. Направление - это ось координат и лево/право.
void save_exchange_data_part_xl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArrayPtr[(def.locNx)*j + (def.locNx) * (def.locNy)*k] = HostBuffer[j + (def.locNy) * k];
			test_nan(HostArrayPtr[(def.locNx)*j + (def.locNx) * (def.locNy)*k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_xr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArrayPtr[(def.locNx) - 1 + (def.locNx)*j + (def.locNx) * (def.locNy)*k] = HostBuffer[j + (def.locNy) * k];
			test_nan(HostArrayPtr[(def.locNx) - 1 + (def.locNx)*j + (def.locNx) * (def.locNy)*k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_yl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArrayPtr[i + (def.locNx) * (def.locNy) * k] = HostBuffer[i + (def.locNx)*k];
			test_nan(HostArrayPtr[i + (def.locNx) * (def.locNy) * k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_yr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArrayPtr[i + (def.locNx) * ((def.locNy) - 1) + (def.locNx) * (def.locNy) * k] = HostBuffer[i + (def.locNx)*k];
			test_nan(HostArrayPtr[i + (def.locNx) * ((def.locNy) - 1) + (def.locNx) * (def.locNy) * k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_zl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostArrayPtr[i + (def.locNx) * j] = HostBuffer[i + (def.locNx)*j];
			test_nan(HostArrayPtr[i + (def.locNx) * j], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_zr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostArrayPtr[i + (def.locNx) * j + (def.locNx) * (def.locNy) * ((def.locNz) - 1)] = HostBuffer[i + (def.locNx)*j];
			test_nan(HostArrayPtr[i + (def.locNx) * j + (def.locNx) * (def.locNy) * ((def.locNz) - 1)], __FILE__, __LINE__);
		}
}

