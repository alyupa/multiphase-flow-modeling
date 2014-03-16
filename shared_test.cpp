#include "defines.h"

// Тестирование

// Функция, вызываемая при ошибке
void print_error(const char *error, const char *file, int line)
{
	printf("Error: %s\nFile: \"%s\"\nLine: %d\n\n", error, file, line);
	fflush(stdout);
#ifdef _WIN32
	printf("\nPress <Enter> to exit...\n");
	fflush(stdout);
	getchar();
#endif
	exit(1);
}

// Функция проверки на выход из допустимого диапазона значений
// во всех точках расчетной области процессора

void test_correct_P_S()
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
			{
				;// TODO: test values
			}
}

// Тест на NaN
// Синтаксис вызова test_nan(x, __FILE__, __LINE__);
void test_nan(double x, const char *file, int line)
{
#ifdef MY_TEST
	if (isnan(x))
	{
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	// Тестовое более жесткое ограничение именно для этой задачи
	/*if (x > 1e30 || x < -1e30)
	{
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}*/
	
#endif
}

// Тест на положительное и не NaN
// Синтаксис вызова test_positive(x, __FILE__, __LINE__);
void test_positive(double x, const char *file, int line)
{
#ifdef MY_TEST
	if (isnan(x))
	{
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (x < 0)
	{
		printf("Error: x<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// Тест на вхождение насыщенностей в [0;1]
// Синтаксис вызова test_S(x, __FILE__, __LINE__);
void test_S(double S, const char *file, int line)
{
#ifdef MY_TEST
	if (isnan(S))
	{
		printf("Error: S=NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (S < 0)
	{
		printf("Error: S<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (S > 1)
	{
		printf("Error: S>1\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// Тест на вхождение скоростей в [-100;100]
// Синтаксис вызова test_u(x, __FILE__, __LINE__);
void test_u(double u, const char *file, int line)
{
#ifdef MY_TEST
	if (isnan(u))
	{
		printf("Error: u=NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (u < -1e8)
	{
		printf("Error: u<-100\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (u > 1e8)
	{
		printf("Error: u>100\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}


// Тест на вхождение плотностей в [-0;1500]
// Синтаксис вызова test_ro(x, __FILE__, __LINE__);
void test_ro(double ro, const char *file, int line)
{
#ifdef MY_TEST
	if (isnan(ro))
	{
		printf("Error: ro = NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (ro < 0)
	{
		printf("Error: ro < 0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (ro > 1e8)
	{
		printf("Error: ro > 1500\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// Функция проверяет, что первый аргумент много больше (по модулю) второго
// Если это не так, печатается предупреждение
void test_arrowhead(double big, double small, const char *file, int line)
{
#ifdef MY_TEST_1
	if (fabs(big / 30) < fabs(small))
	{
		printf("Warning: See task parameters.\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// Тест на корректность параметров задачи
void read_defines_test()
{
#ifdef MY_TEST
	test_positive(def.hx, __FILE__, __LINE__);
	test_positive(def.hy, __FILE__, __LINE__);
	test_positive(def.hz, __FILE__, __LINE__);
	test_positive(def.tau, __FILE__, __LINE__);
	test_positive(def.dt, __FILE__, __LINE__);
	test_positive(def.c_w, __FILE__, __LINE__);
	test_positive(def.c_n, __FILE__, __LINE__);
	test_positive(def.l, __FILE__, __LINE__);
	test_positive(def.beta_w, __FILE__, __LINE__);
	test_positive(def.beta_n, __FILE__, __LINE__);
	test_positive(def.ro0_w, __FILE__, __LINE__);
	test_positive(def.ro0_n, __FILE__, __LINE__);
	test_positive(def.mu_w, __FILE__, __LINE__);
	test_positive(def.mu_n, __FILE__, __LINE__);
	test_positive(def.g_const, __FILE__, __LINE__);
	test_positive(def.P_atm, __FILE__, __LINE__);

	test_positive(def.source, __FILE__, __LINE__);
	test_positive(def.newton_iterations, __FILE__, __LINE__);
	test_positive(def.timeX, __FILE__, __LINE__);
	test_positive(def.save_plots, __FILE__, __LINE__);
	test_positive(strlen(def.plots_dir), __FILE__, __LINE__);
	test_positive(def.print_screen, __FILE__, __LINE__);
	test_positive(def.Nx, __FILE__, __LINE__);
	test_positive(def.Ny, __FILE__, __LINE__);
	test_positive(def.Nz, __FILE__, __LINE__);

	test_positive(def.c_g, __FILE__, __LINE__);
	test_positive(def.beta_g, __FILE__, __LINE__);
	test_positive(def.ro0_g, __FILE__, __LINE__);
	test_positive(def.mu_g, __FILE__, __LINE__);
	test_S(def.S_w_gr, __FILE__, __LINE__);
	test_S(def.S_n_gr, __FILE__, __LINE__);

#ifdef ENERGY
	test_positive(def.T_0, __FILE__, __LINE__);
	test_positive(def.ro_r, __FILE__, __LINE__);
	test_positive(def.lambda0_w, __FILE__, __LINE__);
	test_positive(def.lambda0_n, __FILE__, __LINE__);
	test_positive(def.lambda0_g, __FILE__, __LINE__);
	test_positive(def.lambda0_r, __FILE__, __LINE__);
	test_positive(def.lambdaA_w, __FILE__, __LINE__);
	test_positive(def.lambdaA_n, __FILE__, __LINE__);
	test_positive(def.lambdaA_g, __FILE__, __LINE__);
	test_positive(def.c0_w, __FILE__, __LINE__);
	test_positive(def.c0_n, __FILE__, __LINE__);
	test_positive(def.c0_g, __FILE__, __LINE__);
	test_positive(def.c0_r, __FILE__, __LINE__);
	test_positive(def.C_w, __FILE__, __LINE__);
	test_positive(def.C_w2, __FILE__, __LINE__);
	test_positive(def.C_n, __FILE__, __LINE__);
	test_positive(def.C_g, __FILE__, __LINE__);
	test_positive(def.C_r, __FILE__, __LINE__);
	test_positive(def.alfa_w, __FILE__, __LINE__);
	test_positive(def.alfa_n, __FILE__, __LINE__);
#endif
#endif

}

