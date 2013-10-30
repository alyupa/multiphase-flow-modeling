#include "../../defines.h"

// Вспомогательная функции для метода Ньютона:
// Нахождение обратной матрицы 3*3;
void reverse_matrix(double* a)
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

void compare_matrix(double *a, double *b)
{
	double eps = 1e-10;
	int n = 9;

	for(int i = 0; i < n; i++)
		if (fabs(a[i] - b[i]) >= eps)
		{
			printf("Error: reverse_matrix() doesn't work correctly\n");
			return;
		}
	printf("Success\n");
}

void unit_test_reverse_matrix()
{
	int n = 9;
	double a[9], b[9];

	a[0] = 1; a[1] = 0; a[2] = 0; a[3] = 0; a[4] = 1; a[5] = 0; a[6] = 0; a[7] = 0; a[8] = 1; 
	b[0] = 1; b[1] = 0; b[2] = 0; b[3] = 0; b[4] = 1; b[5] = 0; b[6] = 0; b[7] = 0; b[8] = 1; 
	reverse_matrix(a);
	compare_matrix(a, b);

	a[0] = -4./19; a[1] = 7./38; a[2] = 5./38; a[3] = 3./19; a[4] = 9./38; a[5] = 1./38; a[6] = 2./19; a[7] = -13./38; a[8] = 7./38; 
	b[0] = -2; b[1] = 3; b[2] = 1; b[3] = 1; b[4] = 2; b[5] = -1; b[6] = 3; b[7] = 2; b[8] = 3; 
	reverse_matrix(a);
	compare_matrix(a, b);

	a[0] = 3; a[1] = 0; a[2] = 0; a[3] = 2; a[4] = 2; a[5] = 1; a[6] = 1; a[7] = 0; a[8] = 1; 
	b[0] = 1./3; b[1] = 0; b[2] = 0; b[3] = -1./6; b[4] = 0.5; b[5] = -0.5; b[6] = -1./3; b[7] = 0; b[8] = 1; 
	reverse_matrix(a);
	compare_matrix(a, b);
}

int main(int argc, char* argv[])
{
#ifdef GTEST
	::testing::InitGoogleTest(&argc, argv);
	RUN_ALL_TESTS();
#endif
	unit_test_reverse_matrix();

#ifdef _WIN32
	printf("\nPress <Enter> to exit...\n");
	fflush(stdout);
	getchar();
#endif

	return 0;
}

