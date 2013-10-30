//Нахождение обратной матрицы методом расширенной матрицы (методом Гаусса)

#include "defines.h"

//obraschenie matrici
int reverse_matrix (double *a, int n)
{
	int i, j, k;
	double f, f1;
	double *matr, *row;
	double e = 1E-8;
	int det_flag = 1;
	int n1 = n * 2;

	if (a == NULL || n < 2) 
	{
		printf ("Wrong arguments in reverse_matrix!\n");
		return -2;
	}

	matr = new double[n * n1];
	row = new double[n1];

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			matr [j + n1 * i] = a[j + n * i];

	for (i = 0; i < n; i++)
		for (j = n; j < n1; j++)
			if (j == n + i)
				matr[j + n1 * i] = 1;
			else
				matr[j + n1 * i] = 0;

	for (k = 0; k < n; k++)
	{
		f = matr[k + n1 * k];

		if (fabs(f) < e)
		{
			det_flag = 0;		

			for (i = k + 1; i < n; i++)
			{
				if (matr[k + n1 * i] > e)
				{
					det_flag = 1;

					for (j = 0; j < n1; j++)
					{
						row [j] = matr[j + n1 * k];
						matr[j + n1 * k] = matr[j + n1 * i];
						matr[j + n1 * i] = row [j];
					}

					f = matr[k + n1 * k];
				}
			}

			if (det_flag == 0)
			{
				std::cout << "Determinant = 0 in reverse_matrix!\n";	
				return -1;
			}	
		}

		for (j = 0; j < n1; j++)
		{
			matr[j + n1 * k] /= f;
		}

		for (i = 0; i < n; i++)
			if (i != k)
			{
				f1 = matr[k + n1 * i];

				for (j = 0; j < n1; j++)
				{
					matr[j + n1 * i] -= (matr[j + n1 * k] * f1);
				}
			}	
	}     

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			a[j + n * i] = matr[j + n + n1 * i];

	delete[] matr;
	delete[] row;

	return 0;
}

void mult_matrix_vector (double* result_vect, double* matr, double* vect, int n) 
{
	for (int i = 0; i < n; i++)
	{
		result_vect[i] = 0;

		for (int j = 0; j < n; j++)
			result_vect[i] += matr[j + n * i] * vect[j];

	}
}
