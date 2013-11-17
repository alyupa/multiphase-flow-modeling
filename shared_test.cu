#include "gpu.h"
#include "cuPrintf.cu"

// Проверка ошибок GPU
void checkErrors(const char *label, const char *file, int line)
{
#ifdef MY_TEST
	cudaError_t err;

	err = cudaThreadSynchronize();
	if (err != cudaSuccess)
	{
		char *e = (char*) cudaGetErrorString(err);
		printf("CUDA Error: %s (at %s)\nFile:\"%s\"\nLine:\"%d\"\n\n", e, label, file, line);
	}

	err = cudaGetLastError();
	if (err != cudaSuccess)
	{
		char *e = (char*) cudaGetErrorString(err);
		printf("CUDA Error: %s (at %s)\nFile:\"%s\"\nLine:\"%d\"\n\n", e, label, file, line);
		fflush(stdout);
	}
#endif
}

// Функция, вызываемая при ошибке
__device__ void device_print_error(char *error, const char *file, int line)
{
	CUPRINTF("Error: %s\nFile: \"%s\"\nLine: %d\n\n", error, file, line);
}

// Тест на NaN
// Синтаксис вызова device_test_nan(x, __FILE__, __LINE__);
__device__ void device_test_nan(double x, const char *file, int line)
{
#ifdef MY_TEST
	if ((x > 1e+30) || (x < -1 * 1e+40))
	{
		CUPRINTF("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// Тест на положительное и не NaN
// Синтаксис вызова device_test_positive(x, __FILE__, __LINE__);
__device__ void device_test_positive(double x, const char *file, int line)
{
#ifdef MY_TEST
	if ((x > 1e+30) || (x < 0))
	{
		CUPRINTF("Error: NaN or X<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// Тест на вхождение насыщенностей в [0;1]
// Синтаксис вызова device_test_S(x, __FILE__, __LINE__);
__device__ void device_test_S(double S, const char *file, int line)
{
#ifdef MY_TEST
	if (S < 0)
	{
		CUPRINTF("Error: S<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (S > 1)
	{
		CUPRINTF("Error: S>1\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// Тест на вхождение скоростей в [-100;100]
// Синтаксис вызова test_u(x, __FILE__, __LINE__);
__device__ void device_test_u(double u, const char *file, int line)
{
#ifdef MY_TEST
	if (u < -1e8)
	{
		CUPRINTF("Error: u<-100\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (u > 1e8)
	{
		CUPRINTF("Error: u>100\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// Тест на вхождение плотностей в [0;3000]
// Синтаксис вызова test_ro(x, __FILE__, __LINE__);
__device__ void device_test_ro(double ro, const char *file, int line)
{
#ifdef MY_TEST
	if (ro < 0)
	{
		CUPRINTF("Error: ro < 0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (ro > 3000)
	{
		CUPRINTF("Error: ro > 5000\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}
