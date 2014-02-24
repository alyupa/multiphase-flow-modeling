#include <mpi.h>
#include "defines.h"

extern double *HostBuffer;

// Округление целого в большую сторону
// a - делимое, b - делитель
static int int_ceil(int a, int b)
{
	return (a + b-1) / b;
}
// Время обмена одним double между узлами
static double t_exch (unsigned int N)
{
	return 5.78E-9 * N + 5.54E-5;
}
// Время последовательных загрузки и выгрузки одного double с cpu на gpu
static double t_gpu_load (unsigned int N)
{
	return 3.57E-9 * N + 5.94E-3;
}
// Время расчета одной точки на сpu, с
static double t_cpu_calc (unsigned int N)
{
	return 1.91E-6 * N - 2.96E-2;
}
// Время расчета одной точки на gpu, с
static double t_gpu_calc (unsigned int N)
{
	return 4.7E-8 * N + 5.72E-4;
}

// Недописанная функция оптимального деления сетки по процессорам
void division()
{
	unsigned int Nx, Ny, Nz, size;
	unsigned int N_parameters=20;
	double T=0, T_min=0;
	int flag=0;

	size = def.size;
	Nx = def.Nx;
	Ny = def.Ny;
	Nz = def.Nz;

	unsigned int s_x, s_y, s_z;

	s_x = s_y = s_z = 1;

	for(unsigned int s1=1;s1<=size && s1<Nx;s1++)
		for(unsigned int s2=1;s2<=size/s1 && s2<Ny;s2++)
			for(unsigned int s3=1;s3<=size/(s1*s2) && s3<Nz;s3++)
			{
				double T_calc = t_cpu_calc(int_ceil(Nx, s1) * int_ceil(Ny, s2) * int_ceil(Nz, s3));
				double T_exch = 2. * (my_min(s1-1,1) * t_exch(int_ceil(Ny, s2) * int_ceil(Nz, s3))
					+ my_min(s2-1,1) * t_exch(int_ceil(Nx, s1) * int_ceil(Nz, s3))
					+ my_min(s3-1,1) * t_exch(int_ceil(Ny, s2) * int_ceil(Nx, s1))) * N_parameters;
#ifdef GPU_H
				double T_gpu_cpu = 2. * (my_min(s1-1,1) * t_gpu_load(int_ceil(Ny, s2) * int_ceil(Nz, s3))
					+ min(s2-1,1) * t_gpu_load(int_ceil(Nx, s1) * int_ceil(Nz, s3)) 
					+ min(s3-1,1) * t_gpu_load(int_ceil(Ny, s2) * int_ceil(Nx, s1))) * N_parameters;
#else
				double T_gpu_cpu = 0;
#endif
				T=T_calc + T_exch + T_gpu_cpu;
				if (T < T_min || T_min == 0)
				{
					T_min = T;
					s_x = s1;
					s_y = s2;
					s_z = s3;
				}
			}

		if(!def.rank)
			std::cout<<"s_x="<<s_x<<"  s_y="<<s_y<<"  s_z="<<s_z<<"  T_min="<<T_min<<"  flag="<<flag<<"\n";

		def.sizex = s_x;
		def.sizey = s_y;
		def.sizez = s_z;

}

// Передача и прием данных правой границе
static void right_send_recv(int buffer_size, int destination_rank, int send_recv_id)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer, buffer_size, MPI_DOUBLE, destination_rank, send_recv_id, destination_rank, send_recv_id + 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS)
	{
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\nFile:\"%s\"\nLine:\"%d\"\n\n", __FILE__, __LINE__);
	}
}

// Получение и передача данных на левой границе
static void left_recv_send(int buffer_size, int destination_rank, int send_recv_id)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer, buffer_size, MPI_DOUBLE, destination_rank, send_recv_id + 1, destination_rank, send_recv_id, MPI_COMM_WORLD, &status) == MPI_SUCCESS)
	{
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\nFile:\"%s\"\nLine:\"%d\"\n\n", __FILE__, __LINE__);
	}
}

// Обмен данными на границах между всеми процессорами
// 0. Загружаем данные с усорителя в память хоста
// 1.  Для всех четных процессоров
// 1.1 передаем/получаем правую границу,
// 1.2 получаем/передаем левую границу.
// 2.2 Для нечетных - получаем/передаем левую границу,
// 2.2 передаем/получаем правую.
// Для крайних процессоров соответствующие обмены не требуются
// 3. Загружаем полученные данные в память ускорителя

static void exchange_direct(double* HostArray, double* DevArray, char axis)
{
	if((def.rank) >= (def.sizex) * (def.sizey) * (def.sizez))
		return;

	switch(axis)
	{
	case 'x':
		if(def.sizex > 1) 
		{
			if ((def.rankx) % 2 == 0) // (1)
			{
				if ((def.rankx) != (def.sizex) - 1)
				{
					load_exchange_data_part_xr(HostArray, DevArray); // (0)
					right_send_recv((def.locNy) * (def.locNz), (def.rank) + 1, 500);    // (1.1)
					save_exchange_data_part_xr(HostArray, DevArray); // (3)
				}

				if ((def.rankx) != 0)
				{
					load_exchange_data_part_xl(HostArray, DevArray); // (0)
					left_recv_send((def.locNy) * (def.locNz), (def.rank) - 1, 502);    // (1.2)
					save_exchange_data_part_xl(HostArray, DevArray); // (3)
				}
			}
			else
			{
				if ((def.rankx) != 0) // В принципе, лишняя проверка
				{
					load_exchange_data_part_xl(HostArray, DevArray); // (0)
					left_recv_send((def.locNy) * (def.locNz), (def.rank) - 1, 500);    // (2.1)
					save_exchange_data_part_xl(HostArray, DevArray); // (3)
				}

				if ((def.rankx) != (def.sizex) - 1)
				{
					load_exchange_data_part_xr(HostArray, DevArray); // (0)
					right_send_recv((def.locNy) * (def.locNz), (def.rank) + 1, 502);    // (2.2)
					save_exchange_data_part_xr(HostArray, DevArray); // (3)
				}
			}
		}
		break;
	case 'y':
		if(def.sizey > 1) 
		{
			if ((def.ranky) % 2 == 0) // (1)
			{
				if ((def.ranky) != (def.sizey) - 1)
				{
					load_exchange_data_part_yr(HostArray, DevArray); // (0)
					right_send_recv((def.locNx) * (def.locNz), (def.rank) + (def.sizex), 504);    // (1.1)
					save_exchange_data_part_yr(HostArray, DevArray); // (3)
				}

				if ((def.ranky) != 0)
				{
					load_exchange_data_part_yl(HostArray, DevArray); // (0)
					left_recv_send((def.locNx) * (def.locNz), (def.rank) - (def.sizex), 506);    // (1.2)
					save_exchange_data_part_yl(HostArray, DevArray); // (3)
				}
			}
			else
			{
				if ((def.ranky) != 0) // В принципе, лишняя проверка
				{
					load_exchange_data_part_yl(HostArray, DevArray); // (0)
					left_recv_send((def.locNx) * (def.locNz), (def.rank) - (def.sizex), 504);    // (2.1)
					save_exchange_data_part_yl(HostArray, DevArray); // (3)
				}

				if ((def.ranky) != (def.sizey) - 1)
				{
					load_exchange_data_part_yr(HostArray, DevArray); // (0)
					right_send_recv((def.locNx) * (def.locNz), (def.rank) + (def.sizex), 506);    // (2.2)
					save_exchange_data_part_yr(HostArray, DevArray); // (3)
				}
			}
		}
		break;
	case 'z':
		if(def.sizez > 1) 
		{
			if ((def.rankz) % 2 == 0) // (1)
			{
				if ((def.rankz) != (def.sizez) - 1)
				{
					load_exchange_data_part_zr(HostArray, DevArray); // (0)
					right_send_recv((def.locNx) * (def.locNy), (def.rank) + (def.sizex) * (def.sizey), 508);    // (1.1)
					save_exchange_data_part_zr(HostArray, DevArray); // (3)
				}

				if ((def.rankz) != 0)
				{
					load_exchange_data_part_zl(HostArray, DevArray); // (0)
					left_recv_send((def.locNx) * (def.locNy), (def.rank) - (def.sizex) * (def.sizey), 510);    // (1.2)
					save_exchange_data_part_zl(HostArray, DevArray); // (3)
				}
			}
			else
			{
				if ((def.rankz) != 0) // В принципе, лишняя проверка
				{
					load_exchange_data_part_zl(HostArray, DevArray); // (0)
					left_recv_send((def.locNx) * (def.locNy), (def.rank) - (def.sizex) * (def.sizey), 508);    // (2.1)
					save_exchange_data_part_zl(HostArray, DevArray); // (3)
				}

				if ((def.rankz) != (def.sizez) - 1)
				{
					load_exchange_data_part_zr(HostArray, DevArray); // (0)
					right_send_recv((def.locNx) * (def.locNy), (def.rank) + (def.sizex) * (def.sizey), 510);    // (2.2)
					save_exchange_data_part_zr(HostArray, DevArray); // (3)
				}
			}
		}
		break;
	default:
		break;
	}
}

void exchange(double* HostArray, double* DevArray)
{
	exchange_direct(HostArray, DevArray, 'x');
	exchange_direct(HostArray, DevArray, 'y');
	exchange_direct(HostArray, DevArray, 'z');
}

// Обмен граничными значениями
void exchange_basic_vars()
{
	exchange(HostArraysPtr.P_w, DevArraysPtr.P_w);
	exchange(HostArraysPtr.S_n, DevArraysPtr.S_n);
	exchange(HostArraysPtr.S_w, DevArraysPtr.S_w);
#ifdef ENERGY
	exchange(HostArraysPtr.T, DevArraysPtr.T);
#endif
}

void communication_initialization(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &(def.size)); // The amount of processors
	MPI_Comm_rank(MPI_COMM_WORLD, &(def.rank)); // The number of processor
	//std::cout << "size =" <<defsize<<"  "<<"rank = "<<defrank<<"\n";
}

void communication_finalization(void)
{
	MPI_Finalize();
}

// Реализация фунции Barrier для различных коммуникаций
void barrier(void)
{
	MPI_Barrier(MPI_COMM_WORLD);
}

