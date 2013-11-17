#include <mpi.h>
#include "defines.h"

// Округление целого в большую сторону
// a - делимое, b - делитель
int int_ceil(int a, int b)
{
	return (a + b-1) / b;
}
// Время обмена одним double между узлами
double t_exch (unsigned int N)
{
	return 5.78E-9 * N + 5.54E-5;
}
// Время последовательных загрузки и выгрузки одного double с cpu на gpu
double t_gpu_load (unsigned int N)
{
	return 3.57E-9 * N + 5.94E-3;
}
// Время расчета одной точки на сpu, с
double t_cpu_calc (unsigned int N)
{
	return 1.91E-6 * N - 2.96E-2;
}
// Время расчета одной точки на gpu, с
double t_gpu_calc (unsigned int N)
{
	return 4.7E-8 * N + 5.72E-4;
}

// Недописанная функция оптимального деления сетки по процессорам
void division(consts *def)
{
	unsigned int Nx, Ny, Nz, size;
	unsigned int N_parameters=20;
	double T=0, T_min=0;
	int flag=0;

	size = (*def).size;
	Nx = (*def).Nx;
	Ny = (*def).Ny;
	Nz = (*def).Nz;

	unsigned int s_x, s_y, s_z;

	s_x = s_y = s_z = 1;

	for(unsigned int s1=1;s1<=size && s1<Nx;s1++)
		for(unsigned int s2=1;s2<=size/s1 && s2<Ny;s2++)
			for(unsigned int s3=1;s3<=size/(s1*s2) && s3<Nz;s3++)
			{
				double T_calc = t_cpu_calc(int_ceil(Nx, s1) * int_ceil(Ny, s2) * int_ceil(Nz, s3));
				double T_exch = 2. * (min(s1-1,1) * t_exch(int_ceil(Ny, s2) * int_ceil(Nz, s3)) 
					+ min(s2-1,1) * t_exch(int_ceil(Nx, s1) * int_ceil(Nz, s3)) 
					+ min(s3-1,1) * t_exch(int_ceil(Ny, s2) * int_ceil(Nx, s1))) * N_parameters;
#ifdef GPU_H
				double T_gpu_cpu = 2. * (min(s1-1,1) * t_gpu_load(int_ceil(Ny, s2) * int_ceil(Nz, s3)) 
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

		if(!(*def).rank)
			std::cout<<"s_x="<<s_x<<"  s_y="<<s_y<<"  s_z="<<s_z<<"  T_min="<<T_min<<"  flag="<<flag<<"\n";

		(*def).sizex = s_x;
		(*def).sizey = s_y;
		(*def).sizez = s_z;
}

// Передача и прием данных правой границе
void right_send_recv(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer, buffer_size, MPI_DOUBLE, destination_rank, send_recv_id, destination_rank, send_recv_id + 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS)
	{
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\nFile:\"%s\"\nLine:\"%d\"\n\n", __FILE__, __LINE__);
	}
}

// Получение и передача данных на левой границе
void left_recv_send(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id)
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

void exchange_direct(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, const consts &def, char axis)
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
					load_exchange_data_part_xr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					right_send_recv(HostBuffer, (def.locNy) * (def.locNz), (def.rank) + 1, 500);    // (1.1)
					save_exchange_data_part_xr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}

				if ((def.rankx) != 0)
				{
					load_exchange_data_part_xl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					left_recv_send(HostBuffer, (def.locNy) * (def.locNz), (def.rank) - 1, 502);    // (1.2)
					save_exchange_data_part_xl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}
			}
			else
			{
				if ((def.rankx) != 0) // В принципе, лишняя проверка
				{
					load_exchange_data_part_xl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					left_recv_send(HostBuffer, (def.locNy) * (def.locNz), (def.rank) - 1, 500);    // (2.1)
					save_exchange_data_part_xl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}

				if ((def.rankx) != (def.sizex) - 1)
				{
					load_exchange_data_part_xr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					right_send_recv(HostBuffer, (def.locNy) * (def.locNz), (def.rank) + 1, 502);    // (2.2)
					save_exchange_data_part_xr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
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
					load_exchange_data_part_yr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					right_send_recv(HostBuffer, (def.locNx) * (def.locNz), (def.rank) + (def.sizex), 504);    // (1.1)
					save_exchange_data_part_yr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}

				if ((def.ranky) != 0)
				{
					load_exchange_data_part_yl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					left_recv_send(HostBuffer, (def.locNx) * (def.locNz), (def.rank) - (def.sizex), 506);    // (1.2)
					save_exchange_data_part_yl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}
			}
			else
			{
				if ((def.ranky) != 0) // В принципе, лишняя проверка
				{
					load_exchange_data_part_yl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					left_recv_send(HostBuffer, (def.locNx) * (def.locNz), (def.rank) - (def.sizex), 504);    // (2.1)
					save_exchange_data_part_yl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}

				if ((def.ranky) != (def.sizey) - 1)
				{
					load_exchange_data_part_yr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					right_send_recv(HostBuffer, (def.locNx) * (def.locNz), (def.rank) + (def.sizex), 506);    // (2.2)
					save_exchange_data_part_yr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
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
					load_exchange_data_part_zr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					right_send_recv(HostBuffer, (def.locNx) * (def.locNy), (def.rank) + (def.sizex) * (def.sizey), 508);    // (1.1)
					save_exchange_data_part_zr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}

				if ((def.rankz) != 0)
				{
					load_exchange_data_part_zl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					left_recv_send(HostBuffer, (def.locNx) * (def.locNy), (def.rank) - (def.sizex) * (def.sizey), 510);    // (1.2)
					save_exchange_data_part_zl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}
			}
			else
			{
				if ((def.rankz) != 0) // В принципе, лишняя проверка
				{
					load_exchange_data_part_zl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					left_recv_send(HostBuffer, (def.locNx) * (def.locNy), (def.rank) - (def.sizex) * (def.sizey), 508);    // (2.1)
					save_exchange_data_part_zl(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}

				if ((def.rankz) != (def.sizez) - 1)
				{
					load_exchange_data_part_zr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)
					right_send_recv(HostBuffer, (def.locNx) * (def.locNy), (def.rank) + (def.sizex) * (def.sizey), 510);    // (2.2)
					save_exchange_data_part_zr(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
				}
			}
		}
		break;
	default:
		break;
	}
}

void exchange(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, const consts &def)
{
	exchange_direct(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'x');
	exchange_direct(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'y');
	exchange_direct(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'z');
}

// Обмен граничными значениями давления P2, плотностей ro1 и ro2, Xi между процессорами
void P_ro_Xi_exchange(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def)
{
	exchange(HostArraysPtr.ro_g, DevArraysPtr.ro_g, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.Xi_g, DevArraysPtr.Xi_g, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.P_g, DevArraysPtr.P_g, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.P_n, DevArraysPtr.P_n, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.ro_w, DevArraysPtr.ro_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.ro_n, DevArraysPtr.ro_n, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.Xi_w, DevArraysPtr.Xi_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.Xi_n, DevArraysPtr.Xi_n, HostBuffer, DevBuffer, def);
}


// Обмен граничными значениями скоростей между процессорами
// В случае распределения расчетной области между процессорами по оси X
// передача u1y и u2y не требуется
void u_exchange(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def)
{
	exchange_direct(HostArraysPtr.ux_w, DevArraysPtr.ux_w, HostBuffer, DevBuffer, def, 'x');
	exchange_direct(HostArraysPtr.ux_n, DevArraysPtr.ux_n, HostBuffer, DevBuffer, def, 'x');
	exchange_direct(HostArraysPtr.uy_w, DevArraysPtr.uy_w, HostBuffer, DevBuffer, def, 'y');
	exchange_direct(HostArraysPtr.uy_n, DevArraysPtr.uy_n, HostBuffer, DevBuffer, def, 'y');
	exchange_direct(HostArraysPtr.ux_g, DevArraysPtr.ux_g, HostBuffer, DevBuffer, def, 'x');
	exchange_direct(HostArraysPtr.uy_g, DevArraysPtr.uy_g, HostBuffer, DevBuffer, def, 'y');

	if((def).Nz >= 2)
	{
		exchange_direct(HostArraysPtr.uz_w, DevArraysPtr.uz_w, HostBuffer, DevBuffer, def, 'z');
		exchange_direct(HostArraysPtr.uz_n, DevArraysPtr.uz_n, HostBuffer, DevBuffer, def, 'z');
		exchange_direct(HostArraysPtr.uz_g, DevArraysPtr.uz_g, HostBuffer, DevBuffer, def, 'z');
	}
}

// Обмен граничными значениями давления воды P1 и насыщенности NAPL S2 между процессорами
void P_S_exchange(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def)
{
	exchange(HostArraysPtr.S_w, DevArraysPtr.S_w, HostBuffer, DevBuffer, def);
#ifdef ENERGY
	exchange(HostArraysPtr.T, DevArraysPtr.T, HostBuffer, DevBuffer, def);
#endif
	exchange(HostArraysPtr.P_w, DevArraysPtr.P_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.S_n, DevArraysPtr.S_n, HostBuffer, DevBuffer, def);
}

void communication_initialization(int argc, char* argv[], consts* def)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &((*def).size)); // The amount of processors
	MPI_Comm_rank(MPI_COMM_WORLD, &((*def).rank)); // The number of processor
	//std::cout << "size =" <<(*def).size<<"  "<<"rank = "<<(*def).rank<<"\n";
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

