#include "../../defines.h"

#define MEASURE_COUNT 10

// Буферные массивы для обмена между процессорами
double *HostBuffer;
double *DevBuffer;

clock_t main_part(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, int argc, char* argv[], consts* def);

int main(int argc, char* argv[])
{
#ifdef GTEST
	::testing::InitGoogleTest(&argc, argv);
	RUN_ALL_TESTS();
#endif

	consts def;
	read_defines(argc, argv, &def);

	// Хостовый массив данных расчетной области процессора
	ptr_Arrays HostArraysPtr;
	// GPU-массив данных расчетной области процессора
	ptr_Arrays DevArraysPtr;

	// Счетчик времени исполнения вычислительной части программы
	char fname[] = "3ph_point_times.txt";
	FILE *fp;

	clock_t task_times[MEASURE_COUNT];
	int smesh_sizes[MEASURE_COUNT];
	double task_time_av = 0;

	if (!(fp = fopen(fname, "a")))
		printf("Not open file: %s", fname);
	fprintf(fp, "size\t time, s\n");

	for(int i = 0; i < MEASURE_COUNT; i++)
	{
		def.Nx = 20 + 2 * i;
		def.Ny = 20 + 2 * i;
		def.Nz = 20 + 2 * i;

		task_times[i] = main_part(&HostArraysPtr, &DevArraysPtr, argc, argv, &def);
		smesh_sizes[i] = (def.Nx) * (def.Ny) * (def.Nz);
		fprintf(fp, "%d\t %.5f\n", smesh_sizes[i], (double)task_times[i] / CLOCKS_PER_SEC);
	}

	for(int i = 0; i < MEASURE_COUNT; i++)
		task_time_av += ((double)task_times[i] / smesh_sizes[i]);

	task_time_av /= CLOCKS_PER_SEC;
	task_time_av /= MEASURE_COUNT;
	task_time_av /= ((def.timeX) / (def).dt);

	fprintf(fp, "One point task time in seconds:\t%e\n", task_time_av);
	printf("One point task time in seconds:\t%e\n", task_time_av);

	fclose(fp);
	// При запуске в Windows после работы программы оставлять окно консоли
#ifdef _WIN32
	printf("\nPress <Enter> to exit...\n");
	fflush(stdout);
	getchar();
#endif
	return 0;
}

clock_t main_part(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, int argc, char* argv[], consts* def)
{
	long int time_counter = 0;
	clock_t task_time;

	// Инициализация коммуникаций, перевод глобальных параметров в локальные процессора,
	// выделение памяти, загрузка начальных/сохраненных данных
	initialization(HostArraysPtr, DevArraysPtr, &time_counter, argc, argv, def);

	task_time = clock();

	// Цикл шагов по времени (каждая итерация - новый слой по времени)
	// 1. Проводятся расчеты P1 и S2 на следующем временном слое
	// 2. Каждые (def.print_screen) раз на экран выводится информация о временном слое
	// 3. Каждые save_plots раз данные выгружаются в память хоста и
	//    сохраняются в файлы графиков (**), в файл сохраняется состояние задачи (***)
	for (time_counter++; time_counter <= (*def).timeX / ((*def).dt); time_counter++)
	{
		time_step_function(*HostArraysPtr, *DevArraysPtr, DevBuffer, *def, time_counter * ((*def).dt)); // (1)
	}

	// Вывод информации о времени работы программы в секундах
	task_time = clock() - task_time;

	// Завершение работы и освобождение памяти
	finalization(*HostArraysPtr, *DevArraysPtr, DevBuffer);

	return task_time;
}

// Функция полного цикла расчетов на следующем временном слое
// 1. Расчет плотности жидкостей ro, давлений NAPL P2, переменных Xi
// 2. Обмен между процессорами пограничными значениями P2, ro и Xi
// 3. Расчет скоростей жидкостей
// 4. Обмен между процессорами пограничными значениями скоростей жидкостей
// 5. Расчет переменной roS на следующем временном слое
// 6. Расчет методом Ньютона давления воды P1 и насыщенности DNAPL S2
// 7. Применение граничных условий для P1 и S2
// 8. Обмен между процессорами пограничными значениями P1 и S2
void time_step_function(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer, const consts &def, double t)
{
	P_S_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, def); // (8)
	ro_P_Xi_calculation(HostArraysPtr, DevArraysPtr, def); // (1)
	P_ro_Xi_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, def); // (2)
	u_calculation(HostArraysPtr, DevArraysPtr, def); // (3)
	u_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, def); // (4)
	roS_calculation(HostArraysPtr, DevArraysPtr, t, def); // (5)
	P_S_calculation(HostArraysPtr, DevArraysPtr, def); // (6)
	boundary_conditions(HostArraysPtr, DevArraysPtr, def); // (7)
}

// Преобразование локальных координат процессора к глобальным
// Каждый процессор содержит дополнительную точку в массиве для
// обмена данными, если имеет соседа
// (если 2 соседа с обеих сторон,то +2 точки).
// Глобальные границы хранятся как обычные точки (отсюда и условие на (def.rank)==0)
int local_to_global(int local_index, char axis, const consts &def)
{
	int global_index = local_index;
	switch (axis)
	{
	case 'x':
	{
		global_index += def.rankx * (def.Nx / def.sizex) + min(def.rankx, def.Nx % def.sizex) - min(def.rankx, 1);
		break;
	}
	case 'y':
	{
		global_index += def.ranky * (def.Ny / def.sizey) + min(def.ranky, def.Ny % def.sizey) - min(def.ranky, 1);
		break;
	}
	case 'z':
	{
		global_index += def.rankz * (def.Nz / def.sizez) + min(def.rankz, def.Nz % def.sizez) - min(def.rankz, 1);
		break;
	}
	default:
	{
		print_error("Axis of [local to global] conversation is empty", __FILE__, __LINE__);
	}
	}
	//some_test(global_index);
	return global_index;
}

// Вычисление расчетной области (нагрузки) процессора
// Если поровну не распределяется, то первые (def.NX)%size получают +1 в нагрузку.
// Каждый процессор содержит дополнительную точку в массиве для
// обмена данными, если имеет соседа
// (если 2 соседа с обеих сторон,то +2 точки).
// Глобальные границы хранятся как обычные точки (отсюда и условие на (def.rank)==0)
void global_to_local_vars(consts *def)
{
	(*def).locNx = (*def).Nx / (*def).sizex;

	if ((*def).rankx < (*def).Nx % (*def).sizex)
	{
		((*def).locNx) ++;
	}

	// Крайние процессоры получают по 1 точке для граничных данных,
	// остальные - по 2 на обе границы
	// Если процессор один, то границ у него нет и дополнительные точки не нужны
	if ((*def).sizex > 1)
	{
		if (((*def).rankx == 0) || ((*def).rankx  == (*def).sizex - 1))
		{
			((*def).locNx) ++;
		}
		else
		{
			((*def).locNx) += 2;
		}
	}

	(*def).locNy = (*def).Ny / (*def).sizey;

	if (((*def).ranky < (*def).Ny % (*def).sizey))
	{
		((*def).locNy) ++;
	}

	if ((*def).sizey > 1)
	{
		if (((*def).ranky == 0) || ((*def).ranky == (*def).sizey - 1))
		{
			((*def).locNy) ++;
		}
		else
		{
			((*def).locNy) += 2;
		}
	}

	(*def).locNz = (*def).Nz / (*def).sizez;

	if ((*def).rankz < (*def).Nz % (*def).sizez)
	{
		((*def).locNz) ++;
	}

	if ((*def).sizez > 1)
	{
		if (((*def).rankz == 0) || ((*def).rankz == (*def).sizez - 1))
		{
			((*def).locNz) ++;
		}
		else
		{
			((*def).locNz) += 2;
		}
	}

	if((*def).rank >= (*def).sizex * (*def).sizey * (*def).sizez)
	{
		(*def).locNx = (*def).locNy = (*def).locNz = 0;
	}

	test_positive((*def).locNx, __FILE__, __LINE__);
	test_positive((*def).locNy, __FILE__, __LINE__);
	test_positive((*def).locNz, __FILE__, __LINE__);
}

// Является ли точка активной (т.е. не предназначенной только для обмена на границах)
int is_active_point(int i, int j, int k, const consts &def)
{
	if (((def.rankx) != 0 && i == 0) || ((def.rankx) != (def.sizex) - 1 && i == (def.locNx) - 1)
	    || ((def.ranky) != 0 && j == 0)	|| ((def.ranky) != (def.sizey) - 1 && j == (def.locNy) - 1)
	    || ((((def.rankz) != 0 && k == 0) || ((def.rankz) != (def.sizez) - 1 && k == (def.locNz) - 1)) && (def.Nz) >= 2))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

// Применение начальных данных во всех точках
void sizes_initialization(consts *def)
{
	if((*def).size == 1) 
	{
		(*def).sizex = 1;
		(*def).sizey = 1;
		(*def).sizez = 1;
	} else {
		(*def).sizex = SizeX;
		(*def).sizey = SizeY;
		(*def).sizez = SizeZ;
	}
	(*def).rankx = (*def).rank % (*def).sizex;
	(*def).ranky = ((*def).rank / (*def).sizex) % (*def).sizey;
	(*def).rankz = ((*def).rank / (*def).sizex) / (*def).sizey;
}

void blocks_initialization(consts *def)
{
	(*def).blocksX = 0;
	(*def).blocksY = 0;
	(*def).blocksZ = 0;
}

// Функция вычисления "эффективной" плотности * g * hy
double ro_eff_gdy(ptr_Arrays HostArraysPtr, int local, const consts &def)
{
#ifdef THREE_PHASE
	double ro_g_dy = (HostArraysPtr.ro_g[local] * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local])
					+ HostArraysPtr.ro_w[local] * HostArraysPtr.S_w[local]
					+ HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local]) * (HostArraysPtr.m[local]) * (def.g_const) * (def.hy);

#else
	double ro_g_dy = (HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local]
	                  + HostArraysPtr.ro_w[local] * (1 - HostArraysPtr.S_n[local])) * (HostArraysPtr.m[local]) * (def.g_const) * (def.hy);
#endif
	return ro_g_dy;
}

//----------------------------------------------------------------------------------------------------
// Служебные функции

// Вывод запускаемой задачи
void print_task_name(const consts &def)
{
	// Нулевой процессор выводит название запускаемой задачи
	if (!(def.rank))
	{
#ifdef TWO_PHASE
		char task_name[] = "Two phase filtration";
#endif
#ifdef THREE_PHASE
		char task_name[] = "Three phase filtration";
#endif
#ifdef B_L
		char task_name[] = "Backley-Leverett filtration";
#endif
		std::cout << task_name << " by CAPAZ on " << (def.size) << " node(s).\n";
		read_version();
		fflush(stdout);
	}
}

// Инициализация коммуникаций (1), перевод глобальных параметров в локальные процессора (2),
// инициализация ускорителя (2.5), выделение памяти (3), загрузка начальных/сохраненных данных (4)
// Для задачи Б-Л загрузка проницаемостей из файла.
void initialization(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, long int* time_counter, int argc, char* argv[], consts* def)
{
	FILE *f_save;

	communication_initialization(argc, argv, def); // (1)

	//print_task_name(*def);

	sizes_initialization(def);

	blocks_initialization(def);

	global_to_local_vars(def); // (2)

	device_initialization(def); // (2.5)

	memory_allocation(HostArraysPtr, DevArraysPtr, *def); // (3)

	load_permeability((*HostArraysPtr).K, *def); // (5)
	load_data_to_device((*HostArraysPtr).K, (*DevArraysPtr).K, *def);

	// Если процессор может открыть файл сохраненного состояния,
	// то восстанавливаем состояние, иначе применяем начальные условия
	if (f_save = fopen("save/save.dat", "rb"))
	{
		fclose(f_save);
		restore(*HostArraysPtr, time_counter, *def);
	}
	else
	{
		data_initialization(*HostArraysPtr, time_counter, *def);    // (4)
	}

#ifdef THREE_PHASE
	load_data_to_device((*HostArraysPtr).S_w, (*DevArraysPtr).S_w, *def);
	load_data_to_device((*HostArraysPtr).roS_g_old, (*DevArraysPtr).roS_g_old, *def);
	load_data_to_device((*HostArraysPtr).P_n, (*DevArraysPtr).P_n, *def);
	load_data_to_device((*HostArraysPtr).P_g, (*DevArraysPtr).P_g, *def);
#endif
	load_data_to_device((*HostArraysPtr).P_w, (*DevArraysPtr).P_w, *def);
	load_data_to_device((*HostArraysPtr).S_n, (*DevArraysPtr).S_n, *def);
	load_data_to_device((*HostArraysPtr).roS_w_old, (*DevArraysPtr).roS_w_old, *def);
	load_data_to_device((*HostArraysPtr).roS_n_old, (*DevArraysPtr).roS_n_old, *def);
	load_data_to_device((*HostArraysPtr).m, (*DevArraysPtr).m, *def);
}

// Завершение работы (1), освобождение памяти (2)
void finalization(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer)
{
	memory_free(HostArraysPtr, DevArraysPtr); // (2)
	communication_finalization(); // (1)
	device_finalization(); // (1)
}

// Выделение памяти хоста (1) и ускорителя (2) под массив точек расчетной области
void memory_allocation(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, const consts &def)
{
	host_memory_allocation(HostArraysPtr, def); // (1)
	device_memory_allocation(DevArraysPtr, &DevBuffer, def); // (2)
}

// Освобожение памяти хоста (1) и ускорителя (2) из под массива точек расчетной области
void memory_free(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr)
{
	host_memory_free(HostArraysPtr); // (1)
	device_memory_free(DevArraysPtr, DevBuffer); // (2)
}

// Выделение памяти хоста под массив точек расчетной области
void host_memory_allocation(ptr_Arrays* ArraysPtr, const consts &def)
{
	int buffer_size = 0;

	if(def.sizex > 1)
		buffer_size = (def.locNy) * (def.locNz);
	if(def.sizey > 1 && (def.locNx) * (def.locNz) > buffer_size)
		buffer_size = (def.locNx) * (def.locNz);
	if(def.sizez > 1 && (def.locNx) * (def.locNy) > buffer_size)
		buffer_size = (def.locNx) * (def.locNy);

	if(buffer_size)
	{
		try
		{
			HostBuffer = new double[buffer_size];
		}
		catch(const std::bad_alloc &)
		{
			print_error("Memory for *HostBuffer is not allocated in function host_memory_alloc", __FILE__, __LINE__);
		}
	}

	try
	{
		(*ArraysPtr).S_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).P_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).P_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).ro_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).ro_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).ux_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).uy_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).uz_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).ux_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).uy_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).uz_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).Xi_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).Xi_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).roS_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).roS_w_old = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).roS_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).roS_n_old = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).m = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).K = new double [(def.locNx) * (def.locNy) * (def.locNz)];
#ifdef THREE_PHASE
		(*ArraysPtr).P_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).S_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).ro_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).ux_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).uy_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).uz_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).Xi_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).roS_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).roS_g_old = new double [(def.locNx) * (def.locNy) * (def.locNz)];
#endif
	}
	catch (...)
	{
		print_error("Not enough host memory for ArraysPtr", __FILE__, __LINE__);
	}
}

// Освобожение памяти хоста из под массива точек расчетной области
void host_memory_free(ptr_Arrays ArraysPtr)
{
	delete HostBuffer;
	delete[] ArraysPtr.P_w;
	delete[] ArraysPtr.P_n;
	delete[] ArraysPtr.ro_w;
	delete[] ArraysPtr.ro_n;
	delete[] ArraysPtr.ux_w;
	delete[] ArraysPtr.uy_w;
	delete[] ArraysPtr.uz_w;
	delete[] ArraysPtr.ux_n;
	delete[] ArraysPtr.uy_n;
	delete[] ArraysPtr.uz_n;
	delete[] ArraysPtr.Xi_w;
	delete[] ArraysPtr.Xi_n;
	delete[] ArraysPtr.roS_w;
	delete[] ArraysPtr.roS_w_old;
	delete[] ArraysPtr.roS_n;
	delete[] ArraysPtr.roS_n_old;
	delete[] ArraysPtr.m;
	delete[] ArraysPtr.K;
#ifdef THREE_PHASE
	delete[] ArraysPtr.P_g;
	delete[] ArraysPtr.S_w;
	delete[] ArraysPtr.ro_g;
	delete[] ArraysPtr.ux_g;
	delete[] ArraysPtr.uy_g;
	delete[] ArraysPtr.uz_g;
	delete[] ArraysPtr.Xi_g;
	delete[] ArraysPtr.roS_g;
	delete[] ArraysPtr.roS_g_old;
#endif
	delete[] ArraysPtr.S_n;
}

// Функция сохранения графиков в файлы
void save_data_plots(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, const consts &def)
{
	// Загрузка в память хоста результатов расчета
#ifdef THREE_PHASE
	load_data_to_host(HostArraysPtr.S_w, DevArraysPtr.S_w , def);
	load_data_to_host(HostArraysPtr.ux_w, DevArraysPtr.ux_w , def);
	load_data_to_host(HostArraysPtr.uy_w, DevArraysPtr.uy_w , def);
	load_data_to_host(HostArraysPtr.uz_w, DevArraysPtr.uz_w , def);
	load_data_to_host(HostArraysPtr.ux_g, DevArraysPtr.ux_g , def);
	load_data_to_host(HostArraysPtr.uy_g, DevArraysPtr.uy_g , def);
	load_data_to_host(HostArraysPtr.uz_g, DevArraysPtr.uz_g , def);
#endif
	load_data_to_host(HostArraysPtr.P_w, DevArraysPtr.P_w , def);
	load_data_to_host(HostArraysPtr.S_n, DevArraysPtr.S_n , def);
	load_data_to_host(HostArraysPtr.ux_n, DevArraysPtr.ux_n , def);
	load_data_to_host(HostArraysPtr.uy_n, DevArraysPtr.uy_n , def);
	load_data_to_host(HostArraysPtr.uz_n, DevArraysPtr.uz_n , def);

#ifndef THREE_PHASE
	// Проверка на выход из допустимого диапазона значений P и S
#ifdef MY_TEST
	test_correct_P_S(HostArraysPtr, def);
#endif
#endif

	// Нулевой процессор создает директории, файлы и прописывает заголовки файлов
	if ((def.rank) == 0)
	{
		print_plots_top(t, def);
	}

	if((def.sizey == 1) && (def.sizez == 1))
	{
		// По очереди для каждого из процессоров вызываем функцию вывода на график
		// своей части массива.
		for (int cpu = 0; cpu < ((def.sizex) * (def.sizey) * (def.sizez)); cpu ++)
		{
			// Реализация фунции Barrier для различных коммуникаций
			barrier();
			if ((def.rank) == cpu)
			{
				print_plots(HostArraysPtr, t, def, -1, -1);
			}
		}
	} 
	else 
	{
		// Выстраиваем такую очередь вывода своей части массива каждого из процессоров, 
		// чтобы значения создавали последовательность, корректно отображаемую Tecplot.
		for (int cpux = 0; cpux < (def.sizex); cpux ++)
			for (int i = 0; i < def.Nx / def.sizex + 3; i++)
				for (int cpuy = 0; cpuy < (def.sizey); cpuy ++)
					for (int j = 0; j < def.Ny / def.sizey + 3; j++)
						for (int cpuz = 0; cpuz < (def.sizez); cpuz ++)
						{
							// Реализация фунции Barrier для различных коммуникаций
							barrier();
							if ((def.rank) == cpux + cpuy * (def.sizex) + cpuz * (def.sizex) * (def.sizey) && i < def.locNx && j < def.locNy)
							{
								print_plots(HostArraysPtr, t, def, i, j);
							}
						}
	}

}

// Функция создания директорий, файлов для графиков и сохранения заголовков в них
void print_plots_top(double t, const consts &def)
{
	char fname[30];
	FILE *fp;

	sprintf(fname, "plots/S=%012.4f.dat", t);

#ifdef _WIN32
	_mkdir("plots");
#else
	mkdir("plots", 0000777);
#endif

	// Создание (или перезапись) файла с графиками
	// 1. Для распределения насыщенностей NAPL S_n
	// 2. Для распределения давлений воды P_w
	// 3. Для распределения скоростей {u_x, u_y, u_z}
	// 4. Для распределения типов грунтов
	if (!(fp = fopen(fname, "wt")))
		print_error("Not open file(s) in function SAVE_DATA_PLOTS", __FILE__, __LINE__);

	fprintf(fp, "TITLE =  \"Filtration in time=%5.2f\" \n", t);

	if ((def.Nz) < 2)
	{
#ifdef THREE_PHASE
		//		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"u_x\",\"u_y\",\"porosity\" \n");
		fprintf(fp, "VARIABLES = \"X\",\"Y\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"uw_x\",\"uw_y\",\"un_x\",\"un_y\",\"ug_x\",\"ug_y\",\"porosity\" \n");
#else
		fprintf(fp, "VARIABLES = \"X\",\"Y\",\"S_w\",\"P_w\",\"u_x\", \"u_y\", \"porosity\" \n");
#endif
		fprintf(fp, "ZONE T = \"BIG ZONE\", K=%d,J=%d, F = POINT\n", (def.Nx), (def.Ny));
	}
	else
	{
#ifdef THREE_PHASE
		//		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"u_x\", \"u_y\",\"u_z\",\"porosity\" \n");
		fprintf(fp, "VARIABLES = \"X\",\"Y\",\"Z\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"uw_x\",\"uw_y\",\"uw_z\",\"un_x\",\"un_y\",\"un_z\",\"ug_x\",\"ug_y\",\"ug_z\",\"porosity\" \n");
#else
		fprintf(fp, "VARIABLES = \"X\",\"Y\",\"Z\",\"S_n\",\"P_w\",\"u_x\", \"u_y\", \"u_z\", \"porosity\" \n");
#endif
		fprintf(fp, "ZONE T = \"BIG ZONE\", K=%d,J=%d,I=%d, F = POINT\n", (def.Nx), (def.Ny), (def.Nz));
	}

#ifdef B_L_1
	char fname_xz[30];
	FILE *fp_xz;

	sprintf(fname_xz, "plots/xz=%012.4f.dat", t);
	if (!(fp_xz = fopen(fname_xz, "wt")))
		print_error("Not open file(s) in function SAVE_DATA_PLOTS", __FILE__, __LINE__);

	fprintf(fp_xz, "TITLE =  \"Filtration in time=%5.2f\" \n", t);
	fprintf(fp_xz, "VARIABLES = \"X\",\"Y\",\"S_n\",\"P_w\",\"porosity1\", \"porosity2\",\"porosity3\" \n");
	fprintf(fp_xz, "ZONE T = \"BIG ZONE\", K=%d,J=%d, F = POINT\n", (def.Nx), (def.Nz));
	fclose(fp_xz);
#endif

	fclose(fp);
}

// Функция сохранения данных в файлы графиков
// Если деления между процессорами по y, z нет, то параметры I, J не используются. 
// Если же есть, но хочется писать в файл в простой последовательности, то должно быть I=J=-1
void print_plots(ptr_Arrays HostArraysPtr, double t, const consts &def, int I, int J)
{
	char fname[30];
	FILE *fp;

	sprintf(fname, "plots/S=%012.4f.dat", t);

	// Открытие на дозапись и сохранение графиков
	if (!(fp = fopen(fname, "at")))
		print_error("Not open file(s) in function SAVE_DATA_PLOTS", __FILE__, __LINE__);

	if((def.sizey == 1) && (def.sizez == 1) || (I == -1) && (J == -1))
	{
		for (int i = 0; i < def.locNx; i++)
			for (int j = 0; j < def.locNy; j++)
				for (int k = 0; k < def.locNz; k++)
					if (is_active_point(i, j, k, def))
						print_plot_row(HostArraysPtr, fp, i, j, k, def);
	}
	else
	{
		for (int k = 0; k < def.locNz; k++)
			if (is_active_point(I, J, k, def))
				print_plot_row(HostArraysPtr, fp, I, J, k, def);
	}

	fclose(fp);

#ifdef B_L_1
	char fname_xz[30];
	FILE *fp_xz;

	sprintf(fname_xz, "plots/xz=%012.4f.dat", t);

	if (!(fp_xz = fopen(fname_xz, "at")))
		print_error("Not open file(s) in function SAVE_DATA_PLOTS", __FILE__, __LINE__);

	for (int i = 0; i < def.locNx; i++)
		for (int k = 0; k < def.locNz; k++)
		{
			int j1=def.locNz/3;
			int j2=def.locNz/2;
			int j3=def.locNz/3*2;
			// Если is_active_point(i, j1, k, def) правда, то и для j2, j3 тоже правда
				if (is_active_point(i, j1, k, def))
				{
					local = i + j2 * def.locNx + k * def.locNx * def.locNy;
					int I = local_to_global(i, 'x', def);

					fprintf(fp_xz, "%.2e %.2e %.3e %.3e %.3e %.3e %.3e\n", I * (def.hx), K * (def.hz), HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.K[i + j1 * def.locNx + k * def.locNx * def.locNy], HostArraysPtr.K[i + j2 * def.locNx + k * def.locNx * def.locNy], HostArraysPtr.K[i + j3 * def.locNx + k * def.locNx * def.locNy]); // (1)
				}
		}
	fclose(fp_xz);
#endif
}

// Функция записи строчки значений параметров, которые нужны для построения графиков,  для всех задач.
// Включает в себя:
// 1. Распределение насыщенностей NAPL S_n
// 2. Распределение давлений воды P_w
// 3. Распределение скоростей {u_x, u_y, u_z}
// 4. Распределение типов грунтов
void print_plot_row(ptr_Arrays HostArraysPtr, FILE* fp, int i, int j, int k, const consts &def)
{
	int local = i + j * def.locNx + k * def.locNx * def.locNy;

	// Преобразование локальных координат процессора к глобальным
	int I = local_to_global(i, 'x', def);
	int J = def.Ny - 1 - local_to_global(j, 'y', def);
	int K = local_to_global(k, 'z', def); 

#ifdef THREE_PHASE
	if (def.Nz < 2)
	{
		/*						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", I*(def.hx), J*(def.hy),
									HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local],
									HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.m[local]);
		*/
		fprintf(fp, "%.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", I * (def.hx), J * (def.hy),
				HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local],
				HostArraysPtr.ux_w[local], (-1)*HostArraysPtr.uy_w[local], HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.ux_g[local],
				(-1)*HostArraysPtr.uy_g[local], HostArraysPtr.m[local]);

	}

	else
	{
		/*						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", I*(def.hx), K*(def.hz), J*(def.hy),
									HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local],
									HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.m[local]);
		*/
		fprintf(fp, "%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", I * (def.hx), K * (def.hz), J * (def.hy),
				HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local],
				HostArraysPtr.ux_w[local], HostArraysPtr.uz_w[local], (-1)*HostArraysPtr.uy_w[local],
				HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local],
				HostArraysPtr.ux_g[local], HostArraysPtr.uz_g[local], (-1)*HostArraysPtr.uy_g[local], HostArraysPtr.m[local]);

	}
#endif
#ifdef TWO_PHASE
	if (def.Nz < 2)
	{
		fprintf(fp, "%.2e %.2e %.3e %.3e %.3e %.3e %.3e\n", I * (def.hx), J * (def.hy), HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.m[local]); // (1)

	}
	else
	{
		fprintf(fp, "%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e\n", I * (def.hx), K * (def.hz), J * (def.hy), HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.m[local]); // (1)
	}
#endif
#ifdef B_L
	if (def.Nz < 2)
	{
		fprintf(fp, "%.2e %.2e %.3e %.3e %.3e %.3e %.3e\n", I * (def.hx), J * (def.hy), 1.-HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], HostArraysPtr.uy_n[local], HostArraysPtr.K[local]); // (1)

	}
	else
	{
		fprintf(fp, "%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e\n", I * (def.hx), K * (def.hz), J * (def.hy), 1-HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.K[local]); // (1)
	}
#endif
}

// Функция сохранения данных в файлы графиков формата BjnIO [http://lira.imamod.ru/BjnIO_3D.html]
void print_plots_BjnIO(ptr_Arrays HostArraysPtr, double t, const consts &def)
{
	/*
		char *dname;
		char *targetfuncs="r2d.bjn";
		char *targetgrid="r2dgrid.bjn";
		LPLPLPFLOAT F1;
		double *F1D;

		F1 = alloc_float_mas_n1_n2_n3(def.locNx, def.locNy, def.locNz);
		if(F1 == NULL) exit(0);

		F1D = (double *)malloc(def.locNx, def.locNy, def.locNz*sizeof(double));
		if(F1D == NULL) exit(0);

		// Инициализация файла с данными функции
		int err = WriteBjnGzippedScalar8RecInit(targetfuncs, "r2dgrid.bjn", def.Nx, def.Ny, def.Nz);
		if(err)
			fprintf(stderr, "Can't create file %s.\n", targetfuncs);

		// Запись блока данных значений функции
		err = WriteBjnGzippedScalar8RecFuncByBlock(targetfuncs, "S_n", F1D, def.locNx*(def.rank), 0, 0, def.locNx, def.locNy, def.locNz, 5); // Поправить def.locNx*(def.rank) на точное смещение
		if(err)
			fprintf(stderr, "Can't add func block data err %d\n", err);

		// Запись минимального и максимального значения сеточной функции S_n
		err = WriteBjnGzippedScalar8RecFuncMinMax(targetfuncs, "S_n", 0, 1);
		if(err)
			fprintf(stderr, "Can't add data about block err %d\n", err);

		// Инициализация файла сетки
		err = WriteBjnGzippedScalar8RecInit(targetgrid, "", def.Nx, def.Ny, def.Nz);
		if(err)
			fprintf(stderr, "Can't create file %s.\n", targetgrid);

		// Для каждого из направлений
		for(direct=0;direct<3;direct++)
		{
			for(i1=0; i1<m1; i1++)
				for(i2=0; i2<m2; i2++)
					for(i3=0; i3<m3; i3++)
					{
						for(i=0; i<n1; i++)
							for(j=0; j<n2; j++)
								for(k=0; k<n3; k++)
								{
									float x=(float)i1*(float)n1+(float)i;
									float y=(float)i2*(float)n2+(float)j;
									float z=(float)i3*(float)n3+(float)k;
									switch(direct)
									{
										case 0: F1[i][j][k] = (float)x; dname="x"; break;
										case 1: F1[i][j][k] = (float)y; dname="y"; break;
										case 2: F1[i][j][k] = (float)(ffmin+k*(ffmax-ffmin)/3.f); dname="z"; break;
									}

									fmin=minab(fmin,F1[i][j][k]);
									fmax=maxab(fmax,F1[i][j][k]);
								}

					// Запись блока данных сетки
					err = WriteBjnGzippedScalar8RecFuncByBlock(targetgrid, dname, F1, def.locNx*(def.rank), 0, 0, def.locNx, def.locNy, def.locNz, 9); // Поправить def.locNx*(def.rank) на точное смещение
					if(err)
						fprintf(stderr, "Can't add grid `%s` block data err %d\n", dname, err);
					}

			// Запись максимального и минимального значений сетки
			err = WriteBjnGzippedScalar8RecFuncMinMax(targetgrid, dname, fmin, fmax);
			if(err)
				fprintf(stderr, "Can't add data about block err %d\n", err);
		}
	*/
}

// Функция вывода на экран двумерного массива(для тестирования пересылок)
void print_array_console(double* Arr, const consts &def, char axis)
{
	printf("\n");
	switch(axis)
	{
	case 'x':
		printf("left:\n");
		for (int i = 0; i < def.locNz; i++) 
		{
			for (int j = 0; j < def.locNy; j++)
				printf("%.2f  ", Arr[j * def.locNx + i * def.locNx * def.locNy]);
			printf("\n");
		}
		printf("right:\n");
		for (int i = 0; i < def.locNz; i++) 
		{
			for (int j = 0; j < def.locNy; j++)
				printf("%.2f  ", Arr[def.locNx - 1 + j * def.locNx + i * def.locNx * def.locNy]);
			printf("\n");
		}
		break;
	case 'y':
		printf("left:\n");
		for (int i = 0; i < def.locNz; i++) 
		{
			for (int j = 0; j < def.locNx; j++)
				printf("%.2f  ", Arr[j + i * def.locNx * def.locNy]);
			printf("\n");
		}
		printf("right:\n");
		for (int i = 0; i < def.locNz; i++) 
		{
			for (int j = 0; j < def.locNx; j++)
				printf("%.2f  ", Arr[j + (def.locNy - 1) * def.locNx + i * def.locNx * def.locNy]);
			printf("\n");
		}
		break;
	case 'z':
		printf("left:\n");
		for (int i = 0; i < def.locNy; i++) 
		{
			for (int j = 0; j < def.locNx; j++)
				printf("%.2f  ", Arr[j + i * def.locNx]);
			printf("\n");
		}
		printf("right:\n");
		for (int i = 0; i < def.locNy; i++) 
		{
			for (int j = 0; j < def.locNx; j++)
				printf("%.2f  ", Arr[j + i * def.locNx + (def.locNz - 1) * def.locNx * def.locNy]);
			printf("\n");
		}
		break;
	default:
		break;
	}
	printf("\n");
}

// Функция загрузки файла проницаемостей
void load_permeability(double* K, const consts &def)
{
#ifndef LOAD_K_FROM_FILE
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				K[i + j * def.locNx + k * def.locNx * def.locNy] = def.K[0];
#else
	FILE *input;
	char *file = "../noise.dat";

	if (!(input = fopen(file, "rt")))
	{
		file = "noise.dat";
		if (!(input = fopen(file, "rt")))
		{
			char err[60];
			sprintf(err, "Not open file \"%s\", using value from def.K[0]\n", file);
			std::cout << err;
			//print_error(error, __FILE__, __LINE__);

			for (int i = 0; i < def.locNx; i++)
				for (int j = 0; j < def.locNy; j++)
					for (int k = 0; k < def.locNz; k++)
						K[i + j * def.locNx + k * def.locNx * def.locNy] = def.K[0];
			return;

		}
	}

	int Nx, Ny, Nz;
	if (def.Nz < 2)
	{
		fscanf(input, "%d %d\n", &Nx, &Ny);
		Nz=1;
	}
	else
		fscanf(input, "%d %d %d\n", &Nx, &Ny, &Nz);

	if ((Nx != def.Nx) || (Ny != def.Ny) || (Nz != def.Nz))
	{
		printf("Warning: Nx/Ny/Nz from noise.dat not equal\nError in file \"%s\" at line %d\n", __FILE__, __LINE__);
		fflush(stdout);
	}

	/*
	// Версия для считывания файлов SPE-10 
	double s[6];
	long int row = 0, bigN = 0;
	int index = 0;
	int six = 6;

	while (row * six < def.Nx * (def.Ny) * (def.Nz))
	{
		fscanf(input, "%lf %lf %lf %lf %lf %lf\n", s, s + 1, s + 2, s + 3, s + 4, s + 5);

		for (index = 0; index < six; index++)
		{
			bigN = six * row + index;
			int i = bigN % def.Nx;
			int k = (bigN / def.Nx) % def.Nz;
			int j = bigN / (def.Nz * (def.Nx));

			//K[i + j * def.Nx + k * def.Nx * def.Ny] = 6.64e-11+ s[index] * 10e-13;
			//K[i + j * def.Nx + k * def.Nx * def.Ny] = s[index];
			K[i+j*def.locNx+k*def.locNx*def.locNy]=1e-10 * exp(s[index]);
		}
		row++;
	}

	fclose(input);
	*/
/*
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
			{
				std::cout << "m[" << i << ","<< j << "," << k <<"] = " << m[i + j * def.locNx + k * def.locNx * def.locNy] << "\n";
			}
*/

	// версия для считывания файлов от Антона
	char* str=new char[30*Nx];

	char value[30];
	for(int j=0; j<Ny; j++)
	{
		int n=0;
		fgets(str,30*Nx,input);
		for(int i=0; i<Nx; i++)
		{
			int iter=0;
			if(str[n]==' ')
				n++;
			for (n;str[n]!=' ';n++,iter++)
			{
				value[iter]=str[n];
			}
			value[iter]='\0';
			n++;

			for (int k=0;k<def.locNz;k++)
				K[i+j*def.locNx+k*def.locNx*def.locNy]=1e-10 * exp(atof(value)) * pow(def.upscale_l, 2);
		}
}

fclose(input);
#endif

for (int i = 0; i < def.locNx; i++)
	for (int j = 0; j < def.locNy; j++)
		for (int k = 0; k < def.locNz; k++)
			test_positive(K[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);

}

// Сохранение состояния в файл
void save(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, long int time_counter, const consts &def)
{
	// Загружаем в память хоста данные по roS_old
	// P1 и S2 загружены уже для функции сохранения графиков,
	// x,y и porosity не изменились и загрузки не требуют.
	//load_data_to_host(HostArraysPtr.P1, DevArraysPtr.P1 , localNx);
	//load_data_to_host(HostArraysPtr.S2, DevArraysPtr.S2 , localNx);
	load_data_to_host(HostArraysPtr.roS_w_old, DevArraysPtr.roS_w_old , def);
	load_data_to_host(HostArraysPtr.roS_n_old, DevArraysPtr.roS_n_old , def);
#ifdef THREE_PHASE
	load_data_to_host(HostArraysPtr.roS_g_old, DevArraysPtr.roS_n_old , def);
#endif

	FILE *f_save;

	if ((def.rank) == 0)
	{

#ifdef _WIN32
		_mkdir("save");
#else
		mkdir("save", 0000777);
#endif

		if (!(f_save = fopen("save/save.dat", "wb")))
		{
			printf("\nError: Not open file \"save.dat\"!\n");
			exit(0);
		}
		fclose(f_save);
	}

	for (int cpu = 0; cpu < (def.sizex * (def.sizey) * (def.sizez)); cpu++)
	{
		// Реализация фунции Barrier для различных коммуникаций
		barrier();
		if ((def.rank) == cpu)
		{
			if (!(f_save = fopen("save/save.dat", "ab")))
				print_error("Not open file \"save.dat\"", __FILE__, __LINE__);

			fwrite(&time_counter, sizeof(int), 1, f_save);
#ifdef THREE_PHASE
			fwrite(HostArraysPtr.S_w, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
#endif
			fwrite(HostArraysPtr.P_w, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.S_n, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			//fwrite(HostArraysPtr.x, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			//fwrite(HostArraysPtr.y, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			//fwrite(HostArraysPtr.z, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.roS_w_old, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.roS_n_old, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
#ifdef THREE_PHASE
			fwrite(HostArraysPtr.roS_g_old, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
#endif
			fwrite(HostArraysPtr.m, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fclose(f_save);
		}
	}
}

// Восстановление состояния из файла
void restore(ptr_Arrays HostArraysPtr, long int* time_counter, const consts &def)
{
	FILE *f_save;
	for (int cpu = 0; cpu < (def.sizex * (def.sizey) * (def.sizez)); cpu++)
	{
		// Реализация фунции Barrier для различных коммуникаций
		barrier();

		consts def_tmp;
		def_tmp.locNx = 0;
		def_tmp.locNy = 0;
		def_tmp.locNz = 0;
		def_tmp.Nx = def.Nx;
		def_tmp.Ny = def.Ny;
		def_tmp.Nz = def.Nz;
		def_tmp.size = def.size;
		def_tmp.sizex = def.sizex;
		def_tmp.sizey = def.sizey;
		def_tmp.sizez = def.sizez;

		if ((def.rank) == cpu)
		{
			if (!(f_save = fopen("save/save.dat", "rb")))
				print_error("Not open file \"save.dat\"", __FILE__, __LINE__);

			for (int queue = 0; queue <= (def.rank); queue++)
			{
				def_tmp.rank = queue;
				global_to_local_vars(&def_tmp);
				fread(time_counter, sizeof(int), 1, f_save);
#ifdef THREE_PHASE
				fread(HostArraysPtr.S_w, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
#endif
				fread(HostArraysPtr.P_w, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.S_n, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				//fread(HostArraysPtr.x, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				//fread(HostArraysPtr.y, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				//fread(HostArraysPtr.z, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.roS_w_old, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.roS_n_old, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
#ifdef THREE_PHASE
				fread(HostArraysPtr.roS_g_old, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
#endif
				fread(HostArraysPtr.m, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
			}
			fclose(f_save);
		}
	}
}



//------------------------------------------------------------------------------------------

// Собирает и печатает версию запускаемой программы
void read_version(void)
{

	FILE *rev;
	char str[250] = "";
	int revision;

	if (!(rev = fopen("../.svn/entries", "rt")))
	{
		revision = 0;
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			fgets(str, 250, rev);
		}
		revision = atoi(str);
	}

	printf("Version %s.%d compiled %s %s.\n\n", VERSION, revision, __DATE__, __TIME__);
}

// Масштабирование размерностей входных данных
void resize_defines(consts* def, double l, double t)
{
	// h_i - максимум из hx, hy, hz
	double h_i = max(def->hx, def->hy);
	if (def->Nz > 1)
		h_i = max(h_i, def->hz);

	l=l/h_i;
	t=t/def->dt;
	double m=1.;

	def->upscale_l=l;
	def->upscale_t=t;

	def->dt *= t;
	def->hx *= l;
	def->hy *= l;
	def->hz *= l;
	def->P_atm *= m/(l*l);
	def->Background_Pw *= m/(l*l);
	def->InjWell_Pw *= m/(l*l);
	def->OutWell_Pw *= m/(l*l);
	def->beta_w *= (l*l)/m;
	def->beta_n *= (l*l)/m;
	def->g_const *= l/(t*t);
	def->K[0] *= l*l;
	def->K[1] *= l*l;
	def->Q *= m/(t * l*l*l);
	def->mu_w *= m*t/(l*l);
	def->mu_n *= m*t/(l*l);
	def->ro0_n *= m/(pow(l,3));
	def->ro0_w *= m/(pow(l,3));
	def->l *= l;
	def->c_n *= l/t;
	def->c_w *= l/t;
	def->tau *= t;
	def->timeX *= t;

#ifdef THREE_PHASE
	def->beta_g *= (l*l)/m;
	def->ro0_g *= m/(pow(l,3));
	def->c_g *= l/t;
	def->mu_g *= m*t/(l*l);
#endif
}

#ifdef GTEST
// Тест функции resize_defines
TEST(Main,ResizeDefines)
{
	consts tdef;
	tdef.hx = 1.;
	tdef.hy = 1.;
	tdef.P_atm = 1.;
	tdef.Background_Pw = 1.;
	tdef.OutWell_Pw = 1.;
	tdef.InjWell_Pw = 1.;
	tdef.dt = 1.;
	tdef.beta_n = 1.;
	tdef.beta_w = 2.;
	tdef.g_const = 9.8;
	tdef.K[0] = 1;
	tdef.K[1] = 0.3;
	tdef.Q = 1.;
	tdef.mu_n = 1.;
	tdef.mu_w = 0.5;
	tdef.ro0_w = 1000.;
	tdef.ro0_n = 800.;
	tdef.l = 1e-5;
	tdef.c_n = 1.;
	tdef.c_w = 5.;
	tdef.tau = 1.;
	tdef.timeX = 100;

	resize_defines(&tdef, 0.1, 0.01);
	EXPECT_DOUBLE_EQ(0.1, tdef.hx);
	EXPECT_DOUBLE_EQ(0.1, tdef.hy);
	EXPECT_DOUBLE_EQ(0.01, tdef.dt);
	EXPECT_DOUBLE_EQ(100., tdef.P_atm);
	EXPECT_DOUBLE_EQ(100., tdef.Background_Pw);
	EXPECT_DOUBLE_EQ(100., tdef.OutWell_Pw);
	EXPECT_DOUBLE_EQ(100., tdef.InjWell_Pw);
	EXPECT_DOUBLE_EQ(0.01, tdef.beta_n);
	EXPECT_DOUBLE_EQ(0.02, tdef.beta_w);
	EXPECT_DOUBLE_EQ(9800, tdef.g_const);
	EXPECT_DOUBLE_EQ(1e-2, tdef.K[0]);
	EXPECT_DOUBLE_EQ(3e-3, tdef.K[1]);
	EXPECT_DOUBLE_EQ(1e5, tdef.Q);
	EXPECT_DOUBLE_EQ(1, tdef.mu_n);
	EXPECT_DOUBLE_EQ(0.5, tdef.mu_w);
	EXPECT_DOUBLE_EQ(1e6, tdef.ro0_w);
	EXPECT_DOUBLE_EQ(8e5, tdef.ro0_n);
	EXPECT_DOUBLE_EQ(1e-6, tdef.l);
	EXPECT_DOUBLE_EQ(50, tdef.c_w);
	EXPECT_DOUBLE_EQ(10, tdef.c_n);
	EXPECT_DOUBLE_EQ(0.01, tdef.tau);
	EXPECT_DOUBLE_EQ(1, tdef.timeX);
	EXPECT_DOUBLE_EQ(0.1, tdef.upscale_l);
	EXPECT_DOUBLE_EQ(0.01, tdef.upscale_t);
}
#endif


// Считывание параметров задачи из файла
void read_defines(int argc, char *argv[], consts* def)
{
	FILE *defs;
	char *file;
	char str[250] = "", attr_name[50] = "", attr_value[50] = "";

	file = DEFINES_FILE;

	if (!(defs = fopen(file, "rt")))
	{
		file = "defines.ini";
		if (!(defs = fopen(file, "rt")))
		{
			char error[30];
			sprintf(error, "Not open file \"%s\"", file);
			print_error("Not open file \"defines.ini\"", __FILE__, __LINE__);
		}
	}

	while (!feof(defs))
	{
		unsigned int i, j, a;
		fgets(str, 250, defs);
		if (str[0] == '#')
		{
			continue;
		}

		for (i = 0; str[i] != '='; i++)
		{
			if (i >= strlen(str))
			{
				continue;
			}
			attr_name[i] = str[i];
		}

		attr_name[i] = '\0';
		a = strlen(str);
		for (j = i + 1; str[j] != ' ' && (j < a); j++)
		{
			attr_value[j - i - 1] = str[j];
		}
		attr_value[j - i - 1] = '\0';

		//std::cout << str <<"\n";

		if (!strcmp(attr_name, "HX"))
		{
			def->hx = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "HY"))
		{
			def->hy = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "HZ"))
		{
			def->hz = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "TAU"))
		{
			def->tau = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "DT"))
		{
			def->dt = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_W"))
		{
			def->c_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_N"))
		{
			def->c_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "L"))
		{
			def->l = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "BETA_W"))
		{
			def->beta_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "BETA_N"))
		{
			def->beta_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "RO_W"))
		{
			def->ro0_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "RO_N"))
		{
			def->ro0_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "MU_W"))
		{
			def->mu_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "MU_N"))
		{
			def->mu_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "G_CONST"))
		{
			def->g_const = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "P_ATM"))
		{
			def->P_atm = atof(attr_value);
			continue;
		}

		if (!strcmp(attr_name, "LAMBDA_0"))
		{
			def->lambda[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "LAMBDA_1"))
		{
			def->lambda[1] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "M_0"))
		{
			def->porosity[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "M_1"))
		{
			def->porosity[1] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_WR_0"))
		{
			def->S_wr[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_WR_1"))
		{
			def->S_wr[1] = atof(attr_value);
			continue;
		}

		if (!strcmp(attr_name, "Q"))
		{
			def->Q = atof(attr_value);
			continue;
		}

		if (!strcmp(attr_name, "BACKGROUND_Pw"))
		{
			def->Background_Pw = atof(attr_value);
			continue;
		}

		if (!strcmp(attr_name, "BACKGROUND_Sn"))
		{
			def->Background_Sn = atof(attr_value);
			continue;
		}

		if (!strcmp(attr_name, "INJECTION_WELL_Pw"))
		{
			def->InjWell_Pw = atof(attr_value);
			continue;
		}

		if (!strcmp(attr_name, "INJECTION_WELL_Sn"))
		{
			def->InjWell_Sn = atof(attr_value);
			continue;
		}

		if (!strcmp(attr_name, "OUTPUT_WELL_Pw"))
		{
			def->OutWell_Pw = atof(attr_value);
			continue;
		}

		if (!strcmp(attr_name, "OUTPUT_WELL_Sn"))
		{
			def->OutWell_Sn = atof(attr_value);
			continue;
		}

		if (!strcmp(attr_name, "K_0"))
		{
			def->K[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "K_1"))
		{
			def->K[1] = atof(attr_value);
			continue;
		}

#ifdef TWO_PHASE
		if (!strcmp(attr_name, "P_D_0"))
		{
			def->P_d[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "P_D_1"))
		{
			def->P_d[1] = atof(attr_value);
			continue;
		}
#endif

#ifdef THREE_PHASE
		if (!strcmp(attr_name, "C_G"))
		{
			def->c_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "BETA_G"))
		{
			def->beta_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "RO_G"))
		{
			def->ro0_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "MU_G"))
		{
			def->mu_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "P_D_NW_0"))
		{
			def->P_d_nw[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "P_D_NW_1"))
		{
			def->P_d_nw[1] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "P_D_GN_0"))
		{
			def->P_d_gn[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "P_D_GN_1"))
		{
			def->P_d_gn[1] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_W_GR"))
		{
			def->S_w_gr = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_W_INIT"))
		{
			def->S_w_init = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_N_INIT"))
		{
			def->S_n_init = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_NR_0"))
		{
			def->S_nr[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_NR_1"))
		{
			def->S_nr[1] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_GR_0"))
		{
			def->S_gr[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_GR_1"))
		{
			def->S_gr[1] = atof(attr_value);
			continue;
		}
#endif
		if (!strcmp(attr_name, "S_N_GR"))
		{
			def->S_n_gr = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "SOURCE"))
		{
			def->source = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "ITERATIONS"))
		{
			def->newton_iterations = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "TIMEX"))
		{
			def->timeX = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "SAVE_PLOTS"))
		{
			def->save_plots = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "PRINT_SCREEN"))
		{
			def->print_screen = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "NX"))
		{
			def->Nx = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "NY"))
		{
			def->Ny = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "NZ"))
		{
			def->Nz = atoi(attr_value);
			continue;
		}
	}

	fclose(defs);

	// Если нет масштабирования, то 1
	def->upscale_l=1;
	def->upscale_t=1;

	read_defines_test(*def);

	//resize_defines(def, 0.15, 0.01);
	read_defines_test(*def);
}

// Вывод на экран распределения процессоров в системе и области по процессорам
void print_hosts_configuration(const consts &def) {
	printf("\n  size = %d : (%d, %d, %d)\n  rank = %d : (%d, %d, %d)\n  locN = %d : (%d, %d, %d)\n", 
		 (def).size, (def).sizex, (def).sizey, (def).sizez,
		 (def).rank, (def).rankx, (def).ranky, (def).rankz, 
		 (def).locNx * (def).locNy * (def).locNz, (def).locNx, (def).locNy, (def).locNz);
}

