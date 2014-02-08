#include "defines.h"

// Буферные массивы для обмена между процессорами
double *HostBuffer;
double *DevBuffer;

int main(int argc, char* argv[])
{
	consts def;
	read_defines(argc, argv, &def);

	// Хостовый массив данных расчетной области процессора
	ptr_Arrays HostArraysPtr;
	// GPU-массив данных расчетной области процессора
	ptr_Arrays DevArraysPtr;
	// Счетчик шагов по времени
	long int time_counter = 0;
	// Счетчик времени исполнения вычислительной части программы
	clock_t task_time;

	// Инициализация коммуникаций, перевод глобальных параметров в локальные процессора,
	// выделение памяти, загрузка начальных/сохраненных данных
	initialization(&HostArraysPtr, &DevArraysPtr, &time_counter, argc, argv, &def);

	// Тест
	//print_hosts_configuration(def);
	save_data_plots(HostArraysPtr, DevArraysPtr, 0, def);

	task_time = clock();

	/*!
	* Цикл шагов по времени (каждая итерация - новый слой по времени):
	* -# Проводятся расчеты давлений P и насыщенностей S на следующем временном слое
	* -# Каждые (def.print_screen) раз на экран выводится информация о временном слое
	* -# Каждые save_plots раз данные выгружаются в память хоста и
	*    сохраняются в файлы графиков (**), в файл сохраняется состояние задачи (***) 
	*/
	for (time_counter++; time_counter <= def.timeX / (def.dt); time_counter++)
	{
		if ((time_counter % (def.print_screen) == 0) && (def.rank) == 0) // (2)
		{
			printf("t=%.3f\n", time_counter * (def.dt) /def.upscale_t);
			fflush(stdout);
		}

		time_step_function(HostArraysPtr, DevArraysPtr, DevBuffer, def, time_counter * (def.dt)); // (1)

		if ((time_counter % (def.save_plots)) == 0) // (3)
		{
			// Следующие 2 функции вызываются строго в таком порядке,
			// т.к. save использует данные, загруженные save_data_plots
			save_data_plots(HostArraysPtr, DevArraysPtr, time_counter * (def.dt) /def.upscale_t, def); // (**)
			//save(HostArraysPtr, DevArraysPtr, time_counter, def); // (***)
		}
	}

	// Ждем пока все процессоры закончат расчеты
	barrier();
	// Вывод информации о времени работы программы в секундах
	task_time = clock() - task_time;
	if (!(def.rank))
	{
		printf("Task time in seconds:\t%.2f\n", (double) task_time / CLOCKS_PER_SEC);
	}

	// Завершение работы и освобождение памяти
	finalization(HostArraysPtr, DevArraysPtr, DevBuffer);

	// При запуске в Windows после работы программы оставлять окно консоли
#ifdef _WIN32
	printf("\nPress <Enter> to exit...\n");
	fflush(stdout);
	getchar();
#endif
	return 0;
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
// 9. Вычисление и запись вспомогательной насыщенности(S_w - для двух фаз, S_g - для трех фаз) 
// 10.Вычисление энтальпий фаз и суммарной внутренней энергии на текущем слое
// 11.Вычисление суммарной внутренней энергии на следующем временном слое
void time_step_function(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* DevBuffer, const consts &def, double t)
{
	P_S_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, def); // (8)
	S_calculation(HostArraysPtr, DevArraysPtr, def); // (9)
	ro_P_Xi_calculation(HostArraysPtr, DevArraysPtr, def); // (1)
	P_ro_Xi_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, def); // (2)
	u_calculation(HostArraysPtr, DevArraysPtr, def); // (3)
	u_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, def); // (4)
	roS_calculation(HostArraysPtr, DevArraysPtr, t, def); // (5)
#ifdef ENERGY
	H_E_current_calculation(HostArraysPtr, DevArraysPtr, def); // (10)
	E_calculation(HostArraysPtr, DevArraysPtr, def); // (11)
#endif
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
	if (   ((((def.rankx) != 0 && i == 0) || ((def.rankx) != (def.sizex) - 1 && i == (def.locNx) - 1)) && (def.Nx) >= 2)
	    || ((((def.ranky) != 0 && j == 0) || ((def.ranky) != (def.sizey) - 1 && j == (def.locNy) - 1)) && (def.Ny) >= 2)
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
	division(def);

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
double ro_eff_gdy(const ptr_Arrays &HostArraysPtr, int local, const consts &def)
{
	double ro_g_dy = (HostArraysPtr.ro_g[local] * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local])
					+ HostArraysPtr.ro_w[local] * HostArraysPtr.S_w[local]
					+ HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local]) * (HostArraysPtr.m[local]) * (def.g_const) * (def.hy);
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
		char task_name[] = "Three phase filtration";
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

	print_task_name(*def);

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

	load_data_to_device((*HostArraysPtr).S_w, (*DevArraysPtr).S_w, *def);
	load_data_to_device((*HostArraysPtr).roS_g_old, (*DevArraysPtr).roS_g_old, *def);
	load_data_to_device((*HostArraysPtr).P_n, (*DevArraysPtr).P_n, *def);
	load_data_to_device((*HostArraysPtr).P_g, (*DevArraysPtr).P_g, *def);
	load_data_to_device((*HostArraysPtr).P_w, (*DevArraysPtr).P_w, *def);
	load_data_to_device((*HostArraysPtr).S_n, (*DevArraysPtr).S_n, *def);
	load_data_to_device((*HostArraysPtr).roS_w_old, (*DevArraysPtr).roS_w_old, *def);
	load_data_to_device((*HostArraysPtr).roS_n_old, (*DevArraysPtr).roS_n_old, *def);
	load_data_to_device((*HostArraysPtr).m, (*DevArraysPtr).m, *def);
}

// Завершение работы (1), освобождение памяти (2)
void finalization(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* DevBuffer)
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
void memory_free(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr)
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
		(*ArraysPtr).S_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
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
		(*ArraysPtr).P_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).S_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).ro_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).ux_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).uy_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).uz_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).Xi_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).roS_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).roS_g_old = new double [(def.locNx) * (def.locNy) * (def.locNz)];
#ifdef ENERGY
		(*ArraysPtr).T = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).H_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).H_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).H_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).H_r = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).E = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		(*ArraysPtr).E_new = new double [(def.locNx) * (def.locNy) * (def.locNz)];
#endif
	}
	catch (...)
	{
		print_error("Not enough host memory for ArraysPtr", __FILE__, __LINE__);
	}
}

// Освобожение памяти хоста из под массива точек расчетной области
void host_memory_free(const ptr_Arrays &ArraysPtr)
{
	delete HostBuffer;
	delete[] ArraysPtr.S_n;
	delete[] ArraysPtr.S_w;
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
	delete[] ArraysPtr.P_g;
	delete[] ArraysPtr.S_g;
	delete[] ArraysPtr.ro_g;
	delete[] ArraysPtr.ux_g;
	delete[] ArraysPtr.uy_g;
	delete[] ArraysPtr.uz_g;
	delete[] ArraysPtr.Xi_g;
	delete[] ArraysPtr.roS_g;
	delete[] ArraysPtr.roS_g_old;
#ifdef ENERGY
	delete[] ArraysPtr.T;
	delete[] ArraysPtr.H_w;
	delete[] ArraysPtr.H_n;
	delete[] ArraysPtr.H_g;
	delete[] ArraysPtr.H_r;
	delete[] ArraysPtr.E;
	delete[] ArraysPtr.E_new;
#endif
}

// Функция сохранения графиков в файлы
void save_data_plots(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double t, const consts &def)
{
	// Загрузка в память хоста результатов расчета
	load_data_to_host(HostArraysPtr.S_w, DevArraysPtr.S_w , def);
	load_data_to_host(HostArraysPtr.ux_w, DevArraysPtr.ux_w , def);
	load_data_to_host(HostArraysPtr.uy_w, DevArraysPtr.uy_w , def);
	load_data_to_host(HostArraysPtr.uz_w, DevArraysPtr.uz_w , def);
	load_data_to_host(HostArraysPtr.ux_g, DevArraysPtr.ux_g , def);
	load_data_to_host(HostArraysPtr.uy_g, DevArraysPtr.uy_g , def);
	load_data_to_host(HostArraysPtr.uz_g, DevArraysPtr.uz_g , def);
	load_data_to_host(HostArraysPtr.P_w, DevArraysPtr.P_w , def);
	load_data_to_host(HostArraysPtr.S_n, DevArraysPtr.S_n , def);
	load_data_to_host(HostArraysPtr.ux_n, DevArraysPtr.ux_n , def);
	load_data_to_host(HostArraysPtr.uy_n, DevArraysPtr.uy_n , def);
	load_data_to_host(HostArraysPtr.uz_n, DevArraysPtr.uz_n , def);

	// Проверка на выход из допустимого диапазона значений P и S
#ifdef MY_TEST
	test_correct_P_S(HostArraysPtr, def);
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
	char fname[64];
	FILE *fp;

#ifdef _WIN32
	_mkdir(def.plots_dir);
#else
	mkdir(def.plots_dir, 0000777);
#endif

	sprintf(fname, "%s/F-%012.0f.tec", def.plots_dir, t);

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
		if ((def.Nx) < 2)
		{
#ifdef ENERGY
			fprintf(fp, "VARIABLES = \"Y\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"T\",\"uw_y\",\"un_y\",\"ug_y\",\"porosity\" \n");
#else
			fprintf(fp, "VARIABLES = \"Y\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"uw_y\",\"un_y\",\"ug_y\",\"porosity\" \n");
#endif
			fprintf(fp, "ZONE T = \"BIG ZONE\", K=%d, F = POINT\n", (def.Ny));
		}
		else 
		{

#ifdef ENERGY
			fprintf(fp, "VARIABLES = \"X\",\"Y\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"T\",\"uw_x\",\"uw_y\",\"un_x\",\"un_y\",\"ug_x\",\"ug_y\",\"porosity\" \n");
#else
			fprintf(fp, "VARIABLES = \"X\",\"Y\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"uw_x\",\"uw_y\",\"un_x\",\"un_y\",\"ug_x\",\"ug_y\",\"porosity\" \n");
#endif
			fprintf(fp, "ZONE T = \"BIG ZONE\", K=%d,J=%d, F = POINT\n", (def.Nx), (def.Ny));
		}
	}
	else
	{
#ifdef ENERGY
		fprintf(fp, "VARIABLES = \"X\",\"Y\",\"Z\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"T\",\"uw_x\",\"uw_y\",\"uw_z\",\"un_x\",\"un_y\",\"un_z\",\"ug_x\",\"ug_y\",\"ug_z\",\"porosity\" \n");
#else
		fprintf(fp, "VARIABLES = \"X\",\"Y\",\"Z\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"uw_x\",\"uw_y\",\"uw_z\",\"un_x\",\"un_y\",\"un_z\",\"ug_x\",\"ug_y\",\"ug_z\",\"porosity\" \n");
#endif
		fprintf(fp, "ZONE T = \"BIG ZONE\", K=%d,J=%d,I=%d, F = POINT\n", (def.Nx), (def.Ny), (def.Nz));
	}

	fclose(fp);
}

// Функция сохранения данных в файлы графиков
// Если деления между процессорами по y, z нет, то параметры I, J не используются. 
// Если же есть, но хочется писать в файл в простой последовательности, то должно быть I=J=-1
void print_plots(const ptr_Arrays &HostArraysPtr, double t, const consts &def, int I, int J)
{
	char fname[64];
	FILE *fp;

	sprintf(fname, "%s/F-%012.0f.tec", def.plots_dir, t);

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
}

// Функция записи строчки значений параметров, которые нужны для построения графиков,  для всех задач.
// Включает в себя:
// 1. Распределение насыщенностей NAPL S_n
// 2. Распределение давлений воды P_w
// 3. Распределение скоростей {u_x, u_y, u_z}
// 4. Распределение типов грунтов
void print_plot_row(const ptr_Arrays &HostArraysPtr, FILE* fp, int i, int j, int k, const consts &def)
{
	int local = i + j * def.locNx + k * def.locNx * def.locNy;

	// Преобразование локальных координат процессора к глобальным
	int I = local_to_global(i, 'x', def);
	int J = def.Ny - 1 - local_to_global(j, 'y', def);
	int K = local_to_global(k, 'z', def); 

	if (def.Nz < 2)
	{
		if (def.Nx < 2)	
		{
#ifdef ENERGY
			fprintf(fp, "%.2e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", J * (def.hy),
				HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], HostArraysPtr.S_g[local], HostArraysPtr.P_w[local], HostArraysPtr.T[local], 
				(-1)*HostArraysPtr.uy_w[local], (-1)*HostArraysPtr.uy_n[local], (-1)*HostArraysPtr.uy_g[local], HostArraysPtr.m[local]);
#else
			fprintf(fp, "%.2e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", J * (def.hy),
					HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], HostArraysPtr.S_g[local], HostArraysPtr.P_w[local], 
					(-1)*HostArraysPtr.uy_w[local], (-1)*HostArraysPtr.uy_n[local], (-1)*HostArraysPtr.uy_g[local], HostArraysPtr.m[local]);
#endif
		}
		else
		{
		/*						fprintf(fp,"%.2e %.2e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", I*(def.hx), J*(def.hy),
									HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local],
									HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.m[local]);
		*/
#ifdef ENERGY
			fprintf(fp, "%.2e %.2e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", I * (def.hx), J * (def.hy),
				HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], HostArraysPtr.S_g[local], HostArraysPtr.P_w[local], HostArraysPtr.T[local],
				HostArraysPtr.ux_w[local], (-1)*HostArraysPtr.uy_w[local], HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.ux_g[local],
				(-1)*HostArraysPtr.uy_g[local], HostArraysPtr.m[local]);
#else
			fprintf(fp, "%.2e %.2e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", I * (def.hx), J * (def.hy),
				HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], HostArraysPtr.S_g[local], HostArraysPtr.P_w[local],
				HostArraysPtr.ux_w[local], (-1)*HostArraysPtr.uy_w[local], HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.ux_g[local],
				(-1)*HostArraysPtr.uy_g[local], HostArraysPtr.m[local]);
#endif
		}
	}

	else
	{
		/*						fprintf(fp,"%.2e %.2e %.2e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", I*(def.hx), K*(def.hz), J*(def.hy),
									HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local],
									HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.m[local]);
		*/
#ifdef ENERGY
		fprintf(fp, "%.2e %.2e %.2e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", I * (def.hx), K * (def.hz), J * (def.hy),
			HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], HostArraysPtr.S_g[local], HostArraysPtr.P_w[local], HostArraysPtr.T[local],
			HostArraysPtr.ux_w[local], HostArraysPtr.uz_w[local], (-1)*HostArraysPtr.uy_w[local],
			HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local],
			HostArraysPtr.ux_g[local], HostArraysPtr.uz_g[local], (-1)*HostArraysPtr.uy_g[local], HostArraysPtr.m[local]);
#else
		fprintf(fp, "%.2e %.2e %.2e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", I * (def.hx), K * (def.hz), J * (def.hy),
				HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], HostArraysPtr.S_g[local], HostArraysPtr.P_w[local], 
				HostArraysPtr.ux_w[local], HostArraysPtr.uz_w[local], (-1)*HostArraysPtr.uy_w[local],
				HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local],
				HostArraysPtr.ux_g[local], HostArraysPtr.uz_g[local], (-1)*HostArraysPtr.uy_g[local], HostArraysPtr.m[local]);
#endif
	}
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
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				K[i + j * def.locNx + k * def.locNx * def.locNy] = def.K[0];

	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				test_positive(K[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
}

// Сохранение состояния в файл
void save(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, long int time_counter, const consts &def)
{
	load_data_to_host(HostArraysPtr.roS_w_old, DevArraysPtr.roS_w_old , def);
	load_data_to_host(HostArraysPtr.roS_n_old, DevArraysPtr.roS_n_old , def);
	load_data_to_host(HostArraysPtr.roS_g_old, DevArraysPtr.roS_n_old , def);

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
			fwrite(HostArraysPtr.S_w, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.P_w, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.S_n, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.roS_w_old, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.roS_n_old, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.roS_g_old, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.m, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fclose(f_save);
		}
	}
}

// Восстановление состояния из файла
void restore(const ptr_Arrays &HostArraysPtr, long int* time_counter, const consts &def)
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
				fread(HostArraysPtr.S_w, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.P_w, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.S_n, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.roS_w_old, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.roS_n_old, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.roS_g_old, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
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
	def->beta_w *= (l*l)/m;
	def->beta_n *= (l*l)/m;
	def->g_const *= l/(t*t);
	def->K[0] *= l*l;
	def->K[1] *= l*l;
	def->mu_w *= m*t/(l*l);
	def->mu_n *= m*t/(l*l);
	def->ro0_n *= m/(pow(l,3));
	def->ro0_w *= m/(pow(l,3));
	def->l *= l;
	def->c_n *= l/t;
	def->c_w *= l/t;
	def->tau *= t;
	def->timeX *= t;

	def->beta_g *= (l*l)/m;
	def->ro0_g *= m/(pow(l,3));
	def->c_g *= l/t;
	def->mu_g *= m*t/(l*l);
}

// Считывание параметров задачи из файла
void read_defines(int argc, char *argv[], consts* def)
{
	FILE *defs;
	char str[250] = "", attr_name[50] = "", attr_value[50] = "";
	// Settings file name
	char file[32];

	printf("%d: %s %s \n", argc, argv[0], argv[1]);
	if (argc > 1)
	{
		strcpy(file, argv[1]);
	} else {
		strcpy(file, "defines.ini");
	}

	if (!(defs = fopen(file, "rt")))
	{
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
		// remove \n symbol
		attr_value[j - i - 2] = '\0';

		//std::cout << str <<"\n";

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
		if (!strcmp(attr_name, "PLOTS_DIR"))
		{
			strcpy(def->plots_dir, attr_value);
			continue;
		}
		if (!strcmp(attr_name, "PRINT_SCREEN"))
		{
			def->print_screen = atoi(attr_value);
			continue;
		}
	}

	fclose(defs);

	// Determine hx, hy, hz from Nx, Ny, Nz
	def->hx = 1.0 / (double)def->Nx;
	def->hy = 1.0 / (double)def->Ny;
	def->hz = 1.0 / (double)def->Nz;

	// Если нет масштабирования, то 1
	def->upscale_l=1;
	def->upscale_t=1;

	read_defines_test(*def);

	//resize_defines(def, 0.15, 0.01);
	//read_defines_test(*def);
}

// Вывод на экран распределения процессоров в системе и области по процессорам
void print_hosts_configuration(const consts &def) {
	printf("\n  size = %d : (%d, %d, %d)\n  rank = %d : (%d, %d, %d)\n  locN = %d : (%d, %d, %d)\n", 
		 (def).size, (def).sizex, (def).sizey, (def).sizez,
		 (def).rank, (def).rankx, (def).ranky, (def).rankz, 
		 (def).locNx * (def).locNy * (def).locNz, (def).locNx, (def).locNy, (def).locNz);
}

static void S_local_initialization(const ptr_Arrays &HostArraysPtr, int local, const consts &def)
{
	HostArraysPtr.S_w[local] = def.S_w_init;
	HostArraysPtr.S_n[local] = def.S_n_init;
	HostArraysPtr.S_g[local] = 1. - def.S_w_init - def.S_n_init;
//	HostArraysPtr.S_w[local] = def.S_w_init + 0.1 * cos(0.1 * local) + 0.1 / (local + 1.) + 0.1 * exp(-0.01 * local);
//	HostArraysPtr.S_n[local] = def.S_n_init + 0.1 * sin((double)local) - 0.1 / (local + 1.) - 0.1 * exp(-0.005 * local);;
}

void data_initialization(const ptr_Arrays &HostArraysPtr, long int* t, const consts &def)
{
	*t = 0;
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				if (is_active_point(i, j, k, def))
				{
					int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

					HostArraysPtr.m[local]=def.porosity[0];
					S_local_initialization(HostArraysPtr, local, def);


					/*if ((j == 0) && ((def.source) > 0))
					{
						HostArraysPtr.S_w[local] = def.S_w_gr;
						HostArraysPtr.S_n[local] = def.S_n_gr;
					}
					else
					{
						HostArraysPtr.S_w[local] = def.S_w_init;
						HostArraysPtr.S_n[local] = def.S_n_init;
					}*/

					double ro_g_dy = ((def.ro0_g * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local])
					+ def.ro0_w * HostArraysPtr.S_w[local]
					+ def.ro0_n * HostArraysPtr.S_n[local]) * (HostArraysPtr.m[local]) + (1. - HostArraysPtr.m[local]) * 2000.) * (def.g_const) * (def.hy);

					// Если отдельно задаем значения на границах через градиент
					//if (j == 0)
					{
						HostArraysPtr.P_w[local] = def.P_atm;
						HostArraysPtr.P_n[local] = def.P_atm;
						HostArraysPtr.P_g[local] = def.P_atm;
					}
					/*else
					{
						HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local - (def.locNx)] + ro_g_dy;
						HostArraysPtr.P_n[local] = HostArraysPtr.P_n[local - (def.locNx)] + ro_g_dy;
						HostArraysPtr.P_g[local] = HostArraysPtr.P_g[local - (def.locNx)] + ro_g_dy;
					}*/

					HostArraysPtr.ro_w[local] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm));
					HostArraysPtr.ro_n[local] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[local] - def.P_atm));
					HostArraysPtr.ro_g[local] = def.ro0_g * HostArraysPtr.P_g[local] / def.P_atm;

#ifdef ENERGY
					// !!!! Нужно задать начальные распределения температуры, энтальпии, энергии!
					HostArraysPtr.T[local] = 285;

					test_positive(HostArraysPtr.T[local], __FILE__, __LINE__);
#endif
					test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);
					test_S(HostArraysPtr.S_w[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.P_n[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.P_g[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.ro_w[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.ro_n[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.ro_g[local], __FILE__, __LINE__);
				}
}
