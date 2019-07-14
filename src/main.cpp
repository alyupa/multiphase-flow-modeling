#include "defines.h"

double *HostBuffer;
double *DevBuffer;
consts def;
ptr_Arrays HostArraysPtr;
ptr_Arrays DevArraysPtrLoc[1];

int main(int argc, char* argv[])
{
	read_defines(argc, argv);
	// Счетчик шагов по времени
	long int time_counter = 0;
	// Счетчик времени исполнения вычислительной части программы
	clock_t task_time;

	// Инициализация коммуникаций, перевод глобальных параметров в локальные процессора,
	// выделение памяти, загрузка начальных/сохраненных данных
	initialization(&time_counter, argc, argv);

	// Тест
	//print_hosts_configuration();
	save_data_plots(0);

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

		time_step_function(time_counter * (def.dt)); // (1)

		if ((time_counter % (def.save_plots)) == 0) // (3)
		{
			// Следующие 2 функции вызываются строго в таком порядке,
			// т.к. save использует данные, загруженные save_data_plots
			save_data_plots(time_counter * (def.dt) /def.upscale_t); // (**)
			//save(time_counter); // (***)
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
	finalization();
	printf("\n...THE END...\n");

	// При запуске в Windows после работы программы оставлять окно консоли
#ifdef _WIN32
	printf("\nPress <Enter> to exit...\n");
	fflush(stdout);
	getchar();
#endif
	return 0;
}

// Функция полного цикла расчетов на следующем временном слое
// 1. Расчет плотности жидкостей, давлений, энтальпий фаз и суммарной внутренней энергии на текущем слое
// 3. Расчет скоростей жидкостей
// 5. Расчет переменной roS на следующем временном слое
// 6. Расчет методом Ньютона давления воды P1 и насыщенности DNAPL S2
// 7. Применение граничных условий для P1 и S2
// 8. Обмен между процессорами пограничными значениями P1 и S2
// 9.Вычисление суммарной внутренней энергии на следующем временном слое
void time_step_function(double t)
{
	boundary_conditions(); // (7)
	prepare_all_vars(); // (1)
	u_calculation(); // (3)
	find_values_from_partial_equations(t); // (5)
	solve_nonlinear_system(); // (6)
	exchange_basic_vars(); // (8)
}

// Преобразование локальных координат процессора к глобальным
// Каждый процессор содержит дополнительную точку в массиве для
// обмена данными, если имеет соседа
// (если 2 соседа с обеих сторон,то +2 точки).
// Глобальные границы хранятся как обычные точки (отсюда и условие на (def.rank)==0)
int local_to_global(int local_index, char axis)
{
	int global_index = local_index;
	switch (axis)
	{
	case 'x':
	{
		global_index += def.rankx * (def.Nx / def.sizex) + my_min(def.rankx, def.Nx % def.sizex) - my_min(def.rankx, 1);
		break;
	}
	case 'y':
	{
		global_index += def.ranky * (def.Ny / def.sizey) + my_min(def.ranky, def.Ny % def.sizey) - my_min(def.ranky, 1);
		break;
	}
	case 'z':
	{
		global_index += def.rankz * (def.Nz / def.sizez) + my_min(def.rankz, def.Nz % def.sizez) - my_min(def.rankz, 1);
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
void global_to_local_vars()
{
	def.locNx = def.Nx / def.sizex;

	if (def.rankx < def.Nx % def.sizex)
	{
		(def.locNx) ++;
	}

	// Крайние процессоры получают по 1 точке для граничных данных,
	// остальные - по 2 на обе границы
	// Если процессор один, то границ у него нет и дополнительные точки не нужны
	if (def.sizex > 1)
	{
		if ((def.rankx == 0) || (def.rankx  == def.sizex - 1))
		{
			(def.locNx) ++;
		}
		else
		{
			(def.locNx) += 2;
		}
	}

	def.locNy = def.Ny / def.sizey;

	if ((def.ranky < def.Ny % def.sizey))
	{
		(def.locNy) ++;
	}

	if (def.sizey > 1)
	{
		if ((def.ranky == 0) || (def.ranky == def.sizey - 1))
		{
			(def.locNy) ++;
		}
		else
		{
			(def.locNy) += 2;
		}
	}

	def.locNz = def.Nz / def.sizez;

	if (def.rankz < def.Nz % def.sizez)
	{
		(def.locNz) ++;
	}

	if (def.sizez > 1)
	{
		if ((def.rankz == 0) || (def.rankz == def.sizez - 1))
		{
			(def.locNz) ++;
		}
		else
		{
			(def.locNz) += 2;
		}
	}

	if(def.rank >= def.sizex * def.sizey * def.sizez)
	{
		def.locNx = def.locNy = def.locNz = 0;
	}

	test_positive(def.locNx, __FILE__, __LINE__);
	test_positive(def.locNy, __FILE__, __LINE__);
	test_positive(def.locNz, __FILE__, __LINE__);
}

// Является ли точка активной (т.е. не предназначенной только для обмена на границах)
int is_active_point(int i, int j, int k)
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
void sizes_initialization()
{
	division();

	def.rankx = def.rank % def.sizex;
	def.ranky = (def.rank / def.sizex) % def.sizey;
	def.rankz = (def.rank / def.sizex) / def.sizey;
}

void blocks_initialization()
{
	def.blocksX = 0;
	def.blocksY = 0;
	def.blocksZ = 0;
}

// Функция вычисления "эффективной" плотности * g * hy
double ro_eff_gdy(int local)
{
	double ro_g_dy = (HostArraysPtr.ro_g[local] * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local])
					+ HostArraysPtr.ro_w[local] * HostArraysPtr.S_w[local]
					+ HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local]) * (HostArraysPtr.m[local]) * (def.g_const) * (def.hy);
	return ro_g_dy;
}

//----------------------------------------------------------------------------------------------------
// Служебные функции

// Вывод запускаемой задачи
void print_task_name()
{
	// Нулевой процессор выводит название запускаемой задачи
	if (!(def.rank))
	{
		char task_name[] = "Three phase filtration";
		std::cout << task_name << " by CAPAZ on " << (def.size) << " node(s).\n";
		fflush(stdout);
	}
}

// Инициализация коммуникаций (1), перевод глобальных параметров в локальные процессора (2),
// инициализация ускорителя (2.5), выделение памяти (3), загрузка начальных/сохраненных данных (4)
// Для задачи Б-Л загрузка проницаемостей из файла.
void initialization(long int* time_counter, int argc, char* argv[])
{
	communication_initialization(argc, argv); // (1)

	print_task_name();

	sizes_initialization();

	blocks_initialization();

	global_to_local_vars(); // (2)

	device_initialization(); // (2.5)

	memory_allocation(); // (3)

	load_permeability(HostArraysPtr.K); // (5)

	data_initialization(time_counter);

	load_data_to_device(HostArraysPtr.P_w, DevArraysPtrLoc->P_w);
	load_data_to_device(HostArraysPtr.S_w, DevArraysPtrLoc->S_w);
	load_data_to_device(HostArraysPtr.S_n, DevArraysPtrLoc->S_n);
	load_data_to_device(HostArraysPtr.S_g, DevArraysPtrLoc->S_g);
	load_data_to_device(HostArraysPtr.roS_w_old, DevArraysPtrLoc->roS_w_old);
	load_data_to_device(HostArraysPtr.roS_n_old, DevArraysPtrLoc->roS_n_old);
	load_data_to_device(HostArraysPtr.roS_g_old, DevArraysPtrLoc->roS_g_old);
	load_data_to_device(HostArraysPtr.m, DevArraysPtrLoc->m);
	load_data_to_device(HostArraysPtr.K, DevArraysPtrLoc->K);
#ifdef ENERGY
	load_data_to_device(HostArraysPtr.T, DevArraysPtrLoc->T);
#endif
}

void finalization()
{
	memory_free();
	communication_finalization();
	device_finalization();
}

// Выделение памяти хоста (1) и ускорителя (2) под массив точек расчетной области
void memory_allocation()
{
	host_memory_allocation(); // (1)
	device_memory_allocation(); // (2)
}

// Освобожение памяти хоста (1) и ускорителя (2) из под массива точек расчетной области
void memory_free()
{
	host_memory_free(); // (1)
	device_memory_free(); // (2)
}

// Выделение памяти хоста под массив точек расчетной области
void host_memory_allocation()
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
		HostArraysPtr.S_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.S_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.P_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.P_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.ro_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.ro_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.ux_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.uy_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.uz_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.ux_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.uy_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.uz_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.Xi_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.Xi_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.roS_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.roS_w_old = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.roS_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.roS_n_old = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.m = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.K = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.P_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.S_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.ro_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.ux_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.uy_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.uz_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.Xi_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.roS_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.roS_g_old = new double [(def.locNx) * (def.locNy) * (def.locNz)];
#ifdef ENERGY
		HostArraysPtr.T = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.H_w = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.H_n = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.H_g = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.H_r = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.E = new double [(def.locNx) * (def.locNy) * (def.locNz)];
		HostArraysPtr.E_new = new double [(def.locNx) * (def.locNy) * (def.locNz)];
#endif
	}
	catch (...)
	{
		print_error("Not enough host memory for ArraysPtr", __FILE__, __LINE__);
	}
}

// Освобожение памяти хоста из под массива точек расчетной области
void host_memory_free()
{
	delete HostBuffer;
	delete[] HostArraysPtr.S_n;
	delete[] HostArraysPtr.S_w;
	delete[] HostArraysPtr.P_w;
	delete[] HostArraysPtr.P_n;
	delete[] HostArraysPtr.ro_w;
	delete[] HostArraysPtr.ro_n;
	delete[] HostArraysPtr.ux_w;
	delete[] HostArraysPtr.uy_w;
	delete[] HostArraysPtr.uz_w;
	delete[] HostArraysPtr.ux_n;
	delete[] HostArraysPtr.uy_n;
	delete[] HostArraysPtr.uz_n;
	delete[] HostArraysPtr.Xi_w;
	delete[] HostArraysPtr.Xi_n;
	delete[] HostArraysPtr.roS_w;
	delete[] HostArraysPtr.roS_w_old;
	delete[] HostArraysPtr.roS_n;
	delete[] HostArraysPtr.roS_n_old;
	delete[] HostArraysPtr.m;
	delete[] HostArraysPtr.K;
	delete[] HostArraysPtr.P_g;
	delete[] HostArraysPtr.S_g;
	delete[] HostArraysPtr.ro_g;
	delete[] HostArraysPtr.ux_g;
	delete[] HostArraysPtr.uy_g;
	delete[] HostArraysPtr.uz_g;
	delete[] HostArraysPtr.Xi_g;
	delete[] HostArraysPtr.roS_g;
	delete[] HostArraysPtr.roS_g_old;
#ifdef ENERGY
	delete[] HostArraysPtr.T;
	delete[] HostArraysPtr.H_w;
	delete[] HostArraysPtr.H_n;
	delete[] HostArraysPtr.H_g;
	delete[] HostArraysPtr.H_r;
	delete[] HostArraysPtr.E;
	delete[] HostArraysPtr.E_new;
#endif
}

// Функция сохранения графиков в файлы
void save_data_plots(double t)
{
	// Загрузка в память хоста результатов расчета
	load_data_to_host(HostArraysPtr.S_w, DevArraysPtrLoc->S_w);
	load_data_to_host(HostArraysPtr.ux_w, DevArraysPtrLoc->ux_w);
	load_data_to_host(HostArraysPtr.uy_w, DevArraysPtrLoc->uy_w);
	load_data_to_host(HostArraysPtr.uz_w, DevArraysPtrLoc->uz_w);
	load_data_to_host(HostArraysPtr.ux_g, DevArraysPtrLoc->ux_g);
	load_data_to_host(HostArraysPtr.uy_g, DevArraysPtrLoc->uy_g);
	load_data_to_host(HostArraysPtr.uz_g, DevArraysPtrLoc->uz_g);
	load_data_to_host(HostArraysPtr.P_w, DevArraysPtrLoc->P_w);
	load_data_to_host(HostArraysPtr.S_n, DevArraysPtrLoc->S_n);
	load_data_to_host(HostArraysPtr.S_g, DevArraysPtrLoc->S_g);
	load_data_to_host(HostArraysPtr.ux_n, DevArraysPtrLoc->ux_n);
	load_data_to_host(HostArraysPtr.uy_n, DevArraysPtrLoc->uy_n);
	load_data_to_host(HostArraysPtr.uz_n, DevArraysPtrLoc->uz_n);
#ifdef ENERGY
	load_data_to_host(HostArraysPtr.T, DevArraysPtrLoc->T);
#endif

	// Проверка на выход из допустимого диапазона значений P и S
#ifdef MY_TEST
	test_correct_P_S();
#endif

	// Нулевой процессор создает директории, файлы и прописывает заголовки файлов
	if ((def.rank) == 0)
	{
		print_plots_top(t);
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
				print_plots(t, -1, -1);
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
								print_plots(t, i, j);
							}
						}
	}

}

// Функция создания директорий, файлов для графиков и сохранения заголовков в них
void print_plots_top(double t)
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
void print_plots(double t, int I, int J)
{
	char fname[64];
	FILE *fp;

	sprintf(fname, "%s/F-%012.0f.tec", def.plots_dir, t);

	// Открытие на дозапись и сохранение графиков
	if (!(fp = fopen(fname, "at")))
		print_error("Not open file(s) in function SAVE_DATA_PLOTS", __FILE__, __LINE__);

	if(((def.sizey == 1) && (def.sizez == 1)) || ((I == -1) && (J == -1)))
	{
		for (int i = 0; i < def.locNx; i++)
			for (int j = 0; j < def.locNy; j++)
				for (int k = 0; k < def.locNz; k++)
					if (ACTIVE_POINT)
						print_plot_row(fp, i, j, k);
	}
	else
	{
		for (int k = 0; k < def.locNz; k++)
			if (is_active_point(I, J, k))
				print_plot_row(fp, I, J, k);
	}

	fclose(fp);
}

// Функция записи строчки значений параметров, которые нужны для построения графиков,  для всех задач.
// Включает в себя:
// 1. Распределение насыщенностей NAPL S_n
// 2. Распределение давлений воды P_w
// 3. Распределение скоростей {u_x, u_y, u_z}
// 4. Распределение типов грунтов
void print_plot_row(FILE* fp, int i, int j, int k)
{
	int local = i + j * def.locNx + k * def.locNx * def.locNy;

	// Преобразование локальных координат процессора к глобальным
	int I = local_to_global(i, 'x');
	int J = def.Ny - 1 - local_to_global(j, 'y');
	int K = local_to_global(k, 'z');

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
void print_array_console(double* Arr, char axis)
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
void load_permeability(double* K)
{
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++) {
				K[i + j * def.locNx + k * def.locNx * def.locNy] = def.K[0];
				test_positive(K[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
			}
}

// Считывание параметров задачи из файла
void read_defines(int argc, char *argv[])
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
			def.Nx = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "NY"))
		{
			def.Ny = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "NZ"))
		{
			def.Nz = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "TAU"))
		{
			def.tau = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "DT"))
		{
			def.dt = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_W"))
		{
			def.c_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_N"))
		{
			def.c_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "L"))
		{
			def.l = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "BETA_W"))
		{
			def.beta_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "BETA_N"))
		{
			def.beta_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "RO_W"))
		{
			def.ro0_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "RO_N"))
		{
			def.ro0_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "MU_W"))
		{
			def.mu_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "MU_N"))
		{
			def.mu_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "G_CONST"))
		{
			def.g_const = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "P_ATM"))
		{
			def.P_atm = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "LAMBDA_0"))
		{
			def.lambda[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "M_0"))
		{
			def.porosity[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_WR_0"))
		{
			def.S_wr[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "K_0"))
		{
			def.K[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_G"))
		{
			def.c_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "BETA_G"))
		{
			def.beta_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "RO_G"))
		{
			def.ro0_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "MU_G"))
		{
			def.mu_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "P_D_NW_0"))
		{
			def.P_d_nw[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "P_D_GN_0"))
		{
			def.P_d_gn[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_W_GR"))
		{
			def.S_w_gr = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_W_INIT"))
		{
			def.S_w_init = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_N_INIT"))
		{
			def.S_n_init = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_NR_0"))
		{
			def.S_nr[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_GR_0"))
		{
			def.S_gr[0] = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "S_N_GR"))
		{
			def.S_n_gr = atof(attr_value);
			continue;
		}
#ifdef ENERGY
		if (!strcmp(attr_name, "T_0"))
		{
			def.T_0 = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "RO_R"))
		{
			def.ro_r = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "LAMBDA_0_W"))
		{
			def.lambda0_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "LAMBDA_0_N"))
		{
			def.lambda0_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "LAMBDA_0_G"))
		{
			def.lambda0_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "LAMBDA_0_R"))
		{
			def.lambda0_r = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "LAMBDA_A_W"))
		{
			def.lambdaA_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "LAMBDA_A_N"))
		{
			def.lambdaA_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "LAMBDA_A_G"))
		{
			def.lambdaA_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_0_W"))
		{
			def.c0_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_0_N"))
		{
			def.c0_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_0_G"))
		{
			def.c0_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_0_R"))
		{
			def.c0_r = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_1_W"))
		{
			def.C_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_1_W2"))
		{
			def.C_w2 = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_1_N"))
		{
			def.C_n = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_1_G"))
		{
			def.C_g = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "C_1_R"))
		{
			def.C_r = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "ALFA_W"))
		{
			def.alfa_w = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "ALFA_N"))
		{
			def.alfa_n = atof(attr_value);
			continue;
		}
#endif
		if (!strcmp(attr_name, "SOURCE"))
		{
			def.source = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "ITERATIONS"))
		{
			def.newton_iterations = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "TIMEX"))
		{
			def.timeX = atof(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "SAVE_PLOTS"))
		{
			def.save_plots = atoi(attr_value);
			continue;
		}
		if (!strcmp(attr_name, "PLOTS_DIR"))
		{
			strcpy(def.plots_dir, attr_value);
			continue;
		}
		if (!strcmp(attr_name, "PRINT_SCREEN"))
		{
			def.print_screen = atoi(attr_value);
			continue;
		}
	}

	fclose(defs);

	// Determine hx, hy, hz from Nx, Ny, Nz
	def.hx = 1.0 / (my_max((double)def.Nx - 1, 1));
	def.hy = 1.0 / (my_max((double)def.Ny - 1, 1));
	def.hz = 1.0 / (my_max((double)def.Nz - 1, 1));

	// Если нет масштабирования, то 1
	def.upscale_l=1;
	def.upscale_t=1;

	read_defines_test();
}

// Вывод на экран распределения процессоров в системе и области по процессорам
void print_hosts_configuration() {
	printf("\n  size = %d : (%d, %d, %d)\n  rank = %d : (%d, %d, %d)\n  locN = %d : (%d, %d, %d)\n", 
		 def.size, def.sizex, def.sizey, def.sizez,
		 def.rank, def.rankx, def.ranky, def.rankz,
		 def.locNx * def.locNy * def.locNz, def.locNx, def.locNy, def.locNz);
}

static void S_local_initialization(int local)
{
	HostArraysPtr.S_w[local] = def.S_w_init;
	HostArraysPtr.S_n[local] = def.S_n_init;
	HostArraysPtr.S_g[local] = 1. - def.S_w_init - def.S_n_init;
//	HostArraysPtr.S_w[local] = def.S_w_init + 0.1 * cos(0.1 * local) + 0.1 / (local + 1.) + 0.1 * exp(-0.01 * local);
//	HostArraysPtr.S_n[local] = def.S_n_init + 0.1 * sin((double)local) - 0.1 / (local + 1.) - 0.1 * exp(-0.005 * local);;
}

void data_initialization(long int* t)
{
	*t = 0;
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
			{
				int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

				HostArraysPtr.m[local]=def.porosity[0];
				S_local_initialization(local);

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

				/*double ro_g_dy = ((def.ro0_g * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local])
				+ def.ro0_w * HostArraysPtr.S_w[local]
				+ def.ro0_n * HostArraysPtr.S_n[local]) * (HostArraysPtr.m[local]) + (1. - HostArraysPtr.m[local]) * 2000.) * (def.g_const) * (def.hy);*/

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
				HostArraysPtr.T[local] = 300;

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
