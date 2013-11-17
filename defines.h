#ifndef DEFINES_H
#define DEFINES_H

#define VERSION "1.0"

#define TWO_LAYERS 1
#define NR

// Учитывать ли тепловые процессы
// #define ENERGY

// Количество видеоускорителей на узле кластера
// Для К-100 - 3, для МВС-Экспресс 1 или 2
#define GPU_PER_NODE 3

// Нитей в блоке ускорителя
#define BlockNX 8
#define BlockNY 8
#define BlockNZ 8

// Размеры сетки из процессоров
#define SizeX 1
#define SizeY 1 
#define SizeZ 1

#include <float.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _WIN32
#define isnan _isnan
#include <direct.h>
#else
#include <sys/stat.h>
#endif

// Псевдо-функция минимума/максимума
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

// Вывод графиков BjnIO
//#include "bjnio.h"

/*!
* Структура "параметры расчетной точки"
*/
struct ptr_Arrays_tag
{
	double *S_n, *S_w, *P_w, *P_n, *ro_w, *ro_n, *ux_w, *uy_w, *uz_w, *ux_n, *uy_n, *uz_n, *Xi_w, *Xi_n, *roS_w, *roS_w_old, *roS_n, *roS_n_old;
	double *m, *K;
	double *S_g, *P_g, *ro_g, *ux_g, *uy_g, *uz_g, *Xi_g, *roS_g, *roS_g_old;
#ifdef ENERGY
	double *T, *H_w, *H_n, *H_g, *H_r, *E, *E_new;
#endif
};
typedef struct ptr_Arrays_tag ptr_Arrays;

/*! 
 *  \brief     Структура параметров сред и характерных размеров задачи.
 *  \details   Структура параметров сред и характерных размеров задачи.
 */
struct consts_tag
{
	double upscale_l, upscale_t;
	double lambda[2];
	double S_wr[2];
	double porosity[2];
	double hx, hy, hz, dt, tau, l, c_w, c_n, beta_w, beta_n, P_atm, g_const, mu_w, mu_n, ro0_w, ro0_n, timeX;
	double S_n_gr;
	int Nx, Ny, Nz;
	int source, save_plots, print_screen, newton_iterations;
	double K[2];
	double Q, InjWell_Pw, InjWell_Sn, OutWell_Pw, OutWell_Sn, Background_Sn, Background_Pw; // Дебит скважины
	double S_w_init, S_n_init;
	double S_nr[2];
	double S_gr[2];
	double P_d_nw[2];
	double P_d_gn[2];
	double c_g, beta_g, mu_g, ro0_g, S_w_gr;
	// Локальные размеры
	int locNx, locNy, locNz;
	// Число процессоров и ранг процессора
	int size, rank;
	// Число частей дробления области по измерениям
	int sizex, sizey, sizez;
	// Разложение ранга для трехмерного дробления области между процессорами
	int rankx, ranky, rankz;
	// Количество блоков ускорителя
	int blocksX, blocksY, blocksZ;
	// Plots directory name
	char plots_dir[32];
};
typedef struct consts_tag consts;

extern void time_step_function(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* DevBuffer, const consts &def, double t);
extern void initialization(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, long int* time_counter, int argc, char* argv[], consts* def);
extern void load_permeability(double* K, const consts &def);
extern void finalization(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* DevBuffer);
extern void memory_allocation(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, const consts &def);
extern void host_memory_allocation(ptr_Arrays* ArraysPtr, const consts &def);
extern void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, const consts &def);
extern void memory_free(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr);
extern void host_memory_free(const ptr_Arrays &HostArraysPtr);
extern void device_memory_free(ptr_Arrays DevArraysPtr, double* DevBuffer);
extern void save_data_plots(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double t, const consts &def);
extern void data_initialization(const ptr_Arrays &HostArraysPtr, long int* time_counter, const consts &def);
extern void sizes_initialization(consts* def);
extern void blocks_initialization(consts* def);
extern void communication_initialization(int argc, char* argv[], consts* def);
extern void communication_finalization(void);
extern void global_to_local_vars(consts* def);
extern int local_to_global(int local_index, char axis, const consts &def);
extern int is_active_point(int i, int j, int k, const consts &def);
extern void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, const consts &def);
extern void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, const consts &def);
extern void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, const consts &def);

extern void device_initialization(consts* def);
extern void device_finalization(void);

void division(consts* def);

// Служебные
extern void print_plots_top(double t, const consts &def);
extern void print_plots(const ptr_Arrays &HostArraysPtr, double t, const consts &def, int I, int J);
extern void print_plot_row(const ptr_Arrays &HostArraysPtr, FILE* fp, int i, int j, int k, const consts &def);
extern void print_hosts_configuration(const consts &def);
extern void print_array_console(double* Arr, const consts &def, char axis);
extern void barrier(void);
extern void restore(const ptr_Arrays &HostArraysPtr, long int* time_counter, const consts &def);
extern void save(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, long int time_counter, const consts &def);
extern void read_defines(int argc, char *argv[], consts* def);
extern void read_version(void);
extern void print_error(const char *error, const char *file, int line);

// Unit-тесты
extern void test_correct_P_S(const ptr_Arrays &HostArraysPtr, const consts &def);
extern void test_nan(double x, const char *file, int line);
extern void test_positive(double x, const char *file, int line);
extern void test_S(double S, const char *file, int line);
extern void test_u(double u, const char *file, int line);
extern void test_ro(double ro, const char *file, int line);
extern void test_arrowhead(double big, double small, const char *file, int line);
extern void test_tau(double S_old, double S_now, double S_new, int local, const consts &def, const char *file, int line);
extern void read_defines_test(const consts &def);


// Расчеты в каждой точке
extern double ro_eff_gdy(const ptr_Arrays &HostArraysPtr, int local, const consts &def);
extern void assign_P_Xi(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern void assign_ro(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern void assign_S(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern void assign_u(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern void assign_roS(const ptr_Arrays &HostArraysPtr, double t, int i, int j, int k, const consts &def);
extern void assign_roS_nr(const ptr_Arrays &HostArraysPtr, double t, int i, int j, int k, const consts &def);
extern void Newton(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern void Border_S(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern void Border_P(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern int set_boundary_basic_coordinate(int i, int j, int k, const consts &def);
extern int reverse_matrix (double *a, int n);
extern void mult_matrix_vector (double* result_vect, double* matr, double* vect, int n); 

extern void ro_P_Xi_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
extern void P_ro_Xi_exchange(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def);
extern void u_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
extern void u_exchange(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def);
extern void S_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
extern void roS_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double t, const consts &def);
extern void P_S_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
extern void boundary_conditions(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
extern void P_S_exchange(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def);

extern void load_exchange_data_part_xl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void load_exchange_data_part_xr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void load_exchange_data_part_yl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void load_exchange_data_part_yr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void load_exchange_data_part_zl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void load_exchange_data_part_zr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void save_exchange_data_part_xl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void save_exchange_data_part_xr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void save_exchange_data_part_yl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void save_exchange_data_part_yr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void save_exchange_data_part_zl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
extern void save_exchange_data_part_zr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);

extern int is_injection_well(int i, int j, int k, const consts &def);
extern int is_output_well(int i, int j, int k, const consts &def);
extern void wells_q(const ptr_Arrays &HostArraysPtr, int i, int j, int k, double* q_w, double* q_n, double* q_g, const consts &def);
extern void assing_k(double* k_w, double* k_n, double S_w);

#ifdef ENERGY
extern void assign_H (const ptr_Arrays &HostArraysPtr, int local, const consts &def);
extern double assign_T_flow (const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern double assign_E_flow (const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern void assign_E_current (const ptr_Arrays &HostArraysPtr, int local, const consts &def);
extern void assign_E_new (const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern void H_E_current_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
extern void Border_T(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
extern void E_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);

#endif

#endif

