#ifndef DEFINES_H
#define DEFINES_H

#define VERSION "1.0"

#define TWO_LAYERS 1
#define NR

// Учитывать ли тепловые процессы
// #define ENERGY

// Количество видеоускорителей на узле кластера
// Для К-100 - 3
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

// Point is active
#define ACTIVE_POINT !( ((((def.rankx) != 0 && i == 0) || ((def.rankx) != (def.sizex) - 1 && i == (def.locNx) - 1)) && (def.Nx) >= 2) \
					 || ((((def.ranky) != 0 && j == 0) || ((def.ranky) != (def.sizey) - 1 && j == (def.locNy) - 1)) && (def.Ny) >= 2) \
					 || ((((def.rankz) != 0 && k == 0) || ((def.rankz) != (def.sizez) - 1 && k == (def.locNz) - 1)) && (def.Nz) >= 2))

// Point is internal condition
#define INTERNAL_POINT ((((i != 0) && (i != (def.locNx) - 1)) || ((def.locNx) < 2)) && (j != 0) && (j != (def.locNy) - 1) \
		&& (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))

// Point is on boundary
#define BOUNDARY_POINT ((((i == 0) || (i == (def.locNx) - 1)) && ((def.locNx) >= 2)) || (j == 0) || (j == (def.locNy) - 1) \
		|| (((k == 0) || (k == (def.locNz) - 1)) && ((def.locNz) >= 2)))

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
	double lambda[1];
	double S_wr[1];
	double porosity[1];
	double hx, hy, hz, dt, tau, l, c_w, c_n, beta_w, beta_n, P_atm, g_const, mu_w, mu_n, ro0_w, ro0_n, timeX;
	double S_n_gr;
	int Nx, Ny, Nz;
	int source, save_plots, print_screen, newton_iterations;
	double K[1];
	double S_w_init, S_n_init;
	double S_nr[1];
	double S_gr[1];
	double P_d_nw[1];
	double P_d_gn[1];
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

void time_step_function(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* DevBuffer, const consts &def, double t);
void initialization(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, long int* time_counter, int argc, char* argv[], consts* def);
void load_permeability(double* K, const consts &def);
void finalization(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* DevBuffer);
void memory_allocation(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, const consts &def);
void host_memory_allocation(ptr_Arrays* ArraysPtr, const consts &def);
void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, const consts &def);
void memory_free(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr);
void host_memory_free(const ptr_Arrays &HostArraysPtr);
void device_memory_free(ptr_Arrays DevArraysPtr, double* DevBuffer);
void save_data_plots(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double t, const consts &def);
void data_initialization(const ptr_Arrays &HostArraysPtr, long int* time_counter, const consts &def);
void sizes_initialization(consts* def);
void blocks_initialization(consts* def);
void communication_initialization(int argc, char* argv[], consts* def);
void communication_finalization(void);
void global_to_local_vars(consts* def);
int local_to_global(int local_index, char axis, const consts &def);
int is_active_point(int i, int j, int k, const consts &def);
void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, const consts &def);
void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, const consts &def);
void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, const consts &def);

void device_initialization(consts* def);
void device_finalization(void);

void division(consts* def);

// Служебные
void print_plots_top(double t, const consts &def);
void print_plots(const ptr_Arrays &HostArraysPtr, double t, const consts &def, int I, int J);
void print_plot_row(const ptr_Arrays &HostArraysPtr, FILE* fp, int i, int j, int k, const consts &def);
void print_hosts_configuration(const consts &def);
void print_array_console(double* Arr, const consts &def, char axis);
void barrier(void);
void restore(const ptr_Arrays &HostArraysPtr, long int* time_counter, const consts &def);
void save(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, long int time_counter, const consts &def);
void read_defines(int argc, char *argv[], consts* def);
void read_version(void);
void print_error(const char *error, const char *file, int line);

// Unit-тесты
void test_correct_P_S(const ptr_Arrays &HostArraysPtr, const consts &def);
void test_nan(double x, const char *file, int line);
void test_positive(double x, const char *file, int line);
void test_S(double S, const char *file, int line);
void test_u(double u, const char *file, int line);
void test_ro(double ro, const char *file, int line);
void test_arrowhead(double big, double small, const char *file, int line);
void test_tau(double S_old, double S_now, double S_new, int local, const consts &def, const char *file, int line);
void read_defines_test(const consts &def);


// Расчеты в каждой точке
double ro_eff_gdy(const ptr_Arrays &HostArraysPtr, int local, const consts &def);
void assign_P_Xi(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
void assign_ro(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
void assign_S(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
void assign_u(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
void assign_roS(const ptr_Arrays &HostArraysPtr, double t, int i, int j, int k, const consts &def);
void assign_roS_nr(const ptr_Arrays &HostArraysPtr, double t, int i, int j, int k, const consts &def);
void Newton(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
void Border_S(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
void Border_P(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
int set_boundary_basic_coordinate(int i, int j, int k, const consts &def);
int reverse_matrix (double *a, int n);
void mult_matrix_vector (double* result_vect, double* matr, double* vect, int n);

void ro_P_Xi_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
void P_ro_Xi_exchange(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def);
void u_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
void u_exchange(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def);
void S_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
void roS_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double t, const consts &def);
void P_S_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
void boundary_conditions(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
void P_S_exchange(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def);

void load_exchange_data_part_xl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void load_exchange_data_part_xr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void load_exchange_data_part_yl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void load_exchange_data_part_yr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void load_exchange_data_part_zl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void load_exchange_data_part_zr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void save_exchange_data_part_xl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void save_exchange_data_part_xr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void save_exchange_data_part_yl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void save_exchange_data_part_yr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void save_exchange_data_part_zl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);
void save_exchange_data_part_zr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, const consts &def);

int is_injection_well(int i, int j, int k, const consts &def);
int is_output_well(int i, int j, int k, const consts &def);
void wells_q(const ptr_Arrays &HostArraysPtr, int i, int j, int k, double* q_w, double* q_n, double* q_g, const consts &def);
void assing_k(double* k_w, double* k_n, double S_w);

#ifdef ENERGY
void assign_H (const ptr_Arrays &HostArraysPtr, int local, const consts &def);
double assign_T_flow (const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
double assign_E_flow (const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
void assign_E_current (const ptr_Arrays &HostArraysPtr, int local, const consts &def);
void assign_E_new (const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
void H_E_current_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);
void Border_T(const ptr_Arrays &HostArraysPtr, int i, int j, int k, const consts &def);
void E_calculation(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, const consts &def);

#endif

#endif

