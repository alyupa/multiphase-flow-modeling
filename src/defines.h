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
#define my_min(a,b) (((a) < (b)) ? (a) : (b))
#define my_max(a,b) (((a) > (b)) ? (a) : (b))

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

extern ptr_Arrays HostArraysPtr;
extern ptr_Arrays DevArraysPtrLoc[1];

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
#ifdef ENERGY
	double T_0;
	double ro_r;
	double lambda0_w, lambda0_n, lambda0_g, lambda0_r, lambdaA_w, lambdaA_n, lambdaA_g;
	double c0_w, c0_n, c0_g, c0_r, C_w, C_w2, C_n, C_g, C_r;
	double alfa_w, alfa_n;
#endif
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

extern consts def;

void time_step_function(double t);
void initialization(long int* time_counter, int argc, char* argv[]);
void load_permeability(double* K);
void finalization();
void memory_allocation();
void host_memory_allocation();
void device_memory_allocation();
void memory_free();
void host_memory_free();
void device_memory_free();
void save_data_plots(double t);
void data_initialization(long int* time_counter);
void sizes_initialization();
void blocks_initialization();
void communication_initialization(int argc, char* argv[]);
void communication_finalization(void);
void global_to_local_vars();
int local_to_global(int local_index, char axis);
int is_active_point(int i, int j, int k);
void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr);
void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr);
void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr);

void device_initialization();
void device_finalization(void);

void division();

// Служебные
void print_plots_top(double t);
void print_plots(double t, int I, int J);
void print_plot_row(FILE* fp, int i, int j, int k);
void print_hosts_configuration();
void print_array_console(double* Arr, char axis);
void barrier(void);
void read_defines(int argc, char *argv[]);
void print_error(const char *error, const char *file, int line);

// Unit-тесты
void test_correct_P_S();
void test_nan(double x, const char *file, int line);
void test_positive(double x, const char *file, int line);
void test_S(double S, const char *file, int line);
void test_u(double u, const char *file, int line);
void test_ro(double ro, const char *file, int line);
void test_arrowhead(double big, double small, const char *file, int line);
void read_defines_test();


// Расчеты в каждой точке
double left_difference (double* ptr, char axis);
double right_difference (double* ptr, char axis);
void prepare_local_vars(int i, int j, int k);
double ro_eff_gdy(int local);
void assign_P_Xi(int i, int j, int k);
void assign_ro(int local);
void assign_S(int local);
void assign_u(int i, int j, int k);
void assign_roS(double t, int i, int j, int k);
void assign_roS_nr(double t, int i, int j, int k);
void Newton(int i, int j, int k);
void Border_S(int i, int j, int k);
void Border_P(int i, int j, int k);
int set_boundary_basic_coordinate(int i, int j, int k);
int reverse_matrix (double *a, int n);
void mult_matrix_vector (double* result_vect, double* matr, double* vect, int n);

void prepare_all_vars();
void u_calculation();
void find_values_from_partial_equations(double t);
void solve_nonlinear_system();
void boundary_conditions();
void exchange_basic_vars();

void load_exchange_data_part_xl(double *Array);
void load_exchange_data_part_xr(double *Array);
void load_exchange_data_part_yl(double *Array);
void load_exchange_data_part_yr(double *Array);
void load_exchange_data_part_zl(double *Array);
void load_exchange_data_part_zr(double *Array);
void save_exchange_data_part_xl(double *Array);
void save_exchange_data_part_xr(double *Array);
void save_exchange_data_part_yl(double *Array);
void save_exchange_data_part_yr(double *Array);
void save_exchange_data_part_zl(double *Array);
void save_exchange_data_part_zr(double *Array);

int is_injection_well(int i, int j, int k);
int is_output_well(int i, int j, int k);
void wells_q(int i, int j, int k, double* q_w, double* q_n, double* q_g);
void assing_k(double* k_w, double* k_n, double S_w);

#ifdef ENERGY
void assign_H (int local);
void assign_E_current (int local);
void assign_E_new (int i, int j, int k);
void Border_T(int i, int j, int k);
#endif

#endif

