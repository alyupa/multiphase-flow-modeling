#ifndef MEASURING_H
#define MEASURING_H

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
#include <sys/timeb.h>
#endif

#define MEASURE_COUNT 30
#define MAX_BUFFER_SIZE 1000000
#define MIN_BUFFER_SIZE 100000

//#define USE_GPU
#define GPU_PER_NODE 3

// Нитей в блоке ускорителя
#define BlockNX 8
#define BlockNY 8
#define BlockNZ 8

void exchange(double* HostBuffer, int buffer_size, int size, int rank);
void right_send_recv(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id);
void left_recv_send(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id);

#ifdef USE_GPU
void measuring_host_device_exchange(double *HostBuffer, double *DevBuffer, int rank);
void device_memory_allocation(double** DevBuffer, int buffer_size);
void device_initialization(int rank);
void device_finalization(void);
void load_exchange_data_part(double* HostBuffer, double* DevBuffer, int buffer_size);
void save_exchange_data_part(double* HostBuffer, double* DevBuffer, int buffer_size);
void device_memory_free(double* DevBuffer);
void set_devbuffer_values(double *DevBuffer, int buffer_size);
#endif

#endif
