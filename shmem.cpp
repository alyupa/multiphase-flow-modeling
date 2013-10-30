#define MVSE
#ifdef K100
	#include <mpi.h>
	#include <shmem.h>
	#include <shmem-mpi.h>
#else
	#include <shmem++.h>
#endif
#include "defines.h"

#define MX 20000
static double *pf, *cf[MX];

// Обмен данными на границах между всеми процессорами
// 0. Загружаем данные с усорителя в память хоста
// 1.1 передаем правую границу, 
// 1.2 передаем левую границу.
// Для крайних процессоров соответствующие обмены не требуются
// 3. Загружаем полученные данные в память ускорителя
void exchange(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (0)

	if (def.rank!=def.size-1)
		shmem_double_put (cf[def.rank+1],HostBuffer+(def.Ny)*(def.Nz),sizeof(double)*(def.Ny)*(def.Nz),def.rank+1); // (1.1)
	if (def.rank!=0)
		shmem_double_put (cf[def.rank-1]+(def.Ny)*(def.Nz),HostBuffer,sizeof(double)*(def.Ny)*(def.Nz),def.rank-1); // (1.2)
	shmem_barrier_all();
	HostBuffer=cf[def.rank];

	save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def); // (3)
}

// Обмен граничными значениями давления P2, плотностей ro1 и ro2, Xi между процессорами
void P_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def)
{
#ifdef THREE_PHASE
	exchange(HostArraysPtr.ro_g, DevArraysPtr.ro_g, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.Xi_g, DevArraysPtr.Xi_g, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.P_w, DevArraysPtr.P_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.P_g, DevArraysPtr.P_g, HostBuffer, DevBuffer, def);
#else
	exchange(HostArraysPtr.P_n, DevArraysPtr.P_n, HostBuffer, DevBuffer, def);
#endif
	exchange(HostArraysPtr.ro_w, DevArraysPtr.ro_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.ro_n, DevArraysPtr.ro_n, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.Xi_w, DevArraysPtr.Xi_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.Xi_n, DevArraysPtr.Xi_n, HostBuffer, DevBuffer, def);
}


// Обмен граничными значениями скоростей между процессорами
// В случае распределения расчетной области между процессорами по оси X
// передача u1y и u2y не требуется
void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	exchange(HostArraysPtr.ux_w, DevArraysPtr.ux_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.ux_n, DevArraysPtr.ux_n, HostBuffer, DevBuffer, def);
#ifdef THREE_PHASE
	exchange(HostArraysPtr.ux_g, DevArraysPtr.ux_g, HostBuffer, DevBuffer, def);
#endif
}

// Обмен граничными значениями давления воды P1 и насыщенности NAPL S2 между процессорами
void P_S_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def)
{
#ifdef THREE_PHASE
	exchange(HostArraysPtr.P_n, DevArraysPtr.P_n, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.S_w, DevArraysPtr.S_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.S_g, DevArraysPtr.S_g, HostBuffer, DevBuffer, def);
#else
	exchange(HostArraysPtr.P_w, DevArraysPtr.P_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.S_n, DevArraysPtr.S_n, HostBuffer, DevBuffer, def);
#endif
}

void communication_initialization(int argc, char* argv[], consts* def)
{
	shmem_init(&argc, &argv);
	(*def).size=shmem_n_pes(); // The amount of processors
	(*def).rank=shmem_my_pe(); // The number of processor
	std::cout << "size=" <<(*def).size<<"  "<<"rank= "<<(*def).rank<<"\n";

	pf = (double*)emalloc( sizeof(double)*((*def).Ny)*((*def).Nz));
	shmem_coarray_all( (void*)pf, (long)((*def).Ny*((*def).Nz)*sizeof(*pf)), (void**)cf );
}

void communication_finalization(void)
{
	shmem_finalize();
}

// Реализация фунции Barrier для различных коммуникаций
void barrier(void)
{
	shmem_barrier_all();
}