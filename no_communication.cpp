#include "defines.h"

void division(consts *def)
{
			(*def).sizex = 1;
			(*def).sizey = 1;
			(*def).sizez = 1;
}

void exchange_basic_vars(const ptr_Arrays &HostArraysPtr, const ptr_Arrays &DevArraysPtr, double* HostBuffer, double* DevBuffer, const consts &def)
{
}

void communication_initialization(int argc, char* argv[], consts *def)
{
	(*def).size = 1;
	(*def).rank = 0;
}

void communication_finalization(void)
{
}

// Реализация фунции Barrier для различных коммуникаций
void barrier(void)
{
}

