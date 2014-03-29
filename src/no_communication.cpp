#include "defines.h"

void division()
{
			def.sizex = 1;
			def.sizey = 1;
			def.sizez = 1;
}

void exchange_basic_vars()
{
}

void communication_initialization(int argc, char* argv[])
{
	def.size = 1;
	def.rank = 0;
}

void communication_finalization(void)
{
}

// Реализация фунции Barrier для различных коммуникаций
void barrier(void)
{
}

