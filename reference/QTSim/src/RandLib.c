/*
*  QTLSim 1.0
*
*  copyright 2023
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* GenomeSim is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _WIN32
	#include <windows.h>
#else
	#include <unistd.h>
#endif

#include "RandLib.h"

// Get a process ID
unsigned long	GetProcID(void)
{
	#ifdef _WIN32
		return GetCurrentProcessId();
	#else
		return getpid();
	#endif
}

// Reverse the bits in a long
unsigned long ReverseUSLong(unsigned long x)
{
	unsigned long h;
	int i;

	h = 0;

	for(i = 0; i < 32; i++)
	{
		h = (h << 1) + (x & 1);
		x >>= 1;
	}

	return h;
}

// Get a random seed based on the time and process ID
unsigned long	GetSeed(void)
{
	unsigned long Seed;
	unsigned long Pid;

	Seed = (unsigned)time(NULL);
	Pid = GetProcID();
	Pid = ReverseUSLong(Pid);

	return Seed ^ Pid;
}


