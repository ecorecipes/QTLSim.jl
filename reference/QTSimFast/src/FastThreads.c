/*
*  QTSimFast: a small persistent pool of worker threads (POSIX threads).
*
*  This file is NOT part of the original QTSim. It was written in 2026-11 for
*  QTSimFast, a modified version of QTSim (QTLSim 1.0, copyright 2023 Andrew
*  Meade, University of Reading). See FastThreads.h.
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#include <stdio.h>
#include <stdlib.h>

#include "FastThreads.h"

#ifdef QTSIM_NO_THREADS

struct FAST_POOL { int Unused; };
static struct FAST_POOL Serial;

FAST_POOL*	FastThreadsCreate(int NoThreads) { (void)NoThreads; return &Serial; }
int			FastThreadsNo(FAST_POOL *Pool) { (void)Pool; return 1; }
void		FastThreadsRun(FAST_POOL *Pool, FAST_JOB Fun, void *Ctx) { (void)Pool; Fun(Ctx, 0, 1); }
void		FastThreadsStart(FAST_POOL *Pool, FAST_JOB Fun, void *Ctx) { (void)Pool; Fun(Ctx, 0, 1); }
void		FastThreadsWait(FAST_POOL *Pool) { (void)Pool; }
void		FastThreadsFree(FAST_POOL *Pool) { (void)Pool; }
int			FastThreadsCPUs(void) { return 1; }

#else

#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>

// How long a thread polls before it sleeps on a condition variable.
#define WORKER_SPIN 4000
#define CALLER_SPIN 20000

#if defined(__x86_64__) || defined(__i386__)
	#define CPU_RELAX() __builtin_ia32_pause()
#elif defined(__aarch64__) || defined(__arm__)
	#define CPU_RELAX() __asm__ __volatile__("yield")
#else
	#define CPU_RELAX() ((void)0)
#endif

typedef struct
{
	pthread_t	Thread;
	int			No;
	FAST_POOL	*Pool;
} WORKER;

struct FAST_POOL
{
	int				NoThreads;
	WORKER			*Workers;

	pthread_mutex_t	Lock;
	pthread_cond_t	WorkCond;
	pthread_cond_t	DoneCond;

	atomic_ulong	Epoch;		// incremented for each job
	atomic_int		Pending;	// workers still to finish the current job
	atomic_int		Quit;

	// The current job. Written and read under Lock.
	unsigned long	JobEpoch;
	FAST_JOB		Fun;
	void			*Ctx;
	int				OnlyWorker;	// -1: every worker takes part; else only this one
};

static void* WorkerMain(void *Arg)
{
	WORKER *W;
	FAST_POOL *Pool;
	FAST_JOB Fun;
	void *Ctx;
	unsigned long Seen;
	int Spin, Mine;

	W = (WORKER*)Arg;
	Pool = W->Pool;
	Seen = 0;

	for(;;)
	{
		for(Spin=0;Spin<WORKER_SPIN;Spin++)
		{
			if(atomic_load(&Pool->Epoch) != Seen || atomic_load(&Pool->Quit) != 0)
				break;
			CPU_RELAX();
		}

		pthread_mutex_lock(&Pool->Lock);
		while(Pool->JobEpoch == Seen && atomic_load(&Pool->Quit) == 0)
			pthread_cond_wait(&Pool->WorkCond, &Pool->Lock);

		if(Pool->JobEpoch == Seen)
		{
			pthread_mutex_unlock(&Pool->Lock);
			return NULL;	// Quit
		}

		Seen = Pool->JobEpoch;
		Fun = Pool->Fun;
		Ctx = Pool->Ctx;
		Mine = (Pool->OnlyWorker == -1 || Pool->OnlyWorker == W->No);
		pthread_mutex_unlock(&Pool->Lock);

		if(Mine)
		{
			Fun(Ctx, W->No, Pool->NoThreads);

			if(atomic_fetch_sub(&Pool->Pending, 1) == 1)
			{
				pthread_mutex_lock(&Pool->Lock);
				pthread_cond_signal(&Pool->DoneCond);
				pthread_mutex_unlock(&Pool->Lock);
			}
		}
	}
}

int FastThreadsCPUs(void)
{
	long No;

	No = sysconf(_SC_NPROCESSORS_ONLN);
	if(No < 1)
		return 1;

	return (int)No;
}

FAST_POOL* FastThreadsCreate(int NoThreads)
{
	FAST_POOL *Pool;
	int Index;

	if(NoThreads < 1)
		NoThreads = 1;

	Pool = (FAST_POOL*)calloc(1, sizeof(FAST_POOL));
	if(Pool == NULL)
	{
		printf("Memory allocation error\n");
		exit(1);
	}

	Pool->NoThreads = 1;
	if(NoThreads == 1)
		return Pool;

	pthread_mutex_init(&Pool->Lock, NULL);
	pthread_cond_init(&Pool->WorkCond, NULL);
	pthread_cond_init(&Pool->DoneCond, NULL);
	atomic_init(&Pool->Epoch, 0);
	atomic_init(&Pool->Pending, 0);
	atomic_init(&Pool->Quit, 0);

	Pool->Workers = (WORKER*)calloc(NoThreads, sizeof(WORKER));
	if(Pool->Workers == NULL)
		return Pool;

	// NoThreads is fixed before the first worker starts; if a worker cannot be
	// started, only the workers before it are kept.
	for(Index=1;Index<NoThreads;Index++)
	{
		Pool->Workers[Index].No = Index;
		Pool->Workers[Index].Pool = Pool;
	}

	Pool->NoThreads = NoThreads;
	for(Index=1;Index<NoThreads;Index++)
		if(pthread_create(&Pool->Workers[Index].Thread, NULL, WorkerMain, &Pool->Workers[Index]) != 0)
			break;

	if(Index != NoThreads)
	{
		// Stop the ones that started and run serially.
		pthread_mutex_lock(&Pool->Lock);
		atomic_store(&Pool->Quit, 1);
		pthread_cond_broadcast(&Pool->WorkCond);
		pthread_mutex_unlock(&Pool->Lock);

		while(--Index >= 1)
			pthread_join(Pool->Workers[Index].Thread, NULL);

		Pool->NoThreads = 1;
	}

	return Pool;
}

int FastThreadsNo(FAST_POOL *Pool)
{
	return Pool->NoThreads;
}

static void Post(FAST_POOL *Pool, FAST_JOB Fun, void *Ctx, int OnlyWorker, int Pending)
{
	atomic_store(&Pool->Pending, Pending);

	pthread_mutex_lock(&Pool->Lock);
	Pool->Fun = Fun;
	Pool->Ctx = Ctx;
	Pool->OnlyWorker = OnlyWorker;
	Pool->JobEpoch++;
	atomic_store(&Pool->Epoch, Pool->JobEpoch);
	pthread_cond_broadcast(&Pool->WorkCond);
	pthread_mutex_unlock(&Pool->Lock);
}

void FastThreadsWait(FAST_POOL *Pool)
{
	int Spin;

	if(Pool->NoThreads == 1)
		return;

	for(Spin=0;Spin<CALLER_SPIN;Spin++)
	{
		if(atomic_load(&Pool->Pending) == 0)
			return;
		CPU_RELAX();
	}

	pthread_mutex_lock(&Pool->Lock);
	while(atomic_load(&Pool->Pending) != 0)
		pthread_cond_wait(&Pool->DoneCond, &Pool->Lock);
	pthread_mutex_unlock(&Pool->Lock);
}

void FastThreadsRun(FAST_POOL *Pool, FAST_JOB Fun, void *Ctx)
{
	if(Pool->NoThreads == 1)
	{
		Fun(Ctx, 0, 1);
		return;
	}

	Post(Pool, Fun, Ctx, -1, Pool->NoThreads - 1);
	Fun(Ctx, 0, Pool->NoThreads);
	FastThreadsWait(Pool);
}

void FastThreadsStart(FAST_POOL *Pool, FAST_JOB Fun, void *Ctx)
{
	if(Pool->NoThreads == 1)
	{
		Fun(Ctx, 0, 1);
		return;
	}

	Post(Pool, Fun, Ctx, 1, 1);
}

void FastThreadsFree(FAST_POOL *Pool)
{
	int Index;

	if(Pool->NoThreads > 1)
	{
		pthread_mutex_lock(&Pool->Lock);
		atomic_store(&Pool->Quit, 1);
		pthread_cond_broadcast(&Pool->WorkCond);
		pthread_mutex_unlock(&Pool->Lock);

		for(Index=1;Index<Pool->NoThreads;Index++)
			pthread_join(Pool->Workers[Index].Thread, NULL);
	}

	free(Pool->Workers);
	free(Pool);
}

#endif
