/*
*  QTSimFast: a small persistent pool of worker threads, channels between
*  threads, and a plain thread (POSIX threads, C11 atomics).
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
void		FastThreadsSetSpin(FAST_POOL *Pool, int Spin) { (void)Pool; (void)Spin; }
int			FastThreadsCPUs(void) { return 1; }
int			FastThreadsPerfCPUs(void) { return 0; }
long		FastThreadsPollsPerMicro(void) { return 1; }

// Without threads nothing waits on a channel: it is a number.
struct FAST_CHAN { long Value; };

FAST_CHAN* FastChanCreate(long Value)
{
	FAST_CHAN *Chan;

	Chan = (FAST_CHAN*)calloc(1, sizeof(FAST_CHAN));
	if(Chan == NULL)
	{
		printf("Memory allocation error\n");
		exit(1);
	}
	Chan->Value = Value;

	return Chan;
}

void		FastChanSet(FAST_CHAN *Chan, long Value) { Chan->Value = Value; }
long		FastChanGet(FAST_CHAN *Chan) { return Chan->Value; }
long		FastChanWait(FAST_CHAN *Chan, long Seen, int Spin) { (void)Seen; (void)Spin; return Chan->Value; }
long		FastChanClaim(FAST_CHAN *Chan, long Step) { long Value = Chan->Value; Chan->Value += Step; return Value; }
void		FastChanFree(FAST_CHAN *Chan) { free(Chan); }
FAST_THREAD* FastThreadSpawn(void (*Fun)(void *Ctx), void *Ctx) { (void)Fun; (void)Ctx; return NULL; }
void		FastThreadJoin(FAST_THREAD *Thread) { (void)Thread; }

#else

#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <time.h>

#ifdef __APPLE__
	#include <sys/types.h>
	#include <sys/sysctl.h>
#endif

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

	int				WorkerSpin;	// polls before a worker sleeps (set while no job runs)
	int				CallerSpin;	// polls before the thread that waits for a job sleeps
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
	int Spin, Mine, MaxSpin;

	W = (WORKER*)Arg;
	Pool = W->Pool;
	Seen = 0;
	MaxSpin = WORKER_SPIN;

	for(;;)
	{
		for(Spin=0;Spin<MaxSpin;Spin++)
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
		MaxSpin = Pool->WorkerSpin;
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

long FastThreadsPollsPerMicro(void)
{
	struct timespec A, B;
	atomic_long Value;
	double Nano, Best;
	long Sum;
	int Trial, Index;
	const int Polls = 2000;

	atomic_init(&Value, 0);
	Best = -1;
	Sum = 0;

	// (the fastest of a few trials: a trial that was interrupted counts for nothing)
	for(Trial=0;Trial<5;Trial++)
	{
		if(clock_gettime(CLOCK_MONOTONIC, &A) != 0)
			return 1;

		for(Index=0;Index<Polls;Index++)
		{
			Sum += atomic_load(&Value);
			CPU_RELAX();
		}

		if(clock_gettime(CLOCK_MONOTONIC, &B) != 0)
			return 1;

		Nano = (double)(B.tv_sec - A.tv_sec) * 1.0e9 + (double)(B.tv_nsec - A.tv_nsec);
		if(Nano > 0 && (Best < 0 || Nano < Best))
			Best = Nano;
	}

	if(Best <= 0 || Sum != 0)
		return 1;

	Nano = Polls * 1000.0 / Best;
	if(Nano < 1)
		return 1;
	if(Nano > 1.0e5)
		return 100000;

	return (long)Nano;
}

int FastThreadsPerfCPUs(void)
{
#ifdef __APPLE__
	int No;
	size_t Size;

	No = 0;
	Size = sizeof(No);
	if(sysctlbyname("hw.perflevel0.physicalcpu", &No, &Size, NULL, 0) == 0 && Size == sizeof(No) && No >= 1)
		return No;
#endif

	return 0;
}

// A channel. The value and the count of sleepers are sequentially consistent
// atomics: the setter stores the value and then reads the count, a waiter
// raises the count and then reads the value, so one of them sees the other,
// and a waiter cannot sleep through a Set. The padding keeps two channels,
// written by different threads, off one cache line.
struct FAST_CHAN
{
	char			Pad0[128];
	atomic_long		Value;
	atomic_int		Sleeping;
	pthread_mutex_t	Lock;
	pthread_cond_t	Cond;
	char			Pad1[128];
};

FAST_CHAN* FastChanCreate(long Value)
{
	FAST_CHAN *Chan;

	Chan = (FAST_CHAN*)calloc(1, sizeof(FAST_CHAN));
	if(Chan == NULL)
	{
		printf("Memory allocation error\n");
		exit(1);
	}

	atomic_init(&Chan->Value, Value);
	atomic_init(&Chan->Sleeping, 0);
	pthread_mutex_init(&Chan->Lock, NULL);
	pthread_cond_init(&Chan->Cond, NULL);

	return Chan;
}

void FastChanSet(FAST_CHAN *Chan, long Value)
{
	atomic_store(&Chan->Value, Value);

	if(atomic_load(&Chan->Sleeping) != 0)
	{
		pthread_mutex_lock(&Chan->Lock);
		pthread_cond_broadcast(&Chan->Cond);
		pthread_mutex_unlock(&Chan->Lock);
	}
}

long FastChanGet(FAST_CHAN *Chan)
{
	return atomic_load(&Chan->Value);
}

long FastChanWait(FAST_CHAN *Chan, long Seen, int Spin)
{
	long Value;
	int Index;

	for(Index=0;Index<Spin;Index++)
	{
		Value = atomic_load(&Chan->Value);
		if(Value != Seen)
			return Value;
		CPU_RELAX();
	}

	pthread_mutex_lock(&Chan->Lock);
	atomic_fetch_add(&Chan->Sleeping, 1);
	while((Value = atomic_load(&Chan->Value)) == Seen)
		pthread_cond_wait(&Chan->Cond, &Chan->Lock);
	atomic_fetch_sub(&Chan->Sleeping, 1);
	pthread_mutex_unlock(&Chan->Lock);

	return Value;
}

long FastChanClaim(FAST_CHAN *Chan, long Step)
{
	return atomic_fetch_add(&Chan->Value, Step);
}

void FastChanFree(FAST_CHAN *Chan)
{
	pthread_mutex_destroy(&Chan->Lock);
	pthread_cond_destroy(&Chan->Cond);
	free(Chan);
}

struct FAST_THREAD
{
	pthread_t	Thread;
	void		(*Fun)(void *Ctx);
	void		*Ctx;
};

static void* ThreadMain(void *Arg)
{
	FAST_THREAD *Thread;

	Thread = (FAST_THREAD*)Arg;
	Thread->Fun(Thread->Ctx);

	return NULL;
}

FAST_THREAD* FastThreadSpawn(void (*Fun)(void *Ctx), void *Ctx)
{
	FAST_THREAD *Thread;

	Thread = (FAST_THREAD*)calloc(1, sizeof(FAST_THREAD));
	if(Thread == NULL)
		return NULL;

	Thread->Fun = Fun;
	Thread->Ctx = Ctx;

	if(pthread_create(&Thread->Thread, NULL, ThreadMain, Thread) != 0)
	{
		free(Thread);
		return NULL;
	}

	return Thread;
}

void FastThreadJoin(FAST_THREAD *Thread)
{
	if(Thread == NULL)
		return;

	pthread_join(Thread->Thread, NULL);
	free(Thread);
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
	Pool->WorkerSpin = WORKER_SPIN;
	Pool->CallerSpin = CALLER_SPIN;
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

void FastThreadsSetSpin(FAST_POOL *Pool, int Spin)
{
	if(Pool->NoThreads == 1)
		return;

	// (workers read it with the next job, under the lock)
	pthread_mutex_lock(&Pool->Lock);
	Pool->WorkerSpin = Spin;
	Pool->CallerSpin = (Spin > CALLER_SPIN) ? Spin : CALLER_SPIN;
	pthread_mutex_unlock(&Pool->Lock);
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

	// (CallerSpin is set before the pool's first job and not afterwards)
	for(Spin=0;Spin<Pool->CallerSpin;Spin++)
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
