/*
*  QTSimFast: a small persistent pool of worker threads.
*
*  This file is NOT part of the original QTSim. It was written in 2026-11 for
*  QTSimFast, a modified version of QTSim (QTLSim 1.0, copyright 2023 Andrew
*  Meade, University of Reading) that is faster and whose output is bit-identical
*  to the original's.
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
*
*  A pool runs one job at a time, in one of two ways:
*
*    FastThreadsRun(Pool, Fun, Ctx)    calls Fun(Ctx, t, T) for t = 0 .. T-1, t = 0
*                                      on the calling thread, and returns when all
*                                      have finished;
*    FastThreadsStart(Pool, Fun, Ctx)  calls Fun(Ctx, 1, T) on the first worker and
*    FastThreadsWait(Pool)             returns at once; Wait returns when it is done.
*
*  With one thread (or when built with -DQTSIM_NO_THREADS, which needs no thread
*  library) Run calls Fun(Ctx, 0, 1) and Start runs the job at once on the
*  calling thread.
*
*  Nothing here affects results: the callers give every thread its own part of
*  the output and never reduce across threads.
*/

#ifndef FAST_THREADS_H
#define FAST_THREADS_H

typedef void (*FAST_JOB)(void *Ctx, int ThreadNo, int NoThreads);

typedef struct FAST_POOL FAST_POOL;

// A pool of NoThreads threads, the caller included (NoThreads-1 workers). Fewer
// are used if they cannot be started; FastThreadsNo gives the number in use.
FAST_POOL*	FastThreadsCreate(int NoThreads);
int			FastThreadsNo(FAST_POOL *Pool);
void		FastThreadsRun(FAST_POOL *Pool, FAST_JOB Fun, void *Ctx);
void		FastThreadsStart(FAST_POOL *Pool, FAST_JOB Fun, void *Ctx);
void		FastThreadsWait(FAST_POOL *Pool);
void		FastThreadsFree(FAST_POOL *Pool);

// The number of processors available, 1 if unknown.
int			FastThreadsCPUs(void);

#endif
