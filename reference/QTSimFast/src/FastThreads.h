/*
*  QTSimFast: a small persistent pool of worker threads, and channels between
*  threads.
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
*
*  Added 2026-11-19 for the pipeline of QTSim.c (one thread draws the random
*  numbers, another builds the genomes behind it): a channel, which is a number
*  that one thread sets and other threads read or wait on, and a plain thread.
*
*    FastChanSet(Chan, Value)          publishes Value; everything the setting
*                                      thread wrote before is visible to a thread
*                                      that then reads Value from the channel;
*    FastChanGet(Chan)                 the value last published;
*    FastChanWait(Chan, Seen, Spin)    returns the value once it is not Seen: polls
*                                      Spin times, then sleeps until the next Set.
*
*    FastChanClaim(Chan, Step)         adds Step and returns the value before: of
*                                      several threads that claim, each gets its own.
*
*  FastThreadSpawn returns NULL when built with -DQTSIM_NO_THREADS or when no
*  thread can be started.
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

// How many times a worker polls for the next job before it sleeps, and the
// thread that waits for a job to end at least as many (the defaults are a few
// microseconds' worth). A pool whose jobs follow one another closely answers
// faster with more, at the price of processor time between jobs. To be called
// before the pool's first job.
void		FastThreadsSetSpin(FAST_POOL *Pool, int Spin);

// The number of processors available, 1 if unknown.
int			FastThreadsCPUs(void);

// About how many polls of a waiting thread fit in a microsecond on this
// machine (at least 1), measured when called: for choosing how long to poll.
// The answer varies with the machine's load and affects nothing but waiting.
long		FastThreadsPollsPerMicro(void);

// The number of performance cores where processors are of two kinds and the
// system says so (macOS: hw.perflevel0.physicalcpu), 0 if unknown.
int			FastThreadsPerfCPUs(void);

typedef struct FAST_CHAN FAST_CHAN;
typedef struct FAST_THREAD FAST_THREAD;

FAST_CHAN*	FastChanCreate(long Value);
void		FastChanSet(FAST_CHAN *Chan, long Value);
long		FastChanGet(FAST_CHAN *Chan);
long		FastChanWait(FAST_CHAN *Chan, long Seen, int Spin);
long		FastChanClaim(FAST_CHAN *Chan, long Step);
void		FastChanFree(FAST_CHAN *Chan);

FAST_THREAD* FastThreadSpawn(void (*Fun)(void *Ctx), void *Ctx);
void		FastThreadJoin(FAST_THREAD *Thread);

#endif
