/*
*  QTSimFast: an inlined MT19937 generator and the GSL draws built on it.
*
*  This file is NOT part of the original QTSim. It was written in 2026-11 for
*  QTSimFast, a modified version of QTSim (QTLSim 1.0, copyright 2023 Andrew
*  Meade, University of Reading). See FastRNG.h.
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
#include <string.h>

#include "FastRNG.h"

#define UPPER_MASK 0x80000000U
#define LOWER_MASK 0x7fffffffU
#define MAGIC(y) (((y) & 0x1U) ? 0x9908b0dfU : 0U)

// Advance the state by 624 steps and temper all 624 outputs.
void FastMTRefill(FAST_MT *S)
{
	uint32_t *mt, y;
	int kk;

	mt = S->mt;

	for(kk=0;kk<FMT_N-FMT_M;kk++)
	{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+FMT_M] ^ (y >> 1) ^ MAGIC(y);
	}

	for(;kk<FMT_N-1;kk++)
	{
		y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		mt[kk] = mt[kk+(FMT_M-FMT_N)] ^ (y >> 1) ^ MAGIC(y);
	}

	y = (mt[FMT_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
	mt[FMT_N-1] = mt[FMT_M-1] ^ (y >> 1) ^ MAGIC(y);

	for(kk=0;kk<FMT_N;kk++)
	{
		y = mt[kk];
		y ^= (y >> 11);
		y ^= (y << 7) & 0x9d2c5680U;
		y ^= (y << 15) & 0xefc60000U;
		y ^= (y >> 18);
		S->out[kk] = y;
	}

	S->mti = 0;
}

static unsigned long TypeGet(void *State)
{
	FAST_MT *S;

	S = (FAST_MT*)State;

	if(S->mti >= FMT_N)
		FastMTRefill(S);

	return S->out[S->mti++];
}

static double TypeGetDouble(void *State)
{
	return TypeGet(State) / 4294967296.0;
}

// GSL's mt_set
static void TypeSet(void *State, unsigned long Seed)
{
	FAST_MT *S;
	int i;

	S = (FAST_MT*)State;

	if(Seed == 0)
		Seed = 4357;

	S->mt[0] = (uint32_t)(Seed & 0xffffffffUL);

	for(i=1;i<FMT_N;i++)
		S->mt[i] = 1812433253U * (S->mt[i-1] ^ (S->mt[i-1] >> 30)) + (uint32_t)i;

	S->mti = FMT_N;
}

void FastRngSet(FAST_RNG *R, unsigned long Seed)
{
#ifdef QTSIM_GSL_RNG
	gsl_rng *G;

	G = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(G, Seed);
	R->GRng = *G;
	return;
#endif

	R->GType.name = "qtsimfast_mt19937";
	R->GType.max = 0xffffffffUL;
	R->GType.min = 0;
	R->GType.size = sizeof(FAST_MT);
	R->GType.set = TypeSet;
	R->GType.get = TypeGet;
	R->GType.get_double = TypeGetDouble;

	R->GRng.type = &R->GType;
	R->GRng.state = &R->MT;

	TypeSet(&R->MT, Seed);
}

double FastZigX[FAST_ZIG_C + 1], FastZigRatio[FAST_ZIG_C];

#ifndef QTSIM_GSL_RNG
static uint64_t SplitMix64(uint64_t *x)
{
	uint64_t z;

	z = (*x += 0x9e3779b97f4a7c15ULL);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
	return z ^ (z >> 31);
}

static unsigned long TypeGetX(void *State)
{
	return (unsigned long)(XsNext((FAST_XS*)State) >> 32);
}

static double TypeGetDoubleX(void *State)
{
	return XsUniform((FAST_XS*)State);
}
#endif

// ZIGNOR's tail beyond FAST_ZIG_R (Marsaglia 1964)
double XsZigTail(FAST_XS *X, int Negative)
{
	double x, y;

	do
	{
		x = log(XsOpen(X)) / FAST_ZIG_R;
		y = log(XsOpen(X));
	} while(-2.0 * y < x * x);

	return Negative ? x - FAST_ZIG_R : FAST_ZIG_R - x;
}

// QTSIM_FAST=2: every draw from here on comes from xoshiro256++. Not available
// in the build with GSL's generator (returns 0).
int FastRngSetX(FAST_RNG *R, unsigned long Seed)
{
#ifdef QTSIM_GSL_RNG
	(void)R; (void)Seed;
	return 0;
#else
	uint64_t x;
	double f;
	int i;

	x = (uint64_t)Seed;
	for(i=0;i<4;i++)
		R->Xs.s[i] = SplitMix64(&x);

	f = exp(-0.5 * FAST_ZIG_R * FAST_ZIG_R);
	FastZigX[0] = FAST_ZIG_V / f;
	FastZigX[1] = FAST_ZIG_R;
	FastZigX[FAST_ZIG_C] = 0;
	for(i=2;i<FAST_ZIG_C;i++)
	{
		FastZigX[i] = sqrt(-2.0 * log(FAST_ZIG_V / FastZigX[i - 1] + f));
		f = exp(-0.5 * FastZigX[i] * FastZigX[i]);
	}
	for(i=0;i<FAST_ZIG_C;i++)
		FastZigRatio[i] = FastZigX[i + 1] / FastZigX[i];

	R->X = 1;

	// GSL's own functions read the same stream
	R->GType.name = "qtsimfast_xoshiro256pp";
	R->GType.max = 0xffffffffUL;
	R->GType.min = 0;
	R->GType.get = TypeGetX;
	R->GType.get_double = TypeGetDoubleX;
	R->GRng.type = &R->GType;
	R->GRng.state = &R->Xs;

	return 1;
#endif
}

// Compare the stream, from its current point, with gsl_rng_mt19937 seeded with
// Seed. R itself is not advanced. Returns TRUE (1) if the streams agree.
int FastRngSelfCheck(const FAST_RNG *R, unsigned long Seed)
{
#ifdef QTSIM_GSL_RNG
	(void)R; (void)Seed;
	return 1;
#else
	FAST_MT *Copy;
	gsl_rng *G;
	int Index, Ret;

	Copy = (FAST_MT*)malloc(sizeof(FAST_MT));
	if(Copy == NULL)
		return 0;

	memcpy(Copy, &R->MT, sizeof(FAST_MT));

	G = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(G, Seed);

	Ret = 1;
	if(gsl_rng_max(G) != 0xffffffffUL || gsl_rng_min(G) != 0)
		Ret = 0;

	for(Index=0;Index<3*FMT_N+7 && Ret == 1;Index++)
		if(gsl_rng_get(G) != TypeGet(Copy))
			Ret = 0;

	gsl_rng_free(G);
	free(Copy);

	return Ret;
#endif
}
