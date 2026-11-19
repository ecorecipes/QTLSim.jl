/*
*  QTSimFast: an inlined MT19937 generator and the GSL draws built on it.
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
*  What it does. QTSim draws every random number from gsl_rng_mt19937 through
*  GSL's function pointers. Here the same generator (Matsumoto and Nishimura's
*  MT19937, seeded as GSL's mt_set seeds it) is inlined, with
*
*    FastUniform       = gsl_rng_uniform       (get / 4294967296.0)
*    FastUniformInt    = gsl_rng_uniform_int   (integer arithmetic only)
*    FastPoisson       = gsl_ran_poisson       for mu <= 10: the product-of-uniforms
*                        loop, with exp(-mu) taken from the same libm; for mu > 10
*                        (and NaN is handled as GSL does) it calls gsl_ran_poisson
*
*  and a gsl_rng whose type reads the same state, so that every other draw
*  (gsl_ran_gaussian, gsl_ran_poisson for large means) is made by GSL itself from
*  the same stream. Nothing with floating-point arithmetic that a compiler could
*  contract or reorder is reimplemented: the three functions above use one
*  division by a power of two, one multiplication and comparisons.
*
*  FastRngSelfCheck compares the stream with gsl_rng_mt19937 for the run's seed
*  (three full state regenerations); the program refuses to run if they differ.
*
*  Building with -DQTSIM_GSL_RNG replaces all of this with plain calls to GSL
*  (the original's code path), for comparison.
*/

#ifndef FAST_RNG_H
#define FAST_RNG_H

#include <stdint.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define FMT_N 624
#define FMT_M 397

typedef struct
{
	uint32_t mt[FMT_N];		// the state, as in GSL's mt_state_t
	uint32_t out[FMT_N];	// the tempered outputs of the current state
	int mti;
} FAST_MT;

typedef struct
{
	FAST_MT	MT;
	gsl_rng	GRng;			// a gsl_rng reading MT, for GSL's own distribution functions
	gsl_rng_type GType;
} FAST_RNG;

void	FastRngSet(FAST_RNG *R, unsigned long Seed);
void	FastMTRefill(FAST_MT *S);
int		FastRngSelfCheck(const FAST_RNG *R, unsigned long Seed);

#ifndef QTSIM_GSL_RNG

static inline uint32_t FastGet(FAST_RNG *R)
{
	if(R->MT.mti >= FMT_N)
		FastMTRefill(&R->MT);

	return R->MT.out[R->MT.mti++];
}

// gsl_rng_uniform for mt19937: mt_get_double
static inline double FastUniform(FAST_RNG *R)
{
	return FastGet(R) / 4294967296.0;
}

// gsl_rng_uniform_int, with min = 0 and max = 0xffffffff
static inline unsigned long FastUniformInt(FAST_RNG *R, unsigned long n)
{
	unsigned long scale, k;

	if(n > 0xffffffffUL || n == 0)
		return gsl_rng_uniform_int(&R->GRng, n);	// GSL's error path

	scale = 0xffffffffUL / n;

	do
	{
		k = FastGet(R) / scale;
	} while(k >= n);

	return k;
}

// gsl_ran_poisson(r, mu), where emu = exp(-mu)
static inline unsigned int FastPoisson(FAST_RNG *R, double mu, double emu)
{
	double prod;
	unsigned int k;

	if(mu > 10)
		return gsl_ran_poisson(&R->GRng, mu);

	prod = 1.0;
	k = 0;

	do
	{
		prod *= FastUniform(R);
		k++;
	} while(prod > emu);

	return k - 1;
}

#else

static inline double FastUniform(FAST_RNG *R) { return gsl_rng_uniform(&R->GRng); }
static inline unsigned long FastUniformInt(FAST_RNG *R, unsigned long n) { return gsl_rng_uniform_int(&R->GRng, n); }
static inline unsigned int FastPoisson(FAST_RNG *R, double mu, double emu) { (void)emu; return gsl_ran_poisson(&R->GRng, mu); }

#endif

#endif
