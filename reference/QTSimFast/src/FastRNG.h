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
*  (2026-11-20) With QTSIM_FAST=2, and only then, the draws come from
*  xoshiro256++ by other methods instead (FastRngSetX, below); the output is then
*  no longer QTSim's.
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
	int		X;				// QTSIM_FAST=2: draws from Xs (xoshiro256++) by other methods
	uint64_t Xs[4];
} FAST_RNG;

void	FastRngSet(FAST_RNG *R, unsigned long Seed);
void	FastMTRefill(FAST_MT *S);
int		FastRngSelfCheck(const FAST_RNG *R, unsigned long Seed);
int		FastRngSetX(FAST_RNG *R, unsigned long Seed);	// QTSIM_FAST=2's generator; 0 if unavailable
double	FastZigTail(FAST_RNG *R, int Negative);

/*
*  QTSIM_FAST=2 (not bit-identical to QTSim; see FAST.md): xoshiro256++ (Blackman
*  and Vigna 2018), seeded from the options file's Seed through SplitMix64;
*  uniforms with 53 random bits; integers by Lemire's (2019) multiply-and-reject
*  method; Poisson variates with mean <= 10 by inversion (one uniform each);
*  normal variates by Doornik's (2005) ZIGNOR ziggurat, 128 layers, the layer
*  index and the uniform taken from different bits of one 64-bit draw. GSL's own
*  distribution functions, where still called, read the same xoshiro stream.
*/
static inline uint64_t FastRotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static inline uint64_t FastNextX(FAST_RNG *R)
{
	uint64_t *s = R->Xs;
	const uint64_t result = FastRotl(s[0] + s[3], 23) + s[0];
	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];
	s[2] ^= t;
	s[3] = FastRotl(s[3], 45);

	return result;
}

// (0, 1) and [0, 1), 53 random bits
static inline double FastOpenX(FAST_RNG *R) { return ((double)(FastNextX(R) >> 11) + 0.5) * 0x1.0p-53; }
static inline double FastUniformX(FAST_RNG *R) { return (double)(FastNextX(R) >> 11) * 0x1.0p-53; }

// Lemire's method: 0 .. n-1, 0 < n < 2^32, unbiased
static inline unsigned long FastUniformIntX(FAST_RNG *R, uint32_t n)
{
	uint64_t m;
	uint32_t l, t;

	m = (FastNextX(R) >> 32) * (uint64_t)n;
	l = (uint32_t)m;
	if(l < n)
	{
		t = (uint32_t)(-n) % n;
		while(l < t)
		{
			m = (FastNextX(R) >> 32) * (uint64_t)n;
			l = (uint32_t)m;
		}
	}

	return (unsigned long)(m >> 32);
}

// ZIGNOR; the tables are set by FastRngSetX
#define FAST_ZIG_C 128
#define FAST_ZIG_R 3.442619855899
#define FAST_ZIG_V 9.91256303526217e-3
extern double FastZigX[FAST_ZIG_C + 1], FastZigRatio[FAST_ZIG_C];

static inline double FastNormalX(FAST_RNG *R)
{
	uint64_t b;
	double u, x, f0, f1;
	int i;

	for(;;)
	{
		b = FastNextX(R);
		i = (int)(b & 0x7f);
		u = 2.0 * ((double)(b >> 11) * 0x1.0p-53) - 1.0;

		if(fabs(u) < FastZigRatio[i])
			return u * FastZigX[i];

		if(i == 0)
			return FastZigTail(R, u < 0);

		x = u * FastZigX[i];
		f0 = exp(-0.5 * (FastZigX[i] * FastZigX[i] - x * x));
		f1 = exp(-0.5 * (FastZigX[i + 1] * FastZigX[i + 1] - x * x));
		if(f1 + FastUniformX(R) * (f0 - f1) < 1.0)
			return x;
	}
}

// Poisson by inversion, mu <= 10, emu = exp(-mu)
static inline unsigned int FastPoissonX(FAST_RNG *R, double mu, double emu)
{
	double u, p, F;
	unsigned int k;

	u = FastUniformX(R);
	p = emu;
	F = p;
	k = 0;
	while(u >= F && k < 1000)
	{
		k++;
		p *= mu / k;
		F += p;
	}

	return k;
}

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
	if(R->X)
		return FastUniformX(R);

	return FastGet(R) / 4294967296.0;
}

// gsl_rng_uniform_int, with min = 0 and max = 0xffffffff
static inline unsigned long FastUniformInt(FAST_RNG *R, unsigned long n)
{
	unsigned long scale, k;

	if(R->X && n > 0 && n <= 0xffffffffUL)
		return FastUniformIntX(R, (uint32_t)n);

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

	if(R->X)
		return FastPoissonX(R, mu, emu);

	prod = 1.0;
	k = 0;

	do
	{
		prod *= FastUniform(R);
		k++;
	} while(prod > emu);

	return k - 1;
}

// gsl_ran_gaussian(r, sigma)
static inline double FastGaussian(FAST_RNG *R, double sigma)
{
	if(R->X)
		return sigma * FastNormalX(R);

	return gsl_ran_gaussian(&R->GRng, sigma);
}

#else

static inline double FastGaussian(FAST_RNG *R, double sigma) { return gsl_ran_gaussian(&R->GRng, sigma); }
static inline double FastUniform(FAST_RNG *R) { return gsl_rng_uniform(&R->GRng); }
static inline unsigned long FastUniformInt(FAST_RNG *R, unsigned long n) { return gsl_rng_uniform_int(&R->GRng, n); }
static inline unsigned int FastPoisson(FAST_RNG *R, double mu, double emu) { (void)emu; return gsl_ran_poisson(&R->GRng, mu); }

#endif

#endif
