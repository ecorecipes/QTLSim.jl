/*
*  model_check.c -- an executable cross-check of the Lean model against the C code.
*
*  This file is NOT part of the original QTSim. It was written in 2026-11 for the
*  QTSimFast package (QTSimFast is a modified version of QTSim, QTLSim 1.0,
*  copyright 2023 Andrew Meade, University of Reading).
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version. It is distributed WITHOUT ANY WARRANTY;
*  see the GNU General Public License for more details.
*
*  What it does. It includes one of the two programs' QTSim.c whole (with main
*  renamed), so that it calls the program's own routines, not copies of them:
*
*    -DMODEL_ORIG  the original:  CaclGenotype and RecombinationPoissonChr
*    -DMODEL_FAST  QTSimFast:     SumGenotypes; RecombinationPoissonChr,
*                                 FillGametes and MixBlocks, and the pipeline's
*                                 GameteDraws and FillBatch, which must give the
*                                 same gamete (built with
*                                 -DQTSIM_GSL_RNG so that its draws go through
*                                 a gsl_rng, which is replaced by a script)
*
*  on one tiny instance: 4 diploid individuals, 6 loci (the fast program with
*  blocks of 1, 2, 3 and 6 loci), allele of chromosome g at locus l
*
*    exact:    g*10 + l                      (sums are exact: order cannot matter)
*    inexact:  scale[l] + 1/(g+3), scale = 1, 2^53, 3, -2^53, 0.5, 7
*                                            (sums are rounded: order matters;
*                                            ModelCheck.lean checks that it does)
*
*  It prints every individual's genotype and, for five scripted sets of random
*  decisions, the gamete of individual 1, as the 64-bit patterns of the doubles.
*  The random generator is a script of raw 32-bit outputs chosen so that GSL's
*  gsl_ran_poisson, gsl_rng_uniform_int and gsl_rng_uniform return the decisions
*  listed in CASES below. formal/ModelCheck.lean evaluates the Lean model on the
*  same instance and the same decisions; `make model-check` compares the three
*  outputs byte for byte.
*/

#include <stdint.h>

#define main qtsim_main
#include "QTSim.c"
#undef main

typedef struct { const unsigned long *v; int n, pos; } SCRIPT;

static unsigned long ScriptGet(void *St)
{
	SCRIPT *S = (SCRIPT*)St;
	if(S->pos >= S->n) { fprintf(stderr, "model_check: script exhausted\n"); exit(4); }
	return S->v[S->pos++];
}
static double ScriptGetDouble(void *St) { return ScriptGet(St) / 4294967296.0; }
static void ScriptSet(void *St, unsigned long Seed) { (void)St; (void)Seed; }

static const gsl_rng_type ScriptType = {"script", 0xffffffffUL, 0, sizeof(SCRIPT), ScriptSet, ScriptGet, ScriptGetDouble};

#define NIND 4
#define NLOCI 6

// Raw outputs. HI: a uniform of about 0.988 (keeps the Poisson loop going);
// 0: a uniform of 0 (ends the Poisson loop; as the strand draw, 0 < 0.5 chooses
// the parent's second chromosome); HALF: 0.5, not below 0.5, the first chromosome;
// PT(p): gsl_rng_uniform_int(rng, 5) returns p.
#define HI 0xFD000000UL
#define HALF 0x80000000UL
#define PT(p) ((unsigned long)(p) * (0xffffffffUL / 5))

// The decisions, as the Lean model is given them (formal/ModelCheck.lean):
//   case 1: 2 points: 4, then 4 again (rejected), then 1; second chromosome first
//   case 2: 0 points; first chromosome
//   case 3: 1 point: 3; first chromosome
//   case 4: 3 points: 2, 4, 3; second chromosome first
//   case 5: 6 points drawn (as many as loci: no crossover); second chromosome first
static const unsigned long C1[] = {HI, HI, 0, PT(4), PT(4), PT(1), 0};
static const unsigned long C2[] = {0, HALF};
static const unsigned long C3[] = {HI, 0, PT(3), HALF};
static const unsigned long C4[] = {HI, HI, HI, 0, PT(2), PT(4), PT(3), 0};
static const unsigned long C5[] = {HI, HI, HI, HI, HI, HI, 0, 0};
static const struct { const unsigned long *v; int n; } CASES[] = {
	{C1, 7}, {C2, 2}, {C3, 4}, {C4, 8}, {C5, 8}
};
#define NCASES 5

// magnitudes that differ by 2^53 from locus to locus
static const double Scale[NLOCI] = {1.0, 9007199254740992.0, 3.0, -9007199254740992.0, 0.5, 7.0};

static double Val(int Variant, int G, int L)
{
	if(Variant == 0)
		return (double)(G * 10 + L);

	return Scale[L] + 1.0 / (double)(G + 3);
}

static unsigned long long Bits(double X)
{
	uint64_t B;
	memcpy(&B, &X, sizeof(B));
	return (unsigned long long)B;
}

// One variant's lines, into Out.
static char *Emit(char *Out, int Variant, const double *Geno, double Gam[NCASES][NLOCI])
{
	int I, C, L;

	Out += sprintf(Out, "variant %s\n", Variant == 0 ? "exact" : "inexact");
	for(I=0;I<NIND;I++)
		Out += sprintf(Out, "geno %d %llu\n", I, Bits(Geno[I]));
	for(C=0;C<NCASES;C++)
	{
		Out += sprintf(Out, "gamete %d", C + 1);
		for(L=0;L<NLOCI;L++)
			Out += sprintf(Out, " %llu", Bits(Gam[C][L]));
		Out += sprintf(Out, "\n");
	}
	return Out;
}

static OPTIONS *MakeOptions(void)
{
	OPTIONS *Opt;

	Opt = (OPTIONS*)SMalloc(sizeof(OPTIONS));
	SetDefOptions(Opt);
	Opt->PopSize = NIND;
	Opt->NoQTL = NLOCI;
	Opt->Ploidy = 2;
	Opt->RecombPoissonMean = 3.0;
	Opt->EnvVar = 0;
	Opt->Seed = 1;
	Opt->X = (double*)SMalloc(sizeof(double) * NIND);
	Opt->Y = (double*)SMalloc(sizeof(double) * NIND);
	return Opt;
}

#ifdef MODEL_ORIG

static void Run(int Variant, char *Out)
{
	OPTIONS *Opt;
	INDIVIDUAL *Ind[NIND];
	gsl_rng R;
	SCRIPT Sc;
	double Geno[NIND], Gam[NCASES][NLOCI];
	int I, C, L, K;

	Opt = MakeOptions();
	R.type = &ScriptType;
	R.state = &Sc;
	Opt->rng = &R;

	for(I=0;I<NIND;I++)
	{
		Ind[I] = CreateIndividual(Opt);
		for(K=0;K<2;K++)
			for(L=0;L<NLOCI;L++)
				Ind[I]->Genome[K][L] = Val(Variant, 2 * I + K, L);
	}

	for(I=0;I<NIND;I++)
		Geno[I] = CaclGenotype(Opt, Ind[I]);

	for(C=0;C<NCASES;C++)
	{
		Sc.v = CASES[C].v; Sc.n = CASES[C].n; Sc.pos = 0;
		RecombinationPoissonChr(Opt, Gam[C], Ind[1]->Genome[0], Ind[1]->Genome[1]);
		if(Sc.pos != Sc.n) { fprintf(stderr, "model_check: case %d left draws unused\n", C + 1); exit(4); }
	}

	Emit(Out, Variant, Geno, Gam);
}

int main(void)
{
	static char Out[2][8192];
	int V;

	for(V=0;V<2;V++)
	{
		Run(V, Out[V]);
		fputs(Out[V], stdout);
	}
	return 0;
}

#endif

#ifdef MODEL_FAST

#ifndef QTSIM_GSL_RNG
#error "build the fast variant with -DQTSIM_GSL_RNG"
#endif

static void Run(int Variant, int BlockSize, char *Out)
{
	OPTIONS *Opt;
	SIM *S;
	SCRIPT Sc;
	double Gam[NCASES][NLOCI];
	char Str[32];
	uint32_t Id;
	int G, K, O, C, L;

	sprintf(Str, "%d", BlockSize);
	setenv("QTSIM_BLOCK", Str, 1);
	setenv("QTSIM_THREADS", "1", 1);

	Opt = MakeOptions();
	S = CreateSim(Opt);
	if(S->B != BlockSize) { fprintf(stderr, "model_check: block size not taken\n"); exit(4); }

	S->Rng.GRng.type = &ScriptType;
	S->Rng.GRng.state = &Sc;

	// every chromosome in blocks of its own
	for(G=0;G<(int)S->Row;G++)
		for(K=0;K<S->NB;K++)
		{
			Id = NewBlock(S, 0);
			for(O=0;O<BlockLen(S, K);O++)
				Block(S, Id)[O] = Val(Variant, G, K * S->B + O);
			S->Ids[(size_t)K * S->Row + G] = Id;
		}
	S->NoNew = 0;

	SumGenotypes(S, 0, S->N, NULL);

	for(C=0;C<NCASES;C++)
	{
		// gamete 0 from individual 1 (chromosomes 2 and 3); the other gametes copy chromosome 0
		for(G=0;G<(int)S->Row;G++)
		{
			S->Src1[G] = 0; S->Src2[G] = 0; S->Cut[G] = S->NB;
		}
		Sc.v = CASES[C].v; Sc.n = CASES[C].n; Sc.pos = 0;
		RecombinationPoissonChr(S, 0, 2);
		if(Sc.pos != Sc.n) { fprintf(stderr, "model_check: case %d left draws unused\n", C + 1); exit(4); }
		FillGametes(S, 0, S->NB, S->Cur[0]);
		MixBlocks(S);
		for(L=0;L<NLOCI;L++)
			Gam[C][L] = Block(S, S->NIds[(size_t)(L / S->B) * S->Row + 0])[L % S->B];

		// The same gamete as the pipeline builds it (GameteDraws, FillBatch):
		// the same draws from the same script, the same alleles.
		for(G=0;G<(int)S->Row;G++)
		{
			S->Src1[G] = 0; S->Src2[G] = 0; S->Cut[G] = S->NB; S->Off[G] = -1;
		}
		Sc.pos = 0;
		S->Off[0] = GameteDraws(S, 0, 2);
		if(Sc.pos != Sc.n) { fprintf(stderr, "model_check: case %d left draws unused (pipeline)\n", C + 1); exit(4); }
		FillBatch(S, 0, S->N);
		for(L=0;L<NLOCI;L++)
			if(memcmp(&Gam[C][L], &Block(S, S->NIds[(size_t)(L / S->B) * S->Row + 0])[L % S->B], sizeof(double)) != 0)
			{
				fprintf(stderr, "model_check: case %d, locus %d: the pipeline's gamete differs\n", C + 1, L);
				exit(4);
			}
	}

	Emit(Out, Variant, S->Geno, Gam);
}

int main(void)
{
	static const int Sizes[] = {2, 1, 3, 6};
	static char Out[4][8192];
	int V, B;

	for(V=0;V<2;V++)
	{
		for(B=0;B<4;B++)
		{
			Run(V, Sizes[B], Out[B]);
			if(strcmp(Out[B], Out[0]) != 0)
			{
				fprintf(stderr, "model_check: block size %d gives other values than block size 2\n", Sizes[B]);
				return 1;
			}
		}
		fputs(Out[0], stdout);
	}
	return 0;
}

#endif
