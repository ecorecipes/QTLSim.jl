/*
*  QTSimFast -- a MODIFIED version of QTSim (QTLSim 1.0)
*
*  Modified 2026-11 for speed, for the QTLSim.jl project
*  (https://github.com/ecorecipes/QTLSim.jl), with Claude (Anthropic). The output
*  is bit-identical to the original's for the same options file and seed. The
*  original is by Andrew Meade; its notice follows unchanged.
*
*  What was changed in this file (see ../README.md for the full account):
*   - the simulation core (populations, gametes, mutation, genotypes, output) was
*     rewritten: genomes are stored as shared blocks of loci that are copied only
*     when a crossover or a mutation changes them, instead of one array per
*     individual copied locus by locus; genotypes are summed over loci in the
*     original order, several individuals at a time, optionally on several
*     threads; crossover maps are not allocated or sorted (only the lowest point
*     is ever used, as in the original); statistics that were computed and never
*     used, or computed twice, are computed once or not at all;
*   - the random numbers are those of gsl_rng_mt19937 in the original order, drawn
*     through FastRNG.h;
*   - the options code (defaults, option file, printed header) is unchanged;
*   - functions the original never calls were removed.
*  GenLib.*, PassOptionFile.*, RandLib.* are unmodified copies of the original's.
*
*  ---------------------------------------------------------------------------
*
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "QTSim.h"
#include "RandLib.h"
#include "PassOptionFile.h"
#include "FastRNG.h"
#include "FastThreads.h"

// Draw a random value from a normal distribution, centred at 0 with the specified variance.
double RandNormalVar(gsl_rng* rng, double Var)
{
	// both are the same
//	return gsl_ran_gaussian(rng, 1.0) * sqrt(Var);
	return gsl_ran_gaussian(rng, sqrt(Var));
}

// Allocate and set the mean and SD of the fitness function.
double* SetFitnessNormalParam(double Mean, double SD)
{
	double *Param;

	Param = (double*)SMalloc(sizeof(double) * 2);

	Param[0] = Mean;
	Param[1] = SD;

	return Param;
}

void CaclOptParm(OPTIONS *Opt)
{
//	Opt->OptVarTarget = 0.5;
//	Opt->MuteVar = Opt->EnvVar * Opt->MuteVarScalar;
//	return;

	Opt->OptVarTarget = Opt->Heritabitlity * Opt->PheotypeVar;
	Opt->EnvVar = Opt->PheotypeVar - Opt->OptVarTarget;
	Opt->MuteVar = Opt->EnvVar * Opt->MuteVarScalar;

//	Opt->MuteVar = 1.0 / 1000.0;
}

// Allocate additional options parameters
 void AllocAddtionalOptions(OPTIONS* Opt)
 {
	 CaclOptParm(Opt);

	Opt->rng = gsl_rng_alloc(gsl_rng_mt19937);

	gsl_rng_set(Opt->rng, Opt->Seed);

	Opt->X = (double*)SMalloc(sizeof(double) * Opt->PopSize);
	Opt->Y = (double*)SMalloc(sizeof(double) * Opt->PopSize);
 }


 // Set the default options
void SetDefOptions(OPTIONS* Opt)
{
	Opt->NoGenerations = 20000;
	Opt->Sample = 1;
	Opt->PopSize = 1000;
	Opt->NoQTL = 1;
	Opt->Ploidy = 1;
//	Opt->EnvVar = 3.0;
	Opt->MuteVarScalar = 1.0 / 1000.00;
//	Opt->MuteVar = Opt->EnvVar * Opt->MuteVarScalar;

	Opt->FitnessPar = SetFitnessNormalParam(0.0, 0.1);
	Opt->Seed = GetSeed();
	Opt->Drift = FALSE;

//	Opt->OptVar = TRUE;
	Opt->OptVarTarget = -1.0;

	Opt->ChangeFitnessGen = -1;
	Opt->ChangeFitnessSD = 0;

	Opt->MoveFitnessSDGen = -1;
	Opt->MoveFitnessSD = 1.0;

	Opt->PheotypeVar = 1.0;
	Opt->Heritabitlity = 0.5;

	Opt->StandardisePhenotype = FALSE;

	Opt->MoveType = PHENOTYPE;

	Opt->AssortativeMating = -1;

	Opt->InitFitnessSD = -1;

	Opt->MutationRatePerGamete = MUT_RATE_PER_QTL;

	Opt->RecombPoissonMean = RECOMB_POISSON_DEFAULT;
}


// Create options and read them form the file.
OPTIONS* CreateOptionsFile(int argc, char** argv)
{
	OPTIONS *Opt;

	Opt = (OPTIONS*)SMalloc(sizeof(OPTIONS));

	SetDefOptions(Opt);

	PassOptFile(Opt, argv[1]);

	AllocAddtionalOptions(Opt);

	return Opt;
}


void FreeOptions(OPTIONS* Opt)
{
	gsl_rng_free(Opt->rng);
	free(Opt->FitnessPar);
	free(Opt->X);
	free(Opt->Y);
	free(Opt);
}

// Print the current options.
void PrintOptions(OPTIONS* Opt)
{
	printf("Version:\t%f\n", VERSION);

	printf("Fitness:\t%f\t%f\n", Opt->FitnessPar[0], Opt->FitnessPar[1]);


	printf("Sample:\t%d\n", Opt->Sample);
	printf("PopSize:\t%d\n", Opt->PopSize);
	printf("NoQTL:\t%d\n", Opt->NoQTL);
	printf("Ploidy:\t%d\n", Opt->Ploidy);

	printf("OptimizeVariance:\t%f\n", Opt->OptVarTarget);

	if(Opt->Drift == FALSE)
		printf("Drift:\tFalse\n");
	else
		printf("Drift:\tTrue\n");

	printf("EnvVar:\t%f\n", Opt->EnvVar);
	printf("MutationScalar:\t%f\n", Opt->MuteVarScalar);
	printf("MutationVar:\t%f\n", Opt->MuteVar);
	printf("Generations:\t%d\n", Opt->NoGenerations);


	printf("ChangeFitnessSD:\t%d\t%f\n", Opt->ChangeFitnessGen, Opt->ChangeFitnessSD);

	printf("MoveFitnessSD:\t%d\t%f\n", Opt->MoveFitnessSDGen, Opt->MoveFitnessSD);


	if(Opt->MoveFitnessSDGen != -1)
	{
		printf("MoveType:\t");

		if(Opt->MoveType == FITNESS)
			printf("Fitness\n");

		if(Opt->MoveType == PHENOTYPE)
			printf("Phenotype\n");

		if(Opt->MoveType == CONTINUOUS)
			printf("Continuous\n");
	}

	printf("Seed:\t%ld\n", Opt->Seed);

	printf("PheotypeVar:\t%f\n", Opt->PheotypeVar);
	printf("Heritabitlity:\t%f\n", Opt->Heritabitlity);

	if(Opt->StandardisePhenotype == FALSE)
		printf("StandardisePhenotype:\tFalse\n");
	else
		printf("StandardisePhenotype:\tTrue\n");

/*
	if(Opt->AssortativeMating == -1)
		printf("AssortativeMating:\tFalse\n");
	else
		printf("AssortativeMating:\t%f\n", Opt->AssortativeMating);
*/

	printf("InitFitnessSD:\t%f\n", Opt->InitFitnessSD);

	printf("MutationRatePerGamete\t%f\n", Opt->MutationRatePerGamete);

	printf("RecombinationPoissonMean\t%f\n", Opt->RecombPoissonMean);

	fflush(stdout);
}

/* ---------------------------------------------------------------------------
*  The simulation core of QTSimFast. Everything from here to GetAnalyticalFitSD
*  replaces the original's populations of INDIVIDUALs.
*
*  Bit parity rests on three rules, kept throughout:
*   1. every random draw is made in the original order, with the same
*      distribution function and arguments;
*   2. every floating-point value that is printed or compared is computed by the
*      same operations in the same order as the original computes it (a genotype
*      is ((0 + (a0+b0)) + (a1+b1)) + ..., locus by locus; phenotype, fitness and
*      the statistics are the original's expressions and GSL calls);
*   3. threads and the block size only decide who computes a value and where an
*      allele is stored, never how a value is computed.
* -------------------------------------------------------------------------*/

// The number of individuals whose genotypes are summed side by side. Each sum is
// a chain of additions that must be done in order; chains of different
// individuals are independent, so the processor can work on several at once.
#define LANES 8

// Pairs of doubles added elementwise in one instruction where the compiler has
// vector extensions; elementwise IEEE addition gives the same bits as scalar.
#if (defined(__GNUC__) || defined(__clang__)) && !defined(QTSIM_NO_VECTOR)
	#define FAST_VECTOR
	typedef double V2D __attribute__((vector_size(16)));

	static inline V2D Load2(const double *P)
	{
		V2D V;
		memcpy(&V, P, sizeof(V2D));
		return V;
	}
#endif

typedef struct
{
	int		N, L, P;			// PopSize, NoQTL, Ploidy

	// Genomes. A chromosome is NB block numbers; block b holds the alleles of
	// loci b*B .. b*B+B-1 (the last block LastLen of them). Chromosomes share
	// blocks until a crossover inside a block or a mutation makes a new one.
	int		B, NB, LastLen;
	size_t	IndIds;				// P * NB, block numbers per individual
	double	*Data;				// Cap blocks of B alleles
	uint32_t Cap, HW;			// capacity; blocks HW.. have never been used
	uint32_t *Free, NoFree;		// free blocks below HW
	uint32_t *Owner;			// for a block made this generation, the chromosome
								// that alone holds it (1-based), else 0
	uint32_t *New, NoNew;		// the blocks made this generation
	uint32_t *Ids, *NIds;		// block numbers of the current and the next population

	// Collecting unused blocks: marks are set while the genotypes are summed.
	uint8_t	**Mark;				// one array per thread
	int		GCEvery, SinceGC, Marking;

	// Per individual, current population (QTSim's INDIVIDUAL)
	double	*Geno, *Pheno, *Fit, *Env;
	double	*NEnv;				// environmental deviations of the next population
	int		*Surv;				// QTSim's SPop, as indices
	int		NoSurv;				// QTSim's NoSPop

	FAST_RNG Rng;
	double	EnvSD;				// sqrt(EnvVar)
	double	RecMu, RecEmu;		// crossover number: mean and exp(-mean)
	double	MutMu, MutEmu;		// mutations per chromosome: mean and exp(-mean)
	double	MutSD;				// sqrt of the variance of a mutation's effect
	int		*Points;			// crossover points of one gamete
	uint32_t *Seen, Stamp;		// to find repeated points when there are many

	// The phenotype mean of the current generation, if already computed
	int		HaveMean;
	double	Mean;

	int		UseM;				// gsl_stats_*_m give the bits of gsl_stats_*
	int		Lanes;

	FAST_POOL *Pool;			// for genotypes
	int		ParGeno;			// sum genotypes on the pool's threads
	FAST_POOL *OutPool;			// for output
	int		OutBusy;
} SIM;

// What a line of output is computed from: a copy when it is written by another
// thread while the simulation goes on.
typedef struct
{
	int		Itter, N, NoS;
	double	Theta, Omega, MuteVar;
	int		HaveMean, UseM;
	double	Mean;
	double	*Pheno, *Geno;		// N
	double	*X, *Y;				// NoS survivors: phenotype, fitness
} OUT_JOB;

static OUT_JOB OutJob;

static int EnvInt(const char *Name, int Def, int Min, int Max)
{
	char *Str, *End;
	long Val;

	Str = getenv(Name);
	if(Str == NULL || Str[0] == '\0')
		return Def;

	Val = strtol(Str, &End, 10);
	if(*End != '\0' || Val < Min || Val > Max)
	{
		fprintf(stderr, "QTSimFast: ignoring %s=%s (expecting an integer from %d to %d)\n", Name, Str, Min, Max);
		return Def;
	}

	return (int)Val;
}

static inline double* Block(const SIM *S, uint32_t Id)
{
	return S->Data + (size_t)Id * S->B;
}

static inline int BlockLen(const SIM *S, int BNo)
{
	if(BNo == S->NB - 1)
		return S->LastLen;

	return S->B;
}

// A block for chromosome Tag alone.
static inline uint32_t NewBlock(SIM *S, uint32_t Tag)
{
	uint32_t Id;

	if(S->NoFree > 0)
		Id = S->Free[--S->NoFree];
	else
	{
		if(S->HW >= S->Cap)
		{
			fprintf(stderr, "QTSimFast: internal error, out of blocks\n");
			exit(3);
		}
		Id = S->HW++;
	}

	S->Owner[Id] = Tag;
	S->New[S->NoNew++] = Id;

	return Id;
}

// Does gsl_stats_variance_m(x, mean(x)) give the bits of gsl_stats_variance(x),
// and likewise gsl_stats_sd_m? In GSL's source they are the same code; this
// checks the library in use, and if it ever differed the plain functions would
// be called instead.
static int StatsMSame(void)
{
	static const int Sizes[] = {1, 2, 3, 10, 257, 1000};
	double *V, A, B, Mean;
	gsl_rng *G;
	int SIndex, Index, Size, Ret;

	V = (double*)SMalloc(sizeof(double) * 1000);
	G = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(G, 20261119);

	Ret = TRUE;
	for(SIndex=0;SIndex<(int)(sizeof(Sizes)/sizeof(Sizes[0]));SIndex++)
	{
		Size = Sizes[SIndex];
		for(Index=0;Index<Size;Index++)
			V[Index] = gsl_ran_gaussian(G, 0.37) + 0.11;

		Mean = gsl_stats_mean(V, 1, Size);

		A = gsl_stats_variance(V, 1, Size);
		B = gsl_stats_variance_m(V, 1, Size, Mean);
		if(memcmp(&A, &B, sizeof(double)) != 0)
			Ret = FALSE;

		A = gsl_stats_sd(V, 1, Size);
		B = gsl_stats_sd_m(V, 1, Size, Mean);
		if(memcmp(&A, &B, sizeof(double)) != 0)
			Ret = FALSE;
	}

	gsl_rng_free(G);
	free(V);

	return Ret;
}

static int DefBlockSize(OPTIONS *Opt)
{
	if(Opt->NoQTL <= 16)
		return Opt->NoQTL;

	if(Opt->NoQTL <= 256)
		return 16;

	return 32;
}

// Allocate the populations. Draws the environmental deviations of both of the
// original's populations, as CreatePop does twice.
static SIM* CreateSim(OPTIONS *Opt)
{
	SIM *S;
	size_t NoInd, Cap;
	int Index, NoThreads, MaxThreads;
	char *Str;

	S = (SIM*)SMalloc(sizeof(SIM));
	memset(S, 0, sizeof(SIM));

	S->N = Opt->PopSize;
	S->L = Opt->NoQTL;
	S->P = Opt->Ploidy;

	if(S->N < 1 || S->L < 1)
	{
		fprintf(stderr, "QTSimFast: PopSize and NoQTL must be at least 1.\n");
		exit(1);
	}

	NoInd = (size_t)S->N;

	// Block size: any size gives the same output.
	S->B = EnvInt("QTSIM_BLOCK", 0, 0, 1 << 30);
	if(S->B == 0)
		S->B = DefBlockSize(Opt);
	if(S->B > S->L)
		S->B = S->L;

	for(;;)
	{
		S->NB = (S->L + S->B - 1) / S->B;
		Cap = 2 * NoInd * S->P * S->NB + 2;
		if(Cap < 0xffffffffUL || S->B == S->L)
			break;
		S->B = (S->B * 2 > S->L) ? S->L : S->B * 2;
	}

	S->LastLen = S->L - (S->NB - 1) * S->B;
	S->IndIds = (size_t)S->P * S->NB;
	S->Cap = (uint32_t)Cap;

	// Blocks are taken in order, so memory that is never used is never touched.
	S->Data = (double*)calloc(Cap * S->B, sizeof(double));
	if(S->Data == NULL)
		MallocErr();

	S->Free = (uint32_t*)SMalloc(sizeof(uint32_t) * Cap);
	S->New = (uint32_t*)SMalloc(sizeof(uint32_t) * Cap);
	S->Owner = (uint32_t*)calloc(Cap, sizeof(uint32_t));
	if(S->Owner == NULL)
		MallocErr();

	// Every chromosome starts as block 0, all zeros.
	S->HW = 1;
	S->Ids = (uint32_t*)calloc(NoInd * S->IndIds, sizeof(uint32_t));
	S->NIds = (uint32_t*)calloc(NoInd * S->IndIds, sizeof(uint32_t));
	if(S->Ids == NULL || S->NIds == NULL)
		MallocErr();

	S->Geno = (double*)SMalloc(sizeof(double) * NoInd);
	S->Pheno = (double*)SMalloc(sizeof(double) * NoInd);
	S->Fit = (double*)SMalloc(sizeof(double) * NoInd);
	S->Env = (double*)SMalloc(sizeof(double) * NoInd);
	S->NEnv = (double*)SMalloc(sizeof(double) * NoInd);
	S->Surv = (int*)SMalloc(sizeof(int) * NoInd);

	S->Points = (int*)SMalloc(sizeof(int) * (S->L + 1));
	S->Seen = (uint32_t*)calloc(S->L + 1, sizeof(uint32_t));
	if(S->Seen == NULL)
		MallocErr();

	// Threads. QTSIM_THREADS: a number, or "auto". The default is 1.
	MaxThreads = 256;
	NoThreads = 1;
	Str = getenv("QTSIM_THREADS");
	if(Str != NULL && strcicmp(Str, "auto") == 0)
	{
		// Only large populations gain from threads.
		NoThreads = 1;
		if((double)S->N * S->P * S->L >= 4.0e6)
		{
			NoThreads = FastThreadsCPUs() / 2;
			if(NoThreads > 4)
				NoThreads = 4;
		}
	}
	else
		NoThreads = EnvInt("QTSIM_THREADS", 1, 1, MaxThreads);

	if(NoThreads > S->N / LANES)
		NoThreads = S->N / LANES;
	if(NoThreads < 1)
		NoThreads = 1;

	S->Pool = FastThreadsCreate(NoThreads);
	NoThreads = FastThreadsNo(S->Pool);

	// Small populations are summed faster by one thread than handed out.
	S->ParGeno = FALSE;
	if(NoThreads > 1 && (double)S->N * S->P * S->L >= (double)EnvInt("QTSIM_PAR_MIN", 200000, 0, 1 << 30))
		S->ParGeno = TRUE;

	// With more than one thread the output lines are computed and written by a
	// thread of their own (QTSIM_ASYNC_OUTPUT=0 turns that off).
	S->OutPool = FastThreadsCreate(EnvInt("QTSIM_ASYNC_OUTPUT", NoThreads > 1, 0, 1) ? 2 : 1);

	if(FastThreadsNo(S->OutPool) > 1)
	{
		OutJob.Pheno = (double*)SMalloc(sizeof(double) * NoInd);
		OutJob.Geno = (double*)SMalloc(sizeof(double) * NoInd);
		OutJob.X = (double*)SMalloc(sizeof(double) * NoInd);
		OutJob.Y = (double*)SMalloc(sizeof(double) * NoInd);
	}
	else
	{
		OutJob.X = Opt->X;
		OutJob.Y = Opt->Y;
	}

	S->Mark = (uint8_t**)SMalloc(sizeof(uint8_t*) * NoThreads);
	for(Index=0;Index<NoThreads;Index++)
	{
		S->Mark[Index] = (uint8_t*)calloc(Cap, 1);
		if(S->Mark[Index] == NULL)
			MallocErr();
	}

	S->GCEvery = EnvInt("QTSIM_GC", 1, 1, 1 << 30);
	S->Lanes = EnvInt("QTSIM_LANES", LANES, 1, LANES) == 1 ? 1 : LANES;
	S->UseM = StatsMSame();

	// The random numbers
	FastRngSet(&S->Rng, Opt->Seed);
	if(FastRngSelfCheck(&S->Rng, Opt->Seed) == FALSE)
	{
		fprintf(stderr, "QTSimFast: the built-in generator does not reproduce gsl_rng_mt19937; rebuild with -DQTSIM_GSL_RNG.\n");
		exit(3);
	}

	S->EnvSD = sqrt(Opt->EnvVar);

	S->RecMu = Opt->RecombPoissonMean;
	S->RecEmu = exp(-S->RecMu);

	S->MutMu = Opt->NoQTL * Opt->MutationRatePerGamete;
	S->MutEmu = exp(-S->MutMu);

	// MutateIndividualChrNo
	{
		double MutVar;

		MutVar = Opt->MuteVarScalar / (Opt->Ploidy * Opt->NoQTL * Opt->MutationRatePerGamete);
		MutVar = MutVar * Opt->EnvVar;
		S->MutSD = sqrt(MutVar);
	}

	// CreatePop(CPop), CreatePop(NPop): SetIndividualEnv for each individual.
	// The second population's values are never used, but they are drawn.
	for(Index=0;Index<S->N;Index++)
		S->Env[Index] = (Opt->EnvVar != 0) ? gsl_ran_gaussian(&S->Rng.GRng, S->EnvSD) : 0.0;

	for(Index=0;Index<S->N;Index++)
		S->NEnv[Index] = (Opt->EnvVar != 0) ? gsl_ran_gaussian(&S->Rng.GRng, S->EnvSD) : 0.0;

	return S;
}

// Sum the genotypes of individuals Lo .. Hi-1, W at a time: CaclGenotype.
// Genotype = 0, then Genotype += Genome[0][l] + Genome[1][l] for each locus l in
// order (diploid), or Genotype += Genome[0][l] (haploid).
static inline __attribute__((always_inline)) void SumGenotypesW(const SIM *S, int Lo, int Hi, uint8_t *Mark, const int W, const int Diploid)
{
	const double *A[LANES], *C[LANES];
	const uint32_t *I[LANES];
	double Sum[LANES];
	int Start, J, K, Index, Len, Ind;
	const int NB = S->NB;

	for(Start=Lo;Start<Hi;Start+=W)
	{
		// A short last group is filled up with its last individual.
		for(J=0;J<W;J++)
		{
			Ind = Start + J;
			if(Ind > Hi - 1)
				Ind = Hi - 1;

			I[J] = S->Ids + (size_t)Ind * S->IndIds;
			Sum[J] = 0;
		}

		for(K=0;K<NB;K++)
		{
			Len = (K == NB - 1) ? S->LastLen : S->B;

			for(J=0;J<W;J++)
			{
				A[J] = Block(S, I[J][K]);
				if(Diploid)
					C[J] = Block(S, I[J][NB + K]);
			}

			if(Mark != NULL)
				for(J=0;J<W;J++)
				{
					Mark[I[J][K]] = 1;
					if(Diploid)
						Mark[I[J][NB + K]] = 1;
				}

			Index = 0;

			if(Diploid)
			{
#ifdef FAST_VECTOR
				for(;Index+1<Len;Index+=2)
					for(J=0;J<W;J++)
					{
						V2D Pair = Load2(A[J] + Index) + Load2(C[J] + Index);
						Sum[J] += Pair[0];
						Sum[J] += Pair[1];
					}
#endif
				for(;Index<Len;Index++)
					for(J=0;J<W;J++)
						Sum[J] += A[J][Index] + C[J][Index];
			}
			else
			{
				for(;Index<Len;Index++)
					for(J=0;J<W;J++)
						Sum[J] += A[J][Index];
			}
		}

		for(J=0;J<W && Start+J<Hi;J++)
			S->Geno[Start + J] = Sum[J];
	}
}

static void SumGenotypes(const SIM *S, int Lo, int Hi, uint8_t *Mark)
{
	if(S->P == 2)
	{
		if(S->Lanes == 1)
			SumGenotypesW(S, Lo, Hi, Mark, 1, 1);
		else
			SumGenotypesW(S, Lo, Hi, Mark, LANES, 1);
	}
	else
	{
		if(S->Lanes == 1)
			SumGenotypesW(S, Lo, Hi, Mark, 1, 0);
		else
			SumGenotypesW(S, Lo, Hi, Mark, LANES, 0);
	}
}

static void SumGenotypesJob(void *Ctx, int ThreadNo, int NoThreads)
{
	SIM *S;
	int Lo, Hi, Groups;

	S = (SIM*)Ctx;

	// Whole groups of LANES individuals to each thread.
	Groups = (S->N + LANES - 1) / LANES;
	Lo = (int)((long)Groups * ThreadNo / NoThreads) * LANES;
	Hi = (int)((long)Groups * (ThreadNo + 1) / NoThreads) * LANES;
	if(Hi > S->N)
		Hi = S->N;

	if(Lo < Hi)
		SumGenotypes(S, Lo, Hi, S->Marking ? S->Mark[ThreadNo] : NULL);
}

// Free the blocks no chromosome of the current population uses. Called when
// only the current population is alive and its marks have just been set.
static void CollectBlocks(SIM *S)
{
	uint8_t *M, *O;
	uint32_t Id, HW, NewHW;
	int T, NoT;

	HW = S->HW;
	M = S->Mark[0];

	NoT = S->ParGeno ? FastThreadsNo(S->Pool) : 1;
	for(T=1;T<NoT;T++)
	{
		O = S->Mark[T];
		for(Id=0;Id<HW;Id++)
			M[Id] |= O[Id];
		memset(O, 0, HW);
	}

	// Lower-numbered blocks are reused first.
	S->NoFree = 0;
	NewHW = 0;
	for(Id=HW;Id-- > 0;)
	{
		if(M[Id])
		{
			if(NewHW == 0)
				NewHW = Id + 1;
		}
		else if(NewHW != 0)
			S->Free[S->NoFree++] = Id;
	}

	memset(M, 0, HW);
	S->HW = NewHW;
	S->SinceGC = 0;
}

// Calculate the phenotype for each individual in the population.
static void CaclPopPhenotype(SIM *S)
{
	int Index;

	// Collect unused blocks every GCEvery generations, and whenever a generation
	// that made a new block for every block of every chromosome would run out.
	S->SinceGC++;
	S->Marking = FALSE;
	if(S->SinceGC >= S->GCEvery || (size_t)(S->Cap - S->HW) + S->NoFree < (size_t)S->N * S->IndIds)
		S->Marking = TRUE;

	if(S->ParGeno == TRUE)
		FastThreadsRun(S->Pool, SumGenotypesJob, S);
	else
		SumGenotypes(S, 0, S->N, S->Marking ? S->Mark[0] : NULL);

	if(S->Marking == TRUE)
		CollectBlocks(S);

	for(Index=0;Index<S->N;Index++)
		S->Pheno[Index] = S->Geno[Index] + S->Env[Index];

	S->HaveMean = FALSE;
}

static double PhenotypeMean(SIM *S)
{
	if(S->HaveMean == FALSE)
	{
		S->Mean = gsl_stats_mean(S->Pheno, 1, S->N);
		S->HaveMean = TRUE;
	}

	return S->Mean;
}

// Move the fitness function mean, if appropriate.
// (The original also computes the phenotypes' SD here and does not use it, and
// their mean, which the continuous move does not use.)
static void SetFitnessMean(OPTIONS *Opt, SIM *S)
{
	if(Opt->MoveType == FITNESS)
		Opt->FitnessPar[0] = PhenotypeMean(S) + (Opt->FitnessPar[1] * Opt->MoveFitnessSD);

	if(Opt->MoveType == PHENOTYPE)
		Opt->FitnessPar[0] = PhenotypeMean(S) + (sqrt(Opt->PheotypeVar) * Opt->MoveFitnessSD);

	if(Opt->MoveType == CONTINUOUS)
		Opt->FitnessPar[0] = Opt->FitnessPar[0] + Opt->MoveFitnessSD;
}

// Calculate the fitness of the population
static void PopFitness(OPTIONS *Opt, int Itter, SIM *S)
{
	int Index;
	double PhenotypeSD, Phenotype, Scale;

	CaclPopPhenotype(S);

	PhenotypeSD = 1.0;
	if(Opt->StandardisePhenotype == TRUE)
	{
		if(S->UseM == TRUE)
			PhenotypeSD = gsl_stats_sd_m(S->Pheno, 1, S->N, PhenotypeMean(S));
		else
			PhenotypeSD = gsl_stats_sd(S->Pheno, 1, S->N);
	}

	if(Itter >= Opt->MoveFitnessSDGen && Opt->MoveFitnessSDGen != -1)
		SetFitnessMean(Opt, S);

	// Fitness
	if(Opt->Drift == TRUE)
	{
		for(Index=0;Index<S->N;Index++)
			S->Fit[Index] = 1.0;
		return;
	}

	Scale = (0.4 / Opt->FitnessPar[1]);

	if(Opt->StandardisePhenotype == FALSE)
	{
		for(Index=0;Index<S->N;Index++)
		{
			Phenotype = S->Pheno[Index] - Opt->FitnessPar[0];
			S->Fit[Index] = gsl_ran_gaussian_pdf(Phenotype, Opt->FitnessPar[1]) / Scale;
		}
	}
	else
	{
		for(Index=0;Index<S->N;Index++)
		{
			Phenotype = S->Pheno[Index] - Opt->FitnessPar[0];
			Phenotype = Phenotype / PhenotypeSD;
			S->Fit[Index] = gsl_ran_gaussian_pdf(Phenotype, Opt->FitnessPar[1]) / Scale;
		}
	}
}

// Find the surviving members of the population
static void SetSurviving(SIM *S)
{
	int Index;

	S->NoSurv = 0;
	for(Index=0;Index<S->N;Index++)
		if(FastUniform(&S->Rng) < S->Fit[Index])
			S->Surv[S->NoSurv++] = Index;
}

// The lowest crossover point of a gamete, or NoQTL for none: CreateCombMap.
// The original sorts the points and then switches strand at the first only.
static inline int LowestCrossover(SIM *S)
{
	int NoPoints, Index, Point, Low, PIndex, Found;

	NoPoints = (int)FastPoisson(&S->Rng, S->RecMu, S->RecEmu);

	if(NoPoints == 0 || NoPoints >= S->L)
		return S->L;

	Low = S->L;

	if(NoPoints <= 16)
	{
		for(Index=0;Index<NoPoints;Index++)
		{
			do
			{
				Point = (int)FastUniformInt(&S->Rng, S->L - 1);

				Found = FALSE;
				for(PIndex=0;PIndex<Index;PIndex++)
					if(S->Points[PIndex] == Point)
						Found = TRUE;
			} while(Found == TRUE);

			S->Points[Index] = Point;
			if(Point < Low)
				Low = Point;
		}

		return Low;
	}

	// Many points: remember the points drawn for this gamete in a table.
	if(++S->Stamp == 0)
	{
		memset(S->Seen, 0, sizeof(uint32_t) * (S->L + 1));
		S->Stamp = 1;
	}

	for(Index=0;Index<NoPoints;Index++)
	{
		do
		{
			Point = (int)FastUniformInt(&S->Rng, S->L - 1);
		} while(S->Seen[Point] == S->Stamp);

		S->Seen[Point] = S->Stamp;
		if(Point < Low)
			Low = Point;
	}

	return Low;
}

// Crate a gamete from a parent, using a number of recombination points, drawn
// from a Poisson: RecombinationPoissonChr. Par is the parent's block numbers
// (both chromosomes), Chr the gamete's, Tag the gamete's chromosome number.
static inline void RecombinationPoissonChr(SIM *S, uint32_t *Chr, const uint32_t *Par, uint32_t Tag)
{
	const uint32_t *CP, *OP;
	const double *CD, *OD;
	double *ND;
	uint32_t Id;
	int Low, BNo, Off, Len;
	const int NB = S->NB;

	Low = LowestCrossover(S);

	CP = Par;
	OP = Par + NB;
	if(FastUniform(&S->Rng) < 0.5)
	{
		CP = Par + NB;
		OP = Par;
	}

	if(Low == S->L)
	{
		memcpy(Chr, CP, sizeof(uint32_t) * NB);
		return;
	}

	// Loci 0 .. Low from CP, the rest from OP.
	BNo = Low / S->B;
	Off = Low - BNo * S->B;
	Len = BlockLen(S, BNo);

	memcpy(Chr, CP, sizeof(uint32_t) * BNo);
	memcpy(Chr + BNo + 1, OP + BNo + 1, sizeof(uint32_t) * (NB - BNo - 1));

	if(Off == Len - 1 || CP[BNo] == OP[BNo])
	{
		Chr[BNo] = CP[BNo];
		return;
	}

	Id = NewBlock(S, Tag);
	ND = Block(S, Id);
	CD = Block(S, CP[BNo]);
	OD = Block(S, OP[BNo]);

	memcpy(ND, CD, sizeof(double) * (Off + 1));
	memcpy(ND + Off + 1, OD + Off + 1, sizeof(double) * (Len - Off - 1));

	Chr[BNo] = Id;
}

// Test if the population has died.
int DeadPopulation(OPTIONS *Opt, int PopSize)
{
	if(Opt->Ploidy == 1 && PopSize == 0)
		return TRUE;

	if(Opt->Ploidy == 2 && PopSize < 2)
		return TRUE;

	return FALSE;
}

static void WaitOutput(SIM *S)
{
	if(S->OutBusy == TRUE)
	{
		FastThreadsWait(S->OutPool);
		S->OutBusy = FALSE;
	}
}

// Crate a new population form the old.
static void NewPop(OPTIONS *Opt, SIM *S)
{
	int IIndex, P1, P2;
	uint32_t *Chr;
	const int NB = S->NB;

	SetSurviving(S);

	if(DeadPopulation(Opt, S->NoSurv) == TRUE)
	{
		WaitOutput(S);
		printf("Pop all died.\n");
		exit(0);
	}

	if(Opt->Ploidy == 1)
	{
		for(IIndex=0;IIndex<S->N;IIndex++)
		{
			P1 = S->Surv[FastUniformInt(&S->Rng, S->NoSurv)];
			memcpy(S->NIds + (size_t)IIndex * S->IndIds, S->Ids + (size_t)P1 * S->IndIds, sizeof(uint32_t) * NB);

			S->NEnv[IIndex] = (Opt->EnvVar != 0) ? gsl_ran_gaussian(&S->Rng.GRng, S->EnvSD) : 0.0;
		}

		return;
	}

	for(IIndex=0;IIndex<S->N;IIndex++)
	{
		// Two different parents (ValidDiploidParents; AssortativeMating is not
		// an option, so that is all it tests).
		do
		{
			P1 = (int)FastUniformInt(&S->Rng, S->NoSurv);
			P2 = (int)FastUniformInt(&S->Rng, S->NoSurv);
		} while(P1 == P2);

		Chr = S->NIds + (size_t)IIndex * S->IndIds;

		RecombinationPoissonChr(S, Chr, S->Ids + (size_t)S->Surv[P1] * S->IndIds, (uint32_t)(IIndex * 2 + 1));
		RecombinationPoissonChr(S, Chr + NB, S->Ids + (size_t)S->Surv[P2] * S->IndIds, (uint32_t)(IIndex * 2 + 2));

		S->NEnv[IIndex] = (Opt->EnvVar != 0) ? gsl_ran_gaussian(&S->Rng.GRng, S->EnvSD) : 0.0;
	}
}

// Mutate the next population: MutatePop, MutateIndividual, MutateIndividualChrNo
static void MutatePop(OPTIONS *Opt, SIM *S)
{
	int Index, ChrNo, No, Pos, BNo;
	uint32_t *Chr, Tag, Id, NId;
	const int NB = S->NB;

	for(Index=0;Index<S->N;Index++)
		for(ChrNo=0;ChrNo<S->P;ChrNo++)
		{
			No = (int)FastPoisson(&S->Rng, S->MutMu, S->MutEmu);

			if(No <= 0)
				continue;

			Chr = S->NIds + (size_t)Index * S->IndIds + (size_t)ChrNo * NB;
			Tag = (uint32_t)(Index * S->P + ChrNo + 1);

			while(No > 0)
			{
				Pos = (int)FastUniformInt(&S->Rng, Opt->NoQTL);
				BNo = Pos / S->B;

				// Mutate a block of this chromosome's own, or a copy.
				Id = Chr[BNo];
				if(S->Owner[Id] != Tag)
				{
					NId = NewBlock(S, Tag);
					memcpy(Block(S, NId), Block(S, Id), sizeof(double) * BlockLen(S, BNo));
					Chr[BNo] = NId;
					Id = NId;
				}

				Block(S, Id)[Pos - BNo * S->B] += gsl_ran_gaussian(&S->Rng.GRng, S->MutSD);
				No--;
			}
		}

	// The next generation shares these blocks.
	for(Index=0;Index<(int)S->NoNew;Index++)
		S->Owner[S->New[Index]] = 0;
	S->NoNew = 0;
}

// Write a line of output: Output and OuputReg. (The original also computes the
// variance of the environmental deviations and a correlation, and uses neither.)
static void WriteOutput(void *Ctx, int ThreadNo, int NoThreads)
{
	OUT_JOB *Job;
	double PhenotypeMean, PhenotypeVar, GenotypeVar, PhenotypeFullSD, Slope, Int;
	double *X;
	double a[4];
	int Index;

	(void)ThreadNo; (void)NoThreads;
	Job = (OUT_JOB*)Ctx;

	if(Job->HaveMean == TRUE)
		PhenotypeMean = Job->Mean;
	else
		PhenotypeMean = gsl_stats_mean(Job->Pheno, 1, Job->N);

	if(Job->UseM == TRUE)
		PhenotypeVar = gsl_stats_variance_m(Job->Pheno, 1, Job->N, PhenotypeMean);
	else
		PhenotypeVar = gsl_stats_variance(Job->Pheno, 1, Job->N);

	GenotypeVar = gsl_stats_variance(Job->Geno, 1, Job->N);

	// OuputReg: the regression of fitness on phenotype, phenotypes above the
	// optimum reflected below it
	X = Job->X;
	PhenotypeFullSD = gsl_stats_sd(X, 1, Job->NoS);

	for(Index=0;Index<Job->NoS;Index++)
		if(X[Index] > Job->Theta)
			X[Index] = Job->Theta - (X[Index] - Job->Theta);

	gsl_fit_linear(X, 1, Job->Y, 1, Job->NoS, &Int, &Slope, &a[0], &a[1], &a[2], &a[3]);

	printf("%d\t", Job->Itter);
	printf("%f\t", PhenotypeMean);
	printf("%f\t", PhenotypeMean / sqrt(PhenotypeVar));
	printf("%f\t", GenotypeVar);
	printf("%d\t", Job->NoS);
	printf("%f\t%f\t", Job->Theta, Job->Omega);
	printf("%f\t%f\t", Slope, Slope * PhenotypeFullSD);
	printf("%f\t", sqrt((GenotypeVar*GenotypeVar) / Job->MuteVar));
	printf("\n");

	fflush(stdout);
}

// Write info about current population to stdout
static void Output(OPTIONS *Opt, int Itter, SIM *S)
{
	OUT_JOB *Job;
	int Index, Async;

	WaitOutput(S);

	Job = &OutJob;
	Async = FastThreadsNo(S->OutPool) > 1;

	Job->Itter = Itter;
	Job->N = S->N;
	Job->NoS = S->NoSurv;
	Job->Theta = Opt->FitnessPar[0];
	Job->Omega = Opt->FitnessPar[1];
	Job->MuteVar = Opt->MuteVar;
	Job->HaveMean = S->HaveMean;
	Job->Mean = S->Mean;
	Job->UseM = S->UseM;

	if(Async)
	{
		memcpy(Job->Pheno, S->Pheno, sizeof(double) * S->N);
		memcpy(Job->Geno, S->Geno, sizeof(double) * S->N);
	}
	else
	{
		Job->Pheno = S->Pheno;
		Job->Geno = S->Geno;
	}

	for(Index=0;Index<S->NoSurv;Index++)
	{
		Job->X[Index] = S->Pheno[S->Surv[Index]];
		Job->Y[Index] = S->Fit[S->Surv[Index]];
	}

	if(Async)
	{
		S->OutBusy = TRUE;
		FastThreadsStart(S->OutPool, WriteOutput, Job);
	}
	else
		WriteOutput(Job, 0, 1);
}

static void FreeSim(SIM *S)
{
	int Index, NoThreads;

	WaitOutput(S);

	if(FastThreadsNo(S->OutPool) > 1)
	{
		free(OutJob.Pheno);
		free(OutJob.Geno);
		free(OutJob.X);
		free(OutJob.Y);
	}

	NoThreads = FastThreadsNo(S->Pool);
	FastThreadsFree(S->OutPool);
	FastThreadsFree(S->Pool);

	for(Index=0;Index<NoThreads;Index++)
		free(S->Mark[Index]);
	free(S->Mark);

	free(S->Data);
	free(S->Free);
	free(S->New);
	free(S->Owner);
	free(S->Ids);
	free(S->NIds);
	free(S->Geno);
	free(S->Pheno);
	free(S->Fit);
	free(S->Env);
	free(S->NEnv);
	free(S->Surv);
	free(S->Points);
	free(S->Seen);
	free(S);
}

// Main loop of the simulation
void RunSim(OPTIONS* Opt)
{
	int Itter;
	SIM *S;

	S = CreateSim(Opt);

	printf("Generation\tPhenotype Mean\tPhenotype Ratio\tGenotype Var\tP Survived\tFitness Function Mean\tFitness Function SD\tSurviving Slope\tSurviving Standardised Slope\tw\t");

	printf("\n");

	Opt->Drift = FALSE;
	for(Itter=0;Itter<Opt->NoGenerations;Itter++)
	{
		PopFitness(Opt, Itter, S);

		NewPop(Opt, S);

		if(Itter % Opt->Sample == 0)
			Output(Opt, Itter, S);

		if(Itter == Opt->ChangeFitnessGen)
			Opt->FitnessPar[1] = Opt->ChangeFitnessSD;

		MutatePop(Opt, S);

		swap((void**)&S->Ids, (void**)&S->NIds);
		swap((void**)&S->Env, (void**)&S->NEnv);
	}

	FreeSim(S);
}

// Calculate the analytical fitness SD from Phenotype variance, heritability and mutations scalar.
double GetAnalyticalFitSD(OPTIONS *Opt)
{
	double Num, Dom;

	Num = Opt->PheotypeVar * Opt->Heritabitlity;
	Dom = Opt->PheotypeVar * (1.0 - Opt->Heritabitlity);
	Dom = sqrt(Opt->MuteVarScalar * Dom);

	return Num / Dom;
}


// A function to calculate the analytical fitness SD if not user specified
void NLOptParam(OPTIONS *Opt)
{
	if(Opt->InitFitnessSD != -1)
	{
		Opt->FitnessPar[1] = Opt->InitFitnessSD;
		return;
	}

	Opt->FitnessPar[1] = GetAnalyticalFitSD(Opt);
	return;
}

// Main entry point
int main(int argc, char** argv)
{
	OPTIONS *Opt;

	if(argc != 2)
	{
		printf("%s requires a configuration file.\n", __FILE__);
		return 1;
	}

	Opt = CreateOptionsFile(argc, argv);

	NLOptParam(Opt);

	PrintOptions(Opt);

	RunSim(Opt);

	FreeOptions(Opt);
}
