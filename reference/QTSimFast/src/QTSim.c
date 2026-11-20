/*
*  QTSimFast -- a MODIFIED version of QTSim (QTLSim 1.0)
*
*  Modified 2026-11 for speed, for the QTLSim.jl project
*  (https://github.com/ecorecipes/QTLSim.jl), with Claude (Anthropic). The output
*  is bit-identical to the original's for the same options file and seed, unless
*  the environment variable QTSIM_FAST is set (see below and ../FAST.md). The
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
*   - functions the original never calls were removed;
*   - (2026-11-19) the three means of a line of output are computed in one loop
*     that is gsl_stats_mean's recurrence, checked against the library at start;
*     with two threads or more the draws are made by one thread and the genomes
*     built by another (see "The pipeline" below);
*   - (2026-11-20) QTSIM_FAST=1 or 2, off by default, gives up bit parity: level 1
*     sums genotypes from a cached sum per block (SumGenotypesCached), level 2
*     also draws from another generator by other methods (FastRNG.h). With
*     QTSIM_FAST unset nothing changes.
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
*
*  With one thread a generation is PopFitness, NewPop, Output, MutatePop, as in
*  the original. With two or more (the pipeline) NewPop and MutatePop only make
*  their draws, and a second thread does their work on the genomes.
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
	int		N, L, P;			// PopSize, NoQTL (1 if NoQTL is 0), Ploidy
	int		NoQTL;				// NoQTL as given

	// Genomes. A chromosome is NB block numbers; block b holds the alleles of
	// loci b*B .. b*B+B-1 (the last block LastLen of them). Chromosomes share
	// blocks until a crossover inside a block or a mutation makes a new one.
	int		B, NB, LastLen;
	size_t	Row;				// N * P, the chromosomes of a population
	size_t	IndIds;				// P * NB, block numbers per individual
	double	*Data;				// Cap blocks of B alleles
	uint32_t Cap, HW;			// capacity; blocks HW.. have never been used
	uint32_t *Free, NoFree;		// free blocks below HW
	uint32_t *Owner;			// for a block made this generation, the chromosome
								// that alone holds it (1-based), else 0
	uint32_t *New, NoNew;		// the blocks made this generation
	double	*BSum;				// QTSIM_FAST=1 only: each block's alleles summed (NULL otherwise)
	int		Fast;				// QTSIM_FAST: 0 the exact program, 1 genotypes from block sums,
								// 2 also other random numbers: see "QTSIM_FAST=2" below

	// QTSIM_FAST=2: every individual's draws in a generation come from a stream of
	// its own (XsKey), so they are made on the pool's threads and the output does
	// not depend on the number of threads.
	uint64_t XSeed;
	int		XGen;				// the generation whose draws are made
	int		HasEnv;				// EnvVar != 0
	uint8_t	*Alive;				// survival, before the survivors are listed
	struct XMUT { uint32_t Chr; int Pos; double Eff; } **XMut;	// per thread, in order
	size_t	*XNoMut, *XMutCap;
	uint32_t **XSeen, *XStamp;	// per thread: points drawn for a gamete, when many
	gsl_rng	*XG;				// per thread: a gsl_rng that reads the stream in use
	uint32_t *XMixId;			// the new block of each crossover inside a block, or UINT32_MAX
	double	XTheta, XOmega, XScale, XSD;	// the fitness function, for FitnessJobX
	int		XStd;
	uint32_t *Ids, *NIds;		// block numbers of the current and the next population:
								// block BNo of chromosome ChrNo of individual Ind is
								// Ids[BNo * Row + Ind * P + ChrNo]

	// The gametes of the next population, as drawn: the first Cut blocks come
	// from chromosome Src1 of the current population, the rest from Src2;
	// Mix lists the crossovers inside a block.
	uint32_t *Src1, *Src2;
	int		*Cut;
	struct { uint32_t G; int BNo, Off; } *Mix;
	int		NoMix;
	int		*Head, *Next;		// the gametes whose Cut is a given position, as lists
	uint32_t **Cur;				// one working array per thread

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
	int		UseMean;			// and ThreeMeans those of gsl_stats_mean
	int		Lanes;

	FAST_POOL *Pool;			// for genotypes
	int		ParGeno;			// sum genotypes on the pool's threads
	int		ParFill;			// copy block numbers on the pool's threads
	double	ParMin;				// smallest number of alleles summed on the pool
	FAST_POOL *OutPool;			// for output
	int		OutBusy;

	// What the pool's threads sum (SumGenotypesJob): the individuals
	// JobList[0 .. JobLen-1], or 0 .. JobLen-1 if JobList is NULL.
	const uint32_t *JobIds;
	double	*JobGeno;
	const int *JobList;
	int		JobLen;

	// The pipeline (see "The pipeline" below). The main thread makes the draws
	// and records them; the genome thread builds the next population from the
	// records.
	int		Pipe;				// the pipeline is in use
	int		Spec;				// the genome thread also sums the genotypes
	int		Batch;				// individuals it takes at a time
	FAST_THREAD *PipeThread;
	FAST_CHAN *Prog;			// set by the main thread: how far the draws are
	FAST_CHAN *Done;			// set by the genome thread: the generation it has finished
	FAST_CHAN *Helped;			// set by the main thread: the generation in which it has summed its share
	int		HelpLo, HelpHi;		// its share: List[HelpLo .. HelpHi-1]
	double	HelpMin;			// fewest alleles to sum for which it is asked
	FAST_CHAN *Fill, *Claim;		// within a generation: individuals built, individuals claimed for summing
	long	FBase, FProg;		// the generation's first value of Prog, and the last read
	int		FSum;				// genotypes are summed behind the draws in this generation
	int		ParSpec;			// and on the pool's threads
	long	PGen, Stride;		// Prog is PGen * Stride + progress within the generation
	int		PLast;				// this is the last generation: its offspring are not needed
	int		PStopped;
	int		*Off;				// per gamete: the crossover's place inside block Cut, or -1
	double	*NGeno;				// genotypes of the next population (Spec)
	uint8_t	*MutFlag;			// individuals a mutation has changed
	int		*List;				// individuals to sum
	struct { uint32_t Chr; int Pos; double Eff; } *Mut;	// the mutations, as drawn
	size_t	NoMut, MutCap;
	int		SpinMain, SpinGenome;
} SIM;

// The pipeline's constants. None changes a result.
#define PIPE_PUBLISH 32			// the main thread tells its progress every so many individuals
#define PIPE_BATCH 64			// individuals the genome thread takes at a time
#define PIPE_SPEC_MAX 0.7		// most mutations expected per individual for sums behind the draws
#define PIPE_HELP_MIN 40000		// fewest alleles left to sum for the main thread to take a share
#define PIPE_SPIN_MAIN 300		// microseconds, about, that the main thread polls before it sleeps
#define PIPE_SPIN_GENOME 100	// and the genome thread
#define PIPE_SPIN_POOL 200		// and a thread of the pool, between jobs

static void PipeMain(void *Ctx);

// What a line of output is computed from: a copy when it is written by another
// thread while the simulation goes on.
typedef struct
{
	int		Itter, N, NoS;
	double	Theta, Omega, MuteVar;
	int		HaveMean, UseM, UseMean;
	int		OnePass;			// QTSIM_FAST=2: the statistics in one pass
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

// The means of P[0..N-1], G[0..N-1] and X[0..NoS-1], NoS <= N, in one loop.
// Each is gsl_stats_mean's recurrence, written as GSL writes it:
//
//     long double mean = 0;
//     for(i = 0; i < size; i++) mean += (data[i * stride] - mean) / (i + 1);
//
// with a size_t i. A step is one subtraction, one division and one addition,
// so there is no product for a compiler to fuse into the addition, and each
// operation has one IEEE result. The three chains are independent; advancing
// them side by side lets the processor overlap their divisions, which each
// wait on the step before. StatsMSame compares the result with gsl_stats_mean
// of the library in use, and the plain calls are made if it ever differs.
static void ThreeMeans(const double *P, const double *G, size_t N, const double *X, size_t NoS, double *MeanP, double *MeanG, double *MeanX)
{
	long double MP = 0, MG = 0, MX = 0;
	size_t i;

	for(i=0;i<NoS;i++)
	{
		MP += (P[i] - MP) / (i + 1);
		MG += (G[i] - MG) / (i + 1);
		MX += (X[i] - MX) / (i + 1);
	}

	for(;i<N;i++)
	{
		MP += (P[i] - MP) / (i + 1);
		MG += (G[i] - MG) / (i + 1);
	}

	*MeanP = MP;
	*MeanG = MG;
	*MeanX = MX;
}

// Two checks of the GSL library in use, on test vectors.
//
// Does gsl_stats_variance_m(x, mean(x)) give the bits of gsl_stats_variance(x),
// and likewise gsl_stats_sd_m? In GSL's source they are the same code. If it
// ever differed the plain functions would be called instead (returns FALSE).
//
// Does ThreeMeans give the bits of gsl_stats_mean? If not, *MeanSame is FALSE
// and gsl_stats_mean is called instead.
static int StatsMSame(int *MeanSame)
{
	static const int Sizes[] = {1, 2, 3, 10, 257, 1000};
	double *V, *U, A, B, Mean, M3[3];
	gsl_rng *G;
	int SIndex, Index, Size, Part, Ret;

	V = (double*)SMalloc(sizeof(double) * 2000);
	U = V + 1000;
	G = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(G, 20261119);

	Ret = TRUE;
	*MeanSame = TRUE;
	for(SIndex=0;SIndex<(int)(sizeof(Sizes)/sizeof(Sizes[0]));SIndex++)
	{
		Size = Sizes[SIndex];
		for(Index=0;Index<Size;Index++)
		{
			V[Index] = gsl_ran_gaussian(G, 0.37) + 0.11;
			U[Index] = gsl_ran_gaussian(G, 1.0e-3) * (Index % 7 == 0 ? 1.0e6 : 1.0) - 3.0;
		}

		Mean = gsl_stats_mean(V, 1, Size);

		// (the third vector: a part of the second, as the survivors are a part)
		for(Part=1;Part<=Size;Part+=(Size+1)/2)
		{
			ThreeMeans(V, U, (size_t)Size, U, (size_t)Part, &M3[0], &M3[1], &M3[2]);

			if(memcmp(&M3[0], &Mean, sizeof(double)) != 0)
				*MeanSame = FALSE;

			A = gsl_stats_mean(U, 1, Size);
			if(memcmp(&M3[1], &A, sizeof(double)) != 0)
				*MeanSame = FALSE;

			A = gsl_stats_mean(U, 1, Part);
			if(memcmp(&M3[2], &A, sizeof(double)) != 0)
				*MeanSame = FALSE;
		}

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

// The default block size (QTSIM_BLOCK overrides it). Small enough that the blocks
// in use at one position fit the first-level cache, large enough that the block
// numbers are few beside the alleles.
static int DefBlockSize(OPTIONS *Opt)
{
	if(Opt->NoQTL <= 16)
		return Opt->NoQTL;

	if(Opt->NoQTL < 400)
		return 16;

	if(Opt->NoQTL < 2500)
		return 32;

	return 64;
}

// The default number of generations between collections of unused blocks
// (QTSIM_GC overrides it; any number gives the same output). Collecting less
// often saves the marking and the scan; with many loci the blocks left
// uncollected spread the alleles in use over more memory, which costs more in
// the sums than the collection saves (measured, see REPORT.md).
static int DefGCEvery(const SIM *S)
{
	if(S->L < 2500)
		return 4;

	return 1;
}

// A number of polls that takes about Micro microseconds.
static int SpinPolls(int Micro, long PerMicro)
{
	double Polls;

	Polls = (double)Micro * (double)PerMicro;
	if(Polls > 1.0e9)
		Polls = 1.0e9;

	return (int)Polls;
}

// Allocate the populations. Draws the environmental deviations of both of the
// original's populations, as CreatePop does twice.
static SIM* CreateSim(OPTIONS *Opt)
{
	SIM *S;
	size_t NoInd, Cap;
	int Index, NoThreads, MaxThreads, Requested, Perf;
	long PerMicro;
	char *Str;

	S = (SIM*)SMalloc(sizeof(SIM));
	memset(S, 0, sizeof(SIM));

	S->N = Opt->PopSize;
	S->NoQTL = Opt->NoQTL;
	S->P = Opt->Ploidy;

	if(S->N < 0 || S->NoQTL < 0)
	{
		// (The original fails on these too, each in its own way.)
		fprintf(stderr, "QTSimFast: PopSize and NoQTL must not be negative.\n");
		exit(1);
	}

	// With no loci the original runs, with genotypes of zero. Here that is one
	// locus that no crossover or mutation ever touches: the draws use NoQTL.
	S->L = (S->NoQTL == 0) ? 1 : S->NoQTL;

	NoInd = (size_t)S->N;

	// Block size: any size gives the same output.
	S->B = EnvInt("QTSIM_BLOCK", 0, 0, 1 << 30);
	if(S->B == 0)
		S->B = DefBlockSize(Opt);
	if(S->B > S->L || S->B < 1)
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
	S->Row = NoInd * S->P;
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

	// (block 0, all zeros, sums to 0)
	S->BSum = NULL;
	if(EnvInt("QTSIM_FAST", 0, 0, 2) != 0)
	{
		S->BSum = (double*)calloc(Cap, sizeof(double));
		if(S->BSum == NULL)
			MallocErr();
	}

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

	S->Src1 = (uint32_t*)SMalloc(sizeof(uint32_t) * S->Row);
	S->Src2 = (uint32_t*)SMalloc(sizeof(uint32_t) * S->Row);
	S->Cut = (int*)SMalloc(sizeof(int) * S->Row);
	S->Mix = SMalloc(sizeof(*S->Mix) * S->Row);
	S->Next = (int*)SMalloc(sizeof(int) * S->Row);
	S->Head = (int*)SMalloc(sizeof(int) * (S->NB + 1));
	for(Index=0;Index<=S->NB;Index++)
		S->Head[Index] = -1;

	S->Points = (int*)SMalloc(sizeof(int) * (S->L + 1));
	S->Seen = (uint32_t*)calloc(S->L + 1, sizeof(uint32_t));
	if(S->Seen == NULL)
		MallocErr();

	// Threads. QTSIM_THREADS: a number, or "auto". The default is 1.
	// Perf: the performance cores, where the system tells them apart (0 if not).
	MaxThreads = 256;
	Perf = FastThreadsPerfCPUs();
	Str = getenv("QTSIM_THREADS");
	if(Str != NULL && strcicmp(Str, "auto") == 0)
	{
		// As many as the performance cores, or half the processors, up to 4,
		// when there are enough alleles for the sums to gain from them; else 2,
		// for the pipeline.
		Requested = (Perf > 0) ? Perf : FastThreadsCPUs() / 2;
		if(Requested > 4)
			Requested = 4;
		if((double)S->N * S->P * S->L < 4.0e6 && Requested > 2)
			Requested = 2;
		if(Requested < 1)
			Requested = 1;
	}
	else
		Requested = EnvInt("QTSIM_THREADS", 1, 1, MaxThreads);

	// The pipeline: with two threads or more, unless QTSIM_PIPELINE=0. Then
	// the genome thread also sums the genotypes, behind the draws, when few of
	// them will have to be summed again: when an individual expects at most
	// PIPE_SPEC_MAX mutations (QTSIM_SPEC=0 or 1 decides otherwise).
	S->MutMu = Opt->NoQTL * Opt->MutationRatePerGamete;

	// QTSIM_FAST=1 gives up bit parity with QTSim: see SumGenotypesCached. The
	// pipeline exists to hide the exact sums behind the draws; with cached sums
	// it is not used.
	S->Fast = EnvInt("QTSIM_FAST", 0, 0, 2);
	if(S->Fast != 0)
		fprintf(stderr, "QTSimFast: QTSIM_FAST=%d: %s; the output is NOT bit-identical to QTSim's (see FAST.md)\n", S->Fast,
			S->Fast == 1 ? "genotypes are sums of cached block sums, added in another order" :
			"genotypes from cached block sums, and other random numbers (xoshiro256++, ziggurat normals, inversion Poisson, Lemire integers)");

	S->Pipe = FALSE;
	if(Requested > 1 && S->Fast == 0 && EnvInt("QTSIM_PIPELINE", 1, 0, 1) == 1)
		S->Pipe = TRUE;

	S->Spec = FALSE;
	if(S->Pipe == TRUE)
		S->Spec = EnvInt("QTSIM_SPEC", (S->P * S->MutMu <= PIPE_SPEC_MAX), 0, 1);

	// The threads that sum genotypes: no more than the performance cores
	// (more were slower when measured), unless QTSIM_POOL says how many. When
	// the genome thread sums, it is one of them and the main thread is not.
	NoThreads = Requested;
	if(Perf > 0 && NoThreads > Perf)
		NoThreads = Perf;
	if(S->Spec == TRUE)
		NoThreads--;
	NoThreads = EnvInt("QTSIM_POOL", NoThreads, 1, MaxThreads);

	if(NoThreads > S->N / LANES)
		NoThreads = S->N / LANES;
	if(NoThreads < 1)
		NoThreads = 1;

	S->Pool = FastThreadsCreate(NoThreads);
	NoThreads = FastThreadsNo(S->Pool);

	// Small populations are summed faster by one thread than handed out.
	S->ParMin = (double)EnvInt("QTSIM_PAR_MIN", 2000000, 0, 1 << 30);
	S->ParGeno = FALSE;
	if(NoThreads > 1 && (double)S->N * S->P * S->L >= S->ParMin)
		S->ParGeno = TRUE;

	S->ParFill = FALSE;
	if(S->ParGeno == TRUE && S->NB >= 4 * NoThreads)
		S->ParFill = TRUE;

	// Behind the draws the pool is worth its while for smaller populations:
	// its threads work while the main thread draws.
	S->ParSpec = FALSE;
	if(S->Spec == TRUE && NoThreads > 1 && (double)S->N * S->P * S->L >= S->ParMin / 4)
		S->ParSpec = TRUE;

	// How long waiting threads poll before they sleep (the pipeline's threads
	// wait often and briefly). Without the pipeline the pool is left as it was.
	// (QTSIM_FAST=2 hands the pool several short jobs per generation too)
	PerMicro = (S->Pipe == TRUE || (S->Fast == 2 && FastThreadsNo(S->Pool) > 1)) ? FastThreadsPollsPerMicro() : 1;
	if(S->Pipe == TRUE || (S->Fast == 2 && FastThreadsNo(S->Pool) > 1))
		FastThreadsSetSpin(S->Pool, SpinPolls(EnvInt("QTSIM_SPIN_POOL", PIPE_SPIN_POOL, 0, 1000000), PerMicro));

	// With more than one thread the output lines are computed and written by a
	// thread of their own (QTSIM_ASYNC_OUTPUT=0 turns that off).
	S->OutPool = FastThreadsCreate(EnvInt("QTSIM_ASYNC_OUTPUT", Requested > 1, 0, 1) ? 2 : 1);

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

	S->Cur = (uint32_t**)SMalloc(sizeof(uint32_t*) * NoThreads);
	for(Index=0;Index<NoThreads;Index++)
		S->Cur[Index] = (uint32_t*)SMalloc(sizeof(uint32_t) * S->Row);

	S->Mark = (uint8_t**)SMalloc(sizeof(uint8_t*) * NoThreads);
	for(Index=0;Index<NoThreads;Index++)
	{
		S->Mark[Index] = (uint8_t*)calloc(Cap, 1);
		if(S->Mark[Index] == NULL)
			MallocErr();
	}

	S->GCEvery = EnvInt("QTSIM_GC", DefGCEvery(S), 1, 1 << 30);
	S->Lanes = EnvInt("QTSIM_LANES", LANES, 1, LANES) == 1 ? 1 : LANES;
	S->UseM = StatsMSame(&S->UseMean);
	if(S->UseM == FALSE)
		S->UseMean = FALSE;
	if(EnvInt("QTSIM_PLAIN_STATS", 0, 0, 1) == 1)
		S->UseM = S->UseMean = FALSE;

	// The random numbers
	FastRngSet(&S->Rng, Opt->Seed);
	if(FastRngSelfCheck(&S->Rng, Opt->Seed) == FALSE)
	{
		fprintf(stderr, "QTSimFast: the built-in generator does not reproduce gsl_rng_mt19937; rebuild with -DQTSIM_GSL_RNG.\n");
		exit(3);
	}

	// QTSIM_FAST=2: another generator and other methods, from the first draw on
	if(S->Fast == 2 && FastRngSetX(&S->Rng, Opt->Seed) == 0)
	{
		fprintf(stderr, "QTSimFast: QTSIM_FAST=2 needs the built-in generator (this build uses GSL's).\n");
		exit(3);
	}

	if(S->Fast == 2)
	{
		int T, NT;

		NT = FastThreadsNo(S->Pool);
		S->XSeed = (uint64_t)Opt->Seed;
		S->HasEnv = (Opt->EnvVar != 0);
		S->Alive = (uint8_t*)SMalloc(NoInd + 1);
		S->XMixId = (uint32_t*)SMalloc(sizeof(uint32_t) * (S->Row + 1));
		S->XMut = SMalloc(sizeof(*S->XMut) * NT);
		S->XNoMut = (size_t*)SMalloc(sizeof(size_t) * NT);
		S->XMutCap = (size_t*)SMalloc(sizeof(size_t) * NT);
		S->XSeen = (uint32_t**)SMalloc(sizeof(uint32_t*) * NT);
		S->XStamp = (uint32_t*)SMalloc(sizeof(uint32_t) * NT);
		S->XG = (gsl_rng*)SMalloc(sizeof(gsl_rng) * NT);
		for(T=0;T<NT;T++)
		{
			S->XMutCap[T] = 1024;
			S->XNoMut[T] = 0;
			S->XMut[T] = SMalloc(sizeof(**S->XMut) * S->XMutCap[T]);
			S->XSeen[T] = (uint32_t*)calloc(S->L + 1, sizeof(uint32_t));
			if(S->XSeen[T] == NULL)
				MallocErr();
			S->XStamp[T] = 0;
			S->XG[T].type = S->Rng.GRng.type;	// reads a FAST_XS (FastRngSetX)
			S->XG[T].state = NULL;
		}
	}

	S->EnvSD = sqrt(Opt->EnvVar);

	S->RecMu = Opt->RecombPoissonMean;
	S->RecEmu = exp(-S->RecMu);

	S->MutEmu = exp(-S->MutMu);

	// MutateIndividualChrNo
	{
		double MutVar;

		MutVar = Opt->MuteVarScalar / (Opt->Ploidy * Opt->NoQTL * Opt->MutationRatePerGamete);
		MutVar = MutVar * Opt->EnvVar;
		S->MutSD = sqrt(MutVar);
	}

	// The pipeline's records and its thread.
	S->Off = (int*)SMalloc(sizeof(int) * (S->Row + 1));
	S->NGeno = (double*)SMalloc(sizeof(double) * (NoInd + 1));
	S->List = (int*)SMalloc(sizeof(int) * (NoInd + 1));
	S->MutFlag = (uint8_t*)calloc(NoInd + 1, 1);
	if(S->MutFlag == NULL)
		MallocErr();
	S->MutCap = 1024;
	S->Mut = SMalloc(sizeof(*S->Mut) * S->MutCap);

	S->Stride = (long)S->N + 2;
	S->Prog = FastChanCreate(0);
	S->Done = FastChanCreate(0);
	S->Helped = FastChanCreate(0);
	S->HelpMin = (double)EnvInt("QTSIM_HELP_MIN", PIPE_HELP_MIN, 0, 1 << 30);
	S->Fill = FastChanCreate(0);
	S->Claim = FastChanCreate(0);
	S->SpinMain = SpinPolls(EnvInt("QTSIM_SPIN_MAIN", PIPE_SPIN_MAIN, 0, 1000000), PerMicro);
	S->SpinGenome = SpinPolls(EnvInt("QTSIM_SPIN_GENOME", PIPE_SPIN_GENOME, 0, 1000000), PerMicro);

	S->Batch = EnvInt("QTSIM_BATCH", PIPE_BATCH, 1, 1 << 30);
	S->Batch = (S->Batch + LANES - 1) / LANES * LANES;

	if(S->Pipe == TRUE)
	{
		S->PipeThread = FastThreadSpawn(PipeMain, S);
		if(S->PipeThread == NULL)
			S->Pipe = S->Spec = S->ParSpec = FALSE;
	}

	// CreatePop(CPop), CreatePop(NPop): SetIndividualEnv for each individual.
	// The second population's values are never used, but they are drawn.
	for(Index=0;Index<S->N;Index++)
		S->Env[Index] = (Opt->EnvVar != 0) ? FastGaussian(&S->Rng, S->EnvSD) : 0.0;

	for(Index=0;Index<S->N;Index++)
		S->NEnv[Index] = (Opt->EnvVar != 0) ? FastGaussian(&S->Rng, S->EnvSD) : 0.0;

	return S;
}

// Sum the genotypes of individuals Lo .. Hi-1, W at a time: CaclGenotype.
// Genotype = 0, then Genotype += Genome[0][l] + Genome[1][l] for each locus l in
// order (diploid), or Genotype += Genome[0][l] (haploid).
//
// The sums advance one block of loci at a time over all the individuals, each
// individual's partial sum waiting in Geno, so that the blocks in use at one
// position of the chromosome (a few hundred at most, as a rule) stay in the
// first-level cache while every individual reads them.
//
// The chromosomes are read through the table of block numbers Ids and the
// sums are left in Geno. With a List the individuals are List[Lo .. Hi-1]
// instead of Lo .. Hi-1 (and Mark must be NULL): which individuals are summed
// together, and by which thread, never changes an individual's sum.
static inline __attribute__((always_inline)) void SumGenotypesW(const SIM *S, const uint32_t *Ids, double *Geno, const int *List, int Lo, int Hi, uint8_t *Mark, const int W, const int Diploid)
{
	const double *A[LANES], *C[LANES];
	const uint32_t *Row;
	double Sum[LANES];
	int Start, J, K, Index, Len, Ind;
	const int NB = S->NB;
	const int P = Diploid ? 2 : 1;

	for(Ind=Lo;Ind<Hi;Ind++)
		Geno[(List != NULL) ? List[Ind] : Ind] = 0;

	for(K=0;K<NB;K++)
	{
		Len = (K == NB - 1) ? S->LastLen : S->B;
		Row = Ids + (size_t)K * S->Row;

		if(Mark != NULL)
			for(Index=Lo*P;Index<Hi*P;Index++)
				Mark[Row[Index]] = 1;

		for(Start=Lo;Start<Hi;Start+=W)
		{
			// A short last group is filled up with its last individual.
			for(J=0;J<W;J++)
			{
				Ind = Start + J;
				if(Ind > Hi - 1)
					Ind = Hi - 1;
				if(List != NULL)
					Ind = List[Ind];

				Sum[J] = Geno[Ind];
				A[J] = Block(S, Row[(size_t)Ind * P]);
				if(Diploid)
					C[J] = Block(S, Row[(size_t)Ind * P + 1]);
			}

			Index = 0;

			if(Diploid)
			{
#ifdef FAST_VECTOR
				if(W > 1)
				{
					// Two loci of two individuals at a time: the alleles are
					// added in pairs, a + c, for loci Index and Index+1 of each
					// individual; then both individuals' sums take locus Index,
					// then locus Index+1. Every sum still adds its own terms
					// in order.
					V2D VSum[LANES / 2];

					for(J=0;J<W;J+=2)
						VSum[J / 2] = (V2D){Sum[J], Sum[J + 1]};

					for(;Index+1<Len;Index+=2)
						for(J=0;J<W;J+=2)
						{
							V2D P0 = Load2(A[J] + Index) + Load2(C[J] + Index);
							V2D P1 = Load2(A[J + 1] + Index) + Load2(C[J + 1] + Index);
							VSum[J / 2] += (V2D){P0[0], P1[0]};
							VSum[J / 2] += (V2D){P0[1], P1[1]};
						}

					for(J=0;J<W;J+=2)
					{
						Sum[J] = VSum[J / 2][0];
						Sum[J + 1] = VSum[J / 2][1];
					}
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

			for(J=0;J<W && Start+J<Hi;J++)
				Geno[(List != NULL) ? List[Start + J] : Start + J] = Sum[J];
		}
	}
}

// QTSIM_FAST=1: an individual's genotype is the sum over block positions of its
// chromosomes' cached block sums, (a_K + c_K) for diploids, instead of the sum of
// its alleles locus by locus. The same real number, rounded differently: results
// differ from the original's in the last bits, and runs diverge from there. Cost
// N * P * NB additions instead of N * P * L, plus B per block made.
static void SumGenotypesCached(const SIM *S, const uint32_t *Ids, double *Geno, const int *List, int Lo, int Hi, uint8_t *Mark)
{
	const double *BS = S->BSum;
	const uint32_t *Row;
	int K, I, Ind;
	const int P = S->P;

	for(I=Lo;I<Hi;I++)
		Geno[(List != NULL) ? List[I] : I] = 0;

	for(K=0;K<S->NB;K++)
	{
		Row = Ids + (size_t)K * S->Row;

		if(Mark != NULL)
			for(I=Lo*P;I<Hi*P;I++)
				Mark[Row[I]] = 1;

		if(P == 2)
			for(I=Lo;I<Hi;I++)
			{
				Ind = (List != NULL) ? List[I] : I;
				Geno[Ind] += BS[Row[(size_t)Ind * 2]] + BS[Row[(size_t)Ind * 2 + 1]];
			}
		else
			for(I=Lo;I<Hi;I++)
			{
				Ind = (List != NULL) ? List[I] : I;
				Geno[Ind] += BS[Row[Ind]];
			}
	}
}

// The sum of a block's alleles, kept in BSum (QTSIM_FAST=1).
static inline void SetBlockSum(SIM *S, uint32_t Id, int BNo)
{
	const double *D;
	double Sum;
	int I, Len;

	D = Block(S, Id);
	Len = BlockLen(S, BNo);
	Sum = 0;
	for(I=0;I<Len;I++)
		Sum += D[I];

	S->BSum[Id] = Sum;
}

static void SumGenotypesOf(const SIM *S, const uint32_t *Ids, double *Geno, const int *List, int Lo, int Hi, uint8_t *Mark)
{
	const int W = (S->Lanes == 1) ? 1 : LANES;

	if(S->BSum != NULL)
	{
		SumGenotypesCached(S, Ids, Geno, List, Lo, Hi, Mark);
		return;
	}

	// (one copy of the loops for each case, each with its constants known)
	if(List == NULL)
	{
		if(S->P == 2)
		{
			if(W == 1)
				SumGenotypesW(S, Ids, Geno, NULL, Lo, Hi, Mark, 1, 1);
			else
				SumGenotypesW(S, Ids, Geno, NULL, Lo, Hi, Mark, LANES, 1);
		}
		else
		{
			if(W == 1)
				SumGenotypesW(S, Ids, Geno, NULL, Lo, Hi, Mark, 1, 0);
			else
				SumGenotypesW(S, Ids, Geno, NULL, Lo, Hi, Mark, LANES, 0);
		}

		return;
	}

	if(S->P == 2)
	{
		if(W == 1)
			SumGenotypesW(S, Ids, Geno, List, Lo, Hi, NULL, 1, 1);
		else
			SumGenotypesW(S, Ids, Geno, List, Lo, Hi, NULL, LANES, 1);
	}
	else
	{
		if(W == 1)
			SumGenotypesW(S, Ids, Geno, List, Lo, Hi, NULL, 1, 0);
		else
			SumGenotypesW(S, Ids, Geno, List, Lo, Hi, NULL, LANES, 0);
	}
}

// The genotypes of individuals Lo .. Hi-1 of the current population.
static void SumGenotypes(const SIM *S, int Lo, int Hi, uint8_t *Mark)
{
	SumGenotypesOf(S, S->Ids, S->Geno, NULL, Lo, Hi, Mark);
}

// The part of 0 .. No-1 that thread ThreadNo of NoThreads works on, in whole
// multiples of Unit.
static void ThreadRange(int No, int Unit, int ThreadNo, int NoThreads, int *Lo, int *Hi)
{
	int Groups;

	Groups = (No + Unit - 1) / Unit;
	*Lo = (int)((long)Groups * ThreadNo / NoThreads) * Unit;
	*Hi = (int)((long)Groups * (ThreadNo + 1) / NoThreads) * Unit;
	if(*Hi > No)
		*Hi = No;
}

static void SumGenotypesJob(void *Ctx, int ThreadNo, int NoThreads)
{
	SIM *S;
	int Lo, Hi;

	S = (SIM*)Ctx;

	ThreadRange(S->JobLen, LANES, ThreadNo, NoThreads, &Lo, &Hi);

	if(Lo < Hi)
		SumGenotypesOf(S, S->JobIds, S->JobGeno, S->JobList, Lo, Hi, (S->Marking && S->JobList == NULL) ? S->Mark[ThreadNo] : NULL);
}

// Free the blocks no chromosome of the current population uses. Called when
// only the current population is alive and its marks have just been set, in
// the mark arrays of NoT threads.
static void CollectBlocks(SIM *S, int NoT)
{
	uint8_t *M, *O;
	uint32_t Id, HW, NewHW, NoFree, *Free;
	int T;

	HW = S->HW;
	M = S->Mark[0];

	for(T=1;T<NoT;T++)
	{
		O = S->Mark[T];
		for(Id=0;Id<HW;Id++)
			M[Id] |= O[Id];
		memset(O, 0, HW);
	}

	// The new high-water mark: one above the highest block in use.
	NewHW = 0;
	for(Id=HW;Id-- > 0;)
		if(M[Id])
		{
			NewHW = Id + 1;
			break;
		}

	// The free blocks below it, highest first, so that lower-numbered blocks
	// are reused first. Every number is written and the count advances only
	// past the unmarked ones: the same list as testing each mark with a branch.
	NoFree = 0;
	Free = S->Free;
	for(Id=NewHW;Id-- > 0;)
	{
		Free[NoFree] = Id;
		NoFree += (M[Id] == 0);
	}
	S->NoFree = NoFree;

	memset(M, 0, HW);
	S->HW = NewHW;
	S->SinceGC = 0;
}

// Is this a generation in which unused blocks are collected? Every GCEvery
// generations, and whenever a generation that made a new block for every block
// of every chromosome would run out. Asked once per generation, when the
// population about to be summed is the only one alive.
static int CollectNow(SIM *S)
{
	S->SinceGC++;

	if(S->SinceGC >= S->GCEvery || (size_t)(S->Cap - S->HW) + S->NoFree < (size_t)S->N * S->IndIds)
		return TRUE;

	return FALSE;
}

// Calculate the phenotype for each individual in the population.
static void CaclPopPhenotype(SIM *S, int Itter)
{
	int Index;

	// (When the genome thread sums the genotypes, it has summed these and
	// collected the unused blocks; the first population is summed here.)
	if(S->Spec == FALSE || Itter == 0)
	{
		S->Marking = CollectNow(S);

		if(S->ParGeno == TRUE)
		{
			S->JobIds = S->Ids;
			S->JobGeno = S->Geno;
			S->JobList = NULL;
			S->JobLen = S->N;
			FastThreadsRun(S->Pool, SumGenotypesJob, S);
		}
		else
			SumGenotypes(S, 0, S->N, S->Marking ? S->Mark[0] : NULL);

		if(S->Marking == TRUE)
			CollectBlocks(S, S->ParGeno ? FastThreadsNo(S->Pool) : 1);
	}

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

// QTSIM_FAST=2: the fitness loop of PopFitness on the pool's threads (each value
// computed as in that loop)
static void FitnessJobX(void *Ctx, int ThreadNo, int NoThreads)
{
	SIM *S;
	int Lo, Hi, Index;
	double Phenotype;

	S = (SIM*)Ctx;
	ThreadRange(S->N, 64, ThreadNo, NoThreads, &Lo, &Hi);

	for(Index=Lo;Index<Hi;Index++)
	{
		Phenotype = S->Pheno[Index] - S->XTheta;
		if(S->XStd)
			Phenotype = Phenotype / S->XSD;
		S->Fit[Index] = gsl_ran_gaussian_pdf(Phenotype, S->XOmega) / S->XScale;
	}
}

// Calculate the fitness of the population
static void PopFitness(OPTIONS *Opt, int Itter, SIM *S)
{
	int Index;
	double PhenotypeSD, Phenotype, Scale;

	CaclPopPhenotype(S, Itter);

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

	if(S->Fast == 2)
	{
		S->XTheta = Opt->FitnessPar[0];
		S->XOmega = Opt->FitnessPar[1];
		S->XScale = Scale;
		S->XSD = PhenotypeSD;
		S->XStd = Opt->StandardisePhenotype;
		FastThreadsRun(S->Pool, FitnessJobX, S);
		return;
	}

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
	{
		// if(gsl_rng_uniform(Opt->rng) < In->Fitness) ..., without a branch
		S->Surv[S->NoSurv] = Index;
		S->NoSurv += (FastUniform(&S->Rng) < S->Fit[Index]);
	}
}

// The lowest crossover point of a gamete, or NoQTL for none: CreateCombMap.
// The original sorts the points and then switches strand at the first only.
static inline int LowestCrossover(SIM *S)
{
	int NoPoints, Index, Point, Low, PIndex, Found;

	NoPoints = (int)FastPoisson(&S->Rng, S->RecMu, S->RecEmu);

	if(NoPoints == 0 || NoPoints >= S->NoQTL)
		return S->L;

	Low = S->L;

	if(NoPoints <= 16)
	{
		for(Index=0;Index<NoPoints;Index++)
		{
			do
			{
				Point = (int)FastUniformInt(&S->Rng, S->NoQTL - 1);

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
			Point = (int)FastUniformInt(&S->Rng, S->NoQTL - 1);
		} while(S->Seen[Point] == S->Stamp);

		S->Seen[Point] = S->Stamp;
		if(Point < Low)
			Low = Point;
	}

	return Low;
}

// The draws of a gamete, in the original order: RecombinationPoissonChr. Par
// is the parent's first chromosome and G the gamete, as columns of the tables
// of block numbers. Records the gamete as Src1, Src2 and Cut: its first Cut
// blocks are those of chromosome Src1 and the rest those of Src2. If the
// crossover is inside block Cut, returns its place there (alleles up to it
// from Src1, the others from Src2); else returns -1.
static inline int GameteDraws(SIM *S, uint32_t G, uint32_t Par)
{
	uint32_t CP, OP;
	int Low, BNo, Off, Second;

	Low = LowestCrossover(S);

	// CP = P1; if(gsl_rng_uniform(Opt->rng) < 0.5) CP = P2;
	Second = (FastUniform(&S->Rng) < 0.5);
	CP = Par + (uint32_t)Second;
	OP = Par + 1 - (uint32_t)Second;

	S->Src1[G] = CP;
	S->Src2[G] = OP;

	if(Low == S->L)
	{
		S->Cut[G] = S->NB;
		return -1;
	}

	// Loci 0 .. Low from CP, the rest from OP.
	BNo = Low / S->B;
	Off = Low - BNo * S->B;

	if(Off == BlockLen(S, BNo) - 1)
	{
		// Low is not the last locus, so BNo + 1 is a block.
		S->Cut[G] = BNo + 1;
		return -1;
	}

	// The crossover is inside block BNo.
	S->Cut[G] = BNo;
	return Off;
}

// Crate a gamete from a parent, using a number of recombination points, drawn
// from a Poisson: RecombinationPoissonChr. Only the draws are made here;
// FillGametes copies the block numbers and MixBlocks makes the blocks in which
// a crossover fell, from the lists kept here.
static inline void RecombinationPoissonChr(SIM *S, uint32_t G, uint32_t Par)
{
	int Off, Cut;

	Off = GameteDraws(S, G, Par);
	Cut = S->Cut[G];

	if(Cut == S->NB)
		return;

	S->Next[G] = S->Head[Cut];
	S->Head[Cut] = (int)G;

	if(Off >= 0)
	{
		S->Mix[S->NoMix].G = G;
		S->Mix[S->NoMix].BNo = Cut;
		S->Mix[S->NoMix].Off = Off;
		S->NoMix++;
	}
}

// The block numbers of the gametes at positions Lo .. Hi-1: the first Cut
// from one chromosome of the parent and the rest from the other.
static void FillGametes(const SIM *S, int Lo, int Hi, uint32_t *Cur)
{
	const uint32_t *PRow;
	uint32_t *Row;
	int K, G;
	const int No = (int)S->Row;

	// Cur: the chromosome each gamete is copying at position K
	for(G=0;G<No;G++)
		Cur[G] = (Lo < S->Cut[G]) ? S->Src1[G] : S->Src2[G];

	for(K=Lo;K<Hi;K++)
	{
		for(G=S->Head[K];G!=-1;G=S->Next[G])
			Cur[G] = S->Src2[G];

		PRow = S->Ids + (size_t)K * S->Row;
		Row = S->NIds + (size_t)K * S->Row;

		for(G=0;G<No;G++)
			Row[G] = PRow[Cur[G]];
	}
}

static void FillGametesJob(void *Ctx, int ThreadNo, int NoThreads)
{
	SIM *S;
	int Lo, Hi;

	S = (SIM*)Ctx;

	ThreadRange(S->NB, 1, ThreadNo, NoThreads, &Lo, &Hi);

	if(Lo < Hi)
		FillGametes(S, Lo, Hi, S->Cur[ThreadNo]);
}

// The blocks in which a crossover fell: alleles 0 .. Off from one chromosome
// of the parent, the rest from the other, in a new block unless the parent's
// two blocks are the same block.
static inline void MixBlock(SIM *S, uint32_t G, int BNo, int Off)
{
	const double *CD, *OD;
	double *ND;
	uint32_t CId, OId, Id;
	int I, Len;

	CId = S->Ids[(size_t)BNo * S->Row + S->Src1[G]];
	OId = S->Ids[(size_t)BNo * S->Row + S->Src2[G]];

	if(CId == OId)
		return;	// the block numbers are already right

	Id = NewBlock(S, G + 1);
	ND = Block(S, Id);
	CD = Block(S, CId);
	OD = Block(S, OId);
	Len = BlockLen(S, BNo);

	// (a select per allele, so that the loop has no branch on Off)
	for(I=0;I<Len;I++)
		ND[I] = (I <= Off) ? CD[I] : OD[I];

	if(S->BSum != NULL)
		SetBlockSum(S, Id, BNo);

	S->NIds[(size_t)BNo * S->Row + G] = Id;
}

static void MixBlocks(SIM *S)
{
	int Index;

	for(Index=0;Index<S->NoMix;Index++)
		MixBlock(S, S->Mix[Index].G, S->Mix[Index].BNo, S->Mix[Index].Off);

	S->NoMix = 0;

	for(Index=0;Index<S->NB;Index++)
		S->Head[Index] = -1;
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

// QTSIM_STATS=1: block statistics on standard error at the end of a run.
static void PrintStats(SIM *S, int Itter)
{
	if(EnvInt("QTSIM_STATS", 0, 0, 1) == 1)
		fprintf(stderr, "QTSimFast: generation %d, block size %d, blocks per chromosome %d, blocks in use or free below high-water mark %u (free %u), capacity %u, threads %d, pipeline %s\n", Itter, S->B, S->NB, S->HW, S->NoFree, S->Cap, S->ParGeno ? FastThreadsNo(S->Pool) : 1, S->Pipe ? (S->Spec ? "with sums" : "without sums") : "off");
}

static void WaitOutput(SIM *S)
{
	if(S->OutBusy == TRUE)
	{
		FastThreadsWait(S->OutPool);
		S->OutBusy = FALSE;
	}
}

/* ---------------------------------------------------------------------------
*  The pipeline (QTSIM_THREADS of 2 or more; QTSIM_PIPELINE=0 turns it off).
*
*  After SetSurviving, every draw of a generation (parents, crossover points,
*  strand, environmental deviation; then the mutations' numbers, places and
*  effects) depends on the survivors and the options and on no genome. So the
*  main thread makes all the draws, in the original order, and only records
*  them; a second thread, the genome thread, follows it and does everything
*  that touches genomes: the block numbers of the gametes, the blocks in which
*  a crossover fell, the mutations, and (when few individuals mutate) the
*  genotype sums. Who does the work and when changes no value: each gamete is
*  built from the same record, each mutation is the same addition to the same
*  allele in the original order, and each genotype is the same chain of
*  additions over the loci.
*
*  The two threads talk through two numbers. Prog, set by the main thread, is
*  PGen * Stride + k: in generation PGen the records of individuals 0 .. k-1 are
*  written (k <= N), and k = N + 1 says that the mutations are recorded too.
*  Done, set by the genome thread, is twice the last generation it has finished
*  (one less while it asks the main thread, which has made its draws by then,
*  to sum a share of the genotypes that remain). A
*  thread reads only what the other wrote before the number it has read, and
*  the main thread touches no genome between its first Prog of a generation and
*  the genome thread's Done.
* -------------------------------------------------------------------------*/

// The block numbers of the gametes of individuals Lo .. Hi-1 and the blocks
// in which their crossovers fell: FillGametes and MixBlocks, from the records.
static void FillBatch(SIM *S, int Lo, int Hi)
{
	const uint32_t *PRow, *Src1, *Src2;
	const int *Cut, *Off;
	uint32_t *Row;
	int K, G;
	const int GLo = Lo * S->P, GHi = Hi * S->P, NB = S->NB;

	Src1 = S->Src1;
	Src2 = S->Src2;
	Cut = S->Cut;
	Off = S->Off;

	for(K=0;K<NB;K++)
	{
		PRow = S->Ids + (size_t)K * S->Row;
		Row = S->NIds + (size_t)K * S->Row;

		for(G=GLo;G<GHi;G++)
			Row[G] = PRow[(K < Cut[G]) ? Src1[G] : Src2[G]];
	}

	for(G=GLo;G<GHi;G++)
		if(Off[G] >= 0)
			MixBlock(S, (uint32_t)G, Cut[G], Off[G]);
}

// The mutations of a generation, as the main thread drew them, in that order:
// the second half of MutatePop. If Flag, the individuals changed are flagged.
static void ApplyMutations(SIM *S, int Flag)
{
	uint32_t *Cell, Tag, Id, NId, Chr;
	size_t Index;
	int Pos, BNo;

	for(Index=0;Index<S->NoMut;Index++)
	{
		Chr = S->Mut[Index].Chr;
		Pos = S->Mut[Index].Pos;
		BNo = Pos / S->B;
		Tag = Chr + 1;

		// Mutate a block of this chromosome's own, or a copy.
		Cell = S->NIds + (size_t)BNo * S->Row + Chr;
		Id = *Cell;
		if(S->Owner[Id] != Tag)
		{
			NId = NewBlock(S, Tag);
			memcpy(Block(S, NId), Block(S, Id), sizeof(double) * BlockLen(S, BNo));
			*Cell = NId;
			Id = NId;
		}

		Block(S, Id)[Pos - BNo * S->B] += S->Mut[Index].Eff;

		if(Flag)
			S->MutFlag[Chr / (uint32_t)S->P] = 1;
	}

	// The next generation shares these blocks.
	for(Index=0;Index<S->NoNew;Index++)
		S->Owner[S->New[Index]] = 0;
	S->NoNew = 0;
}

// The first part of a generation behind the draws, as a job of the pool. The
// genome thread (ThreadNo 0) follows the main thread: it builds the gametes of
// the individuals whose draws are recorded, and tells how far it is in Fill.
// If genotypes are summed here (FSum), every thread of the job, the genome
// thread when it has nothing to build, claims PIPE_BATCH built individuals at
// a time from Claim and sums them, until the main thread has recorded the
// mutations; what is not summed by then is summed after them. Fill is the
// number of individuals built, and N + 1 once the claiming is over.
static void PipeFollowJob(void *Ctx, int ThreadNo, int NoThreads)
{
	SIM *S;
	const uint32_t *NIds;
	long V, F, Lo, Hi, PendLo, PendHi;
	int Filled;

	(void)NoThreads;
	S = (SIM*)Ctx;
	NIds = S->NIds;

	const long N = S->N, Base = S->FBase, Batch = S->Batch;

	if(ThreadNo != 0)
	{
		F = FastChanGet(S->Fill);
		for(;;)
		{
			if(F > N)
				return;

			Lo = FastChanGet(S->Claim);
			if(Lo >= N)
				return;

			Hi = (Lo + Batch < N) ? Lo + Batch : N;
			if(F < Hi)
			{
				F = FastChanWait(S->Fill, F, S->SpinGenome);
				continue;
			}

			Lo = FastChanClaim(S->Claim, Batch);
			if(Lo >= N)
				return;

			// (another thread may have claimed in between: wait for these to be built)
			Hi = (Lo + Batch < N) ? Lo + Batch : N;
			while(F < Hi)
				F = FastChanWait(S->Fill, F, S->SpinGenome);

			SumGenotypesOf(S, NIds, S->NGeno, NULL, (int)Lo, (int)Hi, NULL);
			F = FastChanGet(S->Fill);
		}
	}

	V = S->FProg;
	Filled = 0;
	PendLo = PendHi = 0;	// individuals this thread has claimed and not yet summed
	for(;;)
	{
		if(V - Base > N)
			break;	// all the draws are made, the mutations too

		if(V - Base > Filled)
		{
			FillBatch(S, Filled, (int)(V - Base));
			Filled = (int)(V - Base);
			FastChanSet(S->Fill, Filled);
			V = FastChanGet(S->Prog);
			continue;
		}

		if(S->FSum)
		{
			if(PendHi > PendLo)
			{
				if(Filled >= PendHi)
				{
					SumGenotypesOf(S, NIds, S->NGeno, NULL, (int)PendLo, (int)PendHi, NULL);
					PendLo = PendHi = 0;
					V = FastChanGet(S->Prog);
					continue;
				}
			}
			else
			{
				Lo = FastChanGet(S->Claim);
				Hi = (Lo + Batch < N) ? Lo + Batch : N;
				if(Lo < N && Filled >= Hi)
				{
					Lo = FastChanClaim(S->Claim, Batch);
					if(Lo < N)
					{
						PendLo = Lo;
						PendHi = (Lo + Batch < N) ? Lo + Batch : N;
					}
					continue;
				}
			}
		}

		V = FastChanWait(S->Prog, V, S->SpinGenome);
	}

	if(Filled < N)
	{
		FillBatch(S, Filled, (int)N);
		FastChanSet(S->Fill, N);
	}

	if(PendHi > PendLo)
		SumGenotypesOf(S, NIds, S->NGeno, NULL, (int)PendLo, (int)PendHi, NULL);

	FastChanSet(S->Fill, N + 1);
	S->FProg = V;
}

// One generation of the genome thread. V is the value of Prog that told it the
// generation had begun; returns the last value it read.
static long PipeGeneration(SIM *S, long V)
{
	const uint32_t *NIds;
	const int N = S->N;
	size_t Index, No;
	const long Gen = V / S->Stride;
	long Next;
	int Summed, Len, Ind, Collect, Help;

	NIds = S->NIds;

	// The gametes, and genotypes behind the draws. (The offspring of the last
	// generation are never looked at, so they are not summed.)
	S->FBase = Gen * S->Stride;
	S->FProg = V;
	S->FSum = (S->Spec == TRUE && S->PLast == FALSE);
	FastChanSet(S->Fill, 0);
	FastChanSet(S->Claim, 0);

	if(S->FSum && S->ParSpec == TRUE)
		FastThreadsRun(S->Pool, PipeFollowJob, S);
	else
		PipeFollowJob(S, 0, 1);

	V = S->FProg;

	// Every claim was summed: individuals 0 .. Summed-1 have genotypes, of
	// their genomes before mutation.
	Next = FastChanGet(S->Claim);
	Summed = (Next < N) ? (int)Next : N;

	ApplyMutations(S, S->FSum);

	Collect = FALSE;
	if(S->FSum)
	{
		// To sum now: the individuals summed before a mutation changed them,
		// and those not summed yet.
		Len = 0;
		for(Ind=0;Ind<Summed;Ind++)
		{
			S->List[Len] = Ind;
			Len += S->MutFlag[Ind];
		}
		for(Ind=Summed;Ind<N;Ind++)
			S->List[Len++] = Ind;

		for(Index=0;Index<S->NoMut;Index++)
			S->MutFlag[S->Mut[Index].Chr / (uint32_t)S->P] = 0;

		// The main thread has made its draws and waits: if there is enough
		// to sum, it takes a share, the last part of the list.
		Help = 0;
		if((double)Len * S->P * S->L >= S->HelpMin)
		{
			Help = Len / (FastThreadsNo(S->Pool) + 1) / LANES * LANES;
			if(Help > 0)
			{
				Len -= Help;
				S->HelpLo = Len;
				S->HelpHi = Len + Help;
				FastChanSet(S->Done, 2 * Gen - 1);
			}
		}

		if(S->ParSpec == TRUE && (double)Len * S->P * S->L >= S->ParMin / 8)
		{
			S->Marking = FALSE;
			S->JobIds = NIds;
			S->JobGeno = S->NGeno;
			S->JobList = S->List;
			S->JobLen = Len;
			FastThreadsRun(S->Pool, SumGenotypesJob, S);
		}
		else if(Len > 0)
			SumGenotypesOf(S, NIds, S->NGeno, S->List, 0, Len, NULL);

		if(Help > 0)
		{
			long Seen = FastChanGet(S->Helped);

			while(Seen != Gen)
				Seen = FastChanWait(S->Helped, Seen, S->SpinGenome);
		}

		// The population just built is now the only one whose blocks matter:
		// the point at which CaclPopPhenotype asks.
		Collect = CollectNow(S);
	}

	FastChanSet(S->Done, 2 * Gen);

	// Collecting needs nothing from the main thread, which goes on meanwhile
	// (it exchanges the two tables' pointers: NIds was read before).
	if(Collect == TRUE)
	{
		uint8_t *Mark = S->Mark[0];

		No = S->Row * (size_t)S->NB;
		for(Index=0;Index<No;Index++)
			Mark[NIds[Index]] = 1;

		CollectBlocks(S, 1);
	}

	return V;
}

static void PipeMain(void *Ctx)
{
	SIM *S;
	long V, Seen;

	S = (SIM*)Ctx;
	Seen = 0;

	for(;;)
	{
		V = FastChanWait(S->Prog, Seen, S->SpinGenome);
		if(V < 0)
			return;

		Seen = PipeGeneration(S, V);
	}
}

// Stop the genome thread; it finishes what it is doing first.
static void PipeStop(SIM *S)
{
	if(S->Pipe == TRUE && S->PStopped == FALSE)
	{
		FastChanSet(S->Prog, -1);
		FastThreadJoin(S->PipeThread);
		S->PStopped = TRUE;
	}
}

// NewPop's draws with the pipeline: the same draws in the same order, recorded
// for the genome thread.
static void NewPopDraws(OPTIONS *Opt, int Itter, SIM *S)
{
	long Base;
	int IIndex, P1, P2;

	S->PGen++;
	S->PLast = (Itter == Opt->NoGenerations - 1);
	Base = S->PGen * S->Stride;

	if(Opt->Ploidy == 1)
	{
		for(IIndex=0;IIndex<S->N;IIndex++)
		{
			P1 = S->Surv[FastUniformInt(&S->Rng, S->NoSurv)];
			S->Src1[IIndex] = (uint32_t)P1;
			S->Src2[IIndex] = (uint32_t)P1;
			S->Cut[IIndex] = S->NB;
			S->Off[IIndex] = -1;

			S->NEnv[IIndex] = (Opt->EnvVar != 0) ? FastGaussian(&S->Rng, S->EnvSD) : 0.0;

			if((IIndex & (PIPE_PUBLISH - 1)) == PIPE_PUBLISH - 1)
				FastChanSet(S->Prog, Base + IIndex + 1);
		}
	}
	else
	{
		for(IIndex=0;IIndex<S->N;IIndex++)
		{
			do
			{
				P1 = (int)FastUniformInt(&S->Rng, S->NoSurv);
				P2 = (int)FastUniformInt(&S->Rng, S->NoSurv);
			} while(P1 == P2);

			S->Off[IIndex * 2] = GameteDraws(S, (uint32_t)IIndex * 2, (uint32_t)S->Surv[P1] * 2);
			S->Off[IIndex * 2 + 1] = GameteDraws(S, (uint32_t)IIndex * 2 + 1, (uint32_t)S->Surv[P2] * 2);

			S->NEnv[IIndex] = (Opt->EnvVar != 0) ? FastGaussian(&S->Rng, S->EnvSD) : 0.0;

			if((IIndex & (PIPE_PUBLISH - 1)) == PIPE_PUBLISH - 1)
				FastChanSet(S->Prog, Base + IIndex + 1);
		}
	}

	FastChanSet(S->Prog, Base + S->N);
}

// MutatePop's draws with the pipeline: the first half of MutatePop. The number
// of mutations of each chromosome, and for each mutation its locus and then its
// effect, in the original order; the genome thread applies them.
static void MutateDraws(OPTIONS *Opt, SIM *S)
{
	long Seen;
	int Index, No, Pos, Helped;
	const int NoChr = (int)S->Row;

	S->NoMut = 0;
	Helped = FALSE;

	for(Index=0;Index<NoChr;Index++)
	{
		No = (int)FastPoisson(&S->Rng, S->MutMu, S->MutEmu);

		while(No > 0)
		{
			Pos = (int)FastUniformInt(&S->Rng, Opt->NoQTL);

			if(S->NoMut == S->MutCap)
			{
				S->MutCap *= 2;
				S->Mut = realloc(S->Mut, sizeof(*S->Mut) * S->MutCap);
				if(S->Mut == NULL)
					MallocErr();
			}

			S->Mut[S->NoMut].Chr = (uint32_t)Index;
			S->Mut[S->NoMut].Pos = Pos;
			S->Mut[S->NoMut].Eff = FastGaussian(&S->Rng, S->MutSD);
			S->NoMut++;
			No--;
		}
	}

	// The mutations are recorded: the genome thread finishes the generation.
	FastChanSet(S->Prog, S->PGen * S->Stride + S->N + 1);

	// Wait for it, and sum a share of the genotypes if it asks (Done is odd).
	Seen = 2 * S->PGen - 2;
	for(;;)
	{
		Seen = FastChanWait(S->Done, Seen, S->SpinMain);

		if(Seen == 2 * S->PGen)
			break;

		if(Seen == 2 * S->PGen - 1 && Helped == FALSE)
		{
			SumGenotypesOf(S, S->NIds, S->NGeno, S->List, S->HelpLo, S->HelpHi, NULL);
			Helped = TRUE;
			FastChanSet(S->Helped, S->PGen);
		}
	}
}

// Crate a new population form the old.
/* --------------------------------------------------------------------------
*  QTSIM_FAST=2: independent streams. Survival: one uniform per individual from
*  the key (Seed, generation, 1, individual). Offspring i: a xoshiro256++ stream
*  seeded from (Seed, generation, 2, i) makes all of its draws, as QTSim makes
*  them for one offspring (parents, gamete 0, gamete 1, environment), and then the
*  mutations of its chromosomes. The draws are made on the pool's threads, each
*  thread a contiguous range of offspring, and recorded; the genomes are then
*  built from the records in the offspring's order, by one thread.
* -------------------------------------------------------------------------*/

static void SurvivalJobX(void *Ctx, int ThreadNo, int NoThreads)
{
	SIM *S;
	int Lo, Hi, Index;
	double U;

	S = (SIM*)Ctx;
	ThreadRange(S->N, 64, ThreadNo, NoThreads, &Lo, &Hi);

	for(Index=Lo;Index<Hi;Index++)
	{
		U = (double)(XsMix(XsKey(S->XSeed, (uint64_t)S->XGen, 1, (uint64_t)Index)) >> 11) * 0x1.0p-53;
		S->Alive[Index] = (U < S->Fit[Index]);
	}
}

static void SetSurvivingX(SIM *S, int Itter)
{
	int Index;

	S->XGen = Itter;
	FastThreadsRun(S->Pool, SurvivalJobX, S);

	S->NoSurv = 0;
	for(Index=0;Index<S->N;Index++)
	{
		S->Surv[S->NoSurv] = Index;
		S->NoSurv += S->Alive[Index];
	}
}

static inline unsigned int PoissonX(FAST_XS *X, gsl_rng *G, double Mu, double Emu)
{
	if(Mu > 10)
	{
		G->state = X;
		return gsl_ran_poisson(G, Mu);
	}

	return XsPoisson(X, Mu, Emu);
}

// LowestCrossover from stream X, with thread T's table
static inline int LowestCrossoverX(SIM *S, FAST_XS *X, int T)
{
	int NoPoints, Index, Point, Low, PIndex, Found, Points[16];
	uint32_t *Seen, Stamp;

	NoPoints = (int)PoissonX(X, &S->XG[T], S->RecMu, S->RecEmu);

	if(NoPoints == 0 || NoPoints >= S->NoQTL)
		return S->L;

	Low = S->L;

	if(NoPoints <= 16)
	{
		for(Index=0;Index<NoPoints;Index++)
		{
			do
			{
				Point = (int)XsUniformInt(X, (uint32_t)(S->NoQTL - 1));

				Found = FALSE;
				for(PIndex=0;PIndex<Index;PIndex++)
					if(Points[PIndex] == Point)
						Found = TRUE;
			} while(Found == TRUE);

			Points[Index] = Point;
			if(Point < Low)
				Low = Point;
		}

		return Low;
	}

	Seen = S->XSeen[T];
	if(++S->XStamp[T] == 0)
	{
		memset(Seen, 0, sizeof(uint32_t) * (S->L + 1));
		S->XStamp[T] = 1;
	}
	Stamp = S->XStamp[T];

	for(Index=0;Index<NoPoints;Index++)
	{
		do
		{
			Point = (int)XsUniformInt(X, (uint32_t)(S->NoQTL - 1));
		} while(Seen[Point] == Stamp);

		Seen[Point] = Stamp;
		if(Point < Low)
			Low = Point;
	}

	return Low;
}

// GameteDraws from stream X; the crossover's place in block Cut, or -1, in Off[G]
static inline void GameteDrawsX(SIM *S, FAST_XS *X, int T, uint32_t G, uint32_t Par)
{
	int Low, BNo, Off, Second;

	Low = LowestCrossoverX(S, X, T);

	Second = (XsUniform(X) < 0.5);
	S->Src1[G] = Par + (uint32_t)Second;
	S->Src2[G] = Par + 1 - (uint32_t)Second;
	S->Off[G] = -1;

	if(Low == S->L)
	{
		S->Cut[G] = S->NB;
		return;
	}

	BNo = Low / S->B;
	Off = Low - BNo * S->B;

	if(Off == BlockLen(S, BNo) - 1)
	{
		S->Cut[G] = BNo + 1;
		return;
	}

	S->Cut[G] = BNo;
	S->Off[G] = Off;
}

static void OffspringJobX(void *Ctx, int ThreadNo, int NoThreads)
{
	SIM *S;
	FAST_XS X;
	gsl_rng *G;
	struct XMUT *M;
	int Lo, Hi, Ind, P1, P2, C, No;
	uint32_t Chr;

	S = (SIM*)Ctx;
	G = &S->XG[ThreadNo];
	ThreadRange(S->N, 1, ThreadNo, NoThreads, &Lo, &Hi);
	S->XNoMut[ThreadNo] = 0;

	for(Ind=Lo;Ind<Hi;Ind++)
	{
		XsSeed(&X, XsKey(S->XSeed, (uint64_t)S->XGen, 2, (uint64_t)Ind));

		if(S->P == 1)
		{
			P1 = S->Surv[XsUniformInt(&X, (uint32_t)S->NoSurv)];
			S->Src1[Ind] = (uint32_t)P1;
			S->Src2[Ind] = (uint32_t)P1;
			S->Cut[Ind] = S->NB;
			S->Off[Ind] = -1;
		}
		else
		{
			do
			{
				P1 = (int)XsUniformInt(&X, (uint32_t)S->NoSurv);
				P2 = (int)XsUniformInt(&X, (uint32_t)S->NoSurv);
			} while(P1 == P2);

			GameteDrawsX(S, &X, ThreadNo, (uint32_t)Ind * 2, (uint32_t)S->Surv[P1] * 2);
			GameteDrawsX(S, &X, ThreadNo, (uint32_t)Ind * 2 + 1, (uint32_t)S->Surv[P2] * 2);
		}

		S->NEnv[Ind] = S->HasEnv ? S->EnvSD * XsNormal(&X) : 0.0;

		// MutateIndividualChrNo for each of its chromosomes, recorded
		for(C=0;C<S->P;C++)
		{
			Chr = (uint32_t)(Ind * S->P + C);
			No = (int)PoissonX(&X, G, S->MutMu, S->MutEmu);
			for(;No>0;No--)
			{
				if(S->XNoMut[ThreadNo] == S->XMutCap[ThreadNo])
				{
					S->XMutCap[ThreadNo] *= 2;
					S->XMut[ThreadNo] = realloc(S->XMut[ThreadNo], sizeof(**S->XMut) * S->XMutCap[ThreadNo]);
					if(S->XMut[ThreadNo] == NULL)
						MallocErr();
				}
				M = &S->XMut[ThreadNo][S->XNoMut[ThreadNo]++];
				M->Chr = Chr;
				M->Pos = (int)XsUniformInt(&X, (uint32_t)S->NoQTL);
				M->Eff = S->MutSD * XsNormal(&X);
			}
		}
	}
}

static void MixBlocksX(SIM *S);

static void NewPopX(SIM *S, int Itter)
{
	int G;
	const int No = (int)S->Row;

	S->XGen = Itter;
	FastThreadsRun(S->Pool, OffspringJobX, S);

	// the lists of RecombinationPoissonChr, from the records, in order
	for(G=0;G<No;G++)
	{
		if(S->Cut[G] == S->NB)
			continue;

		S->Next[G] = S->Head[S->Cut[G]];
		S->Head[S->Cut[G]] = G;

		if(S->Off[G] >= 0)
		{
			S->Mix[S->NoMix].G = (uint32_t)G;
			S->Mix[S->NoMix].BNo = S->Cut[G];
			S->Mix[S->NoMix].Off = S->Off[G];
			S->NoMix++;
		}
	}

	if(S->ParFill == TRUE)
		FastThreadsRun(S->Pool, FillGametesJob, S);
	else
		FillGametes(S, 0, S->NB, S->Cur[0]);

	MixBlocksX(S);
}

// MixBlocks: the new blocks are numbered here, in order, and filled on the pool
static void MixFillJobX(void *Ctx, int ThreadNo, int NoThreads)
{
	SIM *S;
	int Lo, Hi, Index, I, Len, Off, BNo;
	uint32_t G, Id;
	const double *CD, *OD;
	double *ND;

	S = (SIM*)Ctx;
	ThreadRange(S->NoMix, 16, ThreadNo, NoThreads, &Lo, &Hi);

	for(Index=Lo;Index<Hi;Index++)
	{
		Id = S->XMixId[Index];
		if(Id == UINT32_MAX)
			continue;

		G = S->Mix[Index].G;
		BNo = S->Mix[Index].BNo;
		Off = S->Mix[Index].Off;
		CD = Block(S, S->Ids[(size_t)BNo * S->Row + S->Src1[G]]);
		OD = Block(S, S->Ids[(size_t)BNo * S->Row + S->Src2[G]]);
		ND = Block(S, Id);
		Len = BlockLen(S, BNo);

		for(I=0;I<Len;I++)
			ND[I] = (I <= Off) ? CD[I] : OD[I];

		SetBlockSum(S, Id, BNo);
	}
}

static void MixBlocksX(SIM *S)
{
	int Index, BNo;
	uint32_t G, CId, OId, Id;

	for(Index=0;Index<S->NoMix;Index++)
	{
		G = S->Mix[Index].G;
		BNo = S->Mix[Index].BNo;
		CId = S->Ids[(size_t)BNo * S->Row + S->Src1[G]];
		OId = S->Ids[(size_t)BNo * S->Row + S->Src2[G]];

		Id = UINT32_MAX;	// the block numbers are already right
		if(CId != OId)
		{
			Id = NewBlock(S, G + 1);
			S->NIds[(size_t)BNo * S->Row + G] = Id;
		}
		S->XMixId[Index] = Id;
	}

	FastThreadsRun(S->Pool, MixFillJobX, S);

	S->NoMix = 0;
	for(Index=0;Index<S->NB;Index++)
		S->Head[Index] = -1;
}

// MutatePop from the records, in the offspring's order
static void MutatePopX(SIM *S)
{
	int T, Pos, BNo;
	size_t K;
	uint32_t *Cell, Tag, Id, NId;
	const struct XMUT *M;

	for(T=0;T<FastThreadsNo(S->Pool);T++)
		for(K=0;K<S->XNoMut[T];K++)
		{
			M = &S->XMut[T][K];
			Pos = M->Pos;
			BNo = Pos / S->B;
			Tag = M->Chr + 1;

			Cell = S->NIds + (size_t)BNo * S->Row + M->Chr;
			Id = *Cell;
			if(S->Owner[Id] != Tag)
			{
				NId = NewBlock(S, Tag);
				memcpy(Block(S, NId), Block(S, Id), sizeof(double) * BlockLen(S, BNo));
				*Cell = NId;
				Id = NId;
			}

			Block(S, Id)[Pos - BNo * S->B] += M->Eff;
			SetBlockSum(S, Id, BNo);
		}

	for(K=0;K<S->NoNew;K++)
		S->Owner[S->New[K]] = 0;
	S->NoNew = 0;
}

static void NewPop(OPTIONS *Opt, int Itter, SIM *S)
{
	int IIndex, P1, P2;

	if(S->Fast == 2)
		SetSurvivingX(S, Itter);
	else
		SetSurviving(S);

	if(DeadPopulation(Opt, S->NoSurv) == TRUE)
	{
		PipeStop(S);
		WaitOutput(S);
		PrintStats(S, Itter);
		printf("Pop all died.\n");
		exit(0);
	}

	if(S->Pipe == TRUE)
	{
		NewPopDraws(Opt, Itter, S);
		return;
	}

	if(S->Fast == 2)
	{
		NewPopX(S, Itter);
		return;
	}

	if(Opt->Ploidy == 1)
	{
		// NewIn: a copy of one survivor's chromosome
		for(IIndex=0;IIndex<S->N;IIndex++)
		{
			P1 = S->Surv[FastUniformInt(&S->Rng, S->NoSurv)];
			S->Src1[IIndex] = (uint32_t)P1;
			S->Src2[IIndex] = (uint32_t)P1;
			S->Cut[IIndex] = S->NB;

			S->NEnv[IIndex] = (Opt->EnvVar != 0) ? FastGaussian(&S->Rng, S->EnvSD) : 0.0;
		}
	}
	else
	{
		for(IIndex=0;IIndex<S->N;IIndex++)
		{
			// Two different parents (ValidDiploidParents; AssortativeMating is
			// not an option, so that is all it tests).
			do
			{
				P1 = (int)FastUniformInt(&S->Rng, S->NoSurv);
				P2 = (int)FastUniformInt(&S->Rng, S->NoSurv);
			} while(P1 == P2);

			RecombinationPoissonChr(S, (uint32_t)IIndex * 2, (uint32_t)S->Surv[P1] * 2);
			RecombinationPoissonChr(S, (uint32_t)IIndex * 2 + 1, (uint32_t)S->Surv[P2] * 2);

			S->NEnv[IIndex] = (Opt->EnvVar != 0) ? FastGaussian(&S->Rng, S->EnvSD) : 0.0;
		}
	}

	// The genomes, which no draw depends on.
	if(S->ParFill == TRUE)
		FastThreadsRun(S->Pool, FillGametesJob, S);
	else
		FillGametes(S, 0, S->NB, S->Cur[0]);

	MixBlocks(S);
}

// Mutate the next population: MutatePop, MutateIndividual, MutateIndividualChrNo
static void MutatePop(OPTIONS *Opt, SIM *S)
{
	int Index, No, Pos, BNo;
	uint32_t *Cell, Tag, Id, NId;
	const int NoChr = (int)S->Row;

	// Chromosomes in the original order: individual by individual.
	for(Index=0;Index<NoChr;Index++)
	{
		No = (int)FastPoisson(&S->Rng, S->MutMu, S->MutEmu);

		Tag = (uint32_t)Index + 1;

		while(No > 0)
		{
			Pos = (int)FastUniformInt(&S->Rng, Opt->NoQTL);
			BNo = Pos / S->B;

			// Mutate a block of this chromosome's own, or a copy.
			Cell = S->NIds + (size_t)BNo * S->Row + Index;
			Id = *Cell;
			if(S->Owner[Id] != Tag)
			{
				NId = NewBlock(S, Tag);
				memcpy(Block(S, NId), Block(S, Id), sizeof(double) * BlockLen(S, BNo));
				*Cell = NId;
				Id = NId;
			}

			Block(S, Id)[Pos - BNo * S->B] += FastGaussian(&S->Rng, S->MutSD);
			if(S->BSum != NULL)
				SetBlockSum(S, Id, BNo);
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
	double GenotypeMean, XMean;
	double *X;
	double a[4];
	int Index;

	(void)ThreadNo; (void)NoThreads;
	Job = (OUT_JOB*)Ctx;

	// X: the survivors' phenotypes (OuputReg)
	X = Job->X;

	if(Job->OnePass == TRUE)
	{
		// QTSIM_FAST=2: one loop over the population and one over the survivors,
		// sums about the first value (shifted, for accuracy), instead of GSL's
		// passes; the regression of OuputReg in the same loop, its phenotypes
		// reflected at the optimum as below.
		double P0, G0, SP, SPP, SG, SGG, X0, R0, Y0, SX, SXX, SR, SRR, SY, SRY, D, E, R;
		const double N = (double)Job->N, NS = (double)Job->NoS;

		P0 = Job->Pheno[0];
		G0 = Job->Geno[0];
		SP = SPP = SG = SGG = 0;
		for(Index=0;Index<Job->N;Index++)
		{
			D = Job->Pheno[Index] - P0;
			E = Job->Geno[Index] - G0;
			SP += D;
			SPP += D * D;
			SG += E;
			SGG += E * E;
		}
		PhenotypeMean = P0 + SP / N;
		PhenotypeVar = (SPP - SP * SP / N) / (N - 1);
		GenotypeVar = (SGG - SG * SG / N) / (N - 1);

		X0 = X[0];
		R0 = (X0 > Job->Theta) ? Job->Theta - (X0 - Job->Theta) : X0;
		Y0 = Job->Y[0];
		SX = SXX = SR = SRR = SY = SRY = 0;
		for(Index=0;Index<Job->NoS;Index++)
		{
			D = X[Index] - X0;
			R = ((X[Index] > Job->Theta) ? Job->Theta - (X[Index] - Job->Theta) : X[Index]) - R0;
			E = Job->Y[Index] - Y0;
			SX += D;
			SXX += D * D;
			SR += R;
			SRR += R * R;
			SY += E;
			SRY += R * E;
		}
		PhenotypeFullSD = sqrt((SXX - SX * SX / NS) / (NS - 1));
		Slope = (SRY - SR * SY / NS) / (SRR - SR * SR / NS);

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
		return;
	}

	if(Job->UseMean == TRUE)
	{
		// The three means in one loop, then GSL's variance and SD about a
		// given mean: gsl_stats_variance(x) is gsl_stats_variance_m(x,
		// gsl_stats_mean(x)) in GSL, and StatsMSame has checked both steps.
		ThreeMeans(Job->Pheno, Job->Geno, (size_t)Job->N, X, (size_t)Job->NoS, &PhenotypeMean, &GenotypeMean, &XMean);

		PhenotypeVar = gsl_stats_variance_m(Job->Pheno, 1, Job->N, PhenotypeMean);
		GenotypeVar = gsl_stats_variance_m(Job->Geno, 1, Job->N, GenotypeMean);
		PhenotypeFullSD = gsl_stats_sd_m(X, 1, Job->NoS, XMean);
	}
	else
	{
		if(Job->HaveMean == TRUE)
			PhenotypeMean = Job->Mean;
		else
			PhenotypeMean = gsl_stats_mean(Job->Pheno, 1, Job->N);

		if(Job->UseM == TRUE)
			PhenotypeVar = gsl_stats_variance_m(Job->Pheno, 1, Job->N, PhenotypeMean);
		else
			PhenotypeVar = gsl_stats_variance(Job->Pheno, 1, Job->N);

		GenotypeVar = gsl_stats_variance(Job->Geno, 1, Job->N);
		PhenotypeFullSD = gsl_stats_sd(X, 1, Job->NoS);
	}

	// OuputReg: the regression of fitness on phenotype, phenotypes above the
	// optimum reflected below it

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
	Job->UseMean = S->UseMean;
	Job->OnePass = (S->Fast == 2);

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

	PipeStop(S);
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
	{
		free(S->Mark[Index]);
		free(S->Cur[Index]);
	}
	free(S->Mark);
	free(S->Cur);
	free(S->Head);
	free(S->Next);

	free(S->Data);
	free(S->Free);
	free(S->New);
	free(S->Owner);
	free(S->BSum);
	if(S->Fast == 2)
	{
		for(Index=0;Index<NoThreads;Index++)
		{
			free(S->XMut[Index]);
			free(S->XSeen[Index]);
		}
		free(S->XMut);
		free(S->XNoMut);
		free(S->XMutCap);
		free(S->XSeen);
		free(S->XStamp);
		free(S->XG);
		free(S->Alive);
		free(S->XMixId);
	}
	free(S->Ids);
	free(S->NIds);
	free(S->Geno);
	free(S->Pheno);
	free(S->Fit);
	free(S->Env);
	free(S->NEnv);
	free(S->Surv);
	free(S->Src1);
	free(S->Src2);
	free(S->Cut);
	free(S->Mix);
	free(S->Points);
	free(S->Seen);
	free(S->Off);
	free(S->NGeno);
	free(S->List);
	free(S->MutFlag);
	free(S->Mut);
	FastChanFree(S->Prog);
	FastChanFree(S->Done);
	FastChanFree(S->Helped);
	FastChanFree(S->Fill);
	FastChanFree(S->Claim);
	free(S);
}

// Main loop of the simulation
void RunSim(OPTIONS* Opt)
{
	int Itter;
	SIM *S;

	// A population of no one (the original allocates nothing, draws nothing
	// and finds no survivor in the first generation).
	if(Opt->PopSize == 0)
	{
		printf("Generation\tPhenotype Mean\tPhenotype Ratio\tGenotype Var\tP Survived\tFitness Function Mean\tFitness Function SD\tSurviving Slope\tSurviving Standardised Slope\tw\t");

		printf("\n");

		if(Opt->NoGenerations > 0)
		{
			printf("Pop all died.\n");
			exit(0);
		}

		return;
	}

	S = CreateSim(Opt);

	printf("Generation\tPhenotype Mean\tPhenotype Ratio\tGenotype Var\tP Survived\tFitness Function Mean\tFitness Function SD\tSurviving Slope\tSurviving Standardised Slope\tw\t");

	printf("\n");

	Opt->Drift = FALSE;
	for(Itter=0;Itter<Opt->NoGenerations;Itter++)
	{
		PopFitness(Opt, Itter, S);

		NewPop(Opt, Itter, S);

		if(Itter % Opt->Sample == 0)
			Output(Opt, Itter, S);

		if(Itter == Opt->ChangeFitnessGen)
			Opt->FitnessPar[1] = Opt->ChangeFitnessSD;

		if(S->Pipe == TRUE)
			MutateDraws(Opt, S);
		else if(S->Fast == 2)
			MutatePopX(S);
		else
			MutatePop(Opt, S);

		swap((void**)&S->Ids, (void**)&S->NIds);
		swap((void**)&S->Env, (void**)&S->NEnv);
		if(S->Spec == TRUE)
			swap((void**)&S->Geno, (void**)&S->NGeno);
	}

	PipeStop(S);
	PrintStats(S, Itter);
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
