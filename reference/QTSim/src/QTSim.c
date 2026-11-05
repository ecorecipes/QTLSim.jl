/*
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

#include "QTSim.h"
#include "RandLib.h"
#include "PassOptionFile.h"

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

// Set the environment variance of a newly created individual
void SetIndividualEnv(OPTIONS* Opt, INDIVIDUAL* Ind)
{
#ifdef NO_ENV_VAR
	Ind->Env = 0.0;
	return;
#endif

	if(Opt->EnvVar != 0)
		Ind->Env = RandNormalVar(Opt->rng, Opt->EnvVar);
	else
		Ind->Env = 0.0;
}

// Create an individual
INDIVIDUAL* CreateIndividual(OPTIONS* Opt)
{
	INDIVIDUAL* Ind;
	int Index;

	Ind = (INDIVIDUAL*)SMalloc(sizeof(INDIVIDUAL));

	Ind->Genome = (double**)SMalloc(sizeof(double*) * Opt->Ploidy);
	Ind->Genome[0] = (double*)SMalloc(sizeof(double) * Opt->Ploidy * Opt->NoQTL);


	for(Index=1;Index<Opt->Ploidy;Index++)
		Ind->Genome[Index] = Ind->Genome[0] + Index * Opt->NoQTL;

	for(Index=0;Index<Opt->Ploidy*Opt->NoQTL;Index++)
		Ind->Genome[0][Index] = 0.0;

	SetIndividualEnv(Opt, Ind);

	return Ind;
}

// Free the memory for an individual
void FreeIndividual(INDIVIDUAL* Ind)
{
	free(Ind->Genome[0]);
	free(Ind->Genome);
	free(Ind);
}

// Write an individual to the screen, debugging only
void DumpIndividual(INDIVIDUAL* Ind)
{
	printf("%f\t%f\t%f\t%f\t%f\t%f\n", Ind->Fitness, Ind->Phenotype, Ind->Genotype, Ind->Env, Ind->Genome[0][0], Ind->Genome[0][1]);
}

// Write the population to the screen, debugging only
void DumpPop(OPTIONS* Opt, POP *Pop)
{
	int Index;

	for(Index=0;Index<Opt->PopSize;Index++)
		DumpIndividual(Pop->Pop[Index]);
}

// Allocate and create the initial population from the specified options.
POP* CreatePop(OPTIONS* Opt)
{
	POP *Pop;
	int Index;

	Pop = (POP*)SMalloc(sizeof(POP));

	Pop->Pop = (INDIVIDUAL**)SMalloc(sizeof(INDIVIDUAL*) * Opt->PopSize);
	Pop->SPop = (INDIVIDUAL**)SMalloc(sizeof(INDIVIDUAL*) * Opt->PopSize);

	for(Index=0;Index<Opt->PopSize;Index++)
		Pop->Pop[Index] = CreateIndividual(Opt);


	Pop->TVect = (double*)SMalloc(sizeof(double) * Opt->PopSize);

	Pop->NoSPop = 0;
	return Pop;
}

// Free the memory for a given population.
void FreePop(OPTIONS *Opt, POP *Pop)
{
	int Index;


	free(Pop->TVect);
	free(Pop->SPop);

	for(Index=0;Index<Opt->PopSize;Index++)
		FreeIndividual(Pop->Pop[Index]);

	free(Pop->Pop);
	free(Pop);

}

// Calculate the Genotype for an individual
double  CaclGenotype(OPTIONS* Opt, INDIVIDUAL* Ind)
{
	double Genotype;
	int Index;

	Genotype = 0;
	if(Opt->Ploidy == 1)
	{
		for(Index=0;Index<Opt->NoQTL;Index++)
			Genotype += Ind->Genome[0][Index];

	}

	if(Opt->Ploidy == 2)
	{
		for(Index=0;Index<Opt->NoQTL;Index++)
			Genotype += Ind->Genome[0][Index] + Ind->Genome[1][Index];

	}

	// to controle for QTL
//	Genotype  = Genotype / Opt->NoQTL;

#ifdef QTL_AVE
		Genotype  = Genotype / (Opt->NoQTL * Opt->Ploidy);
#endif

	return Genotype;
}

// Calculate the Phenotype for an individual
void CaclPhenotype(OPTIONS* Opt, INDIVIDUAL* Ind)
{
	Ind->Genotype = CaclGenotype(Opt, Ind);
	Ind->Phenotype = Ind->Genotype + Ind->Env;
}

// Calculate the phenotype for each individual in the population.
void CaclPopPhenotype(OPTIONS* Opt, POP* Pop)
{
	int Index;
	INDIVIDUAL* Ind;

	for(Index=0;Index<Opt->PopSize;Index++)
	{
		Ind = Pop->Pop[Index];
		CaclPhenotype(Opt, Ind);

		Pop->TVect[Index] = Ind->Phenotype;
	}

}

// Not used
void MutateIndividualChrNoOld(OPTIONS* Opt, int ChrNo, INDIVIDUAL* Ind)
{
	int Index;
	double Scale;

	Scale = 1;

	if(Opt->Ploidy == 2)
		Scale = sqrt(2.0);

	for(Index=0;Index<Opt->NoQTL;Index++)
		Ind->Genome[ChrNo][Index] += (RandNormalVar(Opt->rng, Opt->MuteVar) / Scale) * sqrt(Opt->NoQTL);
}

// Mutate a chromosome
void MutateIndividualChrNo(OPTIONS* Opt, int ChrNo, INDIVIDUAL* Ind)
{
	int Pos, No;
	double MutVar;

	No = gsl_ran_poisson(Opt->rng, Opt->NoQTL * Opt->MutationRatePerGamete);


//	if(gsl_rng_uniform(Opt->rng) > Opt->NoQTL * Opt->MutationRatePerGamete)
//		return;

	MutVar = Opt->MuteVarScalar / (Opt->Ploidy * Opt->NoQTL * Opt->MutationRatePerGamete);
	MutVar = MutVar * Opt->EnvVar;

	while(No > 0)
	{
		Pos = gsl_rng_uniform_int(Opt->rng, Opt->NoQTL);
		Ind->Genome[ChrNo][Pos] += RandNormalVar(Opt->rng, MutVar);
		No--;
	}
}

// Mutate an individual
void MutateIndividual(OPTIONS* Opt, INDIVIDUAL* Ind)
{
	int Index;

	for(Index=0;Index<Opt->Ploidy;Index++)
		MutateIndividualChrNo(Opt, Index, Ind);
//		MutateIndividualChrNoPoisson(Opt, Index, Ind);
//		MutateIndividualChrNoOld(Opt, Index, Ind);
}

// Mutate the population
void MutatePop(OPTIONS* Opt, POP* Pop)
{
	int Index;

	for(Index=0;Index<Opt->PopSize;Index++)
		MutateIndividual(Opt, Pop->Pop[Index]);
}


double*	GetChr(OPTIONS* Opt, POP *CPop)
{
	INDIVIDUAL *CIn;

	CIn = CPop->SPop[gsl_rng_uniform_int(Opt->rng, CPop->NoSPop)];

	if(Opt->Ploidy == 1)
		return CIn->Genome[0];

	return CIn->Genome[gsl_rng_uniform_int(Opt->rng, Opt->Ploidy)];
}

// Find the surviving members of the population
void SetSurviving(OPTIONS* Opt, POP* Pop)
{
	INDIVIDUAL *In;
	int Index;

	Pop->NoSPop = 0;
	for(Index=0;Index<Opt->PopSize;Index++)
	{
		In = Pop->Pop[Index];

		if(gsl_rng_uniform(Opt->rng) < In->Fitness)
			Pop->SPop[Pop->NoSPop++] = In;

	}
}

// Not used
int AcceptMatting(double A, double B, double MaxDiff)
{
	double Diff;

	Diff = fabs(A-B);

	if(Diff < MaxDiff)
		return TRUE;


	return FALSE;
}

// Not used
int ValidDiploidParents(OPTIONS* Opt, INDIVIDUAL *P1, INDIVIDUAL *P2, double PhenotypeSD)
{
	if(P1 == P2)
		return FALSE;

	if(Opt->AssortativeMating == -1)
		return TRUE;

//	if(fabs(P1->Fitness - P2->Fitness) > Opt->AssortativeMating)
//		return FALSE;

	return AcceptMatting(P1->Phenotype, P2->Phenotype, PhenotypeSD * Opt->AssortativeMating);
}

// Not used
void Recombination(OPTIONS* Opt, INDIVIDUAL *CIn, INDIVIDUAL *P1, INDIVIDUAL *P2)
{
	int Index, Chr1, Chr2;

	// P1 and P2 are randomly assigned, so it does not matter.

	for(Index=0;Index<Opt->NoQTL;Index++)
	{
		Chr1 = gsl_rng_uniform_int(Opt->rng, Opt->Ploidy);
		CIn->Genome[0][Index] = P1->Genome[Chr1][Index];

		Chr2 = gsl_rng_uniform_int(Opt->rng, Opt->Ploidy);
		CIn->Genome[1][Index] = P2->Genome[Chr2][Index];
	}
}

// Not used
void RecombinationUnEnven(OPTIONS* Opt, INDIVIDUAL *CIn, INDIVIDUAL *P1, INDIVIDUAL *P2)
{
	int Index;
	double *P1D, *P1R, *P2D, *P2R;

	P1D = P1->Genome[0];
	P1R = P1->Genome[1];
	if(gsl_rng_uniform(Opt->rng) < 0.5)
	{
		P1D = P1->Genome[1];
		P1R = P1->Genome[0];
	}

	P2D = P2->Genome[0];
	P2R = P2->Genome[1];
	if(gsl_rng_uniform(Opt->rng) < 0.5)
	{
		P2D = P2->Genome[1];
		P2R = P2->Genome[0];
	}

	// P1 and P2 are randomly assigned, so it does not matter.
	for(Index=0;Index<Opt->NoQTL;Index++)
	{

		if(gsl_rng_uniform(Opt->rng) < RECOMB_BIAS)
			CIn->Genome[0][Index] = P1D[Index];
		else
			CIn->Genome[0][Index] = P1R[Index];


		if(gsl_rng_uniform(Opt->rng) < RECOMB_BIAS)
			CIn->Genome[1][Index] = P2D[Index];
		else
			CIn->Genome[1][Index] = P2R[Index];

//		Chr1 = gsl_rng_uniform_int(Opt->rng, Opt->Ploidy);
//		CIn->Genome[0][Index] = P1->Genome[Chr1][Index];

//		Chr2 = gsl_rng_uniform_int(Opt->rng, Opt->Ploidy);
//		CIn->Genome[1][Index] = P2->Genome[Chr2][Index];
	}
}

// Not used
void RecombinationChr1Point(OPTIONS* Opt, double *In, double *P1, double *P2)
{
	double *CP;
	int Index, Point;

	CP = P1;
	if(gsl_rng_uniform(Opt->rng) < 0.5)
		CP = P2;

	Point = 0;
	if(Opt->NoQTL > 2)
		Point = gsl_rng_uniform_int(Opt->rng, Opt->NoQTL-1);

	for(Index=0;Index<Opt->NoQTL;Index++)
	{
		In[Index] = CP[Index];

		if(Index == Point)
		{
			if(CP == P1)
				CP = P2;
			else
				CP = P1;
		}
	}
}

// Not used
void Recombination1Point(OPTIONS* Opt, INDIVIDUAL *CIn, INDIVIDUAL *P1, INDIVIDUAL *P2)
{
	RecombinationChr1Point(Opt, CIn->Genome[0], P1->Genome[0], P1->Genome[1]);
	RecombinationChr1Point(Opt, CIn->Genome[1], P2->Genome[0], P2->Genome[1]);
}

// Test if a recombination point is already in use
int PointInMap(int *Map, int Point, int Size)
{
	int Index;
	for(Index=0;Index<Size;Index++)
		if(Point == Map[Index])
			return TRUE;

	return FALSE;
}

// Compare two integers for sorting.
int IntComp(const void *A, const void  *B)
{
   return ( *(int*)A - *(int*)B);
}

// Crate a map of recombination locations
int*	CreateCombMap(OPTIONS* Opt, int *NoPoints)
{
	int Index;
	int *Map;

	*NoPoints = gsl_ran_poisson(Opt->rng, Opt->RecombPoissonMean);

	if(*NoPoints == 0 || *NoPoints >= Opt->NoQTL)
	{
		Map = (int*)SMalloc(sizeof(int));
		Map[0] = Opt->NoQTL;
		return Map;
	}

	Map = (int*)SMalloc(sizeof(int) * (*NoPoints+1));
	for(Index=0;Index<*NoPoints;Index++)
	{
		do
		{
			Map[Index] = gsl_rng_uniform_int(Opt->rng, Opt->NoQTL-1);
		} while(PointInMap(Map, Map[Index], Index) == TRUE);
	}

	qsort(Map, *NoPoints, sizeof(int), IntComp);

	Map[*NoPoints] = 0;

	return Map;
}

// Crate a map of recombination locations, using equally spaced points.
int*	CreateEqualCombMap(OPTIONS* Opt, int *NoPoints)
{
	int Index;
	int *Map;
	double ChunkSize;

	*NoPoints = (int)Opt->RecombPoissonMean;

	if(*NoPoints == 0 || *NoPoints >= Opt->NoQTL)
	{
		printf("err:\n");
		exit(1);
	}

	ChunkSize = (double)Opt->NoQTL / (Opt->RecombPoissonMean+1);

	Map = (int*)SMalloc(sizeof(int) * (*NoPoints+1));
	for(Index=0;Index<*NoPoints;Index++)
		Map[Index] = (Index+1) * (int)ChunkSize;

	qsort(Map, *NoPoints, sizeof(int), IntComp);

	Map[*NoPoints] = 0;

	return Map;
}

// Crate a gamete from a parent, using a number of recombination points, drawn from a Poisson
void RecombinationPoissonChr(OPTIONS* Opt, double *Chr, double *P1, double *P2)
{
	int *CombMap;
	int NoCombMap;
	double *CP;
	int Index, PointIndex;

	CombMap = CreateCombMap(Opt, &NoCombMap);
//	CombMap = CreateEqualCombMap(Opt, &NoCombMap);


	CP = P1;
	if(gsl_rng_uniform(Opt->rng) < 0.5)
		CP = P2;

	PointIndex = 0;
	for(Index=0;Index<Opt->NoQTL;Index++)
	{
		Chr[Index] = CP[Index];

		if(Index == CombMap[PointIndex])
		{
			if(CP == P1)
				CP = P2;
			else
				CP = P1;
		}
	}

	free(CombMap);
}

// Crate a gamete based on a number of recombination events
void RecombinationPoisson(OPTIONS* Opt, INDIVIDUAL *CIn, INDIVIDUAL *P1, INDIVIDUAL *P2)
{
	RecombinationPoissonChr(Opt, CIn->Genome[0], P1->Genome[0], P1->Genome[1]);
	RecombinationPoissonChr(Opt, CIn->Genome[1], P2->Genome[0], P2->Genome[1]);
}

// Crate a new induvial
void NewIn(OPTIONS* Opt, POP *CPop, INDIVIDUAL *CIn, double PhenotypeSD)
{
	int CNo;
	INDIVIDUAL *P1, *P2;

	if(Opt->Ploidy == 1)
	{
		P1 = CPop->SPop[gsl_rng_uniform_int(Opt->rng, CPop->NoSPop)];
		memcpy(CIn->Genome[0], P1->Genome[0], sizeof(double)*Opt->NoQTL);
		return;
	}

	do
	{
		P1 = CPop->SPop[gsl_rng_uniform_int(Opt->rng, CPop->NoSPop)];
		P2 = CPop->SPop[gsl_rng_uniform_int(Opt->rng, CPop->NoSPop)];

//		Opt->BreedingT++;
	} while(ValidDiploidParents(Opt, P1, P2, PhenotypeSD) == FALSE);


//	Recombination(Opt, CIn, P1, P2);
//	RecombinationUnEnven(Opt, CIn, P1, P2);
//	Recombination1Point(Opt, CIn, P1, P2);
	RecombinationPoisson(Opt, CIn, P1, P2);
	return;

	CNo = gsl_rng_uniform_int(Opt->rng, Opt->Ploidy);
	memcpy(CIn->Genome[0], P1->Genome[CNo], sizeof(double)*Opt->NoQTL);

	CNo = gsl_rng_uniform_int(Opt->rng, Opt->Ploidy);
	memcpy(CIn->Genome[1], P2->Genome[CNo], sizeof(double)*Opt->NoQTL);
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

// Crate a new population form the old.
void NewPop(OPTIONS* Opt, POP *CPop, POP *NPop)
{
	int IIndex;
	INDIVIDUAL *CIn;
	double PhenotypeSD;

//	Opt->BreedingT = 0;

	SetSurviving(Opt, CPop);

	if(DeadPopulation(Opt, CPop->NoSPop) == TRUE)
	{
		printf("Pop all died.\n");
		exit(0);
	}

	PhenotypeSD = -1;
	if(Opt->AssortativeMating != -1)
	{
		for(IIndex=0;IIndex<Opt->PopSize;IIndex++)
			CPop->TVect[IIndex] = CPop->Pop[IIndex]->Phenotype;

		PhenotypeSD = gsl_stats_sd(CPop->TVect, 1, Opt->PopSize);
	}

	for(IIndex=0;IIndex<Opt->PopSize;IIndex++)
	{
		CIn = NPop->Pop[IIndex];
		NewIn(Opt, CPop, CIn, PhenotypeSD);

		SetIndividualEnv(Opt, CIn);
	}
}

// Calculate the fitness of an individual
void Fitness(OPTIONS *Opt, INDIVIDUAL *Ind, double PhenotypeSD)
{
	double Phenotype;

	if(Opt->Drift == TRUE)
		Ind->Fitness = 1.0;
	else
	{
		Phenotype = Ind->Phenotype - Opt->FitnessPar[0];
		if(Opt->StandardisePhenotype == FALSE)
			Ind->Fitness = gsl_ran_gaussian_pdf(Phenotype, Opt->FitnessPar[1]) / (0.4 / Opt->FitnessPar[1]);
		else
		{
			Phenotype = Phenotype / PhenotypeSD;
			Ind->Fitness = gsl_ran_gaussian_pdf(Phenotype, Opt->FitnessPar[1]) / (0.4 / Opt->FitnessPar[1]);
		}
	}
}

// Move the fitness function mean, if appropriate.
void SetFitnessMean(OPTIONS *Opt,POP *Pop)
{
	double Mean, SD;
	double *TVect;
	int Index;

	TVect = Pop->TVect;

	for(Index=0;Index<Opt->PopSize;Index++)
		TVect[Index] = Pop->Pop[Index]->Phenotype;

	Mean = gsl_stats_mean(TVect, 1, Opt->PopSize);
	SD = gsl_stats_sd(TVect, 1, Opt->PopSize);

	if(Opt->MoveType == FITNESS)
		Opt->FitnessPar[0] = Mean + (Opt->FitnessPar[1] * Opt->MoveFitnessSD);

	if(Opt->MoveType == PHENOTYPE)
		Opt->FitnessPar[0] = Mean + (sqrt(Opt->PheotypeVar) * Opt->MoveFitnessSD);

	if(Opt->MoveType == CONTINUOUS)
		Opt->FitnessPar[0] = Opt->FitnessPar[0] + Opt->MoveFitnessSD;
}

// Calculate the fitness of the population
void PopFitness(OPTIONS *Opt, int Itter, POP *Pop)
{
	int Index;
	double PhenotypeSD;

	CaclPopPhenotype(Opt, Pop);

	PhenotypeSD = 1.0;
	if(Opt->StandardisePhenotype == TRUE)
		PhenotypeSD = gsl_stats_sd(Pop->TVect, 1, Opt->PopSize);

	if(Itter >= Opt->MoveFitnessSDGen && Opt->MoveFitnessSDGen != -1)
		SetFitnessMean(Opt, Pop);


	for(Index=0;Index<Opt->PopSize;Index++)
		Fitness(Opt, Pop->Pop[Index], PhenotypeSD);
}


// Get a vector of Phenotypes
void GetPhenotypeVect(OPTIONS* Opt,POP* Pop,double *Vect)
{
	int Index;

	for(Index=0;Index<Opt->PopSize;Index++)
		Vect[Index] = Pop->Pop[Index]->Phenotype;
}

// Get a vector of Genotypes
void GetGenotypeVect(OPTIONS* Opt,POP* Pop,double *Vect)
{
	int Index;

	for(Index=0;Index<Opt->PopSize;Index++)
		Vect[Index] = Pop->Pop[Index]->Genotype;
}

// Get a vector of fitnesses
void GetFitnessVect(OPTIONS* Opt,POP* Pop,double *Vect)
{
	int Index;

	for(Index=0;Index<Opt->PopSize;Index++)
		Vect[Index] = Pop->Pop[Index]->Fitness;
}

// Get a vector of enviroemntal varainces
void GetEnvVect(OPTIONS* Opt,POP* Pop,double *Vect)
{
	int Index;

	for(Index=0;Index<Opt->PopSize;Index++)
		Vect[Index] = Pop->Pop[Index]->Env;
}

// print the mean of a vector
void PrintVectMean(double* Vect, int Size)
{
	double Mean;
	Mean = gsl_stats_mean(Vect, 1, Size);
	printf("%f\t", Mean);
}

// print the variance of a vector
void PrintVectVar(double* Vect, int Size)
{
	double Var;
	Var = gsl_stats_variance(Vect, 1, Size);
	printf("%f\t", Var);
}

// print the mean and varaince of a vector
void PrintVectMeanVar(double* Vect, int Size)
{
	PrintVectMean(Vect, Size);
	PrintVectVar(Vect, Size);
}

// Print the regression of phenotype (x) on fitness (y), fitness if flipped if > mean of distribution
void OuputReg(OPTIONS* Opt, INDIVIDUAL **Pop, int Size)
{
	double *X, *Y;
	double R2, Slope, Int, PhenotypeFullSD;
	int Index;
	double a[4];


	X = Opt->X;
	Y = Opt->Y;

	for(Index=0;Index<Size;Index++)
		X[Index] = Pop[Index]->Phenotype;
	PhenotypeFullSD = gsl_stats_sd(X, 1, Size);

	for(Index=0;Index<Size;Index++)
	{
		X[Index] = Pop[Index]->Phenotype;
		Y[Index] = Pop[Index]->Fitness;

		if(X[Index] > Opt->FitnessPar[0])
			X[Index] = Opt->FitnessPar[0] - (X[Index] - Opt->FitnessPar[0]);

	}

	gsl_fit_linear(X, 1, Y, 1, Size, &Int, &Slope, &a[0], &a[1], &a[2], &a[3]);
	R2 = gsl_stats_correlation(X, 1, Y, 1, Size);
	R2 = R2 * R2;

    printf("%f\t%f\t", Slope, Slope * PhenotypeFullSD);
}

// Write info about current population to stdout
void	Output(OPTIONS* Opt, int Itter, POP* Pop)
{
	double PhenotypeVar, PhenotypeMean, GenotypeVar, EnvVar;
	double *TVect;

	TVect = Pop->TVect;

	printf("%d\t", Itter);

	GetPhenotypeVect(Opt, Pop, TVect);
	PrintVectMean(TVect, Opt->PopSize);
	PhenotypeMean = gsl_stats_mean(TVect, 1, Opt->PopSize);
	PhenotypeVar = gsl_stats_variance(TVect, 1, Opt->PopSize);
	printf("%f\t", PhenotypeMean / sqrt(PhenotypeVar));

	GetGenotypeVect(Opt, Pop, TVect);

	PrintVectVar(TVect, Opt->PopSize);
	GenotypeVar = gsl_stats_variance(TVect, 1, Opt->PopSize);

	GetFitnessVect(Opt, Pop, TVect);

	GetEnvVect(Opt, Pop, TVect);
	EnvVar = gsl_stats_variance(TVect, 1, Opt->PopSize);

	printf("%d\t", Pop->NoSPop);

	printf("%f\t%f\t", Opt->FitnessPar[0], Opt->FitnessPar[1]);


	OuputReg(Opt, Pop->SPop, Pop->NoSPop);

	printf("%f\t", sqrt((GenotypeVar*GenotypeVar) / Opt->MuteVar));


	printf("\n");

	fflush(stdout);
}

// Main loop of the simulation
void RunSim(OPTIONS* Opt)
{
	int Itter;
	POP *CPop, *NPop;

//	TestModel(Opt);


	CPop = CreatePop(Opt);
	NPop = CreatePop(Opt);

	printf("Generation\tPhenotype Mean\tPhenotype Ratio\tGenotype Var\tP Survived\tFitness Function Mean\tFitness Function SD\tSurviving Slope\tSurviving Standardised Slope\tw\t");

	printf("\n");

	Opt->Drift = FALSE;
	for(Itter=0;Itter<Opt->NoGenerations;Itter++)
	{
		PopFitness(Opt, Itter, CPop);

		NewPop(Opt, CPop, NPop);

		if(Itter % Opt->Sample == 0)
			Output(Opt, Itter, CPop);

		if(Itter == Opt->ChangeFitnessGen)
			Opt->FitnessPar[1] = Opt->ChangeFitnessSD;

		MutatePop(Opt, NPop);

		swap((void**)&CPop, (void**)&NPop);

	}

	FreePop(Opt, CPop);
	FreePop(Opt, NPop);
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

// Code to test a bivariate Gaussian, not used in the project.
void TestBivariateGaussian(OPTIONS *Opt)
{
	int Index;
	double p, x, y;
	gsl_rng* rng;

	p = gsl_ran_bivariate_gaussian_pdf(0.0, 0.0, 1.0, 1.0, 0.9);


	rng = Opt->rng;

	printf("Index\tX\tY\tP\n");
	for(Index=0;Index<10000;Index++)
	{
		x = (gsl_rng_uniform(Opt->rng) * 5) - 2.5;
		y = (gsl_rng_uniform(Opt->rng) * 4) - 2;

		p = gsl_ran_bivariate_gaussian_pdf(x,y, 1.0, 1.0, 0.9);

		printf("%d\t%f\t%f\t%f\n", Index, x, y, p);
	}

	exit(0);
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