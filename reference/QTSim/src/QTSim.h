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

#ifndef GSIM_H
#define GSIM_H

#include "GenLib.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fit.h>

#define	VERSION 0.1

#define GENETIC_TARGET 0.5
#define NO_ITTERS_PER_TEST 10000

#define MUT_RATE_PER_QTL 0.0002 // (2*10^-4)

#define RECOMB_BIAS 0.9

#define RECOMB_POISSON_DEFAULT 1.0


//#define NO_ENV_VAR // Environmental variance 

typedef enum
{
	FITNESS, 
	PHENOTYPE, 
	CONTINUOUS
} OPTIMUM_MOVE_TYPE;

typedef struct
{
	double **Genome;
	double Env;
	double Fitness;
	double Phenotype;
	double Genotype;
} INDIVIDUAL;

typedef struct
{
	INDIVIDUAL** Pop;

	INDIVIDUAL** SPop;
	int NoSPop;

	double *TVect;
} POP;

typedef struct
{
	int NoGenerations;
	int Sample;
	int PopSize;
	int NoQTL;
	int Ploidy;

	double PheotypeVar;
	double Heritabitlity;

	double EnvVar;
	double MuteVar;
	double MuteVarScalar;

	double *FitnessPar;

//	int OptVar;
	double OptVarTarget;

	int ChangeFitnessGen;
	double ChangeFitnessSD;

	int MoveFitnessSDGen;
	double MoveFitnessSD;

	gsl_rng* rng;

	int Drift;

	double *X;
	double *Y;

	unsigned long Seed;

	int StandardisePhenotype;

	OPTIMUM_MOVE_TYPE MoveType;

	double AssortativeMating;

	double InitFitnessSD;

	double MutationRatePerGamete;

	double RecombPoissonMean;

//	int BreedingT;
} OPTIONS;

typedef struct
{
	OPTIONS *Opt;
	POP *CPop;
	POP *NPop;

	int Itters;
	double Taget;

} OPT_STRUCT;


double  CaclGenotype(OPTIONS* Opt, INDIVIDUAL* Ind);
POP* CreatePop(OPTIONS* Opt);
void PopFitness(OPTIONS *Opt, int Itter, POP *Pop);
void NewPop(OPTIONS* Opt, POP *CPop, POP *NPop);
void MutatePop(OPTIONS* Opt, POP* Pop);




#endif
