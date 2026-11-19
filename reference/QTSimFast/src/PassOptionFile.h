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

#ifndef PASS_OPT_FILE_H
#define PASS_OPT_FILE_H

#include "QTSim.h"
#pragma warning(disable:6385)


typedef void (*optfun)(OPTIONS* , int, char**);


void Help(OPTIONS* Opt, int argc, char** argv);
void SetSampleNo(OPTIONS* Opt, int argc, char** argv);
void SetPopSize(OPTIONS* Opt, int argc, char** argv);
void SetNoQTL(OPTIONS* Opt, int argc, char** argv);
void SetPloidy(OPTIONS* Opt, int argc, char** argv);
void SetDrift(OPTIONS* Opt, int argc, char** argv);
void SetEnvVar(OPTIONS* Opt, int argc, char** argv);
void SetMutScalar(OPTIONS* Opt, int argc, char** argv);
void SetGenerations(OPTIONS* Opt, int argc, char** argv);
void SetChangeFitnessSD(OPTIONS* Opt, int argc, char** argv);
void SetMoveFitnessSD(OPTIONS* Opt, int argc, char** argv);
void SetSeed(OPTIONS* Opt, int argc, char** argv);


void SetPheotypeVar(OPTIONS* Opt, int argc, char** argv);
void SetHeritabitlity(OPTIONS* Opt, int argc, char** argv);
void SetStandardisePhenotype(OPTIONS* Opt, int argc, char** argv);
void SetOptimumMoveType(OPTIONS* Opt, int argc, char** argv);
void SetAssortativeMating(OPTIONS* Opt, int argc, char** argv);
void SetInitFitnessSD(OPTIONS* Opt, int argc, char** argv);
void SetMutationRatePerGamete(OPTIONS* Opt, int argc, char** argv);
void SetRecombinationPoissonMean(OPTIONS* Opt, int argc, char** argv);


typedef struct
{
	char *Cmd;
    optfun Fun;
    char *Help;
 } CMD;

static CMD CMD_LIST[] = {
    {"Help", Help, "Print a list of command and their parameters."},
    {"Sample", SetSampleNo, NULL},
    {"PopSize", SetPopSize, NULL},
    {"NoQTL", SetNoQTL, NULL},
    {"Ploidy", SetPloidy, NULL},
    {"Drift", SetDrift, NULL},
//    {"EnvVar", SetEnvVar},
    {"MutationScalar", SetMutScalar, NULL},
    {"Generations", SetGenerations, NULL},
    {"ChangeFitnessSD", SetChangeFitnessSD, NULL},
    {"MoveFitnessSD", SetMoveFitnessSD, NULL},
    {"Seed", SetSeed, NULL},
    {"PhenotypeVar", SetPheotypeVar, NULL},
    {"Heritability", SetHeritabitlity, NULL},
    {"StandardisePhenotype", SetStandardisePhenotype, NULL},
    {"OptimumMoveType", SetOptimumMoveType, NULL},
//    {"AssortativeMating", SetAssortativeMating, NULL},
    {"InitFitnessSD", SetInitFitnessSD, NULL},
    {"MutationRatePerGamete", SetMutationRatePerGamete, NULL},
	{"RecombinationPoissonMean", SetRecombinationPoissonMean, NULL},
    {NULL, NULL, NULL}
};

void PassOptFile(OPTIONS *Opt, char *FName);

#endif

