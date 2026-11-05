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

// General library taken from BayesTraits. 

#ifndef GENLIB
#define GENLIB

#include <stdio.h>

#pragma warning(disable : 4996)

#ifndef TRUE
	#define TRUE 1
#endif

#ifndef FALSE
	#define	FALSE 0
#endif

#define	BUFFERSIZE	655360
#define	MallocErr() MallocErrFull(__FILE__, __LINE__)

void*	smalloc(size_t n, char* FName, unsigned long LineNo);

#define SMalloc(N) smalloc(N, __FILE__, __LINE__)


typedef struct
{
	char	*FileName;
	int		NoOfLines;
	int		MaxLine;
	int		Size;
	char	**Data;
} TEXTFILE;

TEXTFILE*	LoadTextFile(char* Name, char DelComments);
void		FreeTextFile(TEXTFILE* TextFile);
/* Memory */
void		MallocErrFull(char* FileName, int LineNo);
void		swap(void** a, void** b);


/*	Input / Output */
FILE*		OpenWrite(char *FileName);
FILE*		OpenWriteWithExt(char *Base, char *Ext);
FILE*		OpenRead(char *FileName);
FILE*		OpenAppend(char *FileName);

/* Text */
char*		StrMake(char* Str);
int			IsValidInt(char* Str);
int			IsValidDouble(char* Str);
int			MakeArgv(char*	string, char *argv[], int argvsize);
int			strcicmp(char const *a, char const *b);

char*		AppendStr(char *Base, char *Append);

void		ReplaceChar(char Rep, char With, char* String);


#endif