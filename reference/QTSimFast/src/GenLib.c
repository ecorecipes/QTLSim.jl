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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <stdint.h>

#include "GenLib.h"

void		MallocErrFull(char* FileName, int LineNo)
{
	fprintf(stderr, "Memory allocation error in file %s line %d\n", FileName, LineNo);
	fprintf(stdout, "Memory allocation error in file %s line %d\n", FileName, LineNo);

	exit(1);
}

void*	smalloc(size_t n, char* FName, unsigned long LineNo)
{
	void *Ret;

	Ret = malloc(n);

	if(Ret == NULL)
		MallocErrFull(FName, LineNo);


	return Ret;
}

char*		StrMake(char* Str)
{
	char*	Ret;

	Ret = (char*)malloc(sizeof(char) * (strlen(Str) + 1));
	if(Ret == NULL)
		MallocErr();

	strcpy(Ret, Str);

	return Ret;
}

FILE*		OpenRead(char *FileName)
{
	FILE*	Ret;

	Ret = fopen(FileName, "r");
	if(Ret == NULL)
	{
		fprintf(stderr, "Could not open file %s for reading\n", FileName);
		exit(1);
	}

	return Ret;
}

FILE*	OpenAppend(char *FileName)
{
	FILE* Ret;

	Ret = fopen(FileName, "a");

	if(Ret == NULL)
	{
		printf("Could not open file \"%s\" for appending %d::%s\n", FileName, __LINE__, __FILE__);
		exit(1);
	}

	return Ret;
}



FILE*		OpenWrite(char *FileName)
{
	FILE*	Ret;

	Ret = fopen(FileName, "w");
	if(Ret == NULL)
	{
		fprintf(stderr, "Could not open file %s for writting\n", FileName);
		exit(1);
	}

	return Ret;
}

char*		AppendStr(char *Base, char *Append)
{
	char *Buffer;

	Buffer = (char*)SMalloc((strlen(Base) + strlen(Append) + 1) * sizeof(char));

	sprintf(Buffer, "%s%s", Base, Append);

	return Buffer;
}

FILE*		OpenWriteWithExt(char *Base, char *Ext)
{
	char *Buffer;
	FILE *Ret;

	Buffer = AppendStr(Base, Ext);

	Ret = OpenWrite(Buffer);

	free(Buffer);

	return Ret;
}

char	IsNewLine(char *Char)
{
	char	C1;
	char	C2;

	if(*Char == '\0')
		return 'E';


	C1 = *Char;
	Char++;
	C2 = *Char;

	if((C1 == 13) && (C2 == 10))
		return 'P';

	if(C1 == 13)
		return 'M';

	if(C1 == 10)
		return 'U';

	return 'N';
}

int		FindNoOfNL(char* FileBuffer, int FileSize)
{
	int		BIndex;
	char	Line;
	int		Ret;

	Ret=0;
	for(BIndex=0;BIndex<FileSize;BIndex++)
	{
		Line = IsNewLine(&FileBuffer[BIndex]);

		if(Line != 'N')
			Ret++;

		if(Line == 'P')
			BIndex++;
	}

	Ret++;
	return Ret;
}

int		FindLineLen(char *Start)
{
	int	Ret;

	Ret =0;
	while(IsNewLine(Start) == 'N')
	{
		Start++;
		Ret++;
	}

	return Ret;
}

char**	ProcessBinaryFile(char* FileBuffer, int FileSize, int* NoOfLines)
{
	char	*P;
	int		LineIndex;
	int		LineLen;
	char	NL;
	char	**Ret;

	(*NoOfLines) = FindNoOfNL(FileBuffer, FileSize);

	Ret = (char**)malloc(sizeof(char*) * (*NoOfLines));
	if(Ret  == NULL)
		MallocErr();

	P = FileBuffer;
	for(LineIndex=0;LineIndex<*NoOfLines;LineIndex++)
	{
		LineLen = FindLineLen(P);
		Ret[LineIndex] = (char*)malloc(sizeof(char) * (LineLen + 1));
		if(Ret[LineIndex] == NULL)
			MallocErr();

		memcpy(Ret[LineIndex], P, sizeof(char) * LineLen);
		Ret[LineIndex][LineLen] = '\0';

		P += LineLen;

		NL = IsNewLine(P);
		if((NL == 'U') || (NL == 'M'))
			P++;

		if(NL == 'P')
			P+=2;
	}

	return Ret;
}

char*	RemoveComment(char* FileBuffer, int* FileSize, char *FileName)
{
	char*	Block;
	char*	Ret;
	int		BIndex;
	int		DIndex;
	int		InCom;

	Block = (char*)malloc(sizeof(char) * (*FileSize)+1);
	if(Block == NULL)
		MallocErr();

	InCom = FALSE;
	BIndex=0;
	for(DIndex=0;DIndex<(*FileSize)+1;DIndex++)
	{
		if(FileBuffer[DIndex] == '[')
		{
			if(InCom == TRUE)
			{
				printf("File %s has nested comments\n", FileName);
				exit(0);
			}
			InCom = TRUE;
		}

		if(InCom != TRUE)
		{
			Block[BIndex] = FileBuffer[DIndex];
			BIndex++;
		}

		if(FileBuffer[DIndex] == ']')
		{
			if(InCom == FALSE)
			{
				printf("Found close comment with out opening in %s\n", FileName);
				exit(0);
			}
			InCom = FALSE;
		}
	}

	BIndex--;
	*FileSize = BIndex;

	free(FileBuffer);

	Ret = (char*)malloc(sizeof(char) * (*FileSize) + 1);
	if(Ret == NULL)
		MallocErr();

	memcpy(Ret, Block, sizeof(char) * (*FileSize) + 1);
	free(Block);

	return Ret;

}

void		FindTextFileMaxLine(TEXTFILE* TextFile)
{
	int	Index;
	int	Len;

	if(TextFile->NoOfLines == 0)
	{
		TextFile->MaxLine = 0;
		return;
	}

	TextFile->MaxLine = (int)strlen(TextFile->Data[0]);
	for(Index=1;Index<TextFile->NoOfLines;Index++)
	{
		Len = (int)strlen(TextFile->Data[Index]);
		if(Len > TextFile->MaxLine)
			TextFile->MaxLine = Len;
	}
	TextFile->MaxLine++;
}

TEXTFILE*	LoadTextFile(char* Name, char DelComments)
{
	TEXTFILE*	Ret;
	FILE*		InFile;
	char*		Buffer;
	int			NoRead;


	Ret = (TEXTFILE*)malloc(sizeof(TEXTFILE));
	if(Ret == NULL)
		MallocErr();

	Ret->FileName	= StrMake(Name);
	Ret->Size		= 0;
	Ret->NoOfLines	= 0;
	Ret->MaxLine	= 0;
	Ret->Data		= NULL;

	InFile = fopen(Name, "rb");
	if(InFile == NULL)
	{
		fprintf(stdout, "Could not open file %s for reading\n", Name);
		fprintf(stderr, "Could not open file %s for reading\n", Name);
		exit(1);
	}

	Buffer = (char*)malloc(sizeof(char) * BUFFERSIZE);
	if(Buffer == NULL)
		MallocErr();

	do
	{
		NoRead = (int)fread(&Buffer[0], sizeof(char), BUFFERSIZE, InFile);
		Ret->Size += NoRead;
	} while(NoRead == BUFFERSIZE);

	free(Buffer);

	Buffer = (char*)malloc(sizeof(char) * (Ret->Size + 1));
	if(Buffer == NULL)
		MallocErr();

	rewind(InFile);
	fread(&Buffer[0], sizeof(char), Ret->Size, InFile);
	fclose(InFile);
	Buffer[Ret->Size] = '\0';

	if(DelComments == TRUE)
		Buffer = RemoveComment(Buffer, &Ret->Size, Name);

	Ret->Data = ProcessBinaryFile(Buffer, Ret->Size, &Ret->NoOfLines);

	FindTextFileMaxLine(Ret);

	free(Buffer);
	return Ret;
}

void	FreeTextFile(TEXTFILE* TextFile)
{
	int	Index;

	free(TextFile->FileName);

	for(Index=0;Index<TextFile->NoOfLines;Index++)
		free(TextFile->Data[Index]);
	free(TextFile->Data);

	free(TextFile);
}

void swap(void** a, void** b)
{
	void* t;

	t = *a;
	*a = *b;
	*b = t;
}

void	ReplaceChar(char Rep, char With, char* String)
{
	int	Index;
	int	Size;

	Size = (int)strlen(String);

	for(Index=0;Index<Size;Index++)
	{
		if(String[Index] == Rep)
			String[Index] = With;
	}

}


int		IsValidDouble(char* Str)
{
	int	Point=FALSE;

	if(atof(Str) != 0.0)
		return TRUE;
	while(isspace(*Str)) Str++;

	if(*Str == '-')
		Str++;

	while(*Str)
	{
		if(*Str == '.')
		{
			if(Point == FALSE)
				Point = TRUE;
			else
				return FALSE;
		}
		else
		{
			if(!((*Str >= '0') && (*Str <= '9')))
				return FALSE;
		}

		Str++;
	}

	return TRUE;
}

int		IsValidInt(char* Str)
{
	if(*Str == '-')
		Str++;

	while(*Str)
	{
		if(!((*Str >= '0') && (*Str <= '9')))
			return FALSE;
		Str++;
	}

	return TRUE;
}

int	MakeArgv(char*	string, char *argv[], int argvsize)
{
	char*	p = string;
	int		i;
	int		argc = 0;

	for(i=0; i<argvsize; i++)
	{
		while(isspace(*p)) p++;

		if(*p != '\0')
			argv[argc++] = p;
		else
		{
			argv[argc] = 0;
			break;
		}

		while(*p != '\0' && !isspace(*p))
			p++;

		if(*p != '\0' && i < argvsize - 1)
			*p++ = '\0';

	}

	return argc;
}

int strcicmp(char const *a, char const *b)
{
    for (;; a++, b++) {
        int d = tolower(*a) - tolower(*b);
        if (d != 0 || !*a)
            return d;
    }
}
