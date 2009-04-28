/*
Copyright (c) 2007-2008 Michael Specht

This file is part of GPF.

GPF is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GPF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GPF.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "GpfIndexer.h"


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
	if (ai_ArgumentCount < 4)
	{
		printf("Usage: gpfindex [genome fasta file] [GPF index out file] [genome title]\n");
		exit(1);
	}

	// construct index filename
	QString ls_GenomeFilename = ac_Arguments__[1];
	QString ls_IndexFilename = ac_Arguments__[2];
	QString ls_GenomeTitle = ac_Arguments__[3];

	k_GpfIndexer lk_GpfIndexer(ls_GenomeFilename, ls_IndexFilename, ls_GenomeTitle);
	lk_GpfIndexer.compileIndex();
}
