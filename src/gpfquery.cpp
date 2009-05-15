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

#include "GpfIndexFile.h"
#include "GpfQuery.h"
#include "RefPtr.h"
#include "StopWatch.h"


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
	if (ai_ArgumentCount < 3)
	{
		printf("Usage: gpfquery [GPF index file] [peptide]\n");
		exit(1);
	}

	QString ls_IndexFilePath = QString(ac_Arguments__[1]);
	// load index file
	RefPtr<k_GpfIndexFile> lk_pGpfIndexFile(new k_GpfIndexFile(ls_IndexFilePath));
	if (!lk_pGpfIndexFile->isGood())
	{
		printf("Error: Unable to load %s.\n", ac_Arguments__[1]);
		exit(1);
	}
	
	QString ls_Peptide = QString(ac_Arguments__[2]);

	{
		k_StopWatch lk_StopWatch("Processing the query took %1.\n");
		k_GpfQuery lk_Query(*(lk_pGpfIndexFile.get_Pointer()), ls_Peptide);
		lk_Query.execute();
	}
	
	/*
	for (qint64 i = 0; i < lk_pGpfIndexFile->mi_TagCount * 2; ++i)
	{
		qint64 li_HmstCount = lk_pGpfIndexFile->mk_HmstCount[i];
		qint64 li_HmstOffset = lk_pGpfIndexFile->mk_HmstOffset[i];
		qint64 li_MassBits = li_HmstCount * lk_pGpfIndexFile->mi_MassBits;
		qint64 li_GnoBits = li_HmstCount * lk_pGpfIndexFile->mi_OffsetBits;
		qint64 li_OldMass = 0;
		for (qint64 k = 0; k < li_HmstCount; ++k)
		{
			qint64 li_Mass = lk_pGpfIndexFile->readIndexBits((li_HmstOffset + k) * lk_pGpfIndexFile->mi_MassBits, lk_pGpfIndexFile->mi_MassBits);
			if (li_Mass < li_OldMass)
				printf("ERROR!\n");
		}
	}
	*/
}
