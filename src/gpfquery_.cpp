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

#include <stdio.h>
#include <QtCore>
#include "GpfBase.h"
#include "GpfQuery.h"


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
 	if (ai_ArgumentCount < 3)
	{
		fprintf(stderr, "Usage:   gpfquery [GPF index file] [query] [precursor mass (optional)]\n");
		exit(1);
	}

	QString ls_IndexFile = ac_Arguments__[1];
	QString ls_Peptide = ac_Arguments__[2];

	double ld_PrecursorMass = 0.0;
	if (ai_ArgumentCount > 3)
	{
		bool lb_Ok;
		ld_PrecursorMass = QString(ac_Arguments__[3]).toDouble(&lb_Ok);
		if (!lb_Ok)
			ld_PrecursorMass = 0.0;
	}

	k_GpfBase lk_GpfBase(QStringList() << ls_IndexFile);
	k_GpfQuery lk_Query(lk_GpfBase, ls_Peptide, ld_PrecursorMass);
	lk_Query.Execute();
	printf(lk_Query.get_Result().toStdString().c_str());
}
