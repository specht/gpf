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
#include "StopWatch.h"


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
 	if (ai_ArgumentCount < 4)
	{
		fprintf(stderr, "Usage:   gpfbatch [options] [GPF index file] [query file] [result file]\n");
		exit(1);
	}

	QStringList lk_Arguments;
	for (int i = 1; i < ai_ArgumentCount; ++i)
		lk_Arguments << QString(ac_Arguments__[i]);
		
	QMap<QString, QString> lk_Parameters;
	while (lk_Arguments.size() > 3)
	{
		QString ls_Key = lk_Arguments.takeFirst();
		QString ls_Value = lk_Arguments.takeFirst();
		lk_Parameters[ls_Key] = ls_Value;
	}
	
	QString ls_IndexFile = lk_Arguments.takeFirst();
	QString ls_QueryFile = lk_Arguments.takeFirst();
	QString ls_ResultFile = lk_Arguments.takeFirst();

	k_GpfBase lk_GpfBase(QStringList() << ls_IndexFile);
	
	QFile lk_File(ls_QueryFile);
	lk_File.open(QIODevice::ReadOnly);
	unsigned int lui_QueryCount = 0;
	while (!lk_File.atEnd())
	{
		lk_File.readLine();
		++lui_QueryCount;
	}
	lk_File.reset();
	unsigned int lui_ProcessedCount = 0;
	QFile lk_OutFile(ls_ResultFile);
	lk_OutFile.open(QIODevice::WriteOnly);
	k_StopWatch lk_StopWatch(QString("Processing took %1.\n"));
	while (!lk_File.atEnd())
	{
		QString ls_Line = QString(lk_File.readLine()).trimmed();
		QStringList lk_Line = ls_Line.split(";");
		QString ls_Peptide = lk_Line.takeFirst();
		ls_Peptide.replace("peptide=", "");
		double ld_PrecursorMass = 0.0;
		if (!lk_Line.empty())
		{
			QString ls_PrecursorMass = lk_Line.takeFirst();
			ls_PrecursorMass.replace("precursorMass=", "");
			bool lb_Ok = false;
			ld_PrecursorMass = QVariant(ls_PrecursorMass).toDouble(&lb_Ok);
			if (!lb_Ok)
				ld_PrecursorMass = 0.0;
		}
		k_GpfQuery lk_Query(lk_GpfBase, ls_Peptide, ld_PrecursorMass);
		lk_Query.SetParameters(lk_Parameters);
		lk_Query.Execute();
		++lui_ProcessedCount;
		printf("\rProcessed %d of %d queries.", lui_ProcessedCount, lui_QueryCount);
		QString ls_Result = lk_Query.get_Result();
		lk_OutFile.write(QString("\"%1\":\n").arg(ls_Line).toAscii());
		if (!ls_Result.isEmpty())
		{
			ls_Result = "  " + ls_Result.replace("\n-", "\n  -");
			lk_OutFile.write(ls_Result.toAscii());
		}
	}
	printf("\nGPF batch finished.\n");
	lk_OutFile.close();
	lk_File.close();
}
