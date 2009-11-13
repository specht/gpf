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
#include "Hit.h"
#include "GpfQuery.h"
#include "StopWatch.h"


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
 	if (ai_ArgumentCount < 4)
	{
		fprintf(stderr, "Usage:   gpfbatch [options] [GPF index file] [query file]\n");
        fprintf(stderr, "Options:  --yamlResultsPath [path]\n");
        fprintf(stderr, "          --csvResultsPath [path]\n");
		exit(1);
	}

	QStringList lk_Arguments;
	for (int i = 1; i < ai_ArgumentCount; ++i)
		lk_Arguments << QString(ac_Arguments__[i]);

    RefPtr<QTextStream> lk_pYamlStream;
    RefPtr<QTextStream> lk_pCsvStream;
    RefPtr<QFile> lk_pYamlFile;
    RefPtr<QFile> lk_pCsvFile;
    
    if (lk_Arguments.contains("--yamlResultsPath"))
    {
        int li_Index = lk_Arguments.indexOf("--yamlResultsPath");
        QString ls_Path = lk_Arguments[li_Index + 1];
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
        lk_pYamlFile = RefPtr<QFile>(new QFile(ls_Path));
        lk_pYamlFile->open(QIODevice::WriteOnly);
        lk_pYamlStream = RefPtr<QTextStream>(new QTextStream(lk_pYamlFile.get_Pointer()));
    }
		
    if (lk_Arguments.contains("--csvResultsPath"))
    {
        int li_Index = lk_Arguments.indexOf("--csvResultsPath");
        QString ls_Path = lk_Arguments[li_Index + 1];
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
        lk_pCsvFile = RefPtr<QFile>(new QFile(ls_Path));
        lk_pCsvFile->open(QIODevice::WriteOnly);
        lk_pCsvStream = RefPtr<QTextStream>(new QTextStream(lk_pCsvFile.get_Pointer()));
    }
        
	QMap<QString, QString> lk_Parameters;
	while (lk_Arguments.size() > 2)
	{
		QString ls_Key = lk_Arguments.takeFirst();
		QString ls_Value = lk_Arguments.takeFirst();
		lk_Parameters[ls_Key] = ls_Value;
	}
	
	QString ls_IndexFile = lk_Arguments.takeFirst();
	QString ls_QueryFile = lk_Arguments.takeFirst();

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
	k_StopWatch lk_StopWatch(QString("Processing took %1.\n"));
    if (lk_pCsvStream)
    {
        *(lk_pCsvStream.get_Pointer()) << "Query,Peptide,Assembly,Left,Right,Mass,Intron ends,Score,Part scores\n";
    }
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
        if (lk_pYamlStream)
        {
            *(lk_pYamlStream.get_Pointer()) << QString("\"%1\":\n").arg(ls_Line);
            if (!ls_Result.isEmpty())
            {
                ls_Result = "  " + ls_Result.replace("\n-", "\n  -");
                *(lk_pYamlStream.get_Pointer()) << ls_Result;
            }
        }
        if (lk_pCsvStream)
        {
            foreach (QString ls_Key, lk_Query.resultList().uniqueKeys())
            {
                k_Hit* lk_Hit_ = lk_Query.resultList()[ls_Key].get_Pointer();
                tk_StringHash lk_Hash = lk_Hit_->descriptionHash();
                *(lk_pCsvStream.get_Pointer()) << 
                    ls_Peptide << "," <<
                    lk_Hash["peptide"] << "," <<
                    "\"" << lk_Hash["assembly"] << "\"," <<
                    lk_Hash["left"] << "," <<
                    lk_Hash["right"] << "," <<
                    lk_Hash["mass"] << "," <<
                    lk_Hash["intronEnds"] << "," <<
                    lk_Hash["score"] << "," <<
                    "\"" << lk_Hash["partScores"] << "\"" <<
                    "\n";
                //ms_Result += QString("- %1\n").arg(lk_pHit->description());
            }
        }
	}
	printf("\nGPF batch finished.\n");
	lk_File.close();
    
    if (lk_pYamlStream)
        lk_pYamlFile->flush();
    if (lk_pCsvStream)
        lk_pCsvStream->flush();
}
