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
	if (ai_ArgumentCount < 3)
	{
		printf("Usage: gpfindex [options] [Genomic DNA sequence file] [GPF index out file]\n");
        printf("The input file must be in FASTA format.\n\n");
        printf("Options:\n");
        printf("  --title <string>\n");
        printf("    Genome title, is input filename by default.\n");
        printf("  --tagSize <int> (default: 5)\n");
        printf("    Tag size.\n");
        printf("  --enzyme <string> (default: 'RK|')\n");
        printf("    Specify the enzyme for creating peptides. The | symbol denotes the\n");
        printf("    cleavage site. RK| means that the enzyme cleaves after R or K.\n");
        printf("    Inhibitory amino acids do not play a role, because missed cleavages\n");
        printf("    are handled by GPF.\n");
        printf("  --alloc <amount> (default: 512M)\n"); 
        printf("    Amount of RAM to allocate for the index buffer. A bigger index\n");
        printf("    buffer leads to faster execution of this program. The amount may\n");
        printf("    be specified as a number, optionally followed by any of the\n");
        printf("    prefixes K, M, G.\n");
        printf("    Please note that in addition to the index buffer, the whole DNA\n");
        printf("    is loaded into RAM, which uses 3/8 of the DNA input file's size.\n");
        printf("  --massPrecision <int> (default: 10000)\n");
        printf("    Specify the precision of masses here. A mass precision of 10000\n");
        printf("    corresponds to four decimal places.\n");
        printf("  --massBits <int> (default: 27)\n");
        printf("    Specify how many bits should be used to encode masses. Together\n");
        printf("    with the mass precision specified, this defines the highest mass\n");
        printf("    that can be encoded in the index file.\n");
        
		exit(1);
	}
	
    QStringList lk_Arguments;
    for (int i = 1; i < ai_ArgumentCount; ++i)
        lk_Arguments << ac_Arguments__[i];

    QString ls_IndexFilename = lk_Arguments.takeLast();
    QString ls_GenomeFilename = lk_Arguments.takeLast();
    
    QString ls_Title = QFileInfo(ls_GenomeFilename).baseName();
    
    qint32 li_TagSize = 5;
    QString ls_Enzyme = "RK|";
    qint64 li_IndexBufferAllocSize = 512 * 1024 * 1024;
    qint32 li_MassPrecision = 10000;
    qint32 li_MassBits = 27;
    
    while (!lk_Arguments.empty())
    {
        QString ls_Key = lk_Arguments.takeFirst();
        if (ls_Key == "--title")
        {
            QString ls_Title = lk_Arguments.takeFirst();
        }
        else if (ls_Key == "--tagSize")
        {
            bool lb_Ok = false;
            QString ls_Value = lk_Arguments.takeFirst();
            li_TagSize = ls_Value.toInt(&lb_Ok);
            if (!lb_Ok)
            {
                printf("Error: Invalid tag size specified: %s.\n", ls_Value.toStdString().c_str());
                exit(1);
            }
        }
        else if (ls_Key == "--enzyme")
        {
            ls_Enzyme = lk_Arguments.takeFirst().toUpper().trimmed();
        }
        else if (ls_Key == "--alloc")
        {
            QString ls_Value = lk_Arguments.takeFirst();
            ls_Value = ls_Value.trimmed();
            QString ls_Suffix = ls_Value.right(1).toUpper();
            qint64 li_Factor = 1;
            if (ls_Suffix == "K" || ls_Suffix == "M" || ls_Suffix == "G")
            {
                ls_Value = ls_Value.left(ls_Value.length() - 1);
                if (ls_Suffix == "K")
                    li_Factor = 1024;
                else if (ls_Suffix == "M")
                    li_Factor = 1024 * 1024;
                else if (ls_Suffix == "G")
                    li_Factor = 1024 * 1024 * 1024;
            }
            bool lb_Ok = false;
            li_IndexBufferAllocSize = ls_Value.toInt(&lb_Ok);
            if (!lb_Ok)
            {
                printf("Error: Invalid allocation size specified: %s.\n", ls_Value.toStdString().c_str());
                exit(1);
            }
            li_IndexBufferAllocSize *= li_Factor;
        }
        else if (ls_Key == "--massPrecision")
        {
            bool lb_Ok = false;
            QString ls_Value = lk_Arguments.takeFirst();
            li_MassPrecision = ls_Value.toInt(&lb_Ok);
            if (!lb_Ok)
            {
                printf("Error: Invalid mass precision specified: %s.\n", ls_Value.toStdString().c_str());
                exit(1);
            }
        }
        else if (ls_Key == "--massBits")
        {
            bool lb_Ok = false;
            QString ls_Value = lk_Arguments.takeFirst();
            li_MassBits = ls_Value.toInt(&lb_Ok);
            if (!lb_Ok)
            {
                printf("Error: Invalid mass bits specified: %s.\n", ls_Value.toStdString().c_str());
                exit(1);
            }
        }
        else
        {
            printf("Error: Invalid parameter specified: %s.\n", ls_Key.toStdString().c_str());
            exit(1);
        }
    }
    
	k_GpfIndexer lk_GpfIndexer(ls_GenomeFilename, ls_IndexFilename, ls_Title,
                               li_TagSize, ls_Enzyme, li_IndexBufferAllocSize,
                               li_MassPrecision, li_MassBits);
	lk_GpfIndexer.compileIndex();
/*	char* s = "VIHAR";
	printf("%d\n", gk_GpfBase.aminoAcidPolymerCode(s, 5));*/
}
