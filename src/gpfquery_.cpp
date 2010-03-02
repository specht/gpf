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


int parseValue(QString as_Value, QStringList ak_Choices, QString as_Description)
{
    for (int i = 0; i < ak_Choices.size(); ++i)
    {
        if (as_Value == ak_Choices[i])
            return i;
    }
    printf("Error: Invalid %s specified: %s.\n", 
           as_Description.toStdString().c_str(),
           as_Value.toStdString().c_str());
    exit(1);
}


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
    if (ai_ArgumentCount < 3)
    {
        printf("Usage: gpfquery [options] [GPF index file] [peptide 1] [precursor mass 1 (optional)] [peptide 2] ...\n");
        printf("\n");
        printf("For every query peptide, an uncharged precursor mass can be specified, this is useful\n");
        printf("if the peptide sequence comes from de novo prediction. If no precursor mass is\n");
        printf("specified, an appropriate precursor mass is calculated from the peptide sequence.\n");
        printf("\n");
        printf("Options:\n");
        printf("  --massAccuracy <float> (default: 10.0)\n");
        printf("    Specify the mass accuracy in ppm.\n");
        printf("  --similaritySearch <yes|no> (default: yes)\n");
        printf("    Specify whether alignments should be found that are exact in mass,\n");
        printf("    but not necessarily exact in amino acid sequence.\n");
        printf("  --distinguishIL <yes|no> (default: no)\n");
        printf("    Specify whether isoleucine (I) and leucine (L) should be\n");
        printf("    considered the same when checking for sequence equality.\n");
        printf("  --searchImmediate <yes|no> (default: yes)\n");
        printf("    Specify whether immediate (i. e. non-intron split) alignments\n");
        printf("    should be searched for.\n");
        printf("  --searchIntronSplit <yes|no|conditional> (default: conditional)\n");
        printf("    Specify whether non-intron split alignments should be searched\n");
        printf("    for. Conditional means that intron split alignments for a given\n");
        printf("    peptide are only searched for if no immediate alignments have\n");
        printf("    been found.\n");
        printf("  --intronSearchType <exhaustive|quick> (default: exhaustive)\n");
        printf("    Specify whether the search for intron split alignments should be\n");
        printf("    carried out exhaustively (which finds more alignments but takes longer)\n");
        printf("    or quickly (which finds less alignments).\n");
        printf("  --maxIntronLength <int> (default: 2100)\n");
        printf("    Specify the maximum intron length in nucleotides.\n");
        printf("  --intronSpliceSites <string> (default: 'GT|AG,GC|AG')\n");
        printf("    Specify possible splice donor/acceptor site consensus sequences.\n");
        printf("  --tagSize <int>\n");
        printf("    Apart from the tag size the genomic DNA sequence was indexed with,\n");
        printf("    you can specify a higher tag size here to achieve the same effect\n");
        printf("    as if the index file had been created with this tag size, although\n");
        printf("    this is a bit slower compared to using the right index file in the\n");
        printf("    first place.\n");
        // :TODO: 
        // 1. implement this 
        // 2. find out whether this is slower (it should be, right?)
        printf("  --peptidesFile <path>\n");
        printf("    Specify a text file containing the query peptides, one peptide per line.\n");
        printf("    Optionally, each peptide may be followed by the precursor mass, separated\n");
        printf("    by a comma, semicolon, or whitespace.\n");
        printf("  --csvOutputPath <path>\n");
        printf("    Specify a target file for CSV output. By default, the CSV output\n");
        printf("    goes to stdout. Caution: If the file exists, it will be overwritten.\n");
        printf("  --quiet\n");
        printf("    Don't print status messages.\n");
        exit(1);
    }
    
    QStringList lk_Arguments;
    
    for (int i = 1; i < ai_ArgumentCount; ++i)
        lk_Arguments << ac_Arguments__[i];
    
    double ld_MassAccuracy = 10.0;
    bool lb_SimilaritySearch = true;
    bool lb_DistinguishIL = false;
    bool lb_SearchImmediateAlignments = true;
    r_SearchIntronSplitAlignments::Enumeration 
        le_SearchIntronSplitAlignments = r_SearchIntronSplitAlignments::Conditional;
    r_IntronSearchType::Enumeration
        le_IntronSearchType = r_IntronSearchType::Exhaustive;
    int li_MaxIntronLength = 2100;
    QString ls_IntronSpliceSites = "GT|AG,GC|AG";
    bool lb_Quiet = false;

    QFile lk_StdOut;
    lk_StdOut.open(stdout, QIODevice::WriteOnly);
    
    RefPtr<QFile> lk_pCsvOutFile;
    
    QIODevice* lk_CsvDevice_ = &lk_StdOut;
    QStringList lk_PeptideFiles;
    
    while (!lk_Arguments.empty())
    {
        QString ls_Key = lk_Arguments.takeFirst();
        if (ls_Key == "--massAccuracy")
        {
            bool lb_Ok = false;
            QString ls_Value = lk_Arguments.takeFirst();
            ld_MassAccuracy = ls_Value.toDouble(&lb_Ok);
            if (!lb_Ok)
            {
                printf("Error: Invalid mass accuracy specified: %s.\n", ls_Value.toStdString().c_str());
                exit(1);
            }
        }
        else if (ls_Key == "--similaritySearch")
        {
            lb_SimilaritySearch = 
                parseValue(lk_Arguments.takeFirst().trimmed().toLower(),
                           QStringList() << "yes" << "no",
                           "similarity search flag") == 0 ? true : false;
        }
        else if (ls_Key == "--distinguishIL")
        {
            lb_DistinguishIL = 
                parseValue(lk_Arguments.takeFirst().trimmed().toLower(),
                           QStringList() << "yes" << "no",
                           "distinguish I/L flag") == 0 ? true : false;
        }
        else if (ls_Key == "--searchImmediate")
        {
            lb_SearchImmediateAlignments = 
                parseValue(lk_Arguments.takeFirst().trimmed().toLower(),
                           QStringList() << "yes" << "no",
                           "search immediate alignments flag") == 0 ? true : false;
        }
        else if (ls_Key == "--searchIntronSplit")
        {
            int li_Value = 
                parseValue(lk_Arguments.takeFirst().trimmed().toLower(),
                           QStringList() << "yes" << "no" << "conditional",
                           "search intron split alignments flag");
            switch (li_Value)
            {
                case 0: 
                    le_SearchIntronSplitAlignments = r_SearchIntronSplitAlignments::Yes;
                    break;
                case 1: 
                    le_SearchIntronSplitAlignments = r_SearchIntronSplitAlignments::No;
                    break;
                case 2: 
                default: 
                    le_SearchIntronSplitAlignments = r_SearchIntronSplitAlignments::Conditional;
                    break;
            }
        }
        else if (ls_Key == "--intronSearchType")
        {
            int li_Value = 
                parseValue(lk_Arguments.takeFirst().trimmed().toLower(),
                           QStringList() << "exhaustive" << "quick",
                           "intron search type");
            switch (li_Value)
            {
                case 0: 
                    le_IntronSearchType = r_IntronSearchType::Exhaustive;
                    break;
                case 1: 
                default:
                    le_IntronSearchType = r_IntronSearchType::Quick;
                    break;
            }
        }
        else if (ls_Key == "--maxIntronLength")
        {
            bool lb_Ok = false;
            QString ls_Value = lk_Arguments.takeFirst();
            li_MaxIntronLength = ls_Value.toDouble(&lb_Ok);
            if (!lb_Ok)
            {
                printf("Error: Invalid max intron length specified: %s.\n", ls_Value.toStdString().c_str());
                exit(1);
            }
        }
        else if (ls_Key == "--intronSpliceSites")
        {
            ls_IntronSpliceSites = lk_Arguments.takeFirst();
        }
        else if (ls_Key == "--peptidesFile")
        {
            QString ls_Path = lk_Arguments.takeFirst();
            lk_PeptideFiles.append(ls_Path);
        }
        else if (ls_Key == "--csvOutputPath")
        {
            lk_pCsvOutFile = RefPtr<QFile>(new QFile(lk_Arguments.takeFirst()));
            lk_pCsvOutFile->open(QIODevice::WriteOnly);
            lk_CsvDevice_ = lk_pCsvOutFile.get_Pointer();
        }
        else if (ls_Key == "--quiet")
            lb_Quiet = true;
        else
        {
            lk_Arguments.push_front(ls_Key);
            break;
        }
    }

    QString ls_IndexFilePath = lk_Arguments.takeFirst();
    // load index file
    RefPtr<k_GpfIndexFile> lk_pGpfIndexFile(new k_GpfIndexFile(ls_IndexFilePath));
    if (!lk_pGpfIndexFile->isGood())
    {
        printf("Error: Unable to load GPF index file %s.\n", ac_Arguments__[1]);
        exit(1);
    }
    
    QList<tk_StringIntPair> lk_QueryPeptides;
    
    foreach (QString ls_Path, lk_PeptideFiles)
    {
        QFile lk_File(ls_Path);
        if (lk_File.open(QIODevice::ReadOnly))
        {
            QTextStream lk_Stream(&lk_File);
            QRegExp lk_RegExp("[\\s,;]+");
            int li_LineCount = 0;
            while (!lk_Stream.atEnd())
            {
                QString ls_Line = lk_Stream.readLine();
                ++li_LineCount;
                QStringList lk_Line = ls_Line.split(lk_RegExp);
                QString ls_Peptide = lk_Line[0];
                qint64 li_Mass = 0;
                if (lk_Line.size() > 1)
                {
                    bool lb_Ok = false;
                    double ld_Mass = lk_Line[1].toDouble(&lb_Ok);
                    if (!lb_Ok)
                    {
                        printf("Error: Invalid precursor mass in %s, line %d.\n", ls_Path.toStdString().c_str(), li_LineCount);
                        exit(1);
                    }
                    li_Mass = (qint64)(ld_Mass * lk_pGpfIndexFile->mi_MassPrecision);
                }
                else
                    li_Mass = lk_pGpfIndexFile->peptideMass(ls_Peptide);
                lk_QueryPeptides << tk_StringIntPair(ls_Peptide, li_Mass);
            }
            lk_File.close();
        }
        else
        {
            printf("Error: Unable to open peptides file %s.\n", ls_Path.toStdString().c_str());
            exit(1);
        }
    }

    // append remaining command line args to query peptides list
    while (!lk_Arguments.empty())
    {
        QString ls_Peptide = lk_Arguments.takeFirst();
        qint64 li_Mass = 0;
        // try to convert next item to mass if there is more
        if (!lk_Arguments.empty())
        {
            bool lb_Ok = false;
            double ld_PrecursorMass = lk_Arguments.first().toDouble(&lb_Ok);
            if (lb_Ok)
            {
                lk_Arguments.takeFirst();
                li_Mass = (qint64)(ld_PrecursorMass * lk_pGpfIndexFile->mi_MassPrecision);
            }
        }
        if (li_Mass == 0)
            li_Mass = lk_pGpfIndexFile->peptideMass(ls_Peptide);
        lk_QueryPeptides << tk_StringIntPair(ls_Peptide, li_Mass);
    }

    k_GpfQuery lk_Query(*(lk_pGpfIndexFile.get_Pointer()), ld_MassAccuracy,
                        lb_SimilaritySearch, lb_DistinguishIL,
                        lb_SearchImmediateAlignments,
                        le_SearchIntronSplitAlignments,
                        le_IntronSearchType, li_MaxIntronLength,
                        ls_IntronSpliceSites, lb_Quiet, lk_CsvDevice_);
    // open new scope for stop watch
    {
        RefPtr<k_StopWatch> lk_pStopWatch;
        if (!lb_Quiet)
            lk_pStopWatch = RefPtr<k_StopWatch>(new k_StopWatch("GPF search took %1.\n"));

        lk_Query.execute(lk_QueryPeptides);
    }
}
