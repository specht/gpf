/*
Copyright (c) 2007-2010 Michael Specht

This file is part of GPF.

GPF is free software: you can redistribute it and/or modify it under 
the terms of the GNU Lesser General Public License as published by the 
Free Software Foundation, either version 3 of the License, or (at your 
option) any later version.

GPF is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GPF.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "GpfIndexFile.h"
#include "GpfQuery.h"
#include "StopWatch.h"


int parseValue(QString as_Value, QStringList ak_Choices, QString as_Description)
{
    for (int i = 0; i < ak_Choices.size(); ++i)
    {
        if (as_Value == ak_Choices[i])
            return i;
    }
    fprintf(stderr, "Error: Invalid %s specified: %s.\n", 
           as_Description.toStdString().c_str(),
           as_Value.toStdString().c_str());
    exit(1);
}


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
    if (ai_ArgumentCount < 3)
    {
        fprintf(stderr, "Usage: gpfquery [options] [GPF index file] [peptide 1] [precursor mass 1 (optional)] [peptide 2] ...\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "For every query peptide, an uncharged precursor mass can be specified, this is useful\n");
        fprintf(stderr, "if the peptide sequence comes from de novo prediction. If no precursor mass is\n");
        fprintf(stderr, "specified, an appropriate precursor mass is calculated from the peptide sequence.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  --massAccuracy <float> (default: 10.0)\n");
        fprintf(stderr, "    Specify the mass accuracy in ppm.\n");
        fprintf(stderr, "  --similaritySearch <yes|no> (default: yes)\n");
        fprintf(stderr, "    Specify whether alignments should be found that are exact in mass,\n");
        fprintf(stderr, "    but not necessarily exact in amino acid sequence.\n");
        fprintf(stderr, "  --distinguishIL <yes|no> (default: no)\n");
        fprintf(stderr, "    Specify whether isoleucine (I) and leucine (L) should be\n");
        fprintf(stderr, "    considered the same when checking for sequence equality.\n");
        fprintf(stderr, "  --searchImmediate <yes|no> (default: yes)\n");
        fprintf(stderr, "    Specify whether immediate (i. e. non-intron split) alignments\n");
        fprintf(stderr, "    should be searched for.\n");
        fprintf(stderr, "  --searchIntronSplit <yes|no|conditional> (default: conditional)\n");
        fprintf(stderr, "    Specify whether non-intron split alignments should be searched\n");
        fprintf(stderr, "    for. Conditional means that intron split alignments for a given\n");
        fprintf(stderr, "    peptide are only searched for if no immediate alignments have\n");
        fprintf(stderr, "    been found.\n");
        fprintf(stderr, "  --maxIntronLength <int> (default: 2100)\n");
        fprintf(stderr, "    Specify the maximum intron length in nucleotides.\n");
        fprintf(stderr, "  --intronSpliceSites <string> (default: 'GT|AG,GC|AG')\n");
        fprintf(stderr, "    Specify possible splice donor/acceptor site consensus sequences.\n");
        fprintf(stderr, "  --printFlankingResidues <int> (default: 5)\n");
        fprintf(stderr, "    Specify how many amino acids flanking a hit should be printed.\n");
        // :TODO: 
        // 1. implement this 
        // 2. find out whether this is slower (it should be, right?)
/*        fprintf(stderr, "  --tagSize <int>\n");
        fprintf(stderr, "    Apart from the tag size the genomic DNA sequence was indexed with,\n");
        fprintf(stderr, "    you can specify a higher tag size here to achieve the same effect\n");
        fprintf(stderr, "    as if the index file had been created with this tag size, although\n");
        fprintf(stderr, "    this is a bit slower compared to using the right index file in the\n");
        fprintf(stderr, "    first place.\n");*/
        fprintf(stderr, "  --peptidesFile <path>\n");
        fprintf(stderr, "    Specify a text file containing the query peptides, one peptide per line.\n");
        fprintf(stderr, "    Optionally, each peptide may be followed by the precursor mass, separated\n");
        fprintf(stderr, "    by a comma, semicolon, or whitespace.\n");
        fprintf(stderr, "  --csvOutputPath <path> (default: stdout)\n");
        fprintf(stderr, "    Specify a target file for CSV output.\n");
        fprintf(stderr, "  --peptidesOutputPath <path>\n");
        fprintf(stderr, "    Specify a target file which all resulting peptides will be written to.\n");
        fprintf(stderr, "  --quiet\n");
        fprintf(stderr, "    Don't print status messages.\n");
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
    int li_PrintFlankingResidues = 5;
    bool lb_Quiet = false;

    QFile* lk_StdOut_ = new QFile();
    lk_StdOut_->open(stdout, QIODevice::WriteOnly);
    QSharedPointer<QFile> lk_pCsvOutFile(lk_StdOut_);
    QSharedPointer<QFile> lk_pPeptidesOutFile;
    
    QStringList lk_PeptideFiles;
    
    while (!lk_Arguments.empty())
    {
        QString ls_OriginalKey = lk_Arguments.takeFirst();
        QString ls_Key = ls_OriginalKey.trimmed();
        if (ls_Key == "--massAccuracy")
        {
            bool lb_Ok = false;
            QString ls_Value = lk_Arguments.takeFirst();
            ld_MassAccuracy = ls_Value.toDouble(&lb_Ok);
            if (!lb_Ok)
            {
                fprintf(stderr, "Error: Invalid mass accuracy specified: %s.\n", ls_Value.toStdString().c_str());
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
            li_MaxIntronLength = ls_Value.toInt(&lb_Ok);
            if (!lb_Ok)
            {
                fprintf(stderr, "Error: Invalid max intron length specified: %s.\n", ls_Value.toStdString().c_str());
                exit(1);
            }
        }
        else if (ls_Key == "--intronSpliceSites")
        {
            ls_IntronSpliceSites = lk_Arguments.takeFirst();
        }
        else if (ls_Key == "--printFlankingResidues")
        {
            bool lb_Ok = false;
            QString ls_Value = lk_Arguments.takeFirst();
            li_PrintFlankingResidues = ls_Value.toInt(&lb_Ok);
            if (!lb_Ok)
            {
                fprintf(stderr, "Error: Invalid flanking residue count specified: %s.\n", ls_Value.toStdString().c_str());
                exit(1);
            }
        }
        else if (ls_Key == "--peptidesFile")
        {
            QString ls_Path = lk_Arguments.takeFirst();
            lk_PeptideFiles.append(ls_Path);
        }
        else if (ls_Key == "--csvOutputPath")
        {
            lk_pCsvOutFile = QSharedPointer<QFile>(new QFile(lk_Arguments.takeFirst()));
            lk_pCsvOutFile->open(QIODevice::WriteOnly);
        }
        else if (ls_Key == "--peptidesOutputPath")
        {
            lk_pPeptidesOutFile = QSharedPointer<QFile>(new QFile(lk_Arguments.takeFirst()));
            lk_pPeptidesOutFile->open(QIODevice::WriteOnly);
        }
        else if (ls_Key == "--quiet")
            lb_Quiet = true;
        else
        {
            lk_Arguments.push_front(ls_OriginalKey);
            break;
        }
    }

    QString ls_IndexFilePath = lk_Arguments.takeFirst();
    // load index file
    QSharedPointer<k_GpfIndexFile> lk_pGpfIndexFile(new k_GpfIndexFile(ls_IndexFilePath));
    if (!lk_pGpfIndexFile->isGood())
    {
        fprintf(stderr, "Error: Unable to load GPF index file %s.\n", ls_IndexFilePath.toStdString().c_str());
        exit(1);
    }
    
    QList<r_Query> lk_QueryPeptides;
    
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
                if (lk_Line.size() == 5)
                {
                    QString ls_Id = lk_Line[0];
                    bool lb_Ok = false;

                    double ld_NTerminalMass = lk_Line[1].toDouble(&lb_Ok);
                    if (!lb_Ok)
                    {
                        fprintf(stderr, "Error: Invalid N-terminal mass in %s, line %d.\n", ls_Path.toStdString().c_str(), li_LineCount);
                        exit(1);
                    }
                    qint64 li_NTerminalMass = (qint64)(ld_NTerminalMass * lk_pGpfIndexFile->mi_MassPrecision);

                    QString ls_Peptide = lk_Line[2];

                    double ld_CTerminalMass = lk_Line[3].toDouble(&lb_Ok);
                    if (!lb_Ok)
                    {
                        fprintf(stderr, "Error: Invalid C-terminal mass in %s, line %d.\n", ls_Path.toStdString().c_str(), li_LineCount);
                        exit(1);
                    }
                    qint64 li_CTerminalMass = (qint64)(ld_CTerminalMass * lk_pGpfIndexFile->mi_MassPrecision);

                    double ld_PrecursorMass = lk_Line[4].toDouble(&lb_Ok);
                    if (!lb_Ok)
                    {
                        fprintf(stderr, "Error: Invalid precursor mass in %s, line %d.\n", ls_Path.toStdString().c_str(), li_LineCount);
                        exit(1);
                    }
                    qint64 li_PrecursorMass = (qint64)(ld_PrecursorMass * lk_pGpfIndexFile->mi_MassPrecision);
                    lk_QueryPeptides << r_Query(ls_Peptide, li_PrecursorMass, li_NTerminalMass, li_CTerminalMass, ls_Id);
                }
                else
                {
                    QString ls_Peptide = lk_Line[0];
                    qint64 li_Mass = 0;
                    if (lk_Line.size() > 1)
                    {
                        bool lb_Ok = false;
                        double ld_Mass = lk_Line[1].toDouble(&lb_Ok);
                        if (!lb_Ok)
                        {
                            fprintf(stderr, "Error: Invalid precursor mass in %s, line %d.\n", ls_Path.toStdString().c_str(), li_LineCount);
                            exit(1);
                        }
                        li_Mass = (qint64)(ld_Mass * lk_pGpfIndexFile->mi_MassPrecision);
                    }
                    else
                        li_Mass = lk_pGpfIndexFile->peptideMass(ls_Peptide);
                    lk_QueryPeptides << r_Query(ls_Peptide, li_Mass);
                }
            }
            lk_File.close();
        }
        else
        {
            fprintf(stderr, "Error: Unable to open peptides file %s.\n", ls_Path.toStdString().c_str());
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
        lk_QueryPeptides << r_Query(ls_Peptide, li_Mass);
    }

    k_GpfQuery lk_Query(*(lk_pGpfIndexFile.data()), ld_MassAccuracy,
                        lb_SimilaritySearch, lb_DistinguishIL,
                        lb_SearchImmediateAlignments,
                        le_SearchIntronSplitAlignments,
                        le_IntronSearchType, li_MaxIntronLength,
                        ls_IntronSpliceSites, li_PrintFlankingResidues,
                        lb_Quiet, lk_pCsvOutFile.data(), lk_pPeptidesOutFile.data());
    // open new scope for stop watch
    {
        QSharedPointer<k_StopWatch> lk_pStopWatch;
        if (!lb_Quiet)
            lk_pStopWatch = QSharedPointer<k_StopWatch>(new k_StopWatch("GPF search took %1.\n"));

        lk_Query.execute(lk_QueryPeptides);
    }
}
