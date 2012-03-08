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

#include "GpfQuery.h"
#include "GpfIndexFile.h"
#include "GpfBase.h"


k_GpfQuery::k_GpfQuery(k_GpfIndexFile& ak_GpfIndexFile, double ad_MassAccuracy,
                       bool ab_SimilaritySearch, bool ab_DistinguishIL,
                       bool ab_SearchImmediate, r_SearchIntronSplitAlignments::Enumeration
                       ae_SearchIntronSplitAlignments, r_IntronSearchType::Enumeration
                       ae_IntronSearchType, int ai_MaxIntronLength, 
                       QString as_IntronSpliceSites,
                       int ai_PrintFlankingResidues, bool ab_Quiet,
                       QIODevice* ak_CsvOutput_, QIODevice* ak_PeptidesOutput_)
    : mk_GpfIndexFile(ak_GpfIndexFile)
    , md_MassAccuracy(ad_MassAccuracy)
    , mb_SimilaritySearch(ab_SimilaritySearch)
    , mb_DistinguishIL(ab_DistinguishIL)
    , mb_SearchImmediateAlignments(ab_SearchImmediate)
    , me_SearchIntronSplitAlignments(ae_SearchIntronSplitAlignments)
    , me_IntronSearchType(ae_IntronSearchType)
    , mi_MinIntronLength(1)
    , mi_MaxIntronLength(ai_MaxIntronLength)
    , mi_MinExonLength(1)
    , mi_FlankingSequenceLength(ai_PrintFlankingResidues)
    , ms_IntronSpliceSites(as_IntronSpliceSites)
    , mb_Quiet(ab_Quiet)
    , mk_CsvOutput_(ak_CsvOutput_)
    , mk_CsvOutStream(mk_CsvOutput_)
    , mk_PeptidesOutput_(ak_PeptidesOutput_)
    , mk_PeptidesOutStream(mk_PeptidesOutput_)
{
    QStringList lk_IntronSpliceSites = ms_IntronSpliceSites.split(",");
    
    if (mk_CsvOutput_)
    {
        mk_CsvOutStream << "Query";
        if (mi_FlankingSequenceLength > 0)
            mk_CsvOutStream << ",AA N-term";
        mk_CsvOutStream << ",Peptide";
        if (mi_FlankingSequenceLength > 0)
            mk_CsvOutStream << ",AA C-term";
        mk_CsvOutStream << ",Mass,Assembly,Intron length,Splice site\n";
    }
    
    // convert human readable intron splice site consensus sequences 
    // into a form which is useful for GPF
    mi_IntronNTermMaxLength = 0;
    mi_IntronCTermMaxLength = 0;
    
    foreach (QString ls_SpliceSite, lk_IntronSpliceSites)
    {
        QStringList lk_SpliceSites = ls_SpliceSite.split("|");
        if (lk_SpliceSites.size() != 2)
        {
            fprintf(stderr, "Error: Invalid intron splice donor/acceptor site consensus sequence: %s.\n", ls_SpliceSite.toStdString().c_str());
            exit(1);
        }
        for (int i = 0; i < 2; ++i)
        {
            if (lk_SpliceSites[i].length() > 10)
            {
                fprintf(stderr, "Error: A splice site consensus sequence must not be longer than 10 nucleotides.\n");
                exit(1);
            }
        }
        QList<tk_IntPair> lk_SpliceSitesNumbers;
        while (!lk_SpliceSites.empty())
        {
            QString ls_Site = lk_SpliceSites.takeFirst();
            qint32 li_DnaCode = 0;
            for (int i = ls_Site.length() - 1; i >= 0; --i)
            {
                qint32 li_Code = gk_GpfBase.mk_DnaCharToNumber_[(int)ls_Site.at(i).toAscii()];
                if ((li_Code & 4) != 0)
                {
                    fprintf(stderr, "Error: Invalid intron splice donor/acceptor site consensus sequence: %s.\n", ls_SpliceSite.toStdString().c_str());
                    exit(1);
                }
                li_DnaCode <<= 3;
                li_DnaCode |= li_Code;
            }
            lk_SpliceSitesNumbers.append(tk_IntPair(li_DnaCode, ls_Site.length()));
        }
        tk_IntPair lk_NTermCode = lk_SpliceSitesNumbers.first();
        tk_IntPair lk_CTermCode = lk_SpliceSitesNumbers.last();
        tk_IntPair lk_NTermCodeReverse(reverseNucleotides(lk_NTermCode.first, lk_NTermCode.second), lk_NTermCode.second);
        tk_IntPair lk_CTermCodeReverse(reverseNucleotides(lk_CTermCode.first, lk_CTermCode.second), lk_CTermCode.second);
        
        if (!mk_IntronNTerm.contains(lk_NTermCode))
            mk_IntronNTerm[lk_NTermCode] = QList<tk_IntPair>();
        mk_IntronNTerm[lk_NTermCode].append(lk_CTermCode);
        
        if (!mk_IntronCTerm.contains(lk_CTermCode))
            mk_IntronCTerm[lk_CTermCode] = QList<tk_IntPair>();
        mk_IntronCTerm[lk_CTermCode].append(lk_NTermCode);
        
        if (!mk_IntronNTermReverse.contains(lk_NTermCodeReverse))
            mk_IntronNTermReverse[lk_NTermCodeReverse] = QList<tk_IntPair>();
        mk_IntronNTermReverse[lk_NTermCodeReverse].append(lk_CTermCodeReverse);
        
        if (!mk_IntronCTermReverse.contains(lk_CTermCodeReverse))
            mk_IntronCTermReverse[lk_CTermCodeReverse] = QList<tk_IntPair>();
        mk_IntronCTermReverse[lk_CTermCodeReverse].append(lk_NTermCodeReverse);

        // update C- and N-term max length
        if (lk_NTermCode.second > mi_IntronNTermMaxLength)
            mi_IntronNTermMaxLength = lk_NTermCode.second;
        
        if (lk_CTermCode.second > mi_IntronCTermMaxLength)
            mi_IntronCTermMaxLength = lk_CTermCode.second;
    }
    mk_IntronNTermKeys = mk_IntronNTerm.keys();
    qSort(mk_IntronNTermKeys.begin(), mk_IntronNTermKeys.end(), &sortByDecreasingLength);
    mk_IntronCTermKeys = mk_IntronCTerm.keys();
    qSort(mk_IntronCTermKeys.begin(), mk_IntronCTermKeys.end(), &sortByDecreasingLength);
    mk_IntronNTermReverseKeys = mk_IntronNTermReverse.keys();
    qSort(mk_IntronNTermReverseKeys.begin(), mk_IntronNTermReverseKeys.end(), &sortByDecreasingLength);
    mk_IntronCTermReverseKeys = mk_IntronCTermReverse.keys();
    qSort(mk_IntronCTermReverseKeys.begin(), mk_IntronCTermReverseKeys.end(), &sortByDecreasingLength);
    
    foreach (tk_IntPair lk_IntPair, mk_IntronNTerm.keys())
        qSort(mk_IntronNTerm[lk_IntPair].begin(), mk_IntronNTerm[lk_IntPair].end(), &sortByDecreasingLength);
    foreach (tk_IntPair lk_IntPair, mk_IntronCTerm.keys())
        qSort(mk_IntronCTerm[lk_IntPair].begin(), mk_IntronCTerm[lk_IntPair].end(), &sortByDecreasingLength);
    foreach (tk_IntPair lk_IntPair, mk_IntronNTermReverse.keys())
        qSort(mk_IntronNTermReverse[lk_IntPair].begin(), mk_IntronNTermReverse[lk_IntPair].end(), &sortByDecreasingLength);
    foreach (tk_IntPair lk_IntPair, mk_IntronCTermReverse.keys())
        qSort(mk_IntronCTermReverse[lk_IntPair].begin(), mk_IntronCTermReverse[lk_IntPair].end(), &sortByDecreasingLength);
}


k_GpfQuery::~k_GpfQuery()
{
}


void k_GpfQuery::execute(const QString& as_Peptide, qint64 ai_PrecursorMass)
{
    ms_QueryPeptide = as_Peptide;
    ms_QueryPeptideIL = as_Peptide;
    ms_QueryPeptideIL.replace("I", "L");
    
    // check whether this is a valid peptide and determine mass
    qint64 li_Mass = ai_PrecursorMass;
    qint64 li_MassDelta = (qint64)(md_MassAccuracy * li_Mass / 1000000.0);
    mi_AlignmentMinMass = li_Mass - li_MassDelta;
    mi_AlignmentMaxMass = li_Mass + li_MassDelta;
    
    // extract all HMST
    QMultiMap<qint32, qint64> lk_AllHmst;
    
    // extract all left HMST
    qint64 li_HalfMass = 0;
    for (int i = 0; i + mk_GpfIndexFile.mi_TagSize <= as_Peptide.length() && li_HalfMass <= mk_GpfIndexFile.mi_MaxMass; ++i)
    {
        QString ls_Tag = as_Peptide.mid(i, mk_GpfIndexFile.mi_TagSize);
        qint32 li_Tag = gk_GpfBase.aminoAcidPolymerCode(ls_Tag.toStdString().c_str(), mk_GpfIndexFile.mi_TagSize) * 2;
//          fprintf(stderr, "tag: [%1.4f, %s] (%x)\n", (double)li_HalfMass / mk_GpfIndexFile.mi_MassPrecision, ls_Tag.toStdString().c_str(), li_Tag);

        //fprintf(stderr, "%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
        lk_AllHmst.insert(li_Tag, li_HalfMass);
        
        // break loop if no similarity search
        if (!mb_SimilaritySearch)
            break;
        
        li_HalfMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)(as_Peptide.at(i).toAscii())];
    }
    
    // extract all right HMST
    li_HalfMass = 0;
    for (int i = as_Peptide.length() - mk_GpfIndexFile.mi_TagSize; i >= 0 && li_HalfMass <= mk_GpfIndexFile.mi_MaxMass; --i)
    {
        QString ls_Tag = as_Peptide.mid(i, mk_GpfIndexFile.mi_TagSize);
        qint32 li_Tag = gk_GpfBase.aminoAcidPolymerCode(ls_Tag.toStdString().c_str(), mk_GpfIndexFile.mi_TagSize) * 2 + 1;
//         fprintf(stderr, "tag: [%s, %1.4f] (%x)\n", ls_Tag.toStdString().c_str(), (double)li_HalfMass / mk_GpfIndexFile.mi_MassPrecision, li_Tag);
        
//          fprintf(stderr, "%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
        lk_AllHmst.insert(li_Tag, li_HalfMass);

        // break loop if no similarity search
        if (!mb_SimilaritySearch)
            break;
        
        li_HalfMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)(as_Peptide.at(i + mk_GpfIndexFile.mi_TagSize - 1).toAscii())];
    }
    
    // lk_GnoMap: GNO => minimum half mass
    tk_GnoMap lk_GnoMap;
    
    // iterate over all HMST
    for (QMultiMap<qint32, qint64>::const_iterator lk_Iter = lk_AllHmst.begin(); lk_Iter != lk_AllHmst.end(); ++lk_Iter)
    {
        qint32 li_TagDirectionIndex = lk_Iter.key();
        
        qint64 li_HalfMass = lk_Iter.value();
//         fprintf(stderr, "li_HalfMass = %9.4f\n", (double)li_HalfMass / mk_GpfIndexFile.mi_MassPrecision);
        qint64 li_HalfMassDelta = (qint64)(md_MassAccuracy * li_HalfMass / 1000000.0);
        qint64 li_MinMass = li_HalfMass - li_HalfMassDelta;
        qint64 li_MaxMass = li_HalfMass + li_HalfMassDelta;
        
        // determine sub range in HMST list (via min and max masses)
        //fprintf(stderr, "tag/dir %8d: %d entries.\n", li_TagDirectionIndex, (qint32)mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]);
        qint64 li_Start = -1;
        qint64 li_Count = 0;
        QList<qint64> lk_Masses;
        for (qint64 i = 0; i < mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]; ++i)
        {
            qint64 li_Mass = mk_GpfIndexFile.readIndexBits(mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits + i * mk_GpfIndexFile.mi_MassBits, mk_GpfIndexFile.mi_MassBits);
/*          fprintf(stderr, "[%x] %9.4f %9.4f %9.4f\n", 
                (qint32)li_TagDirectionIndex, 
                (double)li_Mass / mk_GpfIndexFile.mi_MassPrecision, 
                (double)li_MinMass / mk_GpfIndexFile.mi_MassPrecision, 
                (double)li_MaxMass / mk_GpfIndexFile.mi_MassPrecision);*/
            if (li_Mass >= li_MinMass && li_Mass <= li_MaxMass && li_Start < 0)
                li_Start = i;
            if (li_Mass >= li_MinMass && li_Mass <= li_MaxMass)
            {
                ++li_Count;
                lk_Masses.append(li_Mass);
            }
        }
        //fprintf(stderr, "%d - %d\n", (qint32)li_Start, (qint32)(li_Count));
        
        // now we have the correct range, read the corresponding GNOs
        for (qint64 i = 0; i < li_Count; ++i)
        {
            qint64 li_Gno = mk_GpfIndexFile.readIndexBits(
                mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits +
                mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex] * mk_GpfIndexFile.mi_MassBits + 
                (i + li_Start) * mk_GpfIndexFile.mi_OffsetBits, mk_GpfIndexFile.mi_OffsetBits);
            tk_GnoMassDirection lk_GnoMassDirection(li_Gno, (li_TagDirectionIndex % 2) == 1);
            
            // insert or update this entry if current half mass is lower
            if (!lk_GnoMap.contains(lk_GnoMassDirection))
                lk_GnoMap.insert(lk_GnoMassDirection, lk_Masses[i]);
            else
                if (lk_Masses[i] < lk_GnoMap[lk_GnoMassDirection])
                    lk_GnoMap.insert(lk_GnoMassDirection, lk_Masses[i]);
        }
    }
    
//     fprintf(stderr, "distinct GNO count (anchors in DNA): %d\n", lk_GnoMap.size());
    
    // now we have determined all interesting places in the genome, take a look at each
    // of them and try to construct alignments with the correct mass

    tk_StringSet lk_FoundAssemblies;
    tk_StringSet lk_AcceptedAssemblies;

    if (mb_SearchImmediateAlignments || (me_SearchIntronSplitAlignments == r_SearchIntronSplitAlignments::Yes))
    {
        this->findAlignments(lk_GnoMap, mb_SearchImmediateAlignments, 
                             me_SearchIntronSplitAlignments == r_SearchIntronSplitAlignments::Yes,
                             lk_FoundAssemblies, lk_AcceptedAssemblies);
    }
    if ((me_SearchIntronSplitAlignments == r_SearchIntronSplitAlignments::Conditional) && lk_AcceptedAssemblies.empty())
    {
        // find intron split alignments because there were no immediate hits
        this->findAlignments(lk_GnoMap, false, true, lk_FoundAssemblies, lk_AcceptedAssemblies);
    }
}


/*
ak_FoundAssemblies contains all assemblies that have been found, but not
necessarily accepted, ak_AcceptedAssemblies contains all found and accepted
assemblies
*/
void k_GpfQuery::findAlignments(const tk_GnoMap& ak_GnoMap,
                                bool ab_SearchImmediate,
                                bool ab_SearchIntronSplit,
                                tk_StringSet& ak_FoundAssemblies,
                                tk_StringSet& ak_AcceptedAssemblies)
{
    for (QMap<tk_GnoMassDirection, qint64>::const_iterator lk_Iter = ak_GnoMap.constBegin(); lk_Iter != ak_GnoMap.constEnd(); ++lk_Iter)
    {
        tk_GnoMassDirection lk_GnoMassDirection = lk_Iter.key();
        qint64 li_HalfMass = lk_Iter.value();
        // initialize numbers
        qint64 li_Gno = lk_GnoMassDirection.first;
        bool lb_RightHalfMass = lk_GnoMassDirection.second;
        bool lb_BackwardsFrame = ((li_Gno & mk_GpfIndexFile.mi_GnoBackwardsBit) != 0);
        // if lb_ProgressIncreasing, the DNA pointer is increased during the search
        bool lb_ProgressIncreasing = !(lb_RightHalfMass ^ lb_BackwardsFrame);
        // li_DnaOffset is the current pointer, anchor is the first exon, hook
        // is the secondary exon, if any
        qint64 li_DnaOffset = li_Gno & (~mk_GpfIndexFile.mi_GnoBackwardsBit);
        // contrary to the GNO offset and length,
        // start and end are in such a way that start <= end
        qint64 li_AnchorExonStart = li_DnaOffset;
        qint64 li_AnchorExonEnd = li_AnchorExonStart;
        
//         fprintf(stderr, "STARTING SEARCH at %d\n", (unsigned int)li_DnaOffset);
        
        // do a binary search to find the right scaffold
        qint32 li_FirstScaffold = 0;
        qint32 li_LastScaffold = mk_GpfIndexFile.mk_ScaffoldStart.size() - 1;
        while (abs(li_FirstScaffold - li_LastScaffold) > 1)
        {
            qint32 li_MidScaffold = li_FirstScaffold + (li_LastScaffold - li_FirstScaffold) / 2;
            if (li_DnaOffset < mk_GpfIndexFile.mk_ScaffoldStart[li_MidScaffold])
                li_LastScaffold = li_MidScaffold;
            else
                li_FirstScaffold = li_MidScaffold;
        }
        // adjust scaffold if necessary
        if (li_DnaOffset >= mk_GpfIndexFile.mk_ScaffoldStart[li_FirstScaffold] + mk_GpfIndexFile.mk_ScaffoldLength[li_FirstScaffold])
            ++li_FirstScaffold;
        
        // li_ScaffoldStart & li_ScaffoldEnd point to the first and last nucleotide
        qint64 li_ScaffoldStart = mk_GpfIndexFile.mk_ScaffoldStart[li_FirstScaffold];
        qint64 li_ScaffoldEnd = li_ScaffoldStart + mk_GpfIndexFile.mk_ScaffoldLength[li_FirstScaffold] - 1;
        
        if (li_DnaOffset < li_ScaffoldStart || li_DnaOffset > li_ScaffoldEnd)
            fprintf(stderr, "Internal error: WRONG scaffold borders: %d, %d (%d)\n", (qint32)li_ScaffoldStart, (qint32)li_ScaffoldEnd, (qint32)li_DnaOffset);
        
//         fprintf(stderr, "Scaffold range is %d - %d.\n", (qint32)li_ScaffoldStart, (qint32)li_ScaffoldEnd);
        
        char* lc_TripletToAminoAcid_ = lb_BackwardsFrame ? 
            gk_GpfBase.mk_TranslationTablesReverse[mk_GpfIndexFile.mi_GeneticCode].data() :
            gk_GpfBase.mk_TranslationTables[mk_GpfIndexFile.mi_GeneticCode].data();
        
        tk_IntPairListHash* lk_IntronStart_ = NULL;
        tk_IntPairList* lk_IntronStartKeys_ = NULL;
        
        if (lb_BackwardsFrame)
        {
            lk_IntronStart_ = lb_RightHalfMass ? &mk_IntronCTermReverse : &mk_IntronNTermReverse;
            lk_IntronStartKeys_ = lb_RightHalfMass ? &mk_IntronCTermReverseKeys : &mk_IntronNTermReverseKeys;
        }
        else
        {
            lk_IntronStart_ = lb_RightHalfMass ? &mk_IntronCTerm : &mk_IntronNTerm;
            lk_IntronStartKeys_ = lb_RightHalfMass ? &mk_IntronCTermKeys : &mk_IntronNTermKeys;
        }
        int li_IntronStartMaxLength = lb_RightHalfMass ? mi_IntronCTermMaxLength : mi_IntronNTermMaxLength;
        int li_IntronEndMaxLength = lb_RightHalfMass ? mi_IntronNTermMaxLength : mi_IntronCTermMaxLength;
        
        qint64 li_Step1 = lb_ProgressIncreasing ? 1 : -1;
        qint64 li_Step2 = li_Step1 * 2;
        qint64 li_Step3 = li_Step1 * 3;
        qint64 li_BackwardsFactor = lb_ProgressIncreasing ? 0 : 1;
        qint64 li_BStep1 = lb_BackwardsFrame ? -1 : 1;
        
//         fprintf(stderr, "progress increasing: %d, step: %d\n", lb_ProgressIncreasing, (qint32)li_Step1);
        // try to assemble an alignment
        
        qint64 li_AssemblyMass = mk_GpfIndexFile.mi_WaterMass;
        
        QString ls_Peptide;
        int li_TagAminoAcidsPassed = 0;

        while (li_AssemblyMass < mi_AlignmentMaxMass)
        {
            // break if out of scaffold
            if (lb_ProgressIncreasing)
            {
                if (li_DnaOffset + 2 > li_ScaffoldEnd)
                    break;
            }
            else
            {
                if (li_DnaOffset - 2 < li_ScaffoldStart)
                    break;
            }
            
            quint16 lui_Triplet = 
                gk_GpfBase.readNucleotideTriplet(
                    mk_GpfIndexFile.muc_pDnaBuffer.data(), 
                    li_DnaOffset + li_BackwardsFactor * li_Step2);
            char lc_AminoAcid = lc_TripletToAminoAcid_[lui_Triplet];
            
            // cancel if invalid amino acid
            if (!gk_GpfBase.mb_IsAminoAcid_[(int)lc_AminoAcid])
                break;

            if (lb_RightHalfMass)
                ls_Peptide.prepend(QChar(lc_AminoAcid));
            else
                ls_Peptide.append(QChar(lc_AminoAcid));
            
            if (li_AssemblyMass >= li_HalfMass)
                ++li_TagAminoAcidsPassed;
            
            li_AssemblyMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)lc_AminoAcid];
            
            li_DnaOffset += li_Step3;
            if (lb_ProgressIncreasing)
                li_AnchorExonEnd += 2;
            else
                li_AnchorExonStart -= 2;

            if (li_TagAminoAcidsPassed >= mk_GpfIndexFile.mi_TagSize)
            {
                // we have passed the half mass and the tag, now we can create alignments!

                // maybe we found an immediate hit already?
                if (ab_SearchImmediate)
                {
                    if (li_AssemblyMass >= mi_AlignmentMinMass && li_AssemblyMass <= mi_AlignmentMaxMass)
                    {
                        qint64 li_PrintAnchorExonStart = li_AnchorExonStart - li_ScaffoldStart;
                        qint64 li_PrintAnchorExonLength = (li_AnchorExonEnd - li_AnchorExonStart) + 1;
                        qint64 li_AnchorExonBegin = li_AnchorExonStart;
                        if (lb_BackwardsFrame)
                        {
                            li_PrintAnchorExonStart += li_PrintAnchorExonLength - 1;
                            li_AnchorExonBegin += li_PrintAnchorExonLength - 1;
                        }
                        QString ls_Assembly = QString("{%1/%2}%3%4:%5")
                            .arg(mk_GpfIndexFile.ms_ShortId)
                            .arg(mk_GpfIndexFile.mk_ScaffoldLabels[li_FirstScaffold])
                            .arg(lb_BackwardsFrame ? "-" : "+")
                            .arg(li_PrintAnchorExonStart)
                            .arg(li_PrintAnchorExonLength);
                        if (!ak_FoundAssemblies.contains(ls_Assembly))
                        {
                            ak_FoundAssemblies << ls_Assembly;
                            bool lb_Accept = true;
                            if (!mb_SimilaritySearch)
                            {
                                if (mb_DistinguishIL)
                                    lb_Accept = ms_QueryPeptide == ls_Peptide;
                                else
                                {
                                    QString ls_Check = ls_Peptide;
                                    ls_Check.replace("I", "L");
                                    lb_Accept = ms_QueryPeptideIL == ls_Check;
                                }
                            }
                            if (lb_Accept)
                            {
                                ak_AcceptedAssemblies << ls_Assembly;
                                if (mk_PeptidesOutput_)
                                {
                                    if (!mk_ResultingPeptides.contains(ls_Peptide))
                                        mk_PeptidesOutStream << ls_Peptide + "\n";
                                    mk_ResultingPeptides << ls_Peptide;
                                }
                                if (mk_CsvOutput_)
                                {
                                    if (mi_FlankingSequenceLength > 0)
                                    {
                                        QString ls_LeftAminoAcids;
                                        QString ls_RightAminoAcids;
                                        readFlankingAminoAcids(li_AnchorExonBegin, 
                                                            li_AnchorExonBegin + li_BStep1 * (li_PrintAnchorExonLength - 1),
                                                            ls_LeftAminoAcids,
                                                            ls_RightAminoAcids, 
                                                            li_BStep1,
                                                            li_ScaffoldStart,
                                                            li_ScaffoldEnd,
                                                            lc_TripletToAminoAcid_);
                                        mk_CsvOutStream << QString("%1,%2,%3,%4,%5,\"%6\",,\n")
                                            .arg(ms_QueryPeptide)
                                            .arg(ls_LeftAminoAcids)
                                            .arg(ls_Peptide)
                                            .arg(ls_RightAminoAcids)
                                            .arg((double)li_AssemblyMass / mk_GpfIndexFile.mi_MassPrecision, 1, 'f', mk_GpfIndexFile.mi_MassDecimalDigits)
                                            .arg(ls_Assembly);
                                    }
                                    else
                                    {
                                        mk_CsvOutStream << QString("%1,%2,%3,\"%4\",,\n")
                                            .arg(ms_QueryPeptide)
                                            .arg(ls_Peptide)
                                            .arg((double)li_AssemblyMass / mk_GpfIndexFile.mi_MassPrecision, 1, 'f', mk_GpfIndexFile.mi_MassDecimalDigits)
                                            .arg(ls_Assembly);
                                    }
                                }
                            }
                        }
                    }
                }
                
                // cancel if any other amino acid would make this assembly too heavy
                if (li_AssemblyMass + mk_GpfIndexFile.mi_MinAminoAcidMass > mi_AlignmentMaxMass)
                    break;
                
                // Now scan the next few nucleotides for intron splice donor/acceptor
                // site consensus sequences. 
/*                qint64 li_SeekStart = li_DnaOffset;
                qint64 li_SeekEnd = li_SeekStart + li_Step1 * (mi_MaxIntronLength - 1);
                if (lb_ProgressIncreasing)
                    li_SeekEnd = std::min(li_SeekEnd, li_ScaffoldEnd - mi_MinExonLength);
                else
                    li_SeekEnd = std::max(li_SeekEnd, li_ScaffoldStart + mi_MinExonLength);*/

                qint64 li_IntronScanPointer = li_DnaOffset;
//                 fprintf(stderr, "starting intron search at %d\n", (qint32)li_IntronScanPointer);

                if (ab_SearchIntronSplit)
                {
                    for (int li_IntronStartOffset = 0; li_IntronStartOffset < 3; ++li_IntronStartOffset)
                    {
                        qint64 li_ReadLength = li_IntronStartMaxLength;
                        if (lb_ProgressIncreasing)
                            li_ReadLength = std::min<qint64>(li_ReadLength, li_ScaffoldEnd - li_IntronScanPointer + 1);
                        else
                            li_ReadLength = std::min<qint64>(li_ReadLength, li_IntronScanPointer - li_ScaffoldStart + 1);
                        if (li_ReadLength > 0)
                        {
    //                         fprintf(stderr, "check %d/%d\n", (qint32)(((li_IntronScanPointer - (li_ReadLength - 1) * li_BackwardsFactor))), (qint32)li_ReadLength);
                            qint32 li_Bit = 
                                readBitsFromBuffer(
                                mk_GpfIndexFile.muc_pDnaBuffer.data(), 
                                ((li_IntronScanPointer - (li_ReadLength - 1) * li_BackwardsFactor)) * 3, 
                                li_ReadLength * 3);
                            qint32 li_BitLength = li_ReadLength;
                            li_Bit &= (1 << (li_BitLength * 3)) - 1;
                            // invert nucleotides if backwards frame
                            if (lb_BackwardsFrame)
                                li_Bit = invertNucleotides(li_Bit, li_BitLength);
                            foreach (tk_IntPair lk_Site, *lk_IntronStartKeys_)
                            {
//                                 fprintf(stderr, "checking [%s\n", gk_GpfBase.nucleotideSequenceForCode(lk_Site.first, lk_Site.second).toStdString().c_str());
                                qint32 li_CutBit = li_Bit;
                                qint32 li_CutBitLength = li_BitLength;
                                if (lk_Site.second < li_CutBitLength)
                                {
                                    // throw away some nucleotides of the bit
                                    if (lb_BackwardsFrame)
                                        li_CutBit >>= (li_CutBitLength - lk_Site.second) * 3;
                                    else
                                        li_CutBit &= ((1 << (lk_Site.second * 3)) - 1);
                                    li_CutBitLength = lk_Site.second;
                                }
                                if ((lk_Site.first == li_CutBit) && (lk_Site.second == li_CutBitLength))
                                {
                                    // we found an intron start site!
    /*                                fprintf(stderr, "[START %s/%s at %d]\n",
                                        gk_GpfBase.nucleotideSequenceForCode(lk_Site.first, lk_Site.second).toStdString().c_str(),
                                        gk_GpfBase.nucleotideSequenceForCode(li_CutBit, li_CutBitLength).toStdString().c_str(),
                                        (qint32)li_IntronAnchorOffset
                                    );
                                    fprintf(stderr, "intron start offset is %d.\n", (qint32)li_IntronStartOffset);*/
                                    // fetch first part of split triplet, if any, from right before 
                                    // the start sequence
                                    qint32 li_SplitTriplet = 0;
                                    if (li_IntronStartOffset > 0)
                                    {
                                        li_SplitTriplet = 
                                            readBitsFromBuffer(
                                            mk_GpfIndexFile.muc_pDnaBuffer.data(), 
                                            ((li_IntronScanPointer - li_IntronStartOffset * li_Step1 - (li_IntronStartOffset - 1) * li_BackwardsFactor)) * 3, 
                                            li_IntronStartOffset * 3);
                                        if (lb_BackwardsFrame)
                                            li_SplitTriplet = transposeNucleotides(li_SplitTriplet, li_IntronStartOffset);
                                    }
    //                                 fprintf(stderr, "start nucleotides are '%s' %08o.\n", gk_GpfBase.nucleotideSequenceForCode(li_SplitTriplet, li_IntronStartOffset).toStdString().c_str(), li_SplitTriplet);
                                    
                                    // now scan the next few nucleotides for an appropriate intron end sequence
                                    qint64 li_SubIntronScanPointer = li_IntronScanPointer;
    //                                 fprintf(stderr, "starting intron search at %d\n", (qint32)li_IntronScanPointer);
                                    qint64 li_SubAnchorExonStart = li_AnchorExonStart;
                                    qint64 li_SubAnchorExonEnd = li_AnchorExonEnd;
                                    
                                    if (lb_ProgressIncreasing)
                                        li_SubAnchorExonEnd += li_IntronStartOffset;
                                    else
                                        li_SubAnchorExonStart -= li_IntronStartOffset;
                                    
                                    for (int li_SubIntronOffset = 0; li_SubIntronOffset < mi_MaxIntronLength - li_IntronStartOffset; ++li_SubIntronOffset)
                                    {
                                        qint64 li_SubReadLength = li_IntronEndMaxLength;
                                        if (lb_ProgressIncreasing)
                                            li_SubReadLength = std::min<qint64>(li_SubReadLength, li_ScaffoldEnd - li_SubIntronScanPointer + 1);
                                        else
                                            li_SubReadLength = std::min<qint64>(li_SubReadLength, li_SubIntronScanPointer - li_ScaffoldStart + 1);
                                        if (li_SubReadLength > 0)
                                        {
                    //                         fprintf(stderr, "check %d/%d\n", (qint32)(((li_IntronScanPointer - (li_ReadLength - 1) * li_BackwardsFactor))), (qint32)li_ReadLength);
                                            qint32 li_SubBit = 
                                                readBitsFromBuffer(
                                                mk_GpfIndexFile.muc_pDnaBuffer.data(), 
                                                ((li_SubIntronScanPointer - (li_SubReadLength - 1) * li_BackwardsFactor)) * 3, 
                                                li_SubReadLength * 3);
                                            qint32 li_SubBitLength = li_SubReadLength;
                                            li_SubBit &= (1 << (li_SubBitLength * 3)) - 1;
                                            // invert nucleotides if backwards frame
                                            if (lb_BackwardsFrame)
                                                li_SubBit = invertNucleotides(li_SubBit, li_SubBitLength);
                                            foreach (tk_IntPair lk_SubSite, (*lk_IntronStart_)[lk_Site])
                                            {
//                                                 fprintf(stderr, "checking %s]\n", gk_GpfBase.nucleotideSequenceForCode(lk_SubSite.first, lk_SubSite.second).toStdString().c_str());
                                                qint32 li_SubCutBit = li_SubBit;
                                                qint32 li_SubCutBitLength = li_SubBitLength;
                                                if (lk_SubSite.second < li_SubCutBitLength)
                                                {
                                                    // throw away some nucleotides of the bit
                                                    if (lb_BackwardsFrame)
                                                        li_SubCutBit >>= (li_SubCutBitLength - lk_SubSite.second) * 3;
                                                    else
                                                        li_SubCutBit &= ((1 << (lk_SubSite.second * 3)) - 1);
                                                    li_SubCutBitLength = lk_SubSite.second;
                                                }
                                                if ((lk_SubSite.first == li_SubCutBit) && (lk_SubSite.second == li_SubCutBitLength))
                                                {
                                                    qint64 li_IntronHookOffset = li_SubIntronScanPointer + li_Step1 * li_SubCutBitLength;
                                                    // we found an intron start site!
    /*                                                fprintf(stderr, "[END %s/%s at %d]\n",
                                                        gk_GpfBase.nucleotideSequenceForCode(lk_SubSite.first, lk_SubSite.second).toStdString().c_str(),
                                                        gk_GpfBase.nucleotideSequenceForCode(li_SubCutBit, li_SubCutBitLength).toStdString().c_str(),
                                                        (qint32)(li_IntronHookOffset - li_Step1)
                                                    );*/
                                                    qint64 li_HookExonEnd = li_IntronHookOffset;
                                                    qint64 li_HookExonStart = li_IntronHookOffset;
                                                    qint64 li_SubDnaOffset = li_IntronHookOffset;
                                                    qint64 li_SubAssemblyMass = li_AssemblyMass;
                                                    
                                                    QString ls_SubPeptide = ls_Peptide;
                                                    bool lb_MustCompleteSplitTriplet = li_IntronStartOffset > 0;
                                                    
                                                    // now continue the peptide
                                                    while (li_SubAssemblyMass < mi_AlignmentMaxMass)
                                                    {
                                                        int li_ReadLength = 3;
                                                        if (lb_MustCompleteSplitTriplet)
                                                        {
                                                            // complete a split triplet instead of reading a whole triplet
                                                            li_ReadLength = 3 - li_IntronStartOffset;
                                                        }
                                                        // break if out of scaffold
                                                        if (lb_ProgressIncreasing)
                                                        {
                                                            if (li_SubDnaOffset + (li_ReadLength - 1) > li_ScaffoldEnd)
                                                                break;
                                                        }
                                                        else
                                                        {
                                                            if (li_SubDnaOffset - (li_ReadLength - 1) < li_ScaffoldStart)
                                                                break;
                                                        }
                                                        
                                                        char lc_SubAminoAcid = 'x';
                                                        if (lb_MustCompleteSplitTriplet)
                                                        {
                                                            // read remaining nucleotides and comlete amino acid
                                                            qint32 li_RemainingNucleotides = 
                                                                readBitsFromBuffer(
                                                                    mk_GpfIndexFile.muc_pDnaBuffer.data(),
                                                                    (li_SubDnaOffset - li_BackwardsFactor * (li_ReadLength - 1)) * 3,
                                                                    li_ReadLength * 3
                                                                );
                                                            if (lb_BackwardsFrame)
                                                                li_RemainingNucleotides = transposeNucleotides(li_RemainingNucleotides, li_ReadLength);
                                                            qint32 li_CombinedTriplet = 0;
                                                            if (lb_RightHalfMass)
                                                                li_CombinedTriplet = concatNucleotides(li_RemainingNucleotides, li_ReadLength, li_SplitTriplet, li_IntronStartOffset);
                                                            else
                                                                li_CombinedTriplet = concatNucleotides(li_SplitTriplet, li_IntronStartOffset, li_RemainingNucleotides, li_ReadLength);
                                                            // make sure it's only the lower 9 bits!
                                                            li_CombinedTriplet &= 511;
                                                            lc_SubAminoAcid = gk_GpfBase.mk_TranslationTables[mk_GpfIndexFile.mi_GeneticCode].data()[li_CombinedTriplet];
    /*                                                        fprintf(stderr, "remaining: %s %08o/ combined: %s / makes %c.\n", 
                                                                gk_GpfBase.nucleotideSequenceForCode(li_RemainingNucleotides, li_ReadLength).toStdString().c_str(), 
                                                                li_RemainingNucleotides,
                                                                gk_GpfBase.nucleotideSequenceForCode(li_CombinedTriplet, 3).toStdString().c_str(), 
                                                                lc_SubAminoAcid);*/
                                                        }
                                                        else
                                                        {
                                                            quint16 lui_SubTriplet = 
                                                                gk_GpfBase.readNucleotideTriplet(
                                                                    mk_GpfIndexFile.muc_pDnaBuffer.data(), 
                                                                    li_SubDnaOffset + li_BackwardsFactor * li_Step2);
                                                            lc_SubAminoAcid = lc_TripletToAminoAcid_[lui_SubTriplet];
                                                        }
                                                        
                                                        // cancel if invalid amino acid
                                                        if (!gk_GpfBase.mb_IsAminoAcid_[(int)lc_SubAminoAcid])
                                                            break;
                                                        
                                                        li_SubDnaOffset += li_Step1 * li_ReadLength;
                                                        if (lb_ProgressIncreasing)
                                                            li_HookExonEnd += (li_ReadLength - 1);
                                                        else
                                                            li_HookExonStart -= (li_ReadLength - 1);
                                                        
                                                        if (lb_RightHalfMass)
                                                            ls_SubPeptide.prepend(QChar(lc_SubAminoAcid));
                                                        else
                                                            ls_SubPeptide.append(QChar(lc_SubAminoAcid));
                                                        
                                                        li_SubAssemblyMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)lc_SubAminoAcid];

                                                        if (li_SubAssemblyMass >= mi_AlignmentMinMass && li_SubAssemblyMass <= mi_AlignmentMaxMass)
                                                        {
                                                            qint64 li_PrintAnchorExonStart = li_SubAnchorExonStart - li_ScaffoldStart;
                                                            qint64 li_PrintAnchorExonLength = (li_SubAnchorExonEnd - li_SubAnchorExonStart) + 1;
                                                            qint64 li_PrintHookExonStart = li_HookExonStart - li_ScaffoldStart;
                                                            qint64 li_PrintHookExonLength = (li_HookExonEnd - li_HookExonStart) + 1;
                                                            qint64 li_FirstExonBegin = li_SubAnchorExonStart;
                                                            qint64 li_FirstExonLength = (li_SubAnchorExonEnd - li_SubAnchorExonStart) + 1;
                                                            qint64 li_SecondExonBegin = li_HookExonStart;
                                                            qint64 li_SecondExonLength = (li_HookExonEnd - li_HookExonStart) + 1;
                                                            if (lb_BackwardsFrame)
                                                            {
                                                                li_FirstExonBegin += li_PrintAnchorExonLength - 1;
                                                                li_SecondExonBegin += li_PrintHookExonLength - 1;
                                                            }
                                                            if (lb_RightHalfMass)
                                                            {
                                                                qint64 li_Temp = li_FirstExonBegin;
                                                                li_FirstExonBegin = li_SecondExonBegin;
                                                                li_SecondExonBegin = li_Temp;
                                                                li_Temp = li_FirstExonLength;
                                                                li_FirstExonLength = li_SecondExonLength;
                                                                li_SecondExonLength = li_Temp;
                                                            }
                                                            if (lb_BackwardsFrame)
                                                            {
                                                                li_PrintAnchorExonStart += li_PrintAnchorExonLength - 1;
                                                                li_PrintHookExonStart += li_PrintHookExonLength - 1;
                                                            }
                                                            if (lb_RightHalfMass)
                                                            {
                                                                // swap exons if right half mass, so that exons are ordered by amino acid sequence
                                                                qint64 li_Temp = li_PrintAnchorExonStart;
                                                                li_PrintAnchorExonStart = li_PrintHookExonStart;
                                                                li_PrintHookExonStart = li_Temp;
                                                                li_Temp = li_PrintAnchorExonLength;
                                                                li_PrintAnchorExonLength = li_PrintHookExonLength;
                                                                li_PrintHookExonLength = li_Temp;
                                                            }
                                                            QString ls_Assembly = QString("{%1/%2}%3%4:%5,%6:%7")
                                                                .arg(mk_GpfIndexFile.ms_ShortId)
                                                                .arg(mk_GpfIndexFile.mk_ScaffoldLabels[li_FirstScaffold])
                                                                .arg(lb_BackwardsFrame ? "-" : "+")
                                                                .arg(li_PrintAnchorExonStart)
                                                                .arg(li_PrintAnchorExonLength)
                                                                .arg(li_PrintHookExonStart)
                                                                .arg(li_PrintHookExonLength);
                                                            if (!ak_FoundAssemblies.contains(ls_Assembly))
                                                            {
                                                                ak_FoundAssemblies << ls_Assembly;
                                                                bool lb_Accept = true;
                                                                if (!mb_SimilaritySearch)
                                                                {
                                                                    if (mb_DistinguishIL)
                                                                        lb_Accept = ms_QueryPeptide == ls_SubPeptide;
                                                                    else
                                                                    {
                                                                        QString ls_Check = ls_SubPeptide;
                                                                        ls_Check.replace("I", "L");
                                                                        lb_Accept = ms_QueryPeptideIL == ls_Check;
                                                                    }
                                                                }
                                                                if (lb_Accept)
                                                                {
                                                                    ak_AcceptedAssemblies << ls_Assembly;
                                                                    if (mk_PeptidesOutput_)
                                                                    {
                                                                        if (!mk_ResultingPeptides.contains(ls_SubPeptide))
                                                                            mk_PeptidesOutStream << ls_SubPeptide + "\n";
                                                                        mk_ResultingPeptides << ls_SubPeptide;
                                                                    }
                                                                    
                                                                    if (mk_CsvOutput_)
                                                                    {
                                                                        QString ls_FirstSite = gk_GpfBase.nucleotideSequenceForCode(lk_Site.first, lk_Site.second);
                                                                        QString ls_SecondSite = gk_GpfBase.nucleotideSequenceForCode(lk_SubSite.first, lk_SubSite.second);
                                                                        if (lb_RightHalfMass)
                                                                        {
                                                                            QString ls_Temp = ls_FirstSite;
                                                                            ls_FirstSite = ls_SecondSite;
                                                                            ls_SecondSite = ls_Temp;
                                                                        }
                                                                        if (lb_BackwardsFrame)
                                                                        {
                                                                            // WTF??! There's no QString::reverse()
                                                                            QString ls_Temp;
                                                                            for (int i = ls_FirstSite.length() - 1; i >= 0; --i)
                                                                                ls_Temp.append(ls_FirstSite.at(i));
                                                                            ls_FirstSite = ls_Temp;
                                                                            ls_Temp = QString();
                                                                            for (int i = ls_SecondSite.length() - 1; i >= 0; --i)
                                                                                ls_Temp.append(ls_SecondSite.at(i));
                                                                            ls_SecondSite = ls_Temp;
                                                                        }
                                                                        QString ls_SpliceSite = QString("%1|%2").arg(ls_FirstSite).arg(ls_SecondSite);

                                                                        int li_IntronLength = abs((li_FirstExonBegin + li_BStep1 * li_FirstExonLength) - 
                                                                            (li_SecondExonBegin));
                                                                        if (mi_FlankingSequenceLength > 0)
                                                                        {
                                                                            QString ls_LeftAminoAcids;
                                                                            QString ls_RightAminoAcids;
                                                                            readFlankingAminoAcids(li_FirstExonBegin, 
                                                                                                li_SecondExonBegin + li_BStep1 * (li_SecondExonLength - 1),
                                                                                                ls_LeftAminoAcids,
                                                                                                ls_RightAminoAcids, 
                                                                                                li_BStep1,
                                                                                                li_ScaffoldStart,
                                                                                                li_ScaffoldEnd,
                                                                                                lc_TripletToAminoAcid_);
                                                                                    
                                                                            mk_CsvOutStream << QString("%1,%2,%3,%4,%5,\"%6\",%7,%8\n")
                                                                                .arg(ms_QueryPeptide)
                                                                                .arg(ls_LeftAminoAcids)
                                                                                .arg(ls_SubPeptide)
                                                                                .arg(ls_RightAminoAcids)
                                                                                .arg((double)li_SubAssemblyMass / mk_GpfIndexFile.mi_MassPrecision, 1, 'f', mk_GpfIndexFile.mi_MassDecimalDigits)
                                                                                .arg(ls_Assembly)
                                                                                .arg(li_IntronLength)
                                                                                .arg(ls_SpliceSite);
                                                                        }
                                                                        else
                                                                        {
                                                                            mk_CsvOutStream << QString("%1,%2,%3,\"%4\",%5,%6\n")
                                                                                .arg(ms_QueryPeptide)
                                                                                .arg(ls_SubPeptide)
                                                                                .arg((double)li_SubAssemblyMass / mk_GpfIndexFile.mi_MassPrecision, 1, 'f', mk_GpfIndexFile.mi_MassDecimalDigits)
                                                                                .arg(ls_Assembly)
                                                                                .arg(li_IntronLength)
                                                                                .arg(ls_SpliceSite);
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        
                                                        if (lb_ProgressIncreasing)
                                                            li_HookExonEnd += 1;
                                                        else
                                                            li_HookExonStart -= 1;
                                                        
                                                        // in any case, read whole triplets from now on
                                                        lb_MustCompleteSplitTriplet = false;
                                                    }
                                                }
                                            }
                                        }
                                        li_SubIntronScanPointer += li_Step1;
                                        // break if we run out of the scaffold
                                        if (li_SubIntronScanPointer < li_ScaffoldStart)
                                            break;
                                        if (li_SubIntronScanPointer > li_ScaffoldEnd)
                                            break;
                                    }
                                    
                                }
                            }
                        }
                        li_IntronScanPointer += li_Step1;
                        // break if we run out of the scaffold
                        if (li_IntronScanPointer < li_ScaffoldStart)
                            break;
                        if (li_IntronScanPointer > li_ScaffoldEnd)
                            break;
                    }
                }
//                 fprintf(stderr, "\n");
            }
            if (lb_ProgressIncreasing)
                li_AnchorExonEnd += 1;
            else
                li_AnchorExonStart -= 1;
        }
    }
}


void k_GpfQuery::execute(const QList<tk_StringIntPair> ak_Peptides)
{
    mk_ResultingPeptides.clear();
    int i = 0;
    foreach (tk_StringIntPair lk_Peptide, ak_Peptides)
    {
        if (!mb_Quiet)
        {
            ++i;
            fprintf(stderr, "\rProcessing query %d of %d... ", i, ak_Peptides.size());
        }
        this->execute(lk_Peptide.first, lk_Peptide.second);
    }
    if (!mb_Quiet)
        fprintf(stderr, " done.\n");
}


int k_GpfQuery::reverseSpliceSequence(int ai_Sequence, int ai_Length)
{
    int li_Result = 0;
    for (int i = 0; i < ai_Length; ++i)
    {
        int li_Nucleotide = ai_Sequence % 7;
        ai_Sequence >>= 3;
        li_Result |= li_Nucleotide << (i * 3);
    }
    return li_Result;
}


void k_GpfQuery::readFlankingAminoAcids(qint64 ai_Start, qint64 ai_Stop,
                                        QString& as_Left, QString& as_Right,
                                        int ai_BStep1, qint64 ai_ScaffoldStart,
                                        qint64 ai_ScaffoldEnd, char* ac_TripletToAminoAcid_)
{
    int li_BBackwardsFactor = ai_BStep1 < 0 ? 1 : 0;
    qint64 li_Pointer = ai_Start - ai_BStep1 * 3;
    QString ls_Left;
    QString ls_Right;
    while (ls_Left.length() < mi_FlankingSequenceLength)
    {
        // break if out of scaffold
        if (li_BBackwardsFactor == 1)
        {
            if (li_Pointer > ai_ScaffoldEnd)
                break;
        }
        else
        {
            if (li_Pointer < ai_ScaffoldStart)
                break;
        }
        // read triplet
        qint32 li_Triplet = gk_GpfBase.
            readNucleotideTriplet(mk_GpfIndexFile.muc_pDnaBuffer.data(), li_Pointer - li_BBackwardsFactor * 2) & 511;
        ls_Left = ac_TripletToAminoAcid_[li_Triplet] + ls_Left;
        li_Pointer -= ai_BStep1 * 3;
    }
    QString ls_RightAminoAcids;
    li_Pointer = ai_Stop + ai_BStep1;
    while (ls_Right.length() < mi_FlankingSequenceLength)
    {
        // break if out of scaffold
        if (li_BBackwardsFactor == 1)
        {
            if (li_Pointer > ai_ScaffoldEnd)
                break;
        }
        else
        {
            if (li_Pointer < ai_ScaffoldStart)
                break;
        }
        // read triplet
        qint32 li_Triplet = gk_GpfBase.
            readNucleotideTriplet(mk_GpfIndexFile.muc_pDnaBuffer.data(), li_Pointer - li_BBackwardsFactor * 2) & 511;
        ls_Right += ac_TripletToAminoAcid_[li_Triplet];
        li_Pointer += ai_BStep1 * 3;
    }
    as_Left = ls_Left;
    as_Right = ls_Right;
}


int sortByDecreasingLength(const tk_IntPair& ak_First, const tk_IntPair& ak_Second)
{
    return ak_First.second > ak_Second.second;
}
