/*
Copyright (c) 2007-2010 Michael Specht

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

#include "GpfQuery.h"
#include "GpfIndexFile.h"
#include "GpfBase.h"


k_GpfQuery::k_GpfQuery(k_GpfIndexFile& ak_GpfIndexFile, QIODevice* ak_Output_)
	: mk_GpfIndexFile(ak_GpfIndexFile)
	, mk_Output_(ak_Output_)
	, md_MassAccuracy(5.0)
	, mi_MinIntronLength(1)
	, mi_MaxIntronLength(2100)
	, ms_IntronSpliceSites("GT|AG,GC|AG")
	, mb_SimilaritySearch(true)
	, mb_ImmediateHitsSufficient(false)
	, mk_pStdOutFile(new QFile())
{
    mk_pStdOutFile->open(stdout, QIODevice::WriteOnly);
    if (!mk_Output_)
        mk_Output_ = mk_pStdOutFile.get_Pointer();
    
	QStringList lk_IntronSpliceSites = ms_IntronSpliceSites.split(",");
	
	// convert human readable intron splice site consensus sequences 
	// into a form which is useful for GPF
	mi_IntronNTermMaxLength = 0;
	mi_IntronCTermMaxLength = 0;
	
	foreach (QString ls_SpliceSite, lk_IntronSpliceSites)
	{
		QStringList lk_SpliceSites = ls_SpliceSite.split("|");
		if (lk_SpliceSites.size() != 2)
		{
			printf("Error: Invalid intron splice donor/acceptor site consensus sequence: %s.\n", ls_SpliceSite.toStdString().c_str());
			exit(1);
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
					printf("Error: Invalid intron splice donor/acceptor site consensus sequence: %s.\n", ls_SpliceSite.toStdString().c_str());
					exit(1);
				}
				li_DnaCode <<= 3;
				li_DnaCode |= li_Code;
			}
			lk_SpliceSitesNumbers.append(tk_IntPair(li_DnaCode, ls_Site.length()));
		}
		tk_IntPair lk_NTermCode = lk_SpliceSitesNumbers.first();
        tk_IntPair lk_CTermCode = lk_SpliceSitesNumbers.last();
        
		if (!mk_IntronNTerm.contains(lk_NTermCode))
			mk_IntronNTerm[lk_NTermCode] = QList<tk_IntPair>();
		mk_IntronNTerm[lk_NTermCode].append(lk_CTermCode);
        
		if (!mk_IntronCTerm.contains(lk_CTermCode))
			mk_IntronCTerm[lk_CTermCode] = QList<tk_IntPair>();
        mk_IntronCTerm[lk_CTermCode].append(lk_NTermCode);
        
		if (lk_NTermCode.second > mi_IntronNTermMaxLength)
			mi_IntronNTermMaxLength = lk_NTermCode.second;
        
		if (lk_CTermCode.second > mi_IntronCTermMaxLength)
			mi_IntronCTermMaxLength = lk_CTermCode.second;
	}
	
/*	foreach (tk_IntPair lk_Key, mk_IntronEnd.keys())
    {
        printf("%x/%d:\n", lk_Key.first, lk_Key.second);
        foreach (tk_IntPair lk_Value, mk_IntronEnd[lk_Key])
            printf("  - %x/%d\n", lk_Value.first, lk_Value.second);
    }*/
}


k_GpfQuery::~k_GpfQuery()
{
}


void k_GpfQuery::execute(const QString& as_Peptide)
{
/*    for (int li_TagDirectionIndex = 0; li_TagDirectionIndex < mk_GpfIndexFile.mi_TagCount * 2; ++li_TagDirectionIndex)
    {
        for (qint64 i = 0; i < mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]; ++i)
        {
            qint64 li_Mass = mk_GpfIndexFile.readIndexBits(mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits + i * mk_GpfIndexFile.mi_MassBits, mk_GpfIndexFile.mi_MassBits);
            printf("%x M %d\n", li_TagDirectionIndex, (unsigned int)li_Mass);
        }
        for (qint64 i = 0; i < mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]; ++i)
        {
            qint64 li_Offset = mk_GpfIndexFile.readIndexBits(mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits + mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex] * mk_GpfIndexFile.mi_MassBits + i * mk_GpfIndexFile.mi_OffsetBits, mk_GpfIndexFile.mi_OffsetBits);
            printf("%x O %d\n", li_TagDirectionIndex, (unsigned int)li_Offset);
        }
    }

    return;*/
    ////////////////////////////////////////////////
    
	printf("Searching for %s...\n", as_Peptide.toStdString().c_str());
    
    QTextStream lk_OutStream(mk_Output_);
	
	// check whether this is a valid peptide and determine mass
	mi_Mass = mk_GpfIndexFile.mi_WaterMass;
	for (int i = 0; i < as_Peptide.length(); ++i)
	{
		if (!gk_GpfBase.mb_IsAminoAcid_[(int)(as_Peptide.at(i).toAscii())])
		{
			printf("Error: Invalid amino acids in %s.\n", as_Peptide.toStdString().c_str());
			return;
		}
		mi_Mass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)as_Peptide.at(i).toAscii()];
	}
	qint64 li_MassDelta = (qint64)(md_MassAccuracy * mi_Mass / 1000000.0);
	mi_MinMass = mi_Mass - li_MassDelta;
	mi_MaxMass = mi_Mass + li_MassDelta;
	
	// determine maximum nucleotide count
	mi_MaxNucleotideSpanLength = (mi_MaxMass / mk_GpfIndexFile.mi_AminoAcidMasses_[(int)'G'] + 1) * 3 + mi_MaxIntronLength;
		
	// extract all HMST
	QMultiMap<qint32, qint64> lk_AllHmst;
	
	// extract all left HMST
	qint64 li_HalfMass = 0;
	for (int i = 0; i + mk_GpfIndexFile.mi_TagSize <= as_Peptide.length() && li_HalfMass <= mk_GpfIndexFile.mi_MaxMass; ++i)
	{
		QString ls_Tag = as_Peptide.mid(i, mk_GpfIndexFile.mi_TagSize);
		qint32 li_Tag = gk_GpfBase.aminoAcidPolymerCode(ls_Tag.toStdString().c_str(), mk_GpfIndexFile.mi_TagSize) * 2;
//          printf("tag: [%1.4f, %s] (%x)\n", (double)li_HalfMass / mk_GpfIndexFile.mi_MassPrecision, ls_Tag.toStdString().c_str(), li_Tag);

		//printf("%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
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
//         printf("tag: [%s, %1.4f] (%x)\n", ls_Tag.toStdString().c_str(), (double)li_HalfMass / mk_GpfIndexFile.mi_MassPrecision, li_Tag);
		
//  		printf("%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
		lk_AllHmst.insert(li_Tag, li_HalfMass);

		// break loop if no similarity search
		if (!mb_SimilaritySearch)
			break;
		
		li_HalfMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)(as_Peptide.at(i + mk_GpfIndexFile.mi_TagSize - 1).toAscii())];
	}
	
	// mass direction: false == left, true == right
	typedef QPair<qint64, bool> tk_GnoMassDirection;
	// lk_GnoHash: GNO => minimum half mass
	QHash<tk_GnoMassDirection, qint64> lk_GnoHash;
	
	// iterate over all HMST
	for (QMultiMap<qint32, qint64>::const_iterator lk_Iter = lk_AllHmst.begin(); lk_Iter != lk_AllHmst.end(); ++lk_Iter)
	{
		qint32 li_TagDirectionIndex = lk_Iter.key();
		
		qint64 li_HalfMass = lk_Iter.value();
//         printf("li_HalfMass = %9.4f\n", (double)li_HalfMass / mk_GpfIndexFile.mi_MassPrecision);
        qint64 li_HalfMassDelta = (qint64)(md_MassAccuracy * li_HalfMass / 1000000.0);
		qint64 li_MinMass = li_HalfMass - li_HalfMassDelta;
        qint64 li_MaxMass = li_HalfMass + li_HalfMassDelta;
		
		// determine sub range in HMST list (via min and max masses)
		//printf("tag/dir %8d: %d entries.\n", li_TagDirectionIndex, (qint32)mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]);
		qint64 li_Start = -1;
		qint64 li_Count = 0;
		QList<qint64> lk_Masses;
		for (qint64 i = 0; i < mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]; ++i)
		{
			qint64 li_Mass = mk_GpfIndexFile.readIndexBits(mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits + i * mk_GpfIndexFile.mi_MassBits, mk_GpfIndexFile.mi_MassBits);
/*			printf("[%x] %9.4f %9.4f %9.4f\n", 
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
		//printf("%d - %d\n", (qint32)li_Start, (qint32)(li_Count));
		
		// now we have the correct range, read the corresponding GNOs
		for (qint64 i = 0; i < li_Count; ++i)
		{
			qint64 li_Gno = mk_GpfIndexFile.readIndexBits(
				mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits +
				mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex] * mk_GpfIndexFile.mi_MassBits + 
				(i + li_Start) * mk_GpfIndexFile.mi_OffsetBits, mk_GpfIndexFile.mi_OffsetBits);
			tk_GnoMassDirection lk_GnoMassDirection(li_Gno, li_TagDirectionIndex % 2 == 1);
			
			// insert or update this entry if current half mass is lower
			if (!lk_GnoHash.contains(lk_GnoMassDirection))
				lk_GnoHash[lk_GnoMassDirection] = lk_Masses[i];
			else
				if (lk_Masses[i] < lk_GnoHash[lk_GnoMassDirection])
					lk_GnoHash[lk_GnoMassDirection] = lk_Masses[i];
		}
	}
//  	printf("distinct GNO count (anchors in DNA): %d\n", lk_GnoHash.size());
	
	// now we have determined all interesting places in the genome, take a look at each
	// of them and try to construct alignments with the correct mass
    
    QSet<QString> lk_FoundAssmeblies;
	
	foreach (tk_GnoMassDirection lk_GnoMassDirection, lk_GnoHash.keys())
	{
        // initialize numbers
		qint64 li_HalfMass = lk_GnoHash[lk_GnoMassDirection];
		qint64 li_Gno = lk_GnoMassDirection.first;
		bool lb_RightHalfMass = lk_GnoMassDirection.second;
		bool lb_BackwardsFrame = ((li_Gno & mk_GpfIndexFile.mi_GnoBackwardsBit) != 0);
        // li_DnaOffset is the current pointer, anchor is the first exon, hook
        // is the secondary exon, if any
		qint64 li_DnaOffset = li_Gno & ~mk_GpfIndexFile.mi_GnoBackwardsBit;
        qint64 li_AnchorExonOffset = li_DnaOffset;
        qint64 li_AnchorExonLength = 0;
        
        printf("STARTING SEARCH at %d\n", (unsigned int)li_DnaOffset);
        
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
			printf("WRONG scaffold borders: %d, %d (%d)\n", (qint32)li_ScaffoldStart, (qint32)li_ScaffoldEnd, (qint32)li_DnaOffset);
		
		char* lc_TripletToAminoAcid_ = lb_BackwardsFrame ? 
			gk_GpfBase.mk_DnaTripletToAminoAcidReverse_ :
			gk_GpfBase.mk_DnaTripletToAminoAcid_;
		
        tk_IntPairListHash* lk_IntronStart_ = lb_RightHalfMass ? &mk_IntronCTerm : &mk_IntronNTerm;
        int li_IntronStartMaxLength = lb_RightHalfMass ? mi_IntronCTermMaxLength : mi_IntronNTermMaxLength;
        int li_IntronEndMaxLength = lb_RightHalfMass ? mi_IntronNTermMaxLength : mi_IntronCTermMaxLength;
        
		if (lb_BackwardsFrame)
			li_DnaOffset -= 2;
		qint64 li_Step1 = (lb_BackwardsFrame ^ lb_RightHalfMass) ? -1 : 1;
        qint64 li_Step2 = li_Step1 * 2;
		qint64 li_Step3 = li_Step1 * 3;
        qint64 li_BackwardsFactor = (lb_BackwardsFrame ^ lb_RightHalfMass) ? 1 : 0;
        
        printf("%s%s\n", lb_BackwardsFrame ? "-" : "+", lb_RightHalfMass ? "R" : "L");
        printf("step is %d\n", (qint32)li_Step1);
		
		// try to assemble an alignment
		
		qint64 li_AssemblyMass = mk_GpfIndexFile.mi_WaterMass;
		
		QString ls_Peptide;
		int li_TagAminoAcidsPassed = 0;

		while (li_AssemblyMass < mi_MaxMass)
		{
			quint16 lui_Triplet = gk_GpfBase.readNucleotideTriplet(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), li_DnaOffset + li_BackwardsFactor * li_Step2);
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
            li_AnchorExonLength += 3;
            if (lb_RightHalfMass)
                li_AnchorExonOffset += li_Step3;

            printf("%c %d\n", lc_AminoAcid, (unsigned int)li_DnaOffset);
            
			if (li_TagAminoAcidsPassed >= mk_GpfIndexFile.mi_TagSize)
			{
				// we have passed the half mass and the tag, now we can create alignments!

				// maybe we found an immediate hit already?
				if (li_AssemblyMass >= mi_MinMass && li_AssemblyMass <= mi_MaxMass)
				{
					//printf("IMMEDIATE HIT: %s %d %d\n", ls_Peptide.toStdString().c_str(), (qint32)(li_AssemblyGno & ~mk_GpfIndexFile.mi_GnoBackwardsBit), (qint32)li_AssemblyLength);
                    QString ls_Assembly = QString("%1:%2;%3%4:%5")
                        .arg(mk_GpfIndexFile.ms_ShortId)
                        .arg(mk_GpfIndexFile.mk_ScaffoldLabels[li_FirstScaffold])
                        .arg(lb_BackwardsFrame ? "-" : "+")
                        .arg(li_AnchorExonOffset - li_ScaffoldStart)
                        .arg(li_AnchorExonLength);
                    if (!lk_FoundAssmeblies.contains(ls_Assembly))
                    {
                        lk_OutStream << QString("%1,%2,%3,\"%4\"\n")
                            .arg(as_Peptide)
                            .arg(ls_Peptide)
                            .arg((double)li_AssemblyMass / mk_GpfIndexFile.mi_MassPrecision)
                            .arg(ls_Assembly);
                    }
                    lk_FoundAssmeblies << ls_Assembly;
				}
				
				// or maybe we find the start of an intron?
				for (int li_IntronStartShift = 0; li_IntronStartShift < 3; ++li_IntronStartShift)
				{
					qint32 li_SplitTriplet = 0;
					if (li_IntronStartShift > 0)
						li_SplitTriplet = readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), (li_DnaOffset + li_BackwardsFactor * (li_IntronStartShift - 1)) * 3, li_IntronStartShift * 3);
					qint32 li_IntronStartNucleotides = readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), (li_DnaOffset + li_Step1 * li_IntronStartShift + li_BackwardsFactor * (li_IntronStartMaxLength - 1)) * 3, li_IntronStartMaxLength * 3);
 					printf("%d %s\n", li_IntronStartShift, gk_GpfBase.nucleotideSequenceForCode(li_IntronStartNucleotides, 2).toStdString().c_str());
					foreach (tk_IntPair lk_IntronStartSite, lk_IntronStart_->keys())
					{
						qint32 li_ThisIntronStartNucleotides = li_IntronStartNucleotides & ((1 << (lk_IntronStartSite.second * 3)) - 1);
						if (li_ThisIntronStartNucleotides == lk_IntronStartSite.first)
						{
							// yay! we found an intron start site
                            qint64 li_IntronStartOffset = li_DnaOffset + li_Step1 * li_IntronStartShift;
                            qint64 li_IntronEndOffset = li_IntronStartOffset + li_Step1 * (mi_MinIntronLength + lk_IntronStartSite.second);
 							printf("INTRON START: %s %d (%s%s) %d @ %d\n", ls_Peptide.toStdString().c_str(), li_IntronStartShift, lb_BackwardsFrame ? "-" : "+", lb_RightHalfMass ? "R": "L", li_ThisIntronStartNucleotides, (unsigned int)li_IntronStartOffset);
                            //qint64 li_FirstExonStart = li_IntronStartOffset;
                            //qint64 li_FirstExonLength = li_AnchorExonLength - 3;
                            // find matching intron end site
                            while (true)
                            {
                                // look for matching intron end sites until max intron length is reached or scaffold ends
                                qint32 li_IntronEndNucleotides = readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), (li_IntronEndOffset + li_Step1 * li_BackwardsFactor * (li_IntronEndMaxLength - 1)) * 3, li_IntronEndMaxLength * 3);
                                foreach (tk_IntPair lk_IntronEndSite, (*lk_IntronStart_)[lk_IntronStartSite])
                                {
                                    qint32 li_ThisIntronEndNucleotides = li_IntronEndNucleotides & ((1 << (lk_IntronEndSite.second * 3)) - 1);
                                    if (li_ThisIntronEndNucleotides == lk_IntronEndSite.first)
                                    {
                                        QString ls_IntronSplitPeptide = ls_Peptide;
                                        qint64 li_IntronSplitAssemblyMass = li_AssemblyMass;
                                        qint64 li_HookExonOffset = li_IntronEndOffset + lk_IntronEndSite.second * li_Step1;
                                        qint64 li_HookExonLength = (3 - li_IntronStartShift) % 3;
                                        qint64 li_SecondExonStart = li_IntronEndOffset + lk_IntronEndSite.second * li_Step1;
                                        qint64 li_SecondExonOffset = li_SecondExonStart;
                                        printf("INTRON END: %d, length: %d (%x)\n", (unsigned int)li_SecondExonStart, (unsigned int)abs(li_SecondExonStart - li_IntronStartOffset), li_ThisIntronEndNucleotides);
                                        
                                        // try to complete a split triplet, if there is any
                                        if (li_IntronStartShift > 0)
                                        {
                                            qint32 li_AssembledSplitTriplet = li_SplitTriplet;
                                            li_AssembledSplitTriplet |= readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), li_SecondExonOffset * li_Step1 * 3, (3 - li_IntronStartShift) * 3) << (li_IntronStartShift * 3);
                                            li_SecondExonOffset += (3 - li_IntronStartShift) * li_Step1;
                                            char lc_AminoAcid = lc_TripletToAminoAcid_[li_AssembledSplitTriplet];
//                                             printf("COMPLETING SPLIT TRIPLET: %c\n", lc_AminoAcid);
                                            if (gk_GpfBase.mb_IsAminoAcid_[(int)lc_AminoAcid])
                                            {
                                                ls_IntronSplitPeptide += QChar(lc_AminoAcid);
                                                li_IntronSplitAssemblyMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)lc_AminoAcid];
                                            }
                                            else
                                                break;
                                        }
                                        
                                        while (li_IntronSplitAssemblyMass <= mi_MaxMass)
                                        {
                                            if (li_IntronSplitAssemblyMass >= mi_MinMass && li_IntronSplitAssemblyMass <= mi_MaxMass)
                                            {
                                                // we found an intron hit!
//                                                 printf("INTRON HIT: %s\n", ls_IntronSplitPeptide.toStdString().c_str());
                                                tk_IntPair lk_StartExon(li_AnchorExonOffset, li_AnchorExonLength);
                                                tk_IntPair lk_EndExon(li_HookExonOffset, li_HookExonLength);
                                                    
                                                if (lb_RightHalfMass)
                                                {
                                                    tk_IntPair lk_Temp = lk_StartExon;
                                                    lk_StartExon = lk_EndExon;
                                                    lk_EndExon = lk_Temp;
                                                }
                                                
                                                QString ls_Assembly = QString("%1:%2;%3%4:%5,%6:%7")
                                                    .arg(mk_GpfIndexFile.ms_ShortId)
                                                    .arg(mk_GpfIndexFile.mk_ScaffoldLabels[li_FirstScaffold])
                                                    .arg(lb_BackwardsFrame ? "-" : "+")
                                                    .arg(lk_StartExon.first)
                                                    .arg(lk_StartExon.second)
                                                    .arg(lk_EndExon.first)
                                                    .arg(lk_EndExon.second);
                                                if (!lk_FoundAssmeblies.contains(ls_Assembly))
                                                {
                                                    lk_OutStream << QString("%1,%2,%3,\"%4\"\n")
                                                        .arg(as_Peptide)
                                                        .arg(ls_IntronSplitPeptide)
                                                        .arg((double)li_IntronSplitAssemblyMass / mk_GpfIndexFile.mi_MassPrecision)
                                                        .arg(ls_Assembly);
                                                }
                                                lk_FoundAssmeblies << ls_Assembly;
                                            }
                                            printf("reading from %d\n", (unsigned int)(li_SecondExonOffset - li_BackwardsFactor * 2));
                                            quint32 li_Triplet = readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), (li_SecondExonOffset - li_BackwardsFactor * 2) * 3, 9);
                                            li_SecondExonOffset += li_Step3;
                                            li_HookExonLength += 3;
                                            if (lb_BackwardsFrame)
                                                li_HookExonOffset += li_Step3;
                                            char lc_AminoAcid = lc_TripletToAminoAcid_[li_Triplet];
                                            if (gk_GpfBase.mb_IsAminoAcid_[(int)lc_AminoAcid])
                                            {
                                                if (lb_RightHalfMass)
                                                    ls_IntronSplitPeptide.prepend(QChar(lc_AminoAcid));
                                                else
                                                    ls_IntronSplitPeptide.append(QChar(lc_AminoAcid));
                                                li_IntronSplitAssemblyMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)lc_AminoAcid];
                                            }
                                            else
                                                break;
                                        }
                                        
                                        // add amino acids while peptide mass is too small
                                    }
                                }
                                li_IntronEndOffset += li_Step1;
                                qint64 li_IntronLength = abs(li_IntronEndOffset - li_IntronStartOffset) + 1;
                                if (li_IntronLength > mi_MaxIntronLength)
                                    break;
                            }
						}
					}
				}
			}
			
			// break if out of scaffold
			if (li_DnaOffset < li_ScaffoldStart || li_DnaOffset + 2 > li_ScaffoldEnd)
				break;
		}
		//printf("\n");
	}
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
