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

#include "GpfQuery.h"
#include "GpfIndexFile.h"
#include "GpfBase.h"

k_GpfQuery::k_GpfQuery(k_GpfIndexFile& ak_GpfIndexFile, const QString& as_Peptide)
	: mk_GpfIndexFile(ak_GpfIndexFile)
	, ms_Peptide(as_Peptide)
	, md_MassAccuracy(5.0)
	, mi_MinIntronLength(0)
	, mi_MaxIntronLength(2100)
	, ms_IntronSpliceSites("GT|AG,GC|AG")
{
	QStringList lk_IntronSpliceSites = ms_IntronSpliceSites.split(",");
	
	// convert human readable intron splice site consensus sequences 
	// into a form which is useful for GPF
	mi_IntronStartMaxLength = 0;
	mi_IntronEndMaxLength = 0;
	
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
		if (!mk_IntronStart.contains(lk_SpliceSitesNumbers.first()))
			mk_IntronStart[lk_SpliceSitesNumbers.first()] = QList<tk_IntPair>();
		mk_IntronStart[lk_SpliceSitesNumbers.first()].append(lk_SpliceSitesNumbers.last());
		if (!mk_IntronEnd.contains(lk_SpliceSitesNumbers.last()))
			mk_IntronEnd[lk_SpliceSitesNumbers.last()] = QList<tk_IntPair>();
		mk_IntronEnd[lk_SpliceSitesNumbers.last()].append(lk_SpliceSitesNumbers.first());
		if (lk_SpliceSitesNumbers.first().second > mi_IntronStartMaxLength)
			mi_IntronStartMaxLength = lk_SpliceSitesNumbers.first().second;
		if (lk_SpliceSitesNumbers.last().second > mi_IntronEndMaxLength)
			mi_IntronEndMaxLength = lk_SpliceSitesNumbers.last().second;
	}
}


k_GpfQuery::~k_GpfQuery()
{
}


void k_GpfQuery::execute()
{
	printf("Searching for %s...\n", ms_Peptide.toStdString().c_str());
	
	// check whether this is a valid peptide and determine mass
	mi_Mass = mk_GpfIndexFile.mi_WaterMass;
	for (int i = 0; i < ms_Peptide.length(); ++i)
	{
		if (!gk_GpfBase.mb_IsAminoAcid_[(int)(ms_Peptide.at(i).toAscii())])
		{
			printf("Error: Invalid amino acids in %s.\n", ms_Peptide.toStdString().c_str());
			return;
		}
		mi_Mass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)ms_Peptide.at(i).toAscii()];
	}
	qint64 li_MassDelta = (qint64)(md_MassAccuracy * mi_Mass / 1000000.0);
	mi_MinMass = mi_Mass - li_MassDelta;
	mi_MaxMass = mi_Mass + li_MassDelta;
	
	// determine maximum nucleotide count
	mi_MaxNucleotideSpanLength = (mi_MaxMass / mk_GpfIndexFile.mi_AminoAcidMasses_[(int)'G'] + 1) * 3 + mi_MaxIntronLength;
		
	// extract all HMST
	QMultiMap<qint32, qint64> lk_AllHmst;
	
	bool lb_SimilaritySearch = true;

	// extract all left HMST
	qint64 li_HalfMass = 0;
	for (int i = 0; i + mk_GpfIndexFile.mi_TagSize <= ms_Peptide.length() && li_HalfMass <= mk_GpfIndexFile.mi_MaxMass; ++i)
	{
		QString ls_Tag = ms_Peptide.mid(i, mk_GpfIndexFile.mi_TagSize);
		qint32 li_Tag = gk_GpfBase.aminoAcidPolymerCode(ls_Tag.toStdString().c_str(), mk_GpfIndexFile.mi_TagSize) * 2;

		//printf("%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
		lk_AllHmst.insert(li_Tag, li_HalfMass);
		
		// break loop if no similarity search
		if (!lb_SimilaritySearch)
			break;
		
		li_HalfMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)(ms_Peptide.at(i).toAscii())];
	}
	
	// extract all right HMST
	li_HalfMass = 0;
	for (int i = ms_Peptide.length() - mk_GpfIndexFile.mi_TagSize; i >= 0 && li_HalfMass <= mk_GpfIndexFile.mi_MaxMass; --i)
	{
		QString ls_Tag = ms_Peptide.mid(i, mk_GpfIndexFile.mi_TagSize);
		qint32 li_Tag = gk_GpfBase.aminoAcidPolymerCode(ls_Tag.toStdString().c_str(), mk_GpfIndexFile.mi_TagSize) * 2 + 1;
		
// 		printf("%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
		lk_AllHmst.insert(li_Tag, li_HalfMass);

		// break loop if no similarity search
		if (!lb_SimilaritySearch)
			break;
		
		li_HalfMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)(ms_Peptide.at(i + mk_GpfIndexFile.mi_TagSize - 1).toAscii())];
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
		qint64 li_MinMass = li_HalfMass;
		qint64 li_MaxMass = li_HalfMass;
		
		// determine sub range in HMST list (via min and max masses)
		//printf("tag/dir %8d: %d entries.\n", li_TagDirectionIndex, (qint32)mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]);
		qint64 li_Start = -1;
		qint64 li_Count = 0;
		QList<qint64> lk_Masses;
		for (qint64 i = 0; i < mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]; ++i)
		{
			qint64 li_Mass = mk_GpfIndexFile.readIndexBits(mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits + i * mk_GpfIndexFile.mi_MassBits, mk_GpfIndexFile.mi_MassBits);
// 			printf("[%d] %d %d %d\n", (qint32)li_TagDirectionIndex, (qint32)li_Mass, (qint32)li_MinMass, (qint32)li_MaxMass);
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
	printf("distinct GNO count (anchors in DNA): %d\n", lk_GnoHash.size());
	
	// now have determined all interesting places in the genome, take a look at each
	// of them and try to construct alignments with the correct mass
	
	foreach (tk_GnoMassDirection lk_GnoMassDirection, lk_GnoHash.keys())
	{
		qint64 li_HalfMass = lk_GnoHash[lk_GnoMassDirection];
		qint64 li_Gno = lk_GnoMassDirection.first;
		bool lb_RightHalfMass = lk_GnoMassDirection.second;
		bool lb_BackwardsFrame = ((li_Gno & mk_GpfIndexFile.mi_GnoBackwardsBit) != 0);
		qint64 li_DnaOffset = li_Gno & ~mk_GpfIndexFile.mi_GnoBackwardsBit;
		//printf("reading GNO 0x%08x (%d): %d/%d\n", (qint32)li_Gno, lb_RightHalfMass, (qint32)li_DnaOffset, lb_BackwardsFrame);
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
		
		qint64 li_ScaffoldStart = mk_GpfIndexFile.mk_ScaffoldStart[li_FirstScaffold];
		qint64 li_ScaffoldEnd = li_ScaffoldStart + mk_GpfIndexFile.mk_ScaffoldLength[li_FirstScaffold] - 1;
		if (li_DnaOffset < li_ScaffoldStart || li_DnaOffset > li_ScaffoldEnd)
			printf("WRONG scaffold borders: %d, %d (%d)\n", (qint32)li_ScaffoldStart, (qint32)li_ScaffoldEnd, (qint32)li_DnaOffset);
		
		char* lc_TripletToAminoAcid_ = lb_BackwardsFrame ? 
			gk_GpfBase.mk_DnaTripletToAminoAcidReverse_ :
			gk_GpfBase.mk_DnaTripletToAminoAcid_;
		
		qint64 li_AssemblyGno = li_Gno;
		if (lb_RightHalfMass)
		{
			li_AssemblyGno += lb_BackwardsFrame ? 2 : -2;
			li_DnaOffset += lb_BackwardsFrame ? 2 : -2;
		}
		
		qint64 li_AssemblyLength = 3;
		
		if (lb_BackwardsFrame)
			li_DnaOffset -= 2;
		qint64 li_Step1 = lb_BackwardsFrame ? -1 : 1;
		if (lb_RightHalfMass)
			li_Step1 = -li_Step1;
		qint64 li_Step3 = li_Step1 * 3;
		
		// try to assemble an alignment
		
		qint64 li_AssemblyMass = mk_GpfIndexFile.mi_WaterMass;
		
		QString ls_Peptide;
		int li_TagAminoAcidsPassed = 0;

		while (li_AssemblyMass < mi_MaxMass)
		{
			quint16 lui_Triplet = gk_GpfBase.readNucleotideTriplet(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), li_DnaOffset);
			char lc_AminoAcid = lc_TripletToAminoAcid_[lui_Triplet];
			
			// cancel if invalid amino acid
			if (!gk_GpfBase.mb_IsAminoAcid_[(int)lc_AminoAcid])
				break;

			if (lb_RightHalfMass)
				ls_Peptide.prepend(QChar(lc_AminoAcid));
			else
				ls_Peptide.append(QChar(lc_AminoAcid));
			//printf("%c", lc_AminoAcid);
			
			if (li_AssemblyMass >= li_HalfMass)
				++li_TagAminoAcidsPassed;
			
			li_AssemblyMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)lc_AminoAcid];

			if (li_TagAminoAcidsPassed >= mk_GpfIndexFile.mi_TagSize)
			{
				// we have passed the half mass and the tag, now we can create alignments!

				// maybe we found an immediate hit already?
				if (li_AssemblyMass >= mi_MinMass && li_AssemblyMass <= mi_MaxMass)
				{
					printf("IMMEDIATE HIT: %s %d %d\n", ls_Peptide.toStdString().c_str(), (qint32)(li_AssemblyGno & ~mk_GpfIndexFile.mi_GnoBackwardsBit), (qint32)li_AssemblyLength);
				}
				
				// or maybe we find the start of an intron?
				for (int li_IntronStartOffset = 0; li_IntronStartOffset < 3; ++li_IntronStartOffset)
				{
					qint32 li_SplitTriplet = 0;
					if (li_IntronStartOffset > 0)
						li_SplitTriplet = readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), (li_DnaOffset + li_Step3) * 3, li_IntronStartOffset * 3);
					qint32 li_IntronStartNucleotides = readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), (li_DnaOffset + li_Step3 + li_Step1 * li_IntronStartOffset) * 3, mi_IntronStartMaxLength * 3);
					printf("%d %d\n", li_IntronStartOffset, li_IntronStartNucleotides);
					foreach (tk_IntPair lk_IntronStartSite, mk_IntronStart.keys())
					{
						qint32 li_ThisIntronStartNucleotides = li_IntronStartNucleotides & ((1 << (lk_IntronStartSite.second * 3)) - 1);
						if (li_ThisIntronStartNucleotides == lk_IntronStartSite.first)
						{
							// yay! we found an intron start site
 							printf("INTRON START: %s %d (%s%s)\n", ls_Peptide.toStdString().c_str(), li_IntronStartOffset, lb_BackwardsFrame ? "-" : "+", lb_RightHalfMass ? "R": "L");
							if (!lb_RightHalfMass)
							{
								qint64 li_IntronEndOffset = li_DnaOffset + li_Step3 + li_Step1 * (li_IntronStartOffset + lk_IntronStartSite.second);
								int li_IntronLength = lk_IntronStartSite.second;
								while (li_IntronLength < mi_MaxIntronLength)
								{
									// :TODO: check how many nucleotides can actually be read (scaffold border?)
									qint32 li_IntronEndNucleotides = readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), li_IntronEndOffset * 3, mi_IntronEndMaxLength * 3);
									foreach (tk_IntPair lk_IntronEndSite, mk_IntronStart[lk_IntronStartSite])
									{
										if (li_IntronEndNucleotides == lk_IntronEndSite.first)
										{
											QString ls_IntronSplitPeptide = ls_Peptide;
											qint64 li_IntronSplitAssemblyLength = li_AssemblyLength;
											qint64 li_IntronSplitAssemblyMass = li_AssemblyMass;
											
											qint32 li_AssembledSplitTriplet = li_SplitTriplet;
											li_AssembledSplitTriplet |= readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), (li_IntronEndOffset + lk_IntronEndSite.second * li_Step1) * 3, (3 - li_IntronStartOffset) * 3) << (li_IntronStartOffset * 3);
											qint64 li_IntronSplitDnaOffset = li_IntronEndOffset + (lk_IntronEndSite.second + (3 - li_IntronStartOffset)) * li_Step1;
											
											char lc_AminoAcid = lc_TripletToAminoAcid_[li_AssembledSplitTriplet];
											
											if (lc_AminoAcid == '$')
												break;
											
											// cancel if invalid amino acid
											if (!gk_GpfBase.mb_IsAminoAcid_[(int)lc_AminoAcid])
												break;

											if (lb_RightHalfMass)
												ls_IntronSplitPeptide.prepend(QChar(lc_AminoAcid));
											else
												ls_IntronSplitPeptide.append(QChar(lc_AminoAcid));
											
											li_IntronSplitAssemblyMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)lc_AminoAcid];
											li_IntronSplitAssemblyLength += 3;
											
											while (li_IntronSplitAssemblyMass < mi_MaxMass)
											{
												quint16 lui_Triplet = gk_GpfBase.readNucleotideTriplet(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), li_IntronSplitDnaOffset);
												char lc_AminoAcid = lc_TripletToAminoAcid_[lui_Triplet];
												
												// cancel if invalid amino acid
												if (!gk_GpfBase.mb_IsAminoAcid_[(int)lc_AminoAcid])
													break;

												if (lb_RightHalfMass)
													ls_IntronSplitPeptide.prepend(QChar(lc_AminoAcid));
												else
													ls_IntronSplitPeptide.append(QChar(lc_AminoAcid));
												
// 												printf("%s\n", ls_IntronSplitPeptide.toStdString().c_str());
												
												li_IntronSplitAssemblyMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)lc_AminoAcid];
												li_IntronSplitAssemblyLength += 3;
											
												// check whether we have found an intron split hit
												if (li_IntronSplitAssemblyMass >= mi_MinMass && li_IntronSplitAssemblyMass <= mi_MaxMass)
												{
													printf("INTRON SPLIT HIT: %s %d\n", ls_IntronSplitPeptide.toStdString().c_str(), (qint32)li_IntronSplitAssemblyLength);
												}
												li_IntronSplitDnaOffset += li_Step3;
												if (li_IntronSplitDnaOffset < li_ScaffoldStart || li_IntronSplitDnaOffset + 2 > li_ScaffoldEnd)
													break;
											}
										}
									}
									li_IntronEndOffset += li_Step1;
									// break if intron end has gone past scaffold end
									if (li_IntronEndOffset > li_ScaffoldEnd)
										break;
								}
/*								int li_IntronEndOffset = li_DnaOffset + li_Step3 + li_Step1 * (li_IntronStartOffset + lk_IntronStartSite.second);
								qint32 li_IntronEndNucleotides = readBitsFromBuffer(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), li_DnaOffset * 3, mi_IntronEndMaxLength * 3);
								printf("%d %d\n", li_IntronEndNucleotides & 7, (li_IntronEndNucleotides >> 3) & 7);*/
								// check every possible intron end consensus sequence
							}
						}
					}
				}
			}
			
			// advance and break if out of scaffold
			li_DnaOffset += li_Step3;
			li_AssemblyLength += 3;
			if (lb_RightHalfMass)
				li_AssemblyGno += li_Step3;
			if (li_DnaOffset < li_ScaffoldStart || li_DnaOffset + 2 > li_ScaffoldEnd)
				break;
		}
		//printf("\n");
	}
}
