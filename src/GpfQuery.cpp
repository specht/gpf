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
{
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
	mi_MinMass = mi_Mass;
	mi_MaxMass = mi_Mass;
		
	// extract all HMST
	QMultiMap<qint32, qint64> lk_AllLeftHmst;
	QMultiMap<qint32, qint64> lk_AllRightHmst;
	
	// extract all left HMST
	qint64 li_HalfMass = 0;
	for (int i = 0; i + mk_GpfIndexFile.mi_TagSize <= ms_Peptide.length() && li_HalfMass <= mk_GpfIndexFile.mi_MaxMass; ++i)
	{
		QString ls_Tag = ms_Peptide.mid(i, mk_GpfIndexFile.mi_TagSize);
		qint32 li_Tag = gk_GpfBase.aminoAcidPolymerCode(ls_Tag.toStdString().c_str(), mk_GpfIndexFile.mi_TagSize) * 2;

		//printf("%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
		lk_AllLeftHmst.insert(li_Tag, li_HalfMass);
		li_HalfMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)(ms_Peptide.at(i).toAscii())];
	}
	
	// extract all right HMST
	li_HalfMass = 0;
	for (int i = ms_Peptide.length() - mk_GpfIndexFile.mi_TagSize; i >= 0 && li_HalfMass <= mk_GpfIndexFile.mi_MaxMass; --i)
	{
		QString ls_Tag = ms_Peptide.mid(i, mk_GpfIndexFile.mi_TagSize);
		qint32 li_Tag = gk_GpfBase.aminoAcidPolymerCode(ls_Tag.toStdString().c_str(), mk_GpfIndexFile.mi_TagSize) * 2 + 1;
		
// 		printf("%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
		lk_AllRightHmst.insert(li_Tag, li_HalfMass);
		li_HalfMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)(ms_Peptide.at(i + mk_GpfIndexFile.mi_TagSize - 1).toAscii())];
	}
	
	// mass direction: false == left, true == right
	typedef QPair<qint64, bool> tk_GnoMassDirection;
	QSet<tk_GnoMassDirection> lk_GnoSet;
	
	bool lb_SimilaritySearch = false;
	// if similarity search is switched off, choose the HMST with the least occurences, because
	// one HMST would be enough => super duper extra speed ON!!
	// no, actually we must choose a smallest HMST for both the left and right half mass HMST
	// because we cannot know whether the input peptide is N- or C-terminal tryptic, or both
	if (!lb_SimilaritySearch)
	{
		qint64 li_MinimumEntryCount = -1;
		QPair<qint32, qint64> lk_BestHmst;
		for (QMultiMap<qint32, qint64>::const_iterator lk_Iter = lk_AllLeftHmst.begin(); lk_Iter != lk_AllLeftHmst.end(); ++lk_Iter)
		{
			qint32 li_TagDirectionIndex = lk_Iter.key();
			if (li_MinimumEntryCount < 0 || mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex] < li_MinimumEntryCount)
			{
				li_MinimumEntryCount = mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex];
				lk_BestHmst = QPair<qint32, qint64>(lk_Iter.key(), lk_Iter.value());
			}
		}
		lk_AllLeftHmst.clear();
		lk_AllLeftHmst.insert(lk_BestHmst.first, lk_BestHmst.second);
		
		// :UGLY: code repeat from above, this time for right HMST
		li_MinimumEntryCount = -1;
		for (QMultiMap<qint32, qint64>::const_iterator lk_Iter = lk_AllRightHmst.begin(); lk_Iter != lk_AllRightHmst.end(); ++lk_Iter)
		{
			qint32 li_TagDirectionIndex = lk_Iter.key();
			if (li_MinimumEntryCount < 0 || mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex] < li_MinimumEntryCount)
			{
				li_MinimumEntryCount = mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex];
				lk_BestHmst = QPair<qint32, qint64>(lk_Iter.key(), lk_Iter.value());
			}
		}
		lk_AllRightHmst.clear();
		lk_AllRightHmst.insert(lk_BestHmst.first, lk_BestHmst.second);
	}
	
	// merge left and right HMST
	QMultiMap<qint32, qint64> lk_AllHmst = lk_AllLeftHmst + lk_AllRightHmst;

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
		for (qint64 i = 0; i < mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]; ++i)
		{
			qint64 li_Mass = mk_GpfIndexFile.readIndexBits(mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits + i * mk_GpfIndexFile.mi_MassBits, mk_GpfIndexFile.mi_MassBits);
// 			printf("[%d] %d %d %d\n", (qint32)li_TagDirectionIndex, (qint32)li_Mass, (qint32)li_MinMass, (qint32)li_MaxMass);
			if (li_Mass >= li_MinMass && li_Mass <= li_MaxMass && li_Start < 0)
				li_Start = i;
			if (li_Mass >= li_MinMass && li_Mass <= li_MaxMass)
				++li_Count;
		}
		//printf("%d - %d\n", (qint32)li_Start, (qint32)(li_Count));
		
		// now we have the correct range, read the corresponding GNOs
		for (qint64 i = 0; i < li_Count; ++i)
		{
			qint64 li_Gno = mk_GpfIndexFile.readIndexBits(
				mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits +
				mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex] * mk_GpfIndexFile.mi_MassBits + 
				(i + li_Start) * mk_GpfIndexFile.mi_OffsetBits, mk_GpfIndexFile.mi_OffsetBits);
			lk_GnoSet.insert(tk_GnoMassDirection(li_Gno, li_TagDirectionIndex % 2 == 1));
		}
	}
	printf("distinct GNO count: %d\n", lk_GnoSet.size());
	
	foreach (tk_GnoMassDirection lk_GnoMassDirection, lk_GnoSet)
	{
		qint64 li_Gno = lk_GnoMassDirection.first;
		bool lb_RightHalfMass = lk_GnoMassDirection.second;
		bool lb_BackwardsFrame = ((li_Gno & mk_GpfIndexFile.mi_GnoBackwardsBit) != 0);
		qint64 li_DnaOffset = li_Gno & ~mk_GpfIndexFile.mi_GnoBackwardsBit;
		printf("reading GNO 0x%08x (%d): %d/%d\n", (qint32)li_Gno, lb_RightHalfMass, (qint32)li_DnaOffset, lb_BackwardsFrame);
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
		// adjust scaffold if neceassary
		if (li_DnaOffset >= mk_GpfIndexFile.mk_ScaffoldStart[li_FirstScaffold] + mk_GpfIndexFile.mk_ScaffoldLength[li_FirstScaffold])
			++li_FirstScaffold;
		
		qint64 li_ScaffoldStart = mk_GpfIndexFile.mk_ScaffoldStart[li_FirstScaffold];
		qint64 li_ScaffoldEnd = li_ScaffoldStart + mk_GpfIndexFile.mk_ScaffoldLength[li_FirstScaffold] - 1;
		if (li_DnaOffset < li_ScaffoldStart || li_DnaOffset > li_ScaffoldEnd)
			printf("WRONG scaffold borders: %d, %d (%d)\n", (qint32)li_ScaffoldStart, (qint32)li_ScaffoldEnd, (qint32)li_DnaOffset);
		char* lc_TripletToAminoAcid_ = lb_BackwardsFrame ? 
			gk_GpfBase.mk_DnaTripletToAminoAcidReverse_ :
			gk_GpfBase.mk_DnaTripletToAminoAcid_;
		qint64 li_DnaAnchor = li_DnaOffset;
		
		qint64 li_AssemblyGno = li_Gno;
		if (lb_RightHalfMass)
		{
			li_AssemblyGno += lb_BackwardsFrame ? 2 : -2;
			li_DnaOffset += lb_BackwardsFrame ? 2 : -2;
		}
		
		qint64 li_AssemblyLength = 3;
		
		if (lb_BackwardsFrame)
			li_DnaOffset -= 2;
		qint64 li_Step = lb_BackwardsFrame ? -3 : 3;
		if (lb_RightHalfMass)
			li_Step = -li_Step;
		
		// try to assemble an immediate hit
		
		qint64 li_AssemblyMass = mk_GpfIndexFile.mi_WaterMass;

		while (li_AssemblyMass < mi_MaxMass)
		{
			quint16 lui_Triplet = gk_GpfBase.readNucleotideTriplet(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), li_DnaOffset);
			char lc_AminoAcid = lc_TripletToAminoAcid_[lui_Triplet];
			
			// cancel if invalid amino acid
			if (!gk_GpfBase.mb_IsAminoAcid_[(int)lc_AminoAcid])
				break;
			
			printf("%c", lc_AminoAcid);
			
			li_AssemblyMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)lc_AminoAcid];
			
			if (li_AssemblyMass >= mi_MinMass && li_AssemblyMass <= mi_MaxMass)
			{
				printf("%d %d\n", (qint32)(li_AssemblyGno & ~mk_GpfIndexFile.mi_GnoBackwardsBit), (qint32)li_AssemblyLength);
			}
			
			// advance and break if out of scaffold
			li_DnaOffset += li_Step;
			li_AssemblyLength += 3;
			if (lb_RightHalfMass)
				li_AssemblyGno += li_Step;
			if (li_DnaOffset < li_ScaffoldStart || li_DnaOffset + 2 > li_ScaffoldEnd)
				break;
		}
		printf("\n");
	}
}
