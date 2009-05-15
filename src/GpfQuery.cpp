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
	// check whether this is a valid peptide
	for (int i = 0; i < ms_Peptide.length(); ++i)
	{
		if (!gk_GpfBase.mb_IsAminoAcid_[(int)(ms_Peptide.at(i).toAscii())])
		{
			printf("Error: Invalid amino acids in %s.\n", ms_Peptide.toStdString().c_str());
			return;
		}
	}
		
	// extract all HMST
	QMultiMap<qint32, qint64> lk_AllHmst;
	
	// extract all left HMST
	qint64 li_HalfMass = 0;
	for (int i = 0; i + mk_GpfIndexFile.mi_TagSize <= ms_Peptide.length() && li_HalfMass <= mk_GpfIndexFile.mi_MaxMass; ++i)
	{
		QString ls_Tag = ms_Peptide.mid(i, mk_GpfIndexFile.mi_TagSize);
		qint32 li_Tag = gk_GpfBase.aminoAcidPolymerCode(ls_Tag.toStdString().c_str(), mk_GpfIndexFile.mi_TagSize) * 2;

		//printf("%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
		lk_AllHmst.insert(li_Tag, li_HalfMass);
		li_HalfMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)(ms_Peptide.at(i).toAscii())];
	}
	
	// extract all right HMST
	li_HalfMass = 0;
	for (int i = ms_Peptide.length() - mk_GpfIndexFile.mi_TagSize; i >= 0 && li_HalfMass <= mk_GpfIndexFile.mi_MaxMass; --i)
	{
		QString ls_Tag = ms_Peptide.mid(i, mk_GpfIndexFile.mi_TagSize);
		qint32 li_Tag = gk_GpfBase.aminoAcidPolymerCode(ls_Tag.toStdString().c_str(), mk_GpfIndexFile.mi_TagSize) * 2 + 1;
		
		//printf("%s %d %d\n", ls_Tag.toStdString().c_str(), li_Tag, (qint32)li_HalfMass);
		lk_AllHmst.insert(li_Tag, li_HalfMass);
		li_HalfMass += mk_GpfIndexFile.mi_AminoAcidMasses_[(int)(ms_Peptide.at(i + mk_GpfIndexFile.mi_TagSize - 1).toAscii())];
	}
	
	qint64 li_GnoCount = 0;
	QSet<qint64> lk_GnoSet;

	// iterate over all HMST
	for (QMultiMap<qint32, qint64>::const_iterator lk_Iter = lk_AllHmst.begin(); lk_Iter != lk_AllHmst.end(); ++lk_Iter)
	{
		qint32 li_TagDirectionIndex = lk_Iter.key();
		
		qint64 li_HalfMass = lk_Iter.value();
		qint64 li_MinMass = li_HalfMass;
		qint64 li_MaxMass = li_HalfMass;
		
		// determine sub range in HMST list (via min and max masses)
		printf("tag/dir %8d: %d entries: ", li_TagDirectionIndex, (qint32)mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]);
		qint64 li_Start = -1;
		qint64 li_Count = 0;
		for (qint64 i = 0; i < mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex]; ++i)
		{
			qint64 li_Mass = mk_GpfIndexFile.readIndexBits(mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits + i * mk_GpfIndexFile.mi_MassBits, mk_GpfIndexFile.mi_MassBits);
			//printf("%d %d %d\n", (qint32)li_Mass, (qint32)li_MinMass, (qint32)li_MaxMass);
			if (li_Mass >= li_MinMass && li_Mass <= li_MaxMass && li_Start < 0)
				li_Start = i;
			if (li_Mass >= li_MinMass && li_Mass <= li_MaxMass)
				++li_Count;
		}
		printf("%d - %d\n", (qint32)li_Start, (qint32)(li_Start + li_Count - 1));
		
		// now we have the correct range, read the corresponding GNOs
		for (qint64 i = 0; i < li_Count; ++i)
		{
			qint64 li_Gno = mk_GpfIndexFile.readIndexBits(
				mk_GpfIndexFile.mk_HmstOffset[li_TagDirectionIndex] * mk_GpfIndexFile.mi_HmstBits +
				mk_GpfIndexFile.mk_HmstCount[li_TagDirectionIndex] * mk_GpfIndexFile.mi_MassBits + 
				(i + li_Start) * mk_GpfIndexFile.mi_OffsetBits, mk_GpfIndexFile.mi_OffsetBits);
			++li_GnoCount;
			lk_GnoSet.insert(li_Gno);
		}
	}
	printf("GNO count: %d\n", (qint32)li_GnoCount);
	printf("distinct GNO count: %d\n", lk_GnoSet.size());
	
	foreach (qint64 li_Gno, lk_GnoSet)
	{
		bool lb_BackwardsFrame = ((li_Gno & mk_GpfIndexFile.mi_GnoBackwardsBit) != 0);
		qint64 li_DnaOffset = li_Gno & ~mk_GpfIndexFile.mi_GnoBackwardsBit;
		printf("reading GNO 0x%08x: %d/%d\n", (qint32)li_Gno, (qint32)li_DnaOffset, lb_BackwardsFrame);
		if (lb_BackwardsFrame)
		{
			for (int i = 0; i < 50; ++i)
			{
				quint16 lui_Triplet = gk_GpfBase.readNucleotideTriplet(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), li_DnaOffset - i * 3 - 2);
				char lc_AminoAcid = gk_GpfBase.mk_DnaTripletToAminoAcidReverse_[lui_Triplet];
				printf("%c", lc_AminoAcid);
			}
		}
		else
		{
			for (int i = 0; i < 50; ++i)
			{
				quint16 lui_Triplet = gk_GpfBase.readNucleotideTriplet(mk_GpfIndexFile.muc_pDnaBuffer.get_Pointer(), li_DnaOffset + i * 3);
				char lc_AminoAcid = gk_GpfBase.mk_DnaTripletToAminoAcid_[lui_Triplet];
				printf("%c", lc_AminoAcid);
			}
		}
		printf("\n");
	}
}
