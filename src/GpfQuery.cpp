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

#include <QReadLocker>
#include <QWriteLocker>
#include "GpfQuery.h"
#include "IndexReader.h"
#include "IntronSplitHit.h"
#include "IntronSplitFixedHit.h"
#include "Sorter.h"
#include "StopWatch.h"


k_GpfQuery::k_GpfQuery(k_GpfBase& ak_GpfBase, QString as_Query, double ad_PrecursorMass, QString as_Label)
	: mk_GpfBase(ak_GpfBase)
	, mk_pIndexFileInfo(RefPtr<k_GpfIndexFileInfo>(NULL))
	, ms_Query(as_Query.trimmed())
	, ms_CollapsedQuery(ak_GpfBase.CollapsePeptide(as_Query))
	, ms_Result(QString())
	, mb_IsFinished(false)
	, md_PrecursorMass(ad_PrecursorMass)
	, md_LeftGapMass(0.0)
	, md_RightGapMass(0.0)
	, ms_Label(as_Label)
{
	#include "generated/GpfParametersInitialize.inc.cpp"

	// set left and right gap masses, if available
	if (ms_Query.startsWith("[") && ms_Query.endsWith("]"))
	{
		ms_Query = ms_Query.mid(1, ms_Query.length() - 2);
		QStringList lk_Parts = ms_Query.split(",");
		if (lk_Parts.size() <= 3)
		{
			bool lb_Ok = false;
			double ld_Value = lk_Parts.first().toDouble(&lb_Ok);
			if (lb_Ok)
			{
				md_LeftGapMass = ld_Value;
				lk_Parts.removeFirst();
			}
			ms_Query = lk_Parts.takeFirst().trimmed();
			if (!lk_Parts.empty())
			{
				lb_Ok = false;
				double ld_Value = lk_Parts.first().toDouble(&lb_Ok);
				if (lb_Ok)
				{
					md_RightGapMass = ld_Value;
					lk_Parts.removeFirst();
				}
			}
		}
	}
	
	// set expected mass of the result peptide
	if (ad_PrecursorMass != 0.0)
		mui_QueryMass = (unsigned int)(ad_PrecursorMass * gui_MassPrecision);
	else
		mui_QueryMass = mk_GpfBase.CalculatePeptideMass(ms_Query) + (unsigned int)((md_LeftGapMass + md_RightGapMass) * gui_MassPrecision);
}


k_GpfQuery::~k_GpfQuery()
{
}


void k_GpfQuery::Execute()
{
	k_StopWatch lk_StopWatch;

	// set min and max result masses
	mi_MassTolerance = ((qint64)mui_QueryMass * GET_INT_PARAMETER(MassError) / 1000000);
	mui_ResultMassMaximum = mui_QueryMass + mi_MassTolerance;
	mui_ResultMassMinimum = mui_QueryMass - mi_MassTolerance;

	bool lb_Error = false;

	// if there was no genome specified in the parameters, pick the default genome
	if (mk_pIndexFileInfo.get_Pointer() == NULL)
		mk_pIndexFileInfo = mk_GpfBase.get_DefaultIndexFileInfo();

	mk_IndexFile.setFileName(mk_pIndexFileInfo->get_Filename());
	mk_IndexFile.open(QFile::ReadOnly | QFile::Unbuffered);

	// check peptide sanity
	bool lb_FoundInvalidAminoAcids = false;
	for (int i = 0; i < ms_Query.length(); ++i) {
		if (mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(i).toAscii()] >= 20)
			lb_FoundInvalidAminoAcids = true;
	}

	if (lb_FoundInvalidAminoAcids)
	{
		// there was an error!
		lb_Error = true;
	}

	mk_Hits.clear();

	if (!lb_Error)
	{
        mk_Hits = mk_Hits.unite(this->QueryImmediateHit());

		if (GET_INT_PARAMETER(SearchIntrons) == r_YesNoOption::Yes)
		{
			//this->QueryIntronHit();
			mk_Hits = mk_Hits.unite(this->QueryCsIntronHit());
		}
	}

	// finish all hits
	tk_ResultList::iterator lk_Iter = mk_Hits.begin();
	for (; lk_Iter != mk_Hits.end(); ++lk_Iter)
		lk_Iter->get_Pointer()->Finish();

	/*
	// discard incomplete hits if appropriate
	if (lb_AnyFullScoreResults && GET_INT_PARAMETER(CropCompleteHits) == r_YesNoOption::Yes)
	{
		lk_Iter = lk_ImmediateHits.begin();
		for (; lk_Iter != lk_ImmediateHits.end(); ++lk_Iter)
			if (!lk_Iter->get_Pointer()->get_HasFullScore())
				lk_Iter->get_Pointer()->MarkDiscarded();
	}
	*/

	// erase all discarded hits
	lk_Iter = mk_Hits.begin();
	while (lk_Iter != mk_Hits.end())
	{
		tk_ResultList::iterator lk_NextIter = lk_Iter + 1;
		if (lk_Iter->get_Pointer()->get_IsDiscarded())
			mk_Hits.erase(lk_Iter);

		lk_Iter = lk_NextIter;
	}

	// write results
	foreach (QString ls_Key, mk_Hits.uniqueKeys())
	{
		RefPtr<k_Hit> lk_pHit = mk_Hits[ls_Key];
		ms_Result += QString("- %1\n").arg(lk_pHit->description());
	}

	// mark finished
	{
		QWriteLocker lk_Locker(&mk_IsFinishedLock);
		mb_IsFinished = true;
	}

	mk_IndexFile.close();

	ms_QueryProcessingTime = lk_StopWatch.getTimeAsString();

	emit finished();
}


QString k_GpfQuery::get_Query() const
{
	return ms_Query;
}


QString k_GpfQuery::get_CollapsedQuery() const
{
	return ms_CollapsedQuery;
}


double k_GpfQuery::get_PrecursorMass() const
{
	return md_PrecursorMass;
}


QString k_GpfQuery::get_Label() const
{
	return ms_Label;
}


QString k_GpfQuery::get_Result() const
{
	return ms_Result;
}


QString k_GpfQuery::get_Parameters()
{
	QString ls_Parameters;
	ls_Parameters = "parameters:\n";
	QMap<r_GpfParameterName::Enumeration, k_GpfParameter*>::iterator lk_ParamIter = mk_Parameters.begin();
	for (; lk_ParamIter != mk_Parameters.end(); ++lk_ParamIter)
	{
		k_GpfParameter* lk_Parameter_ = lk_ParamIter.value();
		ls_Parameters += "  " + lk_Parameter_->getId() + ": '" + lk_Parameter_->getValueAsString() + "'\n";
	}
	return ls_Parameters;
}


QString k_GpfQuery::get_Info()
{
	QString ls_Info;

	ls_Info += "info:\n";
	ls_Info += QString("  query: '%1'\n").arg(ms_Query);
	ls_Info += QString("  duration: %1\n").arg(ms_QueryProcessingTime);
	ls_Info += QString("  precursorMassUsed: %1\n").arg((double)mui_QueryMass / gui_MassPrecision);

	QString ls_Parameters = "  " + this->get_Parameters().replace("\n", "\n  ");
	ls_Info += ls_Parameters;


	return ls_Info;
}


bool k_GpfQuery::get_IsFinished() const
{
	bool lb_IsFinished;
	{
		QReadLocker lk_Locker(&mk_IsFinishedLock);
		lb_IsFinished = mb_IsFinished;
	}
	return lb_IsFinished;
}


k_GpfBase& k_GpfQuery::get_GpfBase()
{
	return mk_GpfBase;
}


QFile& k_GpfQuery::get_IndexFile()
{
	return mk_IndexFile;
}


k_GpfIndexFileInfo* k_GpfQuery::get_GpfIndexFileInfo()
{
	return mk_pIndexFileInfo.get_Pointer();
}


tk_ResultList& k_GpfQuery::resultList()
{
    return mk_Hits;
}


void k_GpfQuery::SetParameters(QMap<QString,QString> ak_Parameters)
{
	for (int i = 0; i < mk_Parameters.size(); ++i)
	{
		r_GpfParameterName::Enumeration le_Parameter = mk_Parameters.keys()[i];
		k_GpfParameter* lk_Parameter_ = mk_Parameters[le_Parameter];
		mk_Parameters[le_Parameter]->reset();
		if (ak_Parameters.contains(lk_Parameter_->getId()))
			mk_Parameters[le_Parameter]->setValue(ak_Parameters[lk_Parameter_->getId()]);
	}

	// apply genome selection
	if (ak_Parameters.contains("genome"))
		mk_pIndexFileInfo = mk_GpfBase.get_IndexFileInfo(ak_Parameters["genome"]);
}


tk_ResultList k_GpfQuery::QueryImmediateHit()
{
	// Immediate Query: Find places where a sequence tag consisting of 
	// [left mass, amino acid trimer, right mass] matches within 
	// a single reading frame of a single scaffold.

	// we need at least three consecutive amino acids for the immediate search
	tk_ResultList lk_Results;

	if (ms_Query.length() < 3)
		return lk_Results;

	QString ls_CollapsedQuery = mk_GpfBase.CollapsePeptide(ms_Query);

	unsigned int lui_LeftMass = (unsigned int)(md_LeftGapMass * gui_MassPrecision);
	unsigned int lui_RightMass = (unsigned int)(md_RightGapMass * gui_MassPrecision);

	for (int i = 3; i < ms_Query.length(); ++i)
	{
		char lc_AminoAcid = ms_Query.at(i).toAscii();
		r_AminoAcid::Enumeration le_AminoAcid = mk_GpfBase.me_CharToAminoAcid_[(int)lc_AminoAcid];
		lui_RightMass += mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)le_AminoAcid];
	}

	QHash<unsigned int, RefPtr<k_Hit> > lk_ImmediateHits;

	unsigned int lui_StartIndex = 0;
	unsigned int lui_LastIndex = ms_Query.length() - 3;

	if (GET_INT_PARAMETER(SearchSimilar) == r_YesNoOption::No)
	{
		// if similarity search is turned off, only search for the center amino acid trimer
		lui_StartIndex = (ms_Query.length() - 3) / 2;
		lui_LastIndex = lui_StartIndex;
	}

	// extract trimers from query
	for (unsigned int lui_AminoAcidIndex = 0; lui_AminoAcidIndex < (unsigned int)ms_Query.length() - 2; ++lui_AminoAcidIndex)
	{
		r_AminoAcid::Enumeration le_First = mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_AminoAcidIndex).toAscii()];
		r_AminoAcid::Enumeration le_Second = mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_AminoAcidIndex + 1).toAscii()];
		r_AminoAcid::Enumeration le_Third = mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_AminoAcidIndex + 2).toAscii()];
		if (lui_AminoAcidIndex >= lui_StartIndex && lui_AminoAcidIndex <= lui_LastIndex)
		{
			// We have a valid trimer, query it.
			// Fix the trimer to accomodate for similar masses
			unsigned int lui_TrimerIndex = 
				(unsigned int)mk_GpfBase.me_AminoAcidMassSimilarityCollapse[le_First] * 400 + 
				(unsigned int)mk_GpfBase.me_AminoAcidMassSimilarityCollapse[le_Second] * 20 + 
				(unsigned int)mk_GpfBase.me_AminoAcidMassSimilarityCollapse[le_Third];

			// query a sequence tag: lui_TrimerIndex, lui_LeftMass, lui_RightMass
			unsigned int lui_EntryCount = mk_pIndexFileInfo->mui_pTrimerCountTable.get_Pointer()[lui_TrimerIndex];

			if (lui_EntryCount > 0)
			{
				// Filter out entries that lie within a certain mass range.
				// The left and right masses lists are sorted by mass, that means that we
				// can determine start and end of the in-range sublist in log(n) time.

				unsigned int lui_Start = mk_pIndexFileInfo->mui_pTrimerOffsetTable.get_Pointer()[lui_TrimerIndex];
				unsigned int lui_LeftStart = lui_Start;
				unsigned int lui_LeftCount = lui_EntryCount;
				unsigned int lui_RightStart = lui_Start;
				unsigned int lui_RightCount = lui_EntryCount;

				CropList(mk_pIndexFileInfo->mi_LeftTagMassesFilePosition, lui_Start, lui_EntryCount, (int)lui_LeftMass - mi_MassTolerance,
					(int)lui_LeftMass + mi_MassTolerance, lui_LeftStart, lui_LeftCount);
				CropList(mk_pIndexFileInfo->mi_RightTagMassesFilePosition, lui_Start, lui_EntryCount, (int)lui_RightMass - mi_MassTolerance,
					(int)lui_RightMass + mi_MassTolerance, lui_RightStart, lui_RightCount);

				RefPtr<unsigned int> lui_pLeftEntries = RefPtr<unsigned int>(new unsigned int[lui_LeftCount]);
				RefPtr<unsigned int> lui_pRightEntries = RefPtr<unsigned int>(new unsigned int[lui_RightCount]);
				mk_IndexFile.seek(mk_pIndexFileInfo->mi_LeftTagGnoFilePosition + lui_LeftStart * sizeof(unsigned int));
				mk_IndexFile.read((char*)lui_pLeftEntries.get_Pointer(), sizeof(unsigned int) * lui_LeftCount);
				mk_IndexFile.seek(mk_pIndexFileInfo->mi_RightTagGnoFilePosition + lui_RightStart * sizeof(unsigned int));
				mk_IndexFile.read((char*)lui_pRightEntries.get_Pointer(), sizeof(unsigned int) * lui_RightCount);
				
				// sort entries by file position
				if (lui_LeftCount > 0)
					SortUnsignedIntGno(lui_pLeftEntries.get_Pointer(), 0, lui_LeftCount - 1);
				if (lui_RightCount > 0)
					SortUnsignedIntGno(lui_pRightEntries.get_Pointer(), 0, lui_RightCount - 1);

				// check if we can find an immediate hit
				unsigned int lui_LeftPosition = 0;
				unsigned int lui_RightPosition = 0;

				// advance pointers in both the left and right sorted mass lists and check for equal GNO
				while (lui_LeftPosition < lui_LeftCount && lui_RightPosition < lui_RightCount)
				{
					if (lui_pLeftEntries.get_Pointer()[lui_LeftPosition] == lui_pRightEntries.get_Pointer()[lui_RightPosition])
					{
						unsigned int lui_Gno = lui_pLeftEntries.get_Pointer()[lui_LeftPosition];
						bool lb_Forward = (lui_Gno & gui_GlobalNucleotideOffsetBackwardFlag) == 0;
						r_Scaffold lr_Scaffold = mk_pIndexFileInfo->mk_Scaffolds[mk_pIndexFileInfo->get_ScaffoldIndexForGno(lui_Gno)];
						unsigned int lui_ReadingFrame = lr_Scaffold.getReadingFrameForGlobalNucleotideOffset(lui_Gno);

						unsigned int lui_PeptideStartGno, lui_PeptideLength;
						mk_pIndexFileInfo->FindPeptide(lui_PeptideStartGno, lui_PeptideLength, lui_Gno, lui_ReadingFrame);

						// insert a new immediate hit, if it hasn't been found already
						if (!lk_ImmediateHits.contains(lui_PeptideStartGno))
						{
							tk_Assembly lk_Assembly = tk_Assembly()
								<< tk_AssemblyPart(lui_PeptideStartGno & ~gui_GlobalNucleotideOffsetBackwardFlag, lui_PeptideLength * 3);
							RefPtr<k_Hit> lk_pHit(new k_Hit(*this, lb_Forward, lk_Assembly));
							lk_ImmediateHits.insert(lui_PeptideStartGno, lk_pHit);
						}
					}

					if (lt_Gno(&lui_pLeftEntries.get_Pointer()[lui_LeftPosition], &lui_pRightEntries.get_Pointer()[lui_RightPosition]))
						++lui_LeftPosition;
					else
						++lui_RightPosition;
				}
			}
		}
		lui_LeftMass += mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)le_First];
		if (lui_AminoAcidIndex + 3 < (unsigned int)ms_Query.length())
			lui_RightMass -= mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_AminoAcidIndex + 3).toAscii()]];
	}

	QList<RefPtr<k_Hit> > lk_Hits = lk_ImmediateHits.values();
	foreach (RefPtr<k_Hit> lk_pHit, lk_Hits)
		lk_Results.insert(lk_pHit->get_Assembly(), lk_pHit);

	return lk_Results;
}

/*
unsigned int k_GpfQuery::QueryIntronHit()
{
	unsigned int lui_ResultCount = 0;
	// Intron Query: Find pairs of places where the left half sequence tag
	// [left mass, amino acid trimer] and the right half sequence tag
	// [another amino acid trimer, right mass] match and are not more than
	// a given number of nucleotides plus the query gap apart. The intron 
	// is in-between. The pair must be located within a single scaffold but 
	// can cross reading frame boundaries.

	// we need at least 6 amino acids for the intron search
	if (ms_Query.length() < 6)
		return 0;

	// determine number of left and right half sequence tags each
	unsigned int lui_ListSize = ms_Query.length() - 5;

	// Semi sequence tag lists for both left and right semi tags (contain arrays of GNOs)
	RefPtr<unsigned int>* lui_pLeftTagList_ = new RefPtr<unsigned int>[lui_ListSize];
	RefPtr<unsigned int>* lui_pRightTagList_ = new RefPtr<unsigned int>[lui_ListSize];

	// These two guys contain the GNO array sizes for all left and right semi tags
	RefPtr<unsigned int> lui_pLeftTagListEntryCount(new unsigned int[lui_ListSize]);
	RefPtr<unsigned int> lui_pRightTagListEntryCount(new unsigned int[lui_ListSize]);

	// The next two structures contain information on how the tag lists are built up,
	// that is, what scaffold and which direction is contained with how many entries
	QList<QPair<unsigned int, unsigned int> >* lk_LeftTagListSpans_ = new QList<QPair<unsigned int, unsigned int> >[lui_ListSize];
	QList<QPair<unsigned int, unsigned int> >* lk_RightTagListSpans_ = new QList<QPair<unsigned int, unsigned int> >[lui_ListSize];

	for (unsigned int i = 0; i < lui_ListSize; ++i)
	{
		lui_pLeftTagListEntryCount.get_Pointer()[i] = 0;
		lui_pRightTagListEntryCount.get_Pointer()[i] = 0;
	}

	unsigned int lui_LeftMass = 0;
	unsigned int lui_RightMass = 0;

	for (unsigned int lui_Index = 3; lui_Index < (unsigned int)ms_Query.length(); ++lui_Index)
		lui_RightMass += mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_Index).toAscii()]];

	int li_MassError = mui_QueryMass * GET_INT_PARAMETER(MassError) / 1000000;

	// build left and right tag lists

	for (int lui_Index = 0; lui_Index <= ms_Query.length() - 3; ++lui_Index)
	{
		unsigned int lui_LeftTagListIndex = lui_Index;
		unsigned int lui_RightTagListIndex = 0;
		if (lui_Index >= 3)
			lui_RightTagListIndex = lui_Index - 3;

		r_AminoAcid::Enumeration le_First = mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_Index).toAscii()];
		r_AminoAcid::Enumeration le_Second = mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_Index + 1).toAscii()];
		r_AminoAcid::Enumeration le_Third = mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_Index + 2).toAscii()];

		if (le_First < 20 && le_Second < 20 && le_Third < 20)
		{
			// collapse the trimer and determine trimer index (0..7999)
			unsigned int lui_TrimerIndex = 
				(unsigned int)mk_GpfBase.me_AminoAcidMassSimilarityCollapse[le_First] * 400 + 
				(unsigned int)mk_GpfBase.me_AminoAcidMassSimilarityCollapse[le_Second] * 20 + 
				(unsigned int)mk_GpfBase.me_AminoAcidMassSimilarityCollapse[le_Third];

			unsigned int lui_EntryCount = mk_pIndexFileInfo->mui_pTrimerCountTable.get_Pointer()[lui_TrimerIndex];

			if (lui_EntryCount > 0)
			{
				if (lui_Index <= ms_Query.length() - 6)
				{
					// skip this if this is not the leftmost half sequence tag and we are not considering similarity
					if (!(GET_INT_PARAMETER(SearchSimilar) == r_YesNoOption::No && lui_Index > 0))
					{
						// add to left tag list
						unsigned int lui_Start = mk_pIndexFileInfo->mui_pTrimerOffsetTable.get_Pointer()[lui_TrimerIndex];
						unsigned int lui_CroppedStart = lui_Start;
						unsigned int lui_CroppedCount = lui_EntryCount;

						CropList(mk_pIndexFileInfo->mi_LeftTagMassesFilePosition, lui_Start, lui_EntryCount, (int)lui_LeftMass - li_MassError,
							(int)lui_LeftMass + li_MassError, lui_CroppedStart, lui_CroppedCount);

						if (lui_CroppedCount > 0)
						{
							lui_pLeftTagListEntryCount.get_Pointer()[lui_LeftTagListIndex] = lui_CroppedCount;

							lui_pLeftTagList_[lui_LeftTagListIndex] = RefPtr<unsigned int>(new unsigned int[lui_CroppedCount]);
							mk_IndexFile.seek(mk_pIndexFileInfo->mi_LeftTagGnoFilePosition + lui_CroppedStart * sizeof(unsigned int));
							mk_IndexFile.read((char*)lui_pLeftTagList_[lui_LeftTagListIndex].get_Pointer(), 
								sizeof(unsigned int) * lui_CroppedCount);

							// sort by plain DNA offset (regardless of direction)
							SortUnsignedIntGno(lui_pLeftTagList_[lui_LeftTagListIndex].get_Pointer(), 0, lui_CroppedCount - 1);

							lk_LeftTagListSpans_[lui_LeftTagListIndex] = mk_pIndexFileInfo->CreateSpansOfEqualScaffoldAndDirection(lui_pLeftTagList_[lui_LeftTagListIndex].get_Pointer(), lui_CroppedCount);
						}
					}
				}
				if (lui_Index >= 3)
				{
					// skip this if this is not the rightmost half sequence tag and we are not considering similarity
					if (!(GET_INT_PARAMETER(SearchSimilar) == r_YesNoOption::No && lui_Index < ms_Query.length() - 3))
					{
						// add to right tag list
						unsigned int lui_Start = mk_pIndexFileInfo->mui_pTrimerOffsetTable.get_Pointer()[lui_TrimerIndex];
						unsigned int lui_CroppedStart = lui_Start;
						unsigned int lui_CroppedCount = lui_EntryCount;

						CropList(mk_pIndexFileInfo->mi_RightTagMassesFilePosition, lui_Start, lui_EntryCount, (int)lui_RightMass - li_MassError,
							(int)lui_RightMass + li_MassError, lui_CroppedStart, lui_CroppedCount);

						if (lui_CroppedCount > 0)
						{
							lui_pRightTagListEntryCount.get_Pointer()[lui_RightTagListIndex] = lui_CroppedCount;

							lui_pRightTagList_[lui_RightTagListIndex] = RefPtr<unsigned int>(new unsigned int[lui_CroppedCount]);
							mk_IndexFile.seek(mk_pIndexFileInfo->mi_RightTagGnoFilePosition + lui_CroppedStart * sizeof(unsigned int));
							mk_IndexFile.read((char*)lui_pRightTagList_[lui_RightTagListIndex].get_Pointer(), 
								sizeof(unsigned int) * lui_CroppedCount);
						
							// sort by plain DNA offset (regardless of direction)
							SortUnsignedIntGno(lui_pRightTagList_[lui_RightTagListIndex].get_Pointer(), 0, lui_CroppedCount - 1);

							lk_RightTagListSpans_[lui_RightTagListIndex] = mk_pIndexFileInfo->CreateSpansOfEqualScaffoldAndDirection(lui_pRightTagList_[lui_RightTagListIndex].get_Pointer(), lui_CroppedCount);
						}
					}
				}
			}
		}

		lui_LeftMass += mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_Index).toAscii()]];
		if (lui_Index + 3 < ms_Query.length())
			lui_RightMass -= mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(lui_Index + 3).toAscii()]];
	}

	// We now have the information where all the half sequence tags occur within the DNA,
	// in the next step, we will try all possible combinations of left and right half
	// sequence tags and see whether we can find a pair which looks promising.
	// If there is a gap in the query between two interesting half sequence tags,
	// we must check whether that gap can actually be fixed in the DNA, which means that
	// most of the intersting hits we will find here will be discarded later.

	QMap<QPair<unsigned int, unsigned int>, k_IntronSplitHit> lk_IntronSplitHits;

	// try all combinations of left and right half sequence tags
	for (unsigned int lui_LeftTagListIndex = 0; lui_LeftTagListIndex < lui_ListSize; ++lui_LeftTagListIndex)
	{
		unsigned int lui_QueryLengthLeft = lui_LeftTagListIndex + 3;
		for (unsigned int lui_RightTagListIndex = lui_LeftTagListIndex; lui_RightTagListIndex < lui_ListSize; ++lui_RightTagListIndex)
		{
			unsigned int lui_QueryLengthRight = ms_Query.length() - 3 - lui_RightTagListIndex;
			unsigned int lui_LeftTagListSize = lui_pLeftTagListEntryCount.get_Pointer()[lui_LeftTagListIndex];
			unsigned int lui_RightTagListSize = lui_pRightTagListEntryCount.get_Pointer()[lui_RightTagListIndex];

			// Skip this pair of half sequence tags if one of them has produced no results.
			if (lui_LeftTagListSize == 0 || lui_RightTagListSize == 0)
				continue;

			QList<QPair<unsigned int, unsigned int> > lk_LeftTagListSpans = lk_LeftTagListSpans_[lui_LeftTagListIndex];
			QList<QPair<unsigned int, unsigned int> > lk_RightTagListSpans = lk_RightTagListSpans_[lui_RightTagListIndex];

			QList<QPair<unsigned int, unsigned int> >::const_iterator lk_LeftSpansIter = lk_LeftTagListSpans.begin();
			QList<QPair<unsigned int, unsigned int> >::const_iterator lk_RightSpansIter = lk_RightTagListSpans.begin();

			// Advance iterators until one is at the end of its list
			while (lk_LeftSpansIter != lk_LeftTagListSpans.end() && lk_RightSpansIter != lk_RightTagListSpans.end())
			{
				if (lk_LeftSpansIter->first == lk_RightSpansIter->first)
				{
					// We now have two lists of half sequence tags hits on the genome which 
					// are on the same scaffold and encoded in the same direction. See if we 
					// can make a few good hit out of this!
					bool lb_ForwardReadingFrame = (lk_LeftSpansIter->first & 1) == 0;
					r_Scaffold lr_Scaffold = mk_pIndexFileInfo->mk_Scaffolds[lk_LeftSpansIter->first / 2];

					unsigned int* lui_LeftEntry_ = &lui_pLeftTagList_[lui_LeftTagListIndex].get_Pointer()[lk_LeftSpansIter->second];
					unsigned int* lui_RightEntry_ = &lui_pRightTagList_[lui_RightTagListIndex].get_Pointer()[lk_RightSpansIter->second];

					unsigned int lui_LeftEntryCount = lui_LeftTagListSize - lk_LeftSpansIter->second;
					QList<QPair<unsigned int, unsigned int> >::const_iterator lk_TempIter = lk_LeftSpansIter;
					if (++lk_TempIter != lk_LeftTagListSpans.end())
						lui_LeftEntryCount = lk_TempIter->second - lk_LeftSpansIter->second;

					unsigned int lui_RightEntryCount = lui_RightTagListSize - lk_RightSpansIter->second;
					lk_TempIter = lk_RightSpansIter;
					if (++lk_TempIter != lk_RightTagListSpans.end())
						lui_RightEntryCount = lk_TempIter->second - lk_RightSpansIter->second;

					unsigned int* lui_LeftEnd_ = lui_LeftEntry_ + lui_LeftEntryCount;
					unsigned int* lui_RightEnd_ = lui_RightEntry_ + lui_RightEntryCount;

					unsigned int* lui_RightStart_ = lui_RightEntry_;
					while (lui_RightStart_ != lui_RightEnd_ && lui_LeftEntry_ != lui_LeftEnd_)
					{
						// Adjust right start so that the right half sequence tag is really right to
						// the left half sequence tag. Also, the gap in between the both must not be
						// of zero length, because we would have found this hit using the immediate
						// search already. Such a hit would later be reported as an intron-split hit
						// but in fact would be an immediate hit
						if (lb_ForwardReadingFrame)
						{
							// 10 means three amino acids (9) plus at least one nucleotide (1)
							while (lui_RightStart_ != lui_RightEnd_ && *lui_RightStart_ < *lui_LeftEntry_ + 10)
								++lui_RightStart_;
						}
						else
						{
							// 10 means three amino acids (9) plus at least one nucleotide (1)
							while (lui_RightStart_ != lui_RightEnd_ && *lui_RightStart_ > *lui_LeftEntry_ - 10)
								++lui_RightStart_;
						}

						lui_RightEntry_ = lui_RightStart_;

						// If there are more right entries, pick all those with a distance below the specified threshold
						while (lui_RightEntry_ != lui_RightEnd_ && 
							(lb_ForwardReadingFrame? 
							(*lui_RightEntry_ < *lui_LeftEntry_ + 9 + GET_INT_PARAMETER(MaxIntronLength)):
							(*lui_RightEntry_ > *lui_LeftEntry_ - 9 - GET_INT_PARAMETER(MaxIntronLength))))
						{
							unsigned int lui_LeftGno = *lui_LeftEntry_;
							unsigned int lui_RightGno = *lui_RightEntry_;

							unsigned int lui_TempStart, lui_TempLength;
							mk_pIndexFileInfo->FindPeptide(lui_TempStart, lui_TempLength, lui_LeftGno, lr_Scaffold.getReadingFrameForGlobalNucleotideOffset(lui_LeftGno));
							unsigned int lui_PeptideStartGno = lui_TempStart;

							mk_pIndexFileInfo->FindPeptide(lui_TempStart, lui_TempLength, lui_RightGno, lr_Scaffold.getReadingFrameForGlobalNucleotideOffset(lui_RightGno));
							unsigned int lui_PeptideEndGno = lb_ForwardReadingFrame? 
								lui_TempStart + lui_TempLength * 3 - 1:
								lui_TempStart - lui_TempLength * 3 + 1;

							QPair<unsigned int, unsigned int> lk_StartAndEndGno = 
								QPair<unsigned int, unsigned int>(lui_PeptideStartGno, lui_PeptideEndGno);

							unsigned int lui_LeftHalfLength = abs((qint64)lui_LeftGno - (qint64)lui_PeptideStartGno) / 3 + 3;
							unsigned int lui_RightHalfLength = (abs((qint64)lui_RightGno - (qint64)lui_PeptideEndGno) + 1) / 3;

							if (!lk_IntronSplitHits.contains(lk_StartAndEndGno))
							{
								k_IntronSplitHit lk_IntronSplitHit(this, lr_Scaffold, ms_Query, 
									lui_PeptideStartGno, lui_LeftHalfLength, 
									lui_PeptideEndGno, lui_RightHalfLength, 
									lui_QueryLengthLeft, lui_QueryLengthRight);

								lk_IntronSplitHits.insert(lk_StartAndEndGno, lk_IntronSplitHit);
							}

							if (lk_IntronSplitHits.contains(lk_StartAndEndGno))
							{
								lk_IntronSplitHits[lk_StartAndEndGno].MarkChainElements(lui_LeftTagListIndex, 3, r_ChainMarker::Left);
								lk_IntronSplitHits[lk_StartAndEndGno].MarkChainElements(lui_RightTagListIndex + 3, 3, r_ChainMarker::Right);
								lk_IntronSplitHits[lk_StartAndEndGno].Update(lui_LeftHalfLength, lui_RightHalfLength, lui_QueryLengthLeft, lui_QueryLengthRight);
							}

							++lui_RightEntry_;
						}

						++lui_LeftEntry_;
					}
				}

				// advance one iterator
				if (lk_LeftSpansIter->first < lk_RightSpansIter->first)
					++lk_LeftSpansIter;
				else
					++lk_RightSpansIter;
			}
		}
	}

	// finish hits
	QList<k_IntronSplitFixedHit> lk_FixedHitList;

	// try to finish all unfixed intron split hits and collect fixed results in the fixed hit list
	QMap<QPair<unsigned int, unsigned int>, k_IntronSplitHit>::iterator lk_Iter = lk_IntronSplitHits.begin();
	for (; lk_Iter != lk_IntronSplitHits.end(); ++lk_Iter)
	{
		lk_Iter->Finish();
		lk_FixedHitList += lk_Iter->get_FixedHitList();
	}

	QList<k_IntronSplitFixedHit>::iterator lk_ListIter;
	bool lb_AnyFullScoreHits = false;

	lk_ListIter = lk_FixedHitList.begin();
	for (; lk_ListIter != lk_FixedHitList.end(); ++lk_ListIter)
		if (lk_ListIter->get_HasFullScore())
			lb_AnyFullScoreHits = true;

	// discard incomplete hits if appropriate
	if (lb_AnyFullScoreHits && GET_INT_PARAMETER(CropCompleteHits) == r_YesNoOption::Yes)
	{
		lk_ListIter = lk_FixedHitList.begin();
		for (; lk_ListIter != lk_FixedHitList.end(); ++lk_ListIter)
			if (!lk_ListIter->get_HasFullScore())
				lk_ListIter->MarkDiscarded();
	}

	// erase all discarded hits in two steps: 
	// 1. determine all list positions that are to be deleted (backwards)
	// 2. delete these items

	QList <int> lk_ClearList;
	for (int i = lk_FixedHitList.size() - 1; i >= 0; --i)
		if (lk_FixedHitList[i].get_IsDiscarded())
			lk_ClearList << i;

	foreach (int i, lk_ClearList)
		lk_FixedHitList.removeAt(i);

	// print results

	foreach (k_IntronSplitFixedHit lk_FixedHit, lk_FixedHitList)
		ms_Result += lk_FixedHit.description();

	delete [] lk_LeftTagListSpans_;
	delete [] lk_RightTagListSpans_;

	delete [] lui_pLeftTagList_;
	delete [] lui_pRightTagList_;

	return lui_ResultCount;
}

*/

tk_ResultList k_GpfQuery::QueryCsIntronHit()
{
	tk_ResultList lk_Results;
	// Consensus Sequence Intron Query: Process left and right halves separately.
	// Find only one half sequence tag and try to construct the peptide via
	// consensus sequences GT-AG and GC-AG. Only one correct amino acid pentamer
	// plus at least one arbitrary amino acid are required for this type of search.

	// we need at least 6 amino acids (one pentamer plus 1) for the CS intron search
	if (ms_Query.length() < 6)
		return lk_Results;

	tk_PeptideLocations lk_Locations;
	// find left pentamers
	lk_Locations = this->FindPentamers(true);
	lk_Results = lk_Results.unite(this->AssembleHalfHits(lk_Locations, true));

	// find right pentamers
	lk_Locations = this->FindPentamers(false);
	lk_Results = lk_Results.unite(this->AssembleHalfHits(lk_Locations, false));

	return lk_Results;
}


tk_PeptideLocations k_GpfQuery::FindPentamers(bool ab_Left)
{
	int li_ListSize = ms_Query.length() - 2;

	RefPtr<unsigned int>* lui_pTagList_ = new RefPtr<unsigned int>[li_ListSize];
	RefPtr<unsigned int> lui_pTagListEntryCount(new unsigned int[li_ListSize]);

	int li_Mass = ab_Left ? 
		(unsigned int)(md_LeftGapMass * gui_MassPrecision) : 
		(unsigned int)(md_RightGapMass * gui_MassPrecision);

	// look up all half mass sequence tags (three amino acids)
	// fill lui_pTagList and lui_pTagListEntryCount with the GNOs for each trimer
	int li_IndexStart = ab_Left? 0: li_ListSize - 1;
	int li_IndexEnd = ab_Left? li_ListSize: -1;
	int li_IndexStep = ab_Left? 1: -1;
	qint64 li_TagMassesFilePosition = ab_Left? mk_pIndexFileInfo->mi_LeftTagMassesFilePosition: mk_pIndexFileInfo->mi_RightTagMassesFilePosition;
	qint64 li_TagGnoFilePosition = ab_Left? mk_pIndexFileInfo->mi_LeftTagGnoFilePosition: mk_pIndexFileInfo->mi_RightTagGnoFilePosition;
	for (int li_Index = li_IndexStart; li_Index != li_IndexEnd; li_Index += li_IndexStep)
	{
		r_AminoAcid::Enumeration le_First = mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(li_Index).toAscii()];
		r_AminoAcid::Enumeration le_Second = mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(li_Index + 1).toAscii()];
		r_AminoAcid::Enumeration le_Third = mk_GpfBase.me_CharToAminoAcid_[(int)ms_Query.at(li_Index + 2).toAscii()];

		// collapse the trimer and determine trimer index (0..7999)
		unsigned int lui_TrimerIndex = 
			(unsigned int)mk_GpfBase.me_AminoAcidMassSimilarityCollapse[le_First] * 400 + 
			(unsigned int)mk_GpfBase.me_AminoAcidMassSimilarityCollapse[le_Second] * 20 + 
			(unsigned int)mk_GpfBase.me_AminoAcidMassSimilarityCollapse[le_Third];

		unsigned int lui_EntryCount = mk_pIndexFileInfo->mui_pTrimerCountTable.get_Pointer()[lui_TrimerIndex];

		unsigned int lui_Start = mk_pIndexFileInfo->mui_pTrimerOffsetTable.get_Pointer()[lui_TrimerIndex];
		unsigned int lui_CroppedStart = lui_Start;
		unsigned int lui_CroppedCount = lui_EntryCount;
		CropList(li_TagMassesFilePosition, lui_Start, lui_EntryCount, li_Mass - mi_MassTolerance,
			li_Mass + mi_MassTolerance, lui_CroppedStart, lui_CroppedCount);

		lui_pTagListEntryCount.get_Pointer()[li_Index] = lui_CroppedCount;
		if (lui_CroppedCount > 0)
		{
			lui_pTagList_[li_Index] = RefPtr<unsigned int>(new unsigned int[lui_CroppedCount]);
	
			mk_IndexFile.seek(li_TagGnoFilePosition + lui_CroppedStart * sizeof(unsigned int));
			mk_IndexFile.read((char*)lui_pTagList_[li_Index].get_Pointer(), sizeof(unsigned int) * lui_CroppedCount);
	
			SortPlainUnsignedIntGno(lui_pTagList_[li_Index].get_Pointer(), 0, lui_CroppedCount - 1);
		}

		li_Mass += mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)(ab_Left? le_First: le_Third)];
	}

	QList<unsigned int> lk_PentamerGnos;

	// find pentamers by searching for pairs of trimers that have a distance of two amino acids
	for (int li_Index = 0; li_Index < li_ListSize - 2; ++li_Index)
	{
		unsigned int* lui_FirstListItem_ = lui_pTagList_[li_Index].get_Pointer();
		unsigned int* lui_FirstListEnd_ = lui_FirstListItem_ + lui_pTagListEntryCount.get_Pointer()[li_Index];
		unsigned int* lui_SecondListItem_ = lui_pTagList_[li_Index + 2].get_Pointer();
		unsigned int* lui_SecondListEnd_ = lui_SecondListItem_ + lui_pTagListEntryCount.get_Pointer()[li_Index + 2];

		if (lui_FirstListItem_ != NULL && lui_SecondListItem_ != NULL)
		{
			while (lui_SecondListItem_ < lui_SecondListEnd_)
			{
				if ((*lui_SecondListItem_ & gui_GlobalNucleotideOffsetBackwardFlag) == 0)
				{
					while (lui_FirstListItem_ < lui_FirstListEnd_ && *lui_FirstListItem_ < *lui_SecondListItem_ - 6)
						++lui_FirstListItem_;
					if (lui_FirstListItem_ < lui_FirstListEnd_ && *lui_FirstListItem_ == *lui_SecondListItem_ - 6)
						lk_PentamerGnos.push_back(*lui_FirstListItem_);
				}
				else
				{
					while (lui_FirstListItem_ < lui_FirstListEnd_ && *lui_FirstListItem_ < *lui_SecondListItem_ + 6)
						++lui_FirstListItem_;
					if (lui_FirstListItem_ < lui_FirstListEnd_ && *lui_FirstListItem_ == *lui_SecondListItem_ + 6)
						lk_PentamerGnos.push_back(*lui_FirstListItem_);
				}
				++lui_SecondListItem_;
			}
		}
	}
	
	// now we all good pentamers in lk_PentamerGnos

	// lk_Peptides contains all peptides that represent a possible half of a hit
	tk_PeptideLocations lk_Peptides;
	foreach (unsigned int lui_Gno, lk_PentamerGnos)
	{
		r_Scaffold& lr_Scaffold = mk_pIndexFileInfo->mk_Scaffolds[mk_pIndexFileInfo->get_ScaffoldIndexForGno(lui_Gno)];
		unsigned int lui_PeptideStart, lui_PeptideLength;
		mk_pIndexFileInfo->FindPeptide(lui_PeptideStart, lui_PeptideLength, lui_Gno, lr_Scaffold.getReadingFrameForGlobalNucleotideOffset(lui_Gno));
		lk_Peptides.insert(tk_PeptideSpan(lui_PeptideStart, lui_PeptideLength));
	}

	delete [] lui_pTagList_;

	return lk_Peptides;
}


tk_ResultList k_GpfQuery::AssembleHalfHits(tk_PeptideLocations ak_Locations, bool ab_Left)
{
	tk_ResultList lk_Results;

	int li_MaxDnaSpanLength = GET_INT_PARAMETER(MaxIntronLength) + (mui_QueryMass / mk_GpfBase.mui_AminoAcidWeight_[r_AminoAcid::Gly] + 1) * 3 + 3;
	foreach (tk_PeptideSpan lk_Peptide, ak_Locations)
	{
		unsigned int lui_Gno = lk_Peptide.first;
		bool lb_Forward = (lui_Gno & gui_GlobalNucleotideOffsetBackwardFlag) == 0;
		unsigned int lui_PeptidePosition = lui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag;
		unsigned int lui_PeptideLength = lk_Peptide.second * 3;

		// read DNA span
		unsigned int lui_DnaSpanLength = 0;
		unsigned int lui_LeftBorder = ab_Left? 0: li_MaxDnaSpanLength - lui_PeptideLength;
		unsigned int lui_RightBorder = ab_Left? li_MaxDnaSpanLength - lui_PeptideLength: 0;
		RefPtr<unsigned char> luc_pNucleotides = mk_pIndexFileInfo->readDnaSpanAsUnsignedChar(lui_Gno, lui_PeptideLength, lui_LeftBorder, lui_RightBorder, lui_DnaSpanLength, &mk_IndexFile);

		// determine maximum length of the fixed part (according to the precursor mass) 
		// and record fragment masses of the fixed part
		QList<unsigned int> lk_FragmentMasses;
		int li_MaxFixedPartLength = 2;
		int li_Mass = mk_GpfBase.mui_WaterMass;
		lk_FragmentMasses.append(li_Mass);
		int li_IndexStart = ab_Left? 0: lui_DnaSpanLength - 3;
		int li_IndexEnd = ab_Left? lui_DnaSpanLength - 2: -1;
		int li_IndexStep = ab_Left? 3: -3;
		for (int li_Index = li_IndexStart; ab_Left? li_Index < li_IndexEnd: li_Index > li_IndexEnd; li_Index += li_IndexStep)
		{
			unsigned int lui_Triplet = 
				((luc_pNucleotides.get_Pointer()[li_Index] & 7) << 6) | 
				((luc_pNucleotides.get_Pointer()[li_Index + 1] & 7) << 3) | 
				(luc_pNucleotides.get_Pointer()[li_Index + 2] & 7);

			r_AminoAcid::Enumeration le_AminoAcid = (r_AminoAcid::Enumeration)(unsigned char)mk_GpfBase.mc_TripletToAminoAcidForward_[lui_Triplet];
			li_Mass += mk_GpfBase.mui_AminoAcidWeight_[le_AminoAcid];

			if (li_Mass > (int)mui_ResultMassMaximum)
				break;

			if (r_AminoAcid::isTrypticCleavageSiteOrStop(le_AminoAcid) || r_AminoAcid::isUnknown(le_AminoAcid))
			{
				if (ab_Left || ((!ab_Left) && (li_Index != li_IndexStart)))
					break;
			}

			lk_FragmentMasses.append((unsigned int)li_Mass);
			li_MaxFixedPartLength = abs(li_Index - li_IndexStart) + 5;
		}

		// determine potential intron start and end locations
		QList<int> lk_FixedPartCuts, lk_VariablePartCuts;
		li_IndexStart = ab_Left? 1: lui_DnaSpanLength - 3;
		li_IndexEnd = ab_Left? lui_DnaSpanLength - 1: -1;
		li_IndexStep = ab_Left? 1: -1;
		for (int li_Index = li_IndexStart; li_Index != li_IndexEnd; li_Index += li_IndexStep)
		{
			unsigned short* luw_Dinucleotide = (unsigned short*)(&luc_pNucleotides.get_Pointer()[li_Index]);
			if (ab_Left)
			{
				if (li_Index <= li_MaxFixedPartLength && (*luw_Dinucleotide == 0x0302 || *luw_Dinucleotide == 0x0102)) // GT or GC
					lk_FixedPartCuts.append(li_Index);
				if (*luw_Dinucleotide == 0x0200) // AG
					lk_VariablePartCuts.append(li_Index + 2);
			}
			else
			{
				if (*luw_Dinucleotide == 0x0302 || *luw_Dinucleotide == 0x0102) // GT or GC
					lk_VariablePartCuts.append(li_Index);
				if (abs(li_Index - li_IndexStart) + 1 <= li_MaxFixedPartLength &&  (*luw_Dinucleotide == 0x0200)) // AG
					lk_FixedPartCuts.append(li_Index + 2);
			}
		}

		// try all combinations of fixed and variable part cut locations and see whether a feasible peptide is assembled
		QList<int>::const_iterator lk_FixedIter = lk_FixedPartCuts.begin();
		for (; lk_FixedIter != lk_FixedPartCuts.end(); ++lk_FixedIter)
		{
			// remove variable part cut locations that have become useless
			while (!lk_VariablePartCuts.empty() && (ab_Left? (lk_VariablePartCuts.first() <= *lk_FixedIter): (lk_VariablePartCuts.first() >= *lk_FixedIter)))
				lk_VariablePartCuts.removeFirst();

			QList<int>::const_iterator lk_VariableIter = lk_VariablePartCuts.begin();
			for (; lk_VariableIter != lk_VariablePartCuts.end(); ++lk_VariableIter)
			{
				// we have found an intron with a valid dinucleotide pair
				int li_FixedPartSplit = *lk_FixedIter;
				int li_VariablePartSplit = *lk_VariableIter;
				int li_FixedPartLength = ab_Left? li_FixedPartSplit: lui_DnaSpanLength - li_FixedPartSplit;
				int li_VariablePartLength = 0;
				int li_MaxVariablePartLength = ab_Left? lui_DnaSpanLength - li_VariablePartSplit: li_VariablePartSplit;

				// we do not have to translate the fixed part because we 
				// already have the masses of every fixed part fragment
				int li_FixedPartAminoAcidCount = li_FixedPartLength / 3;

				// determine fixed part peptide mass
				int li_PeptideMass = lk_FragmentMasses[li_FixedPartAminoAcidCount];
				r_AminoAcid::Enumeration le_AminoAcid = r_AminoAcid::Met;

				// process optional first (split) amino acid
				if (li_FixedPartLength % 3 != 0)
				{
					// construct a new amino acid that is interrupted by an intron
					int li_Split = li_FixedPartLength % 3;
					unsigned int lui_Triplet = 0;

					// collect fixed and variable part nucleotides
					if (!ab_Left)
						li_Split = 3 - li_Split;
					li_VariablePartLength = 3 - li_Split;

					for (int li_SplitIndex = 0; li_SplitIndex < li_Split; ++li_SplitIndex)
						lui_Triplet = (lui_Triplet << 3) | (luc_pNucleotides.get_Pointer()[(ab_Left? li_FixedPartSplit: li_VariablePartSplit) - li_Split + li_SplitIndex] & 7);
					for (int li_SplitIndex = li_Split; li_SplitIndex < 3; ++li_SplitIndex)
						lui_Triplet = (lui_Triplet << 3) | (luc_pNucleotides.get_Pointer()[(ab_Left? li_VariablePartSplit: li_FixedPartSplit) - li_Split + li_SplitIndex] & 7);

					le_AminoAcid = (r_AminoAcid::Enumeration)(unsigned char)mk_GpfBase.mc_TripletToAminoAcidForward_[lui_Triplet];

					li_PeptideMass += mk_GpfBase.mui_AminoAcidWeight_[le_AminoAcid];
				}

				// process remaining amino acids
				while ((!r_AminoAcid::isTrypticCleavageSiteOrStop(le_AminoAcid)) && (!r_AminoAcid::isUnknown(le_AminoAcid)) &&  (li_VariablePartLength <= li_MaxVariablePartLength - 3))
				{
					unsigned int lui_Triplet = 0;
					for (int li_Index = 0; li_Index < 3; ++li_Index)
						lui_Triplet = (lui_Triplet << 3) | (luc_pNucleotides.get_Pointer()[ab_Left? (li_VariablePartSplit + li_VariablePartLength + li_Index): (li_VariablePartSplit - li_VariablePartLength - 3 + li_Index)] & 7);

					le_AminoAcid = (r_AminoAcid::Enumeration)(unsigned char)mk_GpfBase.mc_TripletToAminoAcidForward_[lui_Triplet];

					li_PeptideMass += mk_GpfBase.mui_AminoAcidWeight_[le_AminoAcid];
					li_VariablePartLength += 3;
				}

				// remove trailing stop if there is one
				if (r_AminoAcid::isStop(le_AminoAcid) || !ab_Left)
				{
					li_VariablePartLength -= 3;
					li_PeptideMass -= mk_GpfBase.mui_AminoAcidWeight_[le_AminoAcid];
				}

				// see if something useful came out
				if ((li_PeptideMass >= (int)mui_ResultMassMinimum && li_PeptideMass <= (int)mui_ResultMassMaximum) && r_AminoAcid::isTrypticCleavageSiteOrStop(le_AminoAcid) && li_VariablePartLength > 0)
				{
					tk_Assembly lk_Assembly;

					if (ab_Left)
					{
						lk_Assembly = tk_Assembly() 
							<< tk_AssemblyPart(lui_PeptidePosition, li_FixedPartLength)
							<< tk_AssemblyPart(lb_Forward? lui_PeptidePosition + li_VariablePartSplit: lui_PeptidePosition - li_VariablePartSplit, li_VariablePartLength);
					}
					else
					{
						if (lb_Forward)
							lk_Assembly = tk_Assembly() 
								<< tk_AssemblyPart(lui_PeptidePosition + lui_PeptideLength - lui_DnaSpanLength + li_VariablePartSplit - li_VariablePartLength, li_VariablePartLength)
								<< tk_AssemblyPart(lui_PeptidePosition + lui_PeptideLength - li_FixedPartLength, li_FixedPartLength);
						else
							lk_Assembly = tk_Assembly() 
								<< tk_AssemblyPart(lui_PeptidePosition - lui_PeptideLength + lui_DnaSpanLength - li_VariablePartSplit + li_VariablePartLength, li_VariablePartLength)
								<< tk_AssemblyPart(lui_PeptidePosition - lui_PeptideLength + li_FixedPartLength, li_FixedPartLength);
					}

					RefPtr<k_Hit> lk_pHit = RefPtr<k_Hit>(new k_Hit(*this, lb_Forward, lk_Assembly));
					lk_Results.insert(lk_pHit->get_Assembly(), lk_pHit);

					// add intron ends
					QString ls_IntronEnds;
					ls_IntronEnds += QChar(mk_GpfBase.mc_NucleotideIntToChar_[luc_pNucleotides.get_Pointer()[ab_Left? li_FixedPartSplit: li_VariablePartSplit]]);
					ls_IntronEnds += QChar(mk_GpfBase.mc_NucleotideIntToChar_[luc_pNucleotides.get_Pointer()[ab_Left? (li_FixedPartSplit + 1): (li_VariablePartSplit + 1)]]);
					ls_IntronEnds += QChar(mk_GpfBase.mc_NucleotideIntToChar_[luc_pNucleotides.get_Pointer()[ab_Left? (li_VariablePartSplit - 2): (li_FixedPartSplit - 2)]]);
					ls_IntronEnds += QChar(mk_GpfBase.mc_NucleotideIntToChar_[luc_pNucleotides.get_Pointer()[ab_Left? (li_VariablePartSplit - 1): (li_FixedPartSplit - 1)]]);
					lk_pHit->AddInformation("intronEnds", ls_IntronEnds);
				}
			}
		}
	}
	return lk_Results;
}


void k_GpfQuery::CropList(qint64 ai_FileOffset, unsigned int aui_Start, unsigned int aui_Count,
						  int ai_MinMass, int ai_MaxMass,
						  unsigned int& aui_CropStart, unsigned int& aui_CropCount)
{
	unsigned int lui_CropStart = aui_Start;
	unsigned int lui_CropCount = aui_Count;

	int li_Entry;

	// return immediately if min / max values of this trimer are not within range
	aui_CropStart = aui_Start;
	aui_CropCount = 0;

	int li_LeftBorderMass;
	READ_UINT32_CAST_INT32(li_LeftBorderMass, ai_FileOffset + aui_Start * sizeof(unsigned int));
	if (li_LeftBorderMass > ai_MaxMass)
		return;

	int li_RightBorderMass;
	READ_UINT32_CAST_INT32(li_RightBorderMass, ai_FileOffset + (aui_Start + aui_Count - 1) * sizeof(unsigned int));
	if (li_RightBorderMass < ai_MinMass)
		return;

	if (aui_Count == 0)
		return;

	/*
	int li_ResultStart, li_ResultEnd;
	int high, low, probe, probeValue;

	high = aui_Start + aui_Count;
	low = aui_Start - 1;
	while (high - low > 1)
	{
		probe =  low + (high - low) / 2;
		READ_UINT32_CAST_INT32(probeValue, ai_FileOffset + probe * sizeof(unsigned int));
		if (probeValue <= ai_MinMass)
			low = probe;
		else
			high = probe;
	}
	li_ResultStart = low;

	high = aui_Start + aui_Count;
	low = aui_Start - 1;
	while (high - low > 1)
	{
		probe =  low + (high - low) / 2;
		READ_UINT32_CAST_INT32(probeValue, ai_FileOffset + probe * sizeof(unsigned int));
		if (probeValue >= ai_MaxMass)
			high = probe;
		else
			low = probe;
	}
	li_ResultEnd = high;

	aui_CropStart = 0; aui_CropCount = 0;

	if (li_ResultEnd >= li_ResultStart)
	{
		aui_CropStart = li_ResultStart;
		aui_CropCount = li_ResultEnd - li_ResultStart + 1;
	}

	return;
	*/

	//printf("left border: %7.2f, right border: %7.2f\n", (double)li_LeftBorderMass / gui_MassPrecision, (double)li_RightBorderMass / gui_MassPrecision);

	unsigned int lui_FindLeft = aui_Start;
	unsigned int lui_FindRight = aui_Start + aui_Count - 1;

	if (ai_MinMass > li_LeftBorderMass)
	{
		// adjust left list endpoint
		forever
		{
			if (lui_FindRight - lui_FindLeft > 1)
			{
				unsigned int lui_Mid = lui_FindLeft + (lui_FindRight - lui_FindLeft) / 2;
				// read mid entry
				READ_UINT32_CAST_INT32(li_Entry, ai_FileOffset + lui_Mid * sizeof(unsigned int));
		
				if (li_Entry < ai_MinMass)
					lui_FindLeft = lui_Mid + 1;
				else
					lui_FindRight = lui_Mid;
			}
			else
			{
				// left and right pointer are now next to each other or even equal
				if (lui_FindRight != lui_FindLeft)
				{
					READ_UINT32_CAST_INT32(li_Entry, ai_FileOffset + lui_FindLeft * sizeof(unsigned int));

					if (li_Entry < ai_MinMass)
						++lui_FindLeft;
				}
				lui_CropStart = lui_FindLeft;
				lui_CropCount = aui_Count - (lui_CropStart - aui_Start);
				break;
			}
		} 
	}

	lui_FindRight = aui_Start + aui_Count - 1;

	if (ai_MaxMass < li_RightBorderMass)
	{
		forever
		{
			if (lui_FindRight - lui_FindLeft > 1)
			{
				unsigned int lui_Mid = lui_FindLeft + (lui_FindRight - lui_FindLeft) / 2;
				// read mid entry
				READ_UINT32_CAST_INT32(li_Entry, ai_FileOffset + lui_Mid * sizeof(unsigned int));
		
				if (li_Entry > ai_MaxMass)
					lui_FindRight = lui_Mid - 1;
				else
					lui_FindLeft = lui_Mid;
			}
			else
			{
				// left and right pointer are now next to each other or even equal
				if (lui_FindRight != lui_FindLeft)
				{
					READ_UINT32_CAST_INT32(li_Entry, ai_FileOffset + lui_FindRight * sizeof(unsigned int));
					if (li_Entry > ai_MaxMass)
						--lui_FindRight;
				}
				lui_CropCount = lui_FindRight - lui_CropStart + 1;
				break;
			}
		} 
	}

	aui_CropStart = lui_CropStart;
	aui_CropCount = lui_CropCount;
}
