
#include "IntronSplitHit.h"


k_IntronSplitHit::k_IntronSplitHit()
	: k_Hit(NULL, r_Scaffold(), "")
{
}


k_IntronSplitHit::k_IntronSplitHit(k_GpfQuery* ak_Query_, r_Scaffold ar_Scaffold, QString as_Query,
								   unsigned int aui_PeptideStartGno, unsigned int aui_LeftHalfLength, 
								   unsigned int aui_PeptideEndGno, unsigned int aui_RightHalfLength,
								   unsigned int aui_QueryLengthLeft, unsigned int aui_QueryLengthRight)
	: k_Hit(ak_Query_, ar_Scaffold, as_Query)
	, mui_PeptideStartGno(aui_PeptideStartGno)
	, mui_LeftHalfLength(aui_LeftHalfLength)
	, mui_PeptideEndGno(aui_PeptideEndGno)
	, mui_RightHalfLength(aui_RightHalfLength)
	, mui_QueryLengthLeft(aui_QueryLengthLeft)
	, mui_QueryLengthRight(aui_QueryLengthRight)
	, mui_QueryGapMass(0)
	, mb_ReadUnfixedPeptide(false)
{
}


k_IntronSplitHit::~k_IntronSplitHit()
{
}


void k_IntronSplitHit::Update(unsigned int aui_LeftHalfLength, unsigned int aui_RightHalfLength,
							  unsigned int aui_QueryLengthLeft, unsigned int aui_QueryLengthRight)
{
	mui_LeftHalfLength = std::min<unsigned int>(mui_LeftHalfLength, aui_LeftHalfLength);
	mui_RightHalfLength = std::min<unsigned int>(mui_RightHalfLength, aui_RightHalfLength);
	mui_QueryLengthLeft = std::min<unsigned int>(mui_QueryLengthLeft, aui_QueryLengthLeft);
	mui_QueryLengthRight = std::min<unsigned int>(mui_QueryLengthRight, aui_QueryLengthRight);
}


void k_IntronSplitHit::Finish()
{
	// count matched amino acids
	k_Hit::Finish();

	// try to fix the gap
	FixGap();

	if (mk_FixedHitList.empty())
		MarkDiscarded();
}


const QList<k_IntronSplitFixedHit>& k_IntronSplitHit::get_FixedHitList() const
{
	return mk_FixedHitList;
}


void k_IntronSplitHit::FixGap()
{
	mi_UnfixedIntronLength = (int)(abs((qint64)mui_PeptideStartGno - (qint64)mui_PeptideEndGno)) + 1 - mui_LeftHalfLength * 3 - mui_RightHalfLength * 3;

	// contains the exact gap amino acids, for later comparison if similarity search is disabled
	QString ls_TargetGapFixCollapsed;

	ms_QueryGap = ms_Query.mid(mui_QueryLengthLeft, ms_Query.length() - mui_QueryLengthLeft - mui_QueryLengthRight);

	// if similarity search is disabled, try an early exit and construct the target gap fix (ls_TargetGapFixCollapsed)
	ls_TargetGapFixCollapsed = mk_Query_->get_GpfBase().CollapsePeptide(ms_QueryGap);
	if (mk_Query_->GET_INT_PARAMETER(SearchSimilar) == r_YesNoOption::No && ls_TargetGapFixCollapsed.length() * 3 >= mi_UnfixedIntronLength)
		return;

	mui_QueryGapMass = mk_Query_->get_GpfBase().CalculatePeptideMass(ms_QueryGap) - mk_Query_->get_GpfBase().mui_WaterMass;

	/*
	if (ms_QueryGap.length() == 0)
		mk_GapFixAlternatives[QPair<unsigned short, unsigned short>(0, 0)] = "[yay!]";
		*/

	if (mi_UnfixedIntronLength > 0)
	{
		// examine intron nucleotides and see if we can construct a gap fix that completes the peptide...
		bool lb_ForwardReadingFrame = (mui_PeptideStartGno & gui_GlobalNucleotideOffsetBackwardFlag) == 0;

		// calculate min and max gap fix length
		unsigned int lui_MinWeight = mk_Query_->get_GpfBase().mui_AminoAcidWeight_[(unsigned char)r_AminoAcid::Gly];
		unsigned int lui_MaxWeight = mk_Query_->get_GpfBase().mui_AminoAcidWeight_[(unsigned char)r_AminoAcid::Trp];

		unsigned int lui_MinGapFixLength = mui_QueryGapMass / lui_MaxWeight;
		if (mui_QueryGapMass % lui_MaxWeight != 0)
			++lui_MinGapFixLength;

		unsigned int lui_MaxGapFixLength = std::min<unsigned int>(mui_QueryGapMass / lui_MinWeight, mi_UnfixedIntronLength / 3);

		int li_MassError = ((qint64)mk_Query_->mui_QueryMass * mk_Query_->GET_INT_PARAMETER(MassError) / 1000000);
		unsigned int lui_MinGapFixMass = (unsigned int)std::max<qint64>((qint64)mui_QueryGapMass - li_MassError, 0);
		unsigned int lui_MaxGapFixMass = (unsigned int)std::min<qint64>((qint64)mui_QueryGapMass + li_MassError, 4294967295);

		// read all gap nucleotides 
		RefPtr<char> lc_pGapNucleotides = RefPtr<char>(new char[mi_UnfixedIntronLength]);
		if (lb_ForwardReadingFrame)
		{
			mk_Query_->get_IndexFile().seek(mr_Scaffold.mi_DnaFilePosition + mui_PeptideStartGno - mr_Scaffold.mui_NucleotideOffset + mui_LeftHalfLength * 3);
			mk_Query_->get_IndexFile().read(lc_pGapNucleotides.get_Pointer(), mi_UnfixedIntronLength);
		}
		else
		{
			RefPtr<char> lc_pTemp = RefPtr<char>(new char[mi_UnfixedIntronLength]);
			mk_Query_->get_IndexFile().seek(mr_Scaffold.mi_DnaFilePosition + (mui_PeptideEndGno & ~gui_GlobalNucleotideOffsetBackwardFlag) - mr_Scaffold.mui_NucleotideOffset + mui_RightHalfLength * 3);
			mk_Query_->get_IndexFile().read(lc_pTemp.get_Pointer(), mi_UnfixedIntronLength);
			// reverse and transpose DNA
			for (unsigned int i = 0; i < (unsigned int)mi_UnfixedIntronLength; ++i)
				lc_pGapNucleotides.get_Pointer()[i] = lc_pTemp.get_Pointer()[mi_UnfixedIntronLength - i - 1] ^ 3;
		}

		RefPtr<char> lc_pGapLeftNucleotides = RefPtr<char>(new char[lui_MaxGapFixLength * 3]);
		RefPtr<char> lc_pGapRightNucleotides = RefPtr<char>(new char[lui_MaxGapFixLength * 3]);
		RefPtr<char> lc_pGapFixNucleotides = RefPtr<char>(new char[lui_MaxGapFixLength * 3]);
		RefPtr<char> lc_pGapFixPeptide = RefPtr<char>(new char[lui_MaxGapFixLength]);
		char lc_IntronEnds_[5];
		lc_IntronEnds_[4] = 0;

		for (unsigned int lui_GapFixSize = lui_MinGapFixLength; lui_GapFixSize <= lui_MaxGapFixLength && (int)lui_GapFixSize * 3 < mi_UnfixedIntronLength; ++lui_GapFixSize)
		{
			// try to reconstruct a feasible gap fix with a length of lui_GapFixSize

			memcpy(lc_pGapLeftNucleotides.get_Pointer(), lc_pGapNucleotides.get_Pointer(), lui_GapFixSize * 3);
			memcpy(lc_pGapRightNucleotides.get_Pointer(), lc_pGapNucleotides.get_Pointer() + mi_UnfixedIntronLength - lui_GapFixSize * 3, lui_GapFixSize * 3);

			unsigned int lui_GapNucleotideCount = lui_GapFixSize * 3;
			unsigned int lui_SplitStep = mk_Query_->GET_INT_PARAMETER(SplitTriplets) == r_YesNoOption::Yes? 1: 3;
			// move split from left to right, ignore full left and right positions (these positions don't produce a mix of left and right parts)
			for (unsigned int lui_SplitOffset = 0; lui_SplitOffset <= lui_GapNucleotideCount; lui_SplitOffset += lui_SplitStep)
			{
				// construct gap peptide with split at current split offset
				memcpy(lc_pGapFixNucleotides.get_Pointer(), lc_pGapLeftNucleotides.get_Pointer(), lui_SplitOffset);
				memcpy(lc_pGapFixNucleotides.get_Pointer() + lui_SplitOffset, lc_pGapRightNucleotides.get_Pointer() + lui_SplitOffset, lui_GapNucleotideCount - lui_SplitOffset);

				// construct intron ends
				memset(lc_IntronEnds_, 0, 5);
				if (mi_UnfixedIntronLength - lui_GapFixSize * 3 >= 4)
				{
					memcpy(lc_IntronEnds_, lc_pGapNucleotides.get_Pointer() + lui_SplitOffset, 2);
					memcpy(lc_IntronEnds_ + 2, lc_pGapNucleotides.get_Pointer() + mi_UnfixedIntronLength - lui_GapFixSize * 3 + lui_SplitOffset - 2, 2);
				}
				for (int i = 0; i < 4; ++i)
					lc_IntronEnds_[i] = mk_Query_->get_GpfBase().mc_NucleotideIntToChar_[(unsigned char)lc_IntronEnds_[i]];

				// luc_pGapPeptide contains the current fix proposal, now translate and check it
				unsigned int lui_GapPeptideMass = 0;
				unsigned int lui_UnknownAminoAcidCount = 0;
				unsigned int lui_StopCodonCount = 0;
				for (unsigned int i = 0; i < lui_GapFixSize; ++i)
				{
					unsigned int lui_Triplet = 
						(((int)lc_pGapFixNucleotides.get_Pointer()[i * 3 + 0] & 7) << 6) |
						(((int)lc_pGapFixNucleotides.get_Pointer()[i * 3 + 1] & 7) << 3) |
						(((int)lc_pGapFixNucleotides.get_Pointer()[i * 3 + 2] & 7));
					char lc_AminoAcid = mk_Query_->get_GpfBase().mc_TripletToAminoAcidForward_[lui_Triplet];
					r_AminoAcid::Enumeration le_AminoAcid = (r_AminoAcid::Enumeration)((unsigned char)lc_AminoAcid);
					lc_pGapFixPeptide.get_Pointer()[i] = mk_Query_->get_GpfBase().mc_AminoAcidToChar_[(unsigned char)le_AminoAcid];
					lui_GapPeptideMass += mk_Query_->get_GpfBase().mui_AminoAcidWeight_[(unsigned char)le_AminoAcid];
					if (le_AminoAcid == r_AminoAcid::Stop)
						++lui_StopCodonCount;
					else if (le_AminoAcid == r_AminoAcid::Unknown)
						++lui_UnknownAminoAcidCount;
				}
				if (lui_UnknownAminoAcidCount <= gui_MaxUnknownAminoAcidCount && 
					lui_StopCodonCount == 0 && 
					lui_GapPeptideMass >= lui_MinGapFixMass && 
					lui_GapPeptideMass <= lui_MaxGapFixMass)
				{
					// we found a good gap fix, add it to the list of gap fix alternatives
					QString ls_GapFix, ls_GapFixCollapsed;
					for (unsigned int i = 0; i < lui_GapFixSize; ++i)
					{
						ls_GapFix += QChar(lc_pGapFixPeptide.get_Pointer()[i]);
						ls_GapFixCollapsed += QChar(mk_Query_->get_GpfBase().mc_AminoAcidToChar_[(unsigned int)mk_Query_->get_GpfBase().me_AminoAcidMassSimilarityCollapse[mk_Query_->get_GpfBase().me_CharToAminoAcid_[ls_GapFix[i].toAscii()]]]);
					}

					if ((mk_Query_->GET_INT_PARAMETER(SearchSimilar) == r_YesNoOption::Yes) ||
						(mk_Query_->GET_INT_PARAMETER(SearchSimilar) == r_YesNoOption::No && ls_GapFixCollapsed == ls_TargetGapFixCollapsed))
					{
						unsigned int lui_LeftAdd = lui_SplitOffset;
						unsigned int lui_RightAdd = lui_GapNucleotideCount - lui_SplitOffset;

						ReadUnfixedPeptide();

						k_IntronSplitFixedHit lk_IntronSplitFixedHit(mk_Query_, mr_Scaffold, ms_Query, 
							ms_PeptideLeft, ls_GapFix, ms_PeptideRight,
							mui_PeptideStartGno, mui_PeptideEndGno, 
							mui_LeftHalfLength * 3 + lui_LeftAdd, mui_RightHalfLength * 3 + lui_RightAdd,
							mi_UnfixedIntronLength - lui_LeftAdd - lui_RightAdd,
							ms_LeftSurroundings, ms_RightSurroundings, QString(lc_IntronEnds_));

						lk_IntronSplitFixedHit.MarkChainElements(mui_QueryLengthLeft - 3, 3, r_ChainMarker::Left);
						lk_IntronSplitFixedHit.MarkChainElements(mui_QueryLength - mui_QueryLengthRight, 3, r_ChainMarker::Right);

						// compare ls_GapFixCollapsed to ls_TargetGapFixCollapsed from left and right border
						for (int i = 0; i < std::min<int>(ls_GapFixCollapsed.length(), ls_TargetGapFixCollapsed.length()) && ls_GapFixCollapsed[i] == ls_TargetGapFixCollapsed[i]; ++i)
							lk_IntronSplitFixedHit.MarkChainElements(mui_QueryLengthLeft + i, 1, r_ChainMarker::Left);

						for (int i = 0; i < std::min<int>(ls_GapFixCollapsed.length(), ls_TargetGapFixCollapsed.length()) && ls_GapFixCollapsed[ls_GapFixCollapsed.length() - i - 1] == ls_TargetGapFixCollapsed[ls_TargetGapFixCollapsed.length() - i - 1]; ++i)
							lk_IntronSplitFixedHit.MarkChainElements(mui_QueryLength - mui_QueryLengthRight - 1 - i, 1, r_ChainMarker::Right);

						lk_IntronSplitFixedHit.Finish();

						mk_FixedHitList.append(lk_IntronSplitFixedHit);
					}
				}
			}
		}
	}
}


void k_IntronSplitHit::ReadUnfixedPeptide()
{
	if (mb_ReadUnfixedPeptide)
		return;

	QString ls_Temp;
	bool lb_Forward = (mui_PeptideStartGno & gui_GlobalNucleotideOffsetBackwardFlag) == 0;

	QString ls_Assembly = QString("%1;%2%3:%4")
		.arg(mk_Query_->get_GpfIndexFileInfo()->get_Key())
		.arg(lb_Forward? "+": "-")
		.arg(mui_PeptideStartGno & ~gui_GlobalNucleotideOffsetBackwardFlag)
		.arg(mui_LeftHalfLength * 3);

	mk_Query_->get_GpfIndexFileInfo()->readAssembly(ls_Assembly, ms_PeptideLeft, ms_LeftSurroundings, ls_Temp);

	if (lb_Forward)
	{
		ls_Assembly = QString("%1;%2%3:%4")
			.arg(mk_Query_->get_GpfIndexFileInfo()->get_Key())
			.arg(lb_Forward? "+": "-")
			.arg(mui_PeptideEndGno - mui_RightHalfLength * 3 + 1)
			.arg(mui_RightHalfLength * 3);
	}
	else
	{
		ls_Assembly = QString("%1;%2%3:%4")
			.arg(mk_Query_->get_GpfIndexFileInfo()->get_Key())
			.arg(lb_Forward? "+": "-")
			.arg((mui_PeptideEndGno & ~gui_GlobalNucleotideOffsetBackwardFlag) + mui_RightHalfLength * 3 - 1)
			.arg(mui_RightHalfLength * 3);
	}

	mk_Query_->get_GpfIndexFileInfo()->readAssembly(ls_Assembly, ms_PeptideRight, ls_Temp, ms_RightSurroundings);

	mb_ReadUnfixedPeptide = true;
}


QString k_IntronSplitHit::description()
{
	return QString();
}
