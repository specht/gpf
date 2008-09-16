#include "IntronSplitFixedHit.h"


k_IntronSplitFixedHit::k_IntronSplitFixedHit()
	: k_Hit(NULL, r_Scaffold(), "")
{
}


k_IntronSplitFixedHit::k_IntronSplitFixedHit(k_GpfQuery* ak_Query_, r_Scaffold ar_Scaffold, QString as_Query,
											 QString as_PeptideLeft, QString as_GapFix, QString as_PeptideRight,
											 unsigned int aui_GnoLeft, unsigned int aui_GnoRight,
											 unsigned int aui_NucleotideCountLeft, unsigned int aui_NucleotideCountRight,
											 unsigned int aui_IntronLength,
											 QString as_LeftSurroundings, QString as_RightSurroundings,
											 QString as_IntronEnds)
	: k_Hit(ak_Query_, ar_Scaffold, as_Query)
	, ms_PeptideLeft(as_PeptideLeft)
	, ms_GapFix(as_GapFix)
	, ms_PeptideRight(as_PeptideRight)
	, mui_GnoLeft(aui_GnoLeft)
	, mui_GnoRight(aui_GnoRight)
	, mui_NucleotideCountLeft(aui_NucleotideCountLeft)
	, mui_NucleotideCountRight(aui_NucleotideCountRight)
	, mui_IntronLength(aui_IntronLength)
	, ms_LeftSurroundings(as_LeftSurroundings)
	, ms_RightSurroundings(as_RightSurroundings)
	, ms_IntronEnds(as_IntronEnds)
{
	// +3 / -3 in the next line:
	// the right GNO describes the start of the last amino acid which itself 
	// consists of three nucleotides, so we skip past the whole peptide part 
	// and then go back to the beginning of the peptide
	ms_Assembly = QString("%1;%2%3:%4,%5:%6")
		.arg(mk_Query_->get_GpfIndexFileInfo()->get_Key())
		.arg((mui_GnoLeft & gui_GlobalNucleotideOffsetBackwardFlag) != 0? "-": "+")
		.arg(mui_GnoLeft & ~gui_GlobalNucleotideOffsetBackwardFlag)
		.arg(mui_NucleotideCountLeft)
		.arg((mui_GnoRight & gui_GlobalNucleotideOffsetBackwardFlag) == 0?
		((mui_GnoRight & ~gui_GlobalNucleotideOffsetBackwardFlag) + 1 - mui_NucleotideCountRight):
		((mui_GnoRight & ~gui_GlobalNucleotideOffsetBackwardFlag) - 1 + mui_NucleotideCountRight))
		.arg(mui_NucleotideCountRight);
}


k_IntronSplitFixedHit::~k_IntronSplitFixedHit()
{
}


QString k_IntronSplitFixedHit::description()
{
	if (mk_Query_ == NULL)
		return QString();

	mui_HitMass = mk_Query_->get_GpfBase().CalculatePeptideMass(ms_PeptideLeft + ms_GapFix + ms_PeptideRight);
	if (mui_HitMass < mk_Query_->mui_ResultMassMinimum || mui_HitMass > mk_Query_->mui_ResultMassMaximum)
	{
		// discard this hit because the mass doesnt fit...
		mb_IsDiscarded = true;
		return QString();
	}

	if (mk_Query_->GET_INT_PARAMETER(CropCsIntrons) == r_YesNoOption::Yes)
	{
		if (ms_IntronEnds != "GTAG" && ms_IntronEnds != "GCAG")
		{
			mb_IsDiscarded = true;
			return QString();
		}
	}

	QString ls_Description;

	QString ls_FullDetails;
	if (mk_Query_->GET_INT_PARAMETER(FullDetails) == r_YesNoOption::Yes)
		ls_FullDetails = QString(", details: %1")
			.arg(mk_Query_->get_GpfIndexFileInfo()->get_AssemblyInfoAsYaml(ms_Assembly));

	ls_Description += QString("- { peptide: '%1', score: %2, mass: %3, left: '%4', right: '%5', assembly: '%6', intronEnds: '%7'%8 }\n")
		.arg(ms_PeptideLeft + ms_GapFix + ms_PeptideRight)
		.arg(mui_MaxChainLength)
		.arg((double)mui_HitMass / gui_MassPrecision, 2)
		.arg(ms_LeftSurroundings)
		.arg(ms_RightSurroundings)
		.arg(ms_Assembly)
		.arg(ms_IntronEnds)
		.arg(ls_FullDetails);

	return ls_Description;
}


void k_IntronSplitFixedHit::CalculateScore()
{
	// For the score calculation, determine the score of both the left and 
	// the right half peptide fragment and then return the maximum.
	// If a triplet split occured, count the split amino acid for both halves.
	//return std::max<unsigned int>(CalculateScoreForChainMarker(r_ChainMarker::Left), CalculateScoreForChainMarker(r_ChainMarker::Right));
	mui_Score = CalculateScoreForChainMarker(r_ChainMarker::Left | r_ChainMarker::Right);
	mk_PartScores = QList<unsigned int>();
	mk_PartScores << 0;
	mk_PartScores << 0;
}
