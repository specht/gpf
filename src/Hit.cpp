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

#include "Hit.h"


k_Hit::k_Hit(k_GpfQuery& ak_Query, bool ab_Forward, QList<QPair<unsigned int, unsigned int> > ak_Assembly)
	: mk_Query(ak_Query)
	, mb_Forward(ab_Forward)
	, mk_Assembly(ak_Assembly)
	, mui_QueryLength(ak_Query.get_Query().length())
	, mi_Score(0)
	, mk_PartScores(QList<unsigned int>())
	, mb_IsDiscarded(false)
{
	Initialize();
}


k_Hit::k_Hit(const k_Hit& ak_Other)
	: mk_Query(ak_Other.mk_Query)
	, mb_Forward(ak_Other.mb_Forward)
	, mk_Assembly(ak_Other.mk_Assembly)
	, mui_QueryLength(ak_Other.mui_QueryLength)
	, mi_Score(ak_Other.mi_Score)
	, mk_PartScores(ak_Other.mk_PartScores)
	, mb_IsDiscarded(ak_Other.mb_IsDiscarded)
{
	Initialize();
}


k_Hit::~k_Hit() 
{
}


QString k_Hit::get_Peptide() const
{
	return ms_Peptide;
}


QString k_Hit::get_Assembly() const
{
	return ms_Assembly;
}


unsigned int k_Hit::get_Mass() const
{
	return mui_HitMass;
}


void k_Hit::AddInformation(QString as_Key, QString as_Value)
{
	mk_ResultItems.append(QPair<QString, QString>(as_Key, as_Value));
}


void k_Hit::Finish()
{
	CalculateScore();

	if (mi_Score < mk_Query.GET_INT_PARAMETER(MinChainLength))
	{
		mb_IsDiscarded = true;
		return;
	}

	if (mk_Query.GET_INT_PARAMETER(SearchSimilar) == r_YesNoOption::No && mk_Query.get_GpfBase().CollapsePeptide(ms_Peptide) != mk_Query.get_CollapsedQuery())
	{
		mb_IsDiscarded = true;
		return;
	}


	mk_ResultItems.append(QPair<QString, QString>("score", QString("%1").arg(mi_Score)));
	QStringList lk_PartScoresString;
	foreach (unsigned int li_Score, mk_PartScores)
		lk_PartScoresString << QString("%1").arg(li_Score);
	mk_ResultItems.append(QPair<QString, QString>("partScores", QString("[%1]").arg(lk_PartScoresString.join(", "))));
}


bool k_Hit::get_HasFullScore() const
{
	return mi_Score == mui_QueryLength;
}


bool k_Hit::get_IsDiscarded() const
{
	return mb_IsDiscarded;
}


void k_Hit::MarkDiscarded()
{
	mb_IsDiscarded = true;
}


QString k_Hit::description()
{
	QString ls_Description = "{ ";
	for (int i = 0; i < mk_ResultItems.size(); ++i)
	{
		ls_Description += QString("%1: %2").arg(mk_ResultItems[i].first, mk_ResultItems[i].second);
		if (i < mk_ResultItems.size() - 1)
			ls_Description += ", ";
	}
	ls_Description += "}";

	return ls_Description;
}


void k_Hit::CalculateScore()
{
	// TODO: k_IntronSplitFixedHit is oviously not used anymore since we changed
	// to consensus sequence intron searching... so everything happens within k_Hit,
	// possibly rendering k_IntronSplitHit and k_IntronSplitFixedHit obsolete?
	// ...do some housekeeping!
	// TODO number 2: it might happen that the score is not calculated correctly -
	// maybe we need a better matching algorithm: if a trimer repeats in the peptide
	// on the left and on the right side, then it may happen, that the chunky bacon
	// is missed by the linear matching algorithm which is applied here.
	
	QString ls_CollapsedPeptide = mk_Query.get_GpfBase().CollapsePeptide(ms_Peptide);
	QString ls_CollapsedQuery = mk_Query.get_CollapsedQuery();
	
	int li_TotalScore = 0;	
	mk_PartScores = QList<unsigned int>();
	
	typedef QPair<unsigned int, unsigned int> tk_UIntPair;
	// calculate score for each exon and determine total score
	unsigned int lui_AssemblyOffset = 0;
	foreach (tk_UIntPair lk_AssemblyPart, mk_Assembly)
	{
		QString ls_CollapsedPeptidePartTriple;
		for (int i = 0; i < ls_CollapsedPeptide.size(); ++i)
			ls_CollapsedPeptidePartTriple += "...";
		for (int i = 0; i < lk_AssemblyPart.second; ++i)
		{
			unsigned int lui_Index = i + lui_AssemblyOffset;
			ls_CollapsedPeptidePartTriple[lui_Index] = ls_CollapsedPeptide[lui_Index / 3];
		}
		QString ls_CollapsedPeptidePart;
		// ATTENTION: amino acids that are intron-split are taken into account
		// for both exons! this is what happens in the next few lines.
		for (int i = 0; i < ls_CollapsedPeptide.size(); ++i)
		{
			QString ls_Portion = ls_CollapsedPeptidePartTriple.mid(i * 3, 3);
			if (ls_Portion == "...")
				ls_CollapsedPeptidePart += ".";
			else
			{
				ls_Portion.replace(".", "");
				ls_CollapsedPeptidePart += ls_Portion[0];
			}
		}
		lui_AssemblyOffset += lk_AssemblyPart.second;
		
		int li_Score = 0;
		
		RefPtr<bool> lb_pMarks(new bool[ls_CollapsedQuery.length()]);
		memset(lb_pMarks.get_Pointer(), false, ls_CollapsedQuery.length());
	
		int li_LastIndex = -1;
		for (int i = 0; i < ls_CollapsedQuery.length() - 2; ++i)
		{
			QString ls_QueryTrimer = ls_CollapsedQuery.mid(i, 3);
			int li_Index = -2;
			while (li_Index != -1 && li_Index <= li_LastIndex)
				li_Index = ls_CollapsedPeptidePart.indexOf(ls_QueryTrimer, li_Index == -2? 0: li_Index + 1);
	
			if (li_Index > li_LastIndex)
			{
				li_LastIndex = li_Index;
				for (int k = 0; k < 3; ++k)
					lb_pMarks.get_Pointer()[i + k] = true;
			}
		}
	
		int li_MaxLength = 0;
		int li_CurrentLength = 0;
		bool lb_Inside = false;
		for (int i = 0; i < ls_CollapsedQuery.length(); ++i)
		{
			if (lb_pMarks.get_Pointer()[i])
			{
				if (!lb_Inside)
				{
					lb_Inside = true;
					li_CurrentLength = 1;
				}
				else
					++li_CurrentLength;
			}
			else
			{
				if (lb_Inside)
				{
					li_MaxLength = std::max<int>(li_MaxLength, li_CurrentLength);
					lb_Inside = false;
					li_CurrentLength = 0;
				}
			}
		}
		li_MaxLength = std::max<int>(li_MaxLength, li_CurrentLength);
		li_Score = li_MaxLength;
	
		mk_PartScores << li_Score;
		li_TotalScore = std::max<int>(li_TotalScore, li_Score);
	}
	
	mi_Score = li_TotalScore;
}


void k_Hit::Initialize()
{
	// construct assembly

	QStringList lk_Parts;
	for (int i = 0; i < mk_Assembly.size(); ++i)
		lk_Parts.append(QString("%1:%2").arg(mk_Assembly[i].first).arg(mk_Assembly[i].second));

	ms_Assembly = QString("%1;%2%3")
		.arg(mk_Query.get_GpfIndexFileInfo()->get_Key())
		.arg(mb_Forward? "+": "-")
		.arg(lk_Parts.join(","));

	// read peptide and left and right surrounding amino acids
	mk_Query.get_GpfIndexFileInfo()->readAssembly(ms_Assembly, ms_Peptide, ms_Left, ms_Right);

	mk_ResultItems.append(QPair<QString, QString>("peptide", "'" + ms_Peptide + "'"));
	mk_ResultItems.append(QPair<QString, QString>("assembly", "'" + ms_Assembly + "'"));
	mk_ResultItems.append(QPair<QString, QString>("left", "'" + ms_Left + "'"));
	mk_ResultItems.append(QPair<QString, QString>("right", "'" + ms_Right + "'"));
	mui_HitMass = mk_Query.get_GpfBase().CalculatePeptideMass(ms_Peptide);
	mk_ResultItems.append(QPair<QString, QString>("mass", QString("%1").arg((double)mui_HitMass / gui_MassPrecision)));

	if (mk_Query.GET_INT_PARAMETER(FullDetails) == r_YesNoOption::Yes)
		mk_ResultItems.append(QPair<QString, QString>("details", mk_Query.get_GpfIndexFileInfo()->get_AssemblyInfoAsYaml(ms_Assembly)));

}
