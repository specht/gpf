#include "Hit.h"


k_Hit::k_Hit(k_GpfQuery& ak_Query, bool ab_Forward, QList<QPair<unsigned int, unsigned int> > ak_Assembly)
	: mk_Query(ak_Query)
	, mb_Forward(ab_Forward)
	, mk_Assembly(ak_Assembly)
	, mui_QueryLength(ak_Query.get_Query().length())
	, mui_MaxChainLength(0)
	, mb_IsDiscarded(false)
{
	Initialize();
}


k_Hit::k_Hit(const k_Hit& ak_Other)
	: mk_Query(ak_Other.mk_Query)
	, mb_Forward(ak_Other.mb_Forward)
	, mk_Assembly(ak_Other.mk_Assembly)
	, mui_QueryLength(ak_Other.mui_QueryLength)
	, mui_MaxChainLength(ak_Other.mui_MaxChainLength)
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
	mui_MaxChainLength = CalculateScore();

	if (mui_MaxChainLength < mk_Query.GET_INT_PARAMETER(MinChainLength))
	{
		mb_IsDiscarded = true;
		return;
	}

	if (mk_Query.GET_INT_PARAMETER(SearchSimilar) == r_YesNoOption::No && mk_Query.get_GpfBase().CollapsePeptide(ms_Peptide) != mk_Query.get_CollapsedQuery())
	{
		mb_IsDiscarded = true;
		return;
	}


	mk_ResultItems.append(QPair<QString, QString>("score", QString("%1").arg(mui_MaxChainLength)));
}


bool k_Hit::get_HasFullScore() const
{
	return mui_MaxChainLength == mui_QueryLength;
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


unsigned int k_Hit::CalculateScore() const
{
	QString ls_CollapsedPeptide = mk_Query.get_GpfBase().CollapsePeptide(ms_Peptide);
	QString ls_CollapsedQuery = mk_Query.get_CollapsedQuery();

	QString ls_Marks = ls_CollapsedPeptide;
	RefPtr<bool> lb_pMarks(new bool[ls_CollapsedQuery.length()]);
	memset(lb_pMarks.get_Pointer(), false, ls_CollapsedQuery.length());

	int li_LastIndex = -1;
	for (int i = 0; i < ls_CollapsedQuery.length() - 2; ++i)
	{
		QString ls_QueryTrimer = ls_CollapsedQuery.mid(i, 3);
		int li_Index = -2;
		while (li_Index != -1 && li_Index <= li_LastIndex)
			li_Index = ls_CollapsedPeptide.indexOf(ls_QueryTrimer, li_Index == -2? 0: li_Index + 1);

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

	return li_MaxLength;
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
