#pragma once


// GPF parameter names
struct r_GpfParameterName
{
	enum Enumeration
	{
		Masses,
		Protease,
		MassError,
		SearchSimilar,
		SearchIntrons,
		MaxIntronLength,
		MinChainLength,
		FullDetails
	};
};


// GPF enum definitions
struct r_YesNoOption
{
	enum Enumeration
	{
		No,
		Yes
	};
};


struct r_ProteaseOption
{
	enum Enumeration
	{
		Trypsin
	};
};


struct r_MassesOption
{
	enum Enumeration
	{
		Monoisotopic
	};
};


struct r_SearchIntronOption
{
	enum Enumeration
	{
		Never,
		Always,
		IfNecessary
	};
};


// GPF enum options
extern QHash<QString, int> gk_YesNoOptions;
extern QHash<int, QString> gk_YesNoOptionsReverse;
extern QHash<QString, int> gk_ProteaseOptions;
extern QHash<int, QString> gk_ProteaseOptionsReverse;
extern QHash<QString, int> gk_MassesOptions;
extern QHash<int, QString> gk_MassesOptionsReverse;
extern QHash<QString, int> gk_SearchIntronOptions;
extern QHash<int, QString> gk_SearchIntronOptionsReverse;
