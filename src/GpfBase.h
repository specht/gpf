#pragma once

#include <QtCore>
#include "GpfConfig.h"
#include "RefPtr.h"
#include "generated/GpfParameters.h"


class k_GpfIndexFileInfo;


#pragma pack(push)
#pragma pack(1)

struct r_IndexEntry
{
	unsigned short muw_Trimer;
	unsigned int mui_Gno;
	unsigned int mui_LeftMass;
	unsigned int mui_RightMass;
};


struct r_PeptideSpan
{
	unsigned int mui_Gno;
	unsigned short muw_Length;
};

#pragma pack(pop)


struct r_AminoAcid
{
	enum Enumeration
	{
		Gly = 0, Ala, Ser, Pro, Val, Thr, Cys, Leu, Ile, Asn, Asp, Gln, Lys, Glu, Met, His, Phe, Arg, Tyr, Trp,
		Stop = 0x80, Unknown
	};

	static bool isStop(r_AminoAcid::Enumeration ae_AminoAcid)
	{
		return ae_AminoAcid == r_AminoAcid::Stop;
	}

	static bool isTrypticCleavageSite(r_AminoAcid::Enumeration ae_AminoAcid)
	{
		return ae_AminoAcid == r_AminoAcid::Arg || ae_AminoAcid == r_AminoAcid::Lys;
	}

	static bool isTrypticCleavageSiteOrStop(r_AminoAcid::Enumeration ae_AminoAcid)
	{
		return isStop(ae_AminoAcid) || isTrypticCleavageSite(ae_AminoAcid);
	}

	static bool isUnknown(r_AminoAcid::Enumeration ae_AminoAcid)
	{
		return ae_AminoAcid == r_AminoAcid::Unknown;
	}
};


struct r_PeptideBorder
{
	enum Enumeration
	{
		Start,
		End
	};
};


struct r_ReadingFrame
{
	unsigned int mui_AminoAcidCount;
	qint64 mi_ReadingFrameFilePosition;
};


struct r_Scaffold
{
	r_Scaffold()
		: ms_Id(QString())
		, mi_DnaFilePosition(0)
		, mui_NucleotideCount(0)
		, mui_NucleotideOffset(0)
		, mui_TotalAminoAcidCount(0)
	{
	}

	r_Scaffold(QString as_Id, qint64 ai_FilePosition)
		: ms_Id(as_Id)
		, mi_DnaFilePosition(ai_FilePosition)
		, mui_NucleotideCount(0)
		, mui_NucleotideOffset(0)
		, mui_TotalAminoAcidCount(0)
	{
	}

	// Note: The global nucleotide offset is the offset of the first nucleotide
	// encoding the specified amino acid in the specified reading frame
	inline unsigned int getGlobalNucleotideOffset(unsigned int aui_ReadingFrame, unsigned int aui_AminoAcidPosition)
	{
		/*
		if (aui_AminoAcidPosition >= mui_NucleotideCount)
			printf("ATTENTION: getGlobalNucleotideOffset called with %d, count is %d.\n", aui_AminoAcidPosition, mui_NucleotideCount);
			*/

		unsigned int lui_Offset;

		if (aui_ReadingFrame < 3)
			// forward reading frame
			lui_Offset = aui_AminoAcidPosition * 3 + aui_ReadingFrame;
		else
			// backward reading frame
			lui_Offset = (mui_NucleotideCount - 3 - aui_AminoAcidPosition * 3 - (aui_ReadingFrame - 3) + 2) | gui_GlobalNucleotideOffsetBackwardFlag;

		return lui_Offset + mui_NucleotideOffset;
	}

	inline qint64 getAminoAcidFilePositionForGlobalNucleotideOffset(unsigned int aui_Gno)
	{
		bool lb_Forward = (aui_Gno & gui_GlobalNucleotideOffsetBackwardFlag) == 0;
		aui_Gno &= ~gui_GlobalNucleotideOffsetBackwardFlag;
		aui_Gno -= mui_NucleotideOffset;

		if (lb_Forward)
		{
			unsigned int lui_FrameNumber = aui_Gno % 3;
			unsigned int lui_AminoAcidOffset = aui_Gno / 3;
			return mr_ReadingFrames[lui_FrameNumber].mi_ReadingFrameFilePosition + lui_AminoAcidOffset;
		}
		else
		{
			unsigned int lui_FrameNumber = (mui_NucleotideCount - aui_Gno - 1) % 3 + 3;
			unsigned int lui_AminoAcidOffset = (mui_NucleotideCount - aui_Gno - 1) / 3;
			return mr_ReadingFrames[lui_FrameNumber].mi_ReadingFrameFilePosition + lui_AminoAcidOffset;
		}
	}

	inline unsigned int getReadingFrameForGlobalNucleotideOffset(unsigned int aui_Gno)
	{
		bool lb_Forward = (aui_Gno & gui_GlobalNucleotideOffsetBackwardFlag) == 0;
		aui_Gno &= ~gui_GlobalNucleotideOffsetBackwardFlag;
		aui_Gno -= mui_NucleotideOffset;

		if (lb_Forward)
			return aui_Gno % 3;
		else
			return (mui_NucleotideCount - aui_Gno - 1) % 3 + 3;
	}

	QString ms_Id;
	qint64 mi_DnaFilePosition;
	unsigned int mui_NucleotideCount;
	unsigned int mui_NucleotideOffset;
	r_ReadingFrame mr_ReadingFrames[6];
	qint64 mi_AminoAcidScaffoldFilePosition;
	unsigned int mui_TotalAminoAcidCount;
};


class k_GpfParameter
{
public:
	k_GpfParameter()
	{
	}

	k_GpfParameter(QString as_Id, QString as_Label, QString as_Description)
		: ms_Id(as_Id)
		, ms_Label(as_Label)
		, ms_Description(as_Description)
	{
	}

	virtual ~k_GpfParameter() {}

	QString getId() const
	{
		return ms_Id;
	}

	QString getLabel() const
	{
		return ms_Label;
	}

	QString getDescription() const
	{
		return ms_Description;
	}

	virtual qint64 getValueAsInt() const
	{
		return 0;
	}

	virtual double getValueAsDouble() const
	{
		return 0;
	}

	virtual QString getValueAsString() const
	{
		return "";
	}

	virtual void reset() {}

	virtual void setValue(QString as_Value) = 0;

protected:
	QString ms_Id;
	QString ms_Label;
	QString ms_Description;
};


template <typename T> class k_GpfParameterT: public k_GpfParameter
{
public:
	k_GpfParameterT()
		: k_GpfParameter()
	{
	}

	k_GpfParameterT(QString as_Id, QString as_Label, QString as_Description, T at_DefaultValue)
		: k_GpfParameter(as_Id, as_Label, as_Description)
		, mt_DefaultValue(at_DefaultValue)
	{
		this->reset();
	}

	virtual ~k_GpfParameterT() {}

	virtual void reset()
	{
		mt_Value = mt_DefaultValue;
	}

	T getValue() const
	{
		return mt_Value;
	}

	void setValue(T at_Value)
	{
		mt_Value = at_Value;
	}

protected:
	T mt_DefaultValue;
	T mt_Value;
};


class k_GpfParameterInt: public k_GpfParameterT<qint64> 
{
public:
	k_GpfParameterInt()
		: k_GpfParameterT<qint64>()
	{
	}

	k_GpfParameterInt(QString as_Id, QString as_Label, QString as_Description, qint64 ai_DefaultValue)
		: k_GpfParameterT<qint64>(as_Id, as_Label, as_Description, ai_DefaultValue)
	{
	}

	virtual qint64 getValueAsInt() const
	{
		return mt_Value;
	}

	virtual double getValueAsDouble() const
	{
		return (double)mt_Value;
	}

	virtual QString getValueAsString() const
	{
		return QString("%1").arg(mt_Value);
	}

	virtual void setValue(QString as_Value)
	{
		mt_Value = as_Value.toInt();
	}
};


class k_GpfParameterEnum: public k_GpfParameterInt
{
public:
	k_GpfParameterEnum()
		: k_GpfParameterInt()
		, mk_Options_(NULL)
		, mk_OptionsReverse_(NULL)
	{
	}

	k_GpfParameterEnum(QString as_Id, QString as_Label, QString as_Description, qint64 ai_DefaultValue, 
		               QHash<QString, int>& ak_Options,
		               QHash<int, QString>& ak_OptionsReverse)
		: k_GpfParameterInt(as_Id, as_Label, as_Description, ai_DefaultValue)
		, mk_Options_(&ak_Options)
		, mk_OptionsReverse_(&ak_OptionsReverse)
	{
	}

	virtual qint64 getValueAsInt() const
	{
		return mt_Value;
	}

	virtual double getValueAsDouble() const
	{
		return (double)mt_Value;
	}

	virtual QString getValueAsString() const
	{
		if (mk_OptionsReverse_ != NULL && mk_OptionsReverse_->contains(mt_Value))
			return (*mk_OptionsReverse_)[mt_Value];

		return QString("%1").arg(mt_Value);
	}

	virtual void setValue(QString as_Value)
	{
		if (mk_Options_ == NULL)
			return;

		if (mk_Options_->contains(as_Value))
			mt_Value = mk_Options_->value(as_Value);
	}

	QHash<QString, int>* mk_Options_;
	QHash<int, QString>* mk_OptionsReverse_;
};


class k_GpfParameterString: public k_GpfParameterT<QString> 
{
public:
	k_GpfParameterString()
		: k_GpfParameterT<QString>()
	{
	}

	k_GpfParameterString(QString as_Id, QString as_Label, QString as_Description, QString as_DefaultValue)
		: k_GpfParameterT<QString>(as_Id, as_Label, as_Description, as_DefaultValue)
	{
	}

	virtual QString getValueAsString() const
	{
		return mt_Value;
	}

	virtual void setValue(QString as_Value)
	{
		mt_Value = as_Value;
	}
};


class k_GpfParameterDouble: public k_GpfParameterT<double> 
{
public:
	k_GpfParameterDouble()
		: k_GpfParameterT<double>()
	{
	}

	k_GpfParameterDouble(QString as_Id, QString as_Label, QString as_Description, double ad_DefaultValue)
		: k_GpfParameterT<double>(as_Id, as_Label, as_Description, ad_DefaultValue)
	{
	}

	virtual qint64 getValueAsInt() const
	{
		return (qint64)mt_Value;
	}

	virtual double getValueAsDouble() const
	{
		return mt_Value;
	}

	virtual QString getValueAsString() const
	{
		return QString("%1").arg(mt_Value);
	}

	virtual void setValue(QString as_Value)
	{
		mt_Value = as_Value.toDouble();
	}
};


#define GET_INT_PARAMETER(as_Name) mk_Parameters[r_GpfParameterName::as_Name]->getValueAsInt()
#define GET_STRING_PARAMETER(as_Name) mk_Parameters[r_GpfParameterName::as_Name]->getValueAsString()
#define GET_DOUBLE_PARAMETER(as_Name) mk_Parameters[r_GpfParameterName::as_Name]->getValueAsDouble()


// various comparison functions
bool lt_Gno(unsigned int* aui_First_, unsigned int* aui_Second_);
bool lt_Direction_Gno(unsigned int* aui_First_, unsigned int* aui_Second_);
bool lt_IndexEntry_Trimer_LeftMass_Gno(r_IndexEntry* ar_First_, r_IndexEntry* ar_Second_);
bool lt_IndexEntry_Trimer_RightMass_Gno(r_IndexEntry* ar_First_, r_IndexEntry* ar_Second_);
bool lt_PeptideSpan_Gno(r_PeptideSpan* ar_First_, r_PeptideSpan* ar_Second_);


class k_GpfBase
{
public:
	// use the following ctor if you want to compile a new index
	k_GpfBase();
	// use the following constructor if you already have index files and want them to be loaded
	k_GpfBase(QStringList ak_IndexFiles);

	virtual ~k_GpfBase();

	RefPtr<k_GpfIndexFileInfo> get_IndexFileInfo(QString as_Id);
	RefPtr<k_GpfIndexFileInfo> get_IndexFileInfoFromAssembly(QString as_Assembly);
	RefPtr<k_GpfIndexFileInfo> get_DefaultIndexFileInfo();
	QString get_GpfParameters() const;
	unsigned int CalculatePeptideMass(QString as_Peptide);
	QString CollapsePeptide(QString as_Peptide);

	// GPF-wide constants and tables, without accessor methods to increase
	// execution speed - please don't change these values in other classes!
	r_AminoAcid::Enumeration me_CharToAminoAcid_[256];
	char mc_AminoAcidToChar_[256];
	char mc_NucleotideCharToInt_[256];
	char mc_NucleotideIntToChar_[256];
	double md_AminoAcidWeight_[256];
	unsigned int mui_AminoAcidWeight_[256];
	r_AminoAcid::Enumeration me_AminoAcidMassSimilarityCollapse[256];
	char mc_TripletToAminoAcidForward_[512];
	char mc_TripletToAminoAcidBackward_[512];
	unsigned int mui_WaterMass;

private:
	void FillTables();

	QMap<QString, RefPtr<k_GpfIndexFileInfo> > mk_IndexFiles;

	QString ms_GpfParameters;
};
