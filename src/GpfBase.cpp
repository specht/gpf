#include "GpfBase.h"
#include "GpfIndexFileInfo.h"
#include "generated/GpfOptionsDeclare.inc.cpp"


k_GpfBase::k_GpfBase()
{
	FillTables();
}


k_GpfBase::k_GpfBase(QStringList ak_IndexFiles)
{
	// load index files
	foreach (QString ls_Filename, ak_IndexFiles)
	{
		RefPtr<k_GpfIndexFileInfo> lk_pIndexFile(new k_GpfIndexFileInfo(*this, ls_Filename));
		mk_IndexFiles[lk_pIndexFile->get_Key()] = lk_pIndexFile;
	}

	if (mk_IndexFiles.empty())
	{
		printf(QString("Error: No GPF index files loaded.\n").toStdString().c_str());
		exit(1);
	}

	FillTables();
}


k_GpfBase::~k_GpfBase()
{
}


RefPtr<k_GpfIndexFileInfo> k_GpfBase::get_IndexFileInfo(QString as_Id)
{
	if (mk_IndexFiles.contains(as_Id))
		return mk_IndexFiles[as_Id];
	else
		return RefPtr<k_GpfIndexFileInfo>(NULL);
}


RefPtr<k_GpfIndexFileInfo> k_GpfBase::get_IndexFileInfoFromAssembly(QString as_Assembly)
{
	QString ls_Genome = as_Assembly.split(QChar(';')).first();
	return get_IndexFileInfo(ls_Genome);
}


RefPtr<k_GpfIndexFileInfo> k_GpfBase::get_DefaultIndexFileInfo()
{
	return mk_IndexFiles[mk_IndexFiles.keys().first()];
}


QString k_GpfBase::get_GpfParameters() const
{
	return ms_GpfParameters;
}


unsigned int k_GpfBase::CalculatePeptideMass(QString as_Peptide)
{
	unsigned int lui_Mass = mui_WaterMass;
	for (int i = 0; i < as_Peptide.length(); ++i)
		lui_Mass += mui_AminoAcidWeight_[(unsigned char)me_CharToAminoAcid_[(int)as_Peptide.at(i).toAscii()]];

	return lui_Mass;
}


QString k_GpfBase::CollapsePeptide(QString as_Peptide)
{
	QString ls_Collapsed;
	for (int i = 0; i < as_Peptide.length(); ++i)
		ls_Collapsed += QChar(mc_AminoAcidToChar_[(unsigned char)me_AminoAcidMassSimilarityCollapse[me_CharToAminoAcid_[(unsigned char)as_Peptide[i].toAscii()]]]);

	return ls_Collapsed;
}


void k_GpfBase::FillTables()
{
	Q_INIT_RESOURCE(libgpf);

	// load GPF parameters and inject available genomes
	QFile lk_File(":res/GpfParameters.yaml");
	lk_File.open(QFile::ReadOnly);
	ms_GpfParameters = QString(lk_File.readAll());
	lk_File.close();

	QString ls_Inject = "  GenomeOption: &Genome\n";
	foreach (QString ls_Key, mk_IndexFiles.keys())
	{
		QString ls_GenomeTitle = mk_IndexFiles[ls_Key]->get_Title();
		if (ls_GenomeTitle.size() == 0)
			ls_Inject += "    - " + ls_Key + "\n";
		else
			ls_Inject += "    - " + ls_Key + ": '" + ls_GenomeTitle + "'\n";
	}
	ms_GpfParameters.insert(ms_GpfParameters.indexOf("enums:\n") + 7, ls_Inject);

	if (!mk_IndexFiles.empty())
	{
		ls_Inject = QString("  - id: genome\n    type: *Genome\n    label: Genome\n    default: %1\n    description: The genome that should be used for the search\n").
			arg(mk_IndexFiles.keys().first());
		ms_GpfParameters.insert(ms_GpfParameters.indexOf("parameters:\n") + 12, ls_Inject);
	}

	// clear tables
	for (int i = 0; i < 256; ++i)
		me_CharToAminoAcid_[i] = r_AminoAcid::Unknown;
	memset(mc_AminoAcidToChar_, 'X', 256);
	memset(md_AminoAcidWeight_, 0, sizeof(double) * 256);
	memset(mc_TripletToAminoAcidForward_, r_AminoAcid::Unknown, 512);
	memset(mc_TripletToAminoAcidBackward_, r_AminoAcid::Unknown, 512);

	// fill nucleotide table
	memset(mc_NucleotideCharToInt_, -1, 256);
	mc_NucleotideCharToInt_['A'] = 0;
	mc_NucleotideCharToInt_['C'] = 1;
	mc_NucleotideCharToInt_['G'] = 2;
	mc_NucleotideCharToInt_['T'] = 3;

	memset(mc_NucleotideIntToChar_, '.', 256);
	mc_NucleotideIntToChar_[0] = 'A';
	mc_NucleotideIntToChar_[1] = 'C';
	mc_NucleotideIntToChar_[2] = 'G';
	mc_NucleotideIntToChar_[3] = 'T';

	// read amino acid information
	char lc_Buffer_[1024];

	lk_File.setFileName(":res/AminoAcids.csv");
	lk_File.open(QFile::ReadOnly);

	forever
	{
		qint64 li_Length = lk_File.readLine(lc_Buffer_, sizeof(lc_Buffer_));
		if (li_Length == -1)
			break;
		if (li_Length == 0 || lc_Buffer_[0] == '#')
			continue;

		QString ls_Line(lc_Buffer_);
		ls_Line = ls_Line.trimmed();
		QStringList lk_List = ls_Line.split(QChar(';'));
		r_AminoAcid::Enumeration le_AminoAcid = (r_AminoAcid::Enumeration)lk_List[0].toInt();
		mc_AminoAcidToChar_[le_AminoAcid] = lk_List[3][0].toAscii();
		me_CharToAminoAcid_[(unsigned char)lk_List[3][0].toAscii()] = le_AminoAcid;
		md_AminoAcidWeight_[le_AminoAcid] = lk_List[4].toDouble();
	}
	lk_File.close();

	// determine unsigned int amino acid masses
	for (int i = 0; i < 256; ++i)
		mui_AminoAcidWeight_[i] = (unsigned int)(md_AminoAcidWeight_[i] * gui_MassPrecision);

	// DNA triplet to amino acid translation
	QFile lk_TripletsFile(":res/DnaToAminoAcid.csv");
	lk_TripletsFile.open(QFile::ReadOnly);

	forever
	{
		qint64 li_Length = lk_TripletsFile.readLine(lc_Buffer_, sizeof(lc_Buffer_));
		if (li_Length == -1)
			break;
		if (li_Length == 0)
			continue;

		QString ls_Line(lc_Buffer_);
		ls_Line = ls_Line.trimmed();
		ls_Line.replace('U', 'T');
		QStringList lk_List = ls_Line.split(QChar(';'));
		QString ls_Triplet = lk_List[0];
		QString ls_StrictAminoAcid = lk_List.last();
		QString ls_ProbableAminoAcid = "";
		if (lk_List.size() > 2)
			ls_ProbableAminoAcid = lk_List[1];

		unsigned int lui_Triplet = 0;
		for (int i = 0; i < 3; ++i)
			lui_Triplet = (lui_Triplet << 3) | (mc_NucleotideCharToInt_[(unsigned char)ls_Triplet[i].toAscii()] & 7);
		mc_TripletToAminoAcidForward_[lui_Triplet] = me_CharToAminoAcid_[(unsigned char)ls_StrictAminoAcid[0].toAscii()];

		lui_Triplet = 0;
		for (int i = 0; i < 3; ++i)
			lui_Triplet = (lui_Triplet << 3) | (mc_NucleotideCharToInt_[(unsigned char)ls_Triplet[2 - i].toAscii()] & 7);
		mc_TripletToAminoAcidBackward_[lui_Triplet ^ gui_InvertTripletMask] = me_CharToAminoAcid_[(unsigned char)ls_StrictAminoAcid[0].toAscii()];
	}

	// fill mass similarity collapse table: I becomes L, and K becomes Q
	for (int i = 0; i < 256; ++i)
		me_AminoAcidMassSimilarityCollapse[i] = (r_AminoAcid::Enumeration)i;
	me_AminoAcidMassSimilarityCollapse[r_AminoAcid::Ile] = r_AminoAcid::Leu;
	me_AminoAcidMassSimilarityCollapse[r_AminoAcid::Lys] = r_AminoAcid::Gln;

	mui_WaterMass = (unsigned int)(gd_WaterMass * gui_MassPrecision);
}


bool lt_Gno(unsigned int* aui_First_, unsigned int* aui_Second_)
{
	return (*aui_First_ & ~gui_GlobalNucleotideOffsetBackwardFlag) < (*aui_Second_ & ~gui_GlobalNucleotideOffsetBackwardFlag);
}


bool lt_Direction_Gno(unsigned int* aui_First_, unsigned int* aui_Second_)
{
	// forward frames come first
	if ((*aui_First_ & gui_GlobalNucleotideOffsetBackwardFlag) != 0 && (*aui_Second_ & gui_GlobalNucleotideOffsetBackwardFlag) == 0)
		return false;

	if ((*aui_First_ & gui_GlobalNucleotideOffsetBackwardFlag) == 0 && (*aui_Second_ & gui_GlobalNucleotideOffsetBackwardFlag) != 0)
		return true;

	if ((*aui_First_ & gui_GlobalNucleotideOffsetBackwardFlag) == 0)
		return (*aui_First_ & ~gui_GlobalNucleotideOffsetBackwardFlag) < (*aui_Second_ & ~gui_GlobalNucleotideOffsetBackwardFlag);
	else
		return (*aui_Second_ & ~gui_GlobalNucleotideOffsetBackwardFlag) < (*aui_First_ & ~gui_GlobalNucleotideOffsetBackwardFlag);
}


bool lt_IndexEntry_Trimer_LeftMass_Gno(r_IndexEntry* ar_First_, r_IndexEntry* ar_Second_)
{
	if (ar_First_->muw_Trimer > ar_Second_->muw_Trimer)
		return false;

	if (ar_First_->muw_Trimer < ar_Second_->muw_Trimer)
		return true;

	if (ar_First_->mui_LeftMass > ar_Second_->mui_LeftMass)
		return false;

	if (ar_First_->mui_LeftMass < ar_Second_->mui_LeftMass)
		return true;

	return ((ar_First_->mui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag) < (ar_Second_->mui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag));
}


bool lt_IndexEntry_Trimer_RightMass_Gno(r_IndexEntry* ar_First_, r_IndexEntry* ar_Second_)
{
	if (ar_First_->muw_Trimer > ar_Second_->muw_Trimer)
		return false;

	if (ar_First_->muw_Trimer < ar_Second_->muw_Trimer)
		return true;

	if (ar_First_->mui_RightMass > ar_Second_->mui_RightMass)
		return false;

	if (ar_First_->mui_RightMass < ar_Second_->mui_RightMass)
		return true;

	return ((ar_First_->mui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag) < (ar_Second_->mui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag));
}


bool lt_PeptideSpan_Gno(r_PeptideSpan* ar_First_, r_PeptideSpan* ar_Second_)
{
	return ar_First_->mui_Gno < ar_Second_->mui_Gno;
}
