#pragma once

#include <QtCore>


struct r_Marker
{
	enum Enumeration
	{
		Undefined = 0,
		Dna = 1,
		DnaScaffold = 2,
		ScaffoldLabels = 3,
		AminoAcids = 4,
		AminoAcidsScaffold = 5,
		AminoAcidsReadingFrame = 6,
		Index = 7,
		IndexOffsets = 8,
		IndexLeftTagMasses = 9,
		IndexLeftTagGno = 10,
		IndexRightTagMasses = 11,
		IndexRightTagGno = 12,
		PeptideSpanList = 13,
		PeptideSpanListReadingFrame = 14,
		GenomeTitle = 15
	};
};


#pragma pack(push)
#pragma pack(1)


struct r_Header
{
	r_Header(QFile& ak_File, r_Marker::Enumeration ae_Type)
		: mk_File(ak_File)
		, me_Type(ae_Type)
		, mi_NextEntry(0)
	{
	}

	virtual void WritePreliminary()
	{
		mk_File.write((const char*)(&me_Type), sizeof(me_Type));
		mi_HeaderPosition = mk_File.pos();
		mk_File.write((const char*)(&mi_NextEntry), sizeof(mi_NextEntry));
	}

	virtual void WriteFinal()
	{
		mi_NextEntry = mk_File.pos();
		mk_File.seek(mi_HeaderPosition);
		mk_File.write((const char*)(&mi_NextEntry), sizeof(mi_NextEntry));
	}

	virtual void WriteFinal(qint64 ai_NextEntry)
	{
		mi_NextEntry = mk_File.pos();
		mk_File.seek(mi_HeaderPosition);
		mk_File.write((const char*)(&ai_NextEntry), sizeof(ai_NextEntry));
	}

	virtual void Read()
	{
		mk_File.read((char*)(&me_Type), sizeof(me_Type));
		mk_File.read((char*)(&mi_NextEntry), sizeof(mi_NextEntry));
	}

	static int getSize()
	{
		return sizeof(r_Marker::Enumeration) + sizeof(qint64);
	}

	QFile& mk_File;
	qint64 mi_HeaderPosition;
	r_Marker::Enumeration me_Type;
	qint64 mi_NextEntry;
};


#define HEADER_0(classname, Marker) \
struct classname: public r_Header \
{ \
	classname(QFile& ak_File) \
		: r_Header(ak_File, Marker) \
	{ \
	} \
 \
	virtual void WritePreliminary() \
	{ \
		r_Header::WritePreliminary(); \
	} \
 \
	virtual void WriteFinal() \
	{ \
		r_Header::WriteFinal(); \
		mk_File.seek(mi_NextEntry); \
	} \
 \
	virtual void WriteFinal(qint64 ai_NextEntry) \
	{ \
		r_Header::WriteFinal(ai_NextEntry); \
		mk_File.seek(mi_NextEntry); \
	} \
\
	virtual void Read() \
	{ \
		r_Header::Read(); \
	} \
\
	static int getSize() \
	{ \
		return sizeof(r_Marker::Enumeration) + sizeof(qint64); \
	} \
};


#define HEADER_1(classname, Marker, type1, name1, default1) \
struct classname: public r_Header \
{ \
	classname(QFile& ak_File) \
		: r_Header(ak_File, Marker) \
		, name1(default1) \
	{ \
	} \
 \
	virtual void WritePreliminary() \
	{ \
		r_Header::WritePreliminary(); \
		mk_File.write((const char*)(&name1), sizeof(name1)); \
	} \
 \
	virtual void WriteFinal() \
	{ \
		r_Header::WriteFinal(); \
		mk_File.write((const char*)(&name1), sizeof(name1)); \
		mk_File.seek(mi_NextEntry); \
	} \
\
	virtual void WriteFinal(qint64 ai_NextEntry) \
	{ \
		r_Header::WriteFinal(ai_NextEntry); \
		mk_File.write((const char*)(&name1), sizeof(name1)); \
		mk_File.seek(mi_NextEntry); \
	} \
\
	virtual void Read() \
	{ \
		r_Header::Read(); \
		mk_File.read((char*)(&name1), sizeof(name1)); \
	} \
\
	static int getSize() \
	{ \
		return sizeof(r_Marker::Enumeration) + sizeof(qint64) + sizeof(type1); \
	} \
\
	type1 name1; \
};


#define HEADER_2(classname, Marker, type1, name1, default1, type2, name2, default2) \
struct classname: public r_Header \
{ \
	classname(QFile& ak_File) \
		: r_Header(ak_File, Marker) \
		, name1(default1) \
		, name2(default2) \
	{ \
	} \
 \
	virtual void WritePreliminary() \
	{ \
		r_Header::WritePreliminary(); \
		mk_File.write((const char*)(&name1), sizeof(name1)); \
		mk_File.write((const char*)(&name2), sizeof(name2)); \
	} \
 \
	virtual void WriteFinal() \
	{ \
		r_Header::WriteFinal(); \
		mk_File.write((const char*)(&name1), sizeof(name1)); \
		mk_File.write((const char*)(&name2), sizeof(name2)); \
		mk_File.seek(mi_NextEntry); \
	} \
 \
	virtual void WriteFinal(qint64 ai_NextEntry) \
	{ \
		r_Header::WriteFinal(ai_NextEntry); \
		mk_File.write((const char*)(&name1), sizeof(name1)); \
		mk_File.write((const char*)(&name2), sizeof(name2)); \
		mk_File.seek(mi_NextEntry); \
	} \
\
	virtual void Read() \
	{ \
		r_Header::Read(); \
		mk_File.read((char*)(&name1), sizeof(name1)); \
		mk_File.read((char*)(&name2), sizeof(name2)); \
	} \
\
	static int getSize() \
	{ \
		return sizeof(r_Marker::Enumeration) + sizeof(qint64) + sizeof(type1) + sizeof(type2); \
	} \
\
	type1 name1; \
	type2 name2; \
};


HEADER_1(r_DnaHeader, r_Marker::Dna, unsigned int, mui_ScaffoldCount, 0);
HEADER_1(r_DnaScaffoldHeader, r_Marker::DnaScaffold, unsigned int, mui_NucleotideCount, 0);
HEADER_0(r_ScaffoldLabelsHeader, r_Marker::ScaffoldLabels);
HEADER_0(r_AminoAcidsHeader, r_Marker::AminoAcids);
HEADER_0(r_AminoAcidsScaffoldHeader, r_Marker::AminoAcidsScaffold);
HEADER_2(r_AminoAcidsReadingFrameHeader, r_Marker::AminoAcidsReadingFrame, unsigned int, mui_FrameNumber, 0, unsigned int, mui_AminoAcidCount, 0);
HEADER_1(r_IndexHeader, r_Marker::Index, unsigned int, mui_TotalIndexedTrimerCount, 0);
HEADER_0(r_IndexOffsetsHeader, r_Marker::IndexOffsets);
HEADER_0(r_IndexLeftTagMassesHeader, r_Marker::IndexLeftTagMasses);
HEADER_0(r_IndexLeftTagGnoHeader, r_Marker::IndexLeftTagGno);
HEADER_0(r_IndexRightTagMassesHeader, r_Marker::IndexRightTagMasses);
HEADER_0(r_IndexRightTagGnoHeader, r_Marker::IndexRightTagGno);
HEADER_0(r_PeptideSpanListHeader, r_Marker::PeptideSpanList);
HEADER_1(r_PeptideSpanListReadingFrameHeader, r_Marker::PeptideSpanListReadingFrame, unsigned int, mui_PeptideCount, 0);
HEADER_0(r_GenomeTitleHeader, r_Marker::GenomeTitle);


#pragma pack(pop)
