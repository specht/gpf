#pragma once

#include <QtCore>
#include "GpfBase.h"


class k_GpfIndexFileInfo
{
public:
	// use the following constructor for loading a GPF index file
	k_GpfIndexFileInfo(k_GpfBase& ak_GpfBase, QString as_Filename);
	// use the following ctor for compiling a new index
	k_GpfIndexFileInfo(k_GpfBase& ak_GpfBase);

	virtual ~k_GpfIndexFileInfo();

	QString get_Filename() const;
	QString get_Key() const;
	QString get_Title() const;
	int get_ScaffoldIndexForGno(unsigned int aui_Gno) const;
	void FindPeptide(unsigned int& aui_PeptideStartGno, unsigned int& aui_PeptideLength, unsigned int aui_Gno, unsigned int aui_ReadingFrame) const;
	QList<QPair<unsigned int, unsigned int> > CreateSpansOfEqualScaffoldAndDirection(unsigned int* aui_Entries_, unsigned int aui_Size);
	QString DecodeGno(unsigned int aui_Gno);
	QString Browse(unsigned int aui_StartPosition, unsigned int aui_Length);
	QString get_AssemblyInfoAsYaml(QString as_Assembly);
	RefPtr<unsigned char> readDnaSpanAsUnsignedChar(unsigned int aui_Gno, unsigned int aui_Length, unsigned int &aui_LeftBorder, unsigned int &aui_RightBorder, unsigned int &aui_ResultLength, QFile* ak_File_ = NULL);
	QString readDnaSpan(unsigned int aui_Gno, unsigned int aui_Length, unsigned int &aui_LeftBorder, unsigned int &aui_RightBorder, QFile* ak_File_ = NULL);
	void readAssembly(QString as_Assembly, QString& as_Peptide, QString& as_AminoAcidsLeft, QString& as_AminoAcidsRight);
	QString readAssemblyAsYaml(QString as_Assembly);

	// index file-wide constants and tables, without accessor methods to increase
	// execution speed - please don't change these values in other classes!
	QVector<r_Scaffold> mk_Scaffolds;
	unsigned int mui_TotalNucleotideCount;
	unsigned int mui_TotalIndexedTrimerCount;
	qint64 mi_LeftTagMassesFilePosition;
	qint64 mi_LeftTagGnoFilePosition;
	qint64 mi_RightTagMassesFilePosition;
	qint64 mi_RightTagGnoFilePosition;
	RefPtr<unsigned int> mui_pTrimerOffsetTable;
	RefPtr<unsigned int> mui_pTrimerCountTable;
	RefPtr<r_PeptideSpan> mr_PeptideSpanTables_[6];
	unsigned int mui_PeptideSpanTablesSize_[6];
	unsigned int mui_MaxScaffoldIdLength;

private:
	k_GpfBase& mk_GpfBase;
	QString ms_Filename;
	QString ms_Key;
	QString ms_Title;
};
