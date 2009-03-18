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

#include "GpfIndexFileInfo.h"
#include "GpfIndexHeader.h"
#include "StopWatch.h"
#include "Sorter.h"


k_GpfIndexFileInfo::k_GpfIndexFileInfo(k_GpfBase& ak_GpfBase, QString as_Filename)
	: mui_TotalNucleotideCount(0)
	, mui_TotalIndexedTrimerCount(0)
	, mi_LeftTagMassesFilePosition(0)
	, mi_LeftTagGnoFilePosition(0)
	, mi_RightTagMassesFilePosition(0)
	, mi_RightTagGnoFilePosition(0)
	, mui_pTrimerOffsetTable(RefPtr<unsigned int>(new unsigned int[8000]))
	, mui_pTrimerCountTable(RefPtr<unsigned int>(new unsigned int[8000]))
	, mui_MaxScaffoldIdLength(0)
	, mk_GpfBase(ak_GpfBase)
	, ms_Filename(as_Filename)
	, ms_Key("n/a")
{
	ms_Key = ms_Filename.left(ms_Filename.indexOf(QChar('.')));
	if (ms_Key.contains(QChar('/')))
		ms_Key = ms_Key.right(ms_Key.length() - ms_Key.lastIndexOf(QChar('/')) - 1);

	// load index file
	//k_StopWatch lk_StopWatch(QString("Parsing %1 took ").arg(ms_Key) + "%1.\n");

	QFile lk_IndexFile(as_Filename);
	if (!lk_IndexFile.open(QFile::ReadOnly | QFile::Unbuffered))
	{
		printf("Unable to open input file.\n");
		exit(1);
	}

	char lc_GpfDnaTag_[7];
	lc_GpfDnaTag_[6] = 0;
	lk_IndexFile.read(lc_GpfDnaTag_, 6);

	if (QString(lc_GpfDnaTag_) != "gpfdna")
	{
		printf("Invalid GPFDNA tag in input file.\n");
		exit(1);
	}

	quint16 luw_VersionNumberHigh;
	quint16 luw_VersionNumberLow;
	lk_IndexFile.read((char*)&luw_VersionNumberHigh, sizeof(luw_VersionNumberHigh));
	lk_IndexFile.read((char*)&luw_VersionNumberLow, sizeof(luw_VersionNumberLow));
	if (luw_VersionNumberHigh != 0 || luw_VersionNumberLow != 1)
	{
		printf("Invalid version number in input file.\n");
		exit(1);
	}

	while (!lk_IndexFile.atEnd())
	{
		qint64 li_MarkerPosition = lk_IndexFile.pos();
		r_Header lr_Header(lk_IndexFile, r_Marker::Undefined);
		lr_Header.Read();
		lk_IndexFile.seek(li_MarkerPosition);
		if (lr_Header.me_Type == r_Marker::GenomeTitle)
		{
			r_GenomeTitleHeader lr_GenomeTitleHeader(lk_IndexFile);
			lr_GenomeTitleHeader.Read();
			unsigned int lui_MaxLabelLength = 1024;
			RefPtr<char> lr_pLabel(new char[lui_MaxLabelLength]);
			quint16 lui_Length;
			lk_IndexFile.read((char*)&lui_Length, sizeof(lui_Length));
			lk_IndexFile.read(lr_pLabel.get_Pointer(), lui_Length);
			lr_pLabel.get_Pointer()[lui_Length] = 0;
			ms_Title = QString(lr_pLabel.get_Pointer());
			//lk_StopWatch.setExitMessage(QString("Parsing %1 (%2) took ").arg(ms_Title).arg(ms_Key) + "%1.\n");
		} 
		else if (lr_Header.me_Type == r_Marker::Dna)
		{
			r_DnaHeader lr_DnaHeader(lk_IndexFile);
			lr_DnaHeader.Read();

			for (unsigned int lui_ScaffoldIndex = 0; lui_ScaffoldIndex < lr_DnaHeader.mui_ScaffoldCount; ++lui_ScaffoldIndex)
			{
				qint64 li_MarkerPosition = lk_IndexFile.pos();
				r_Header lr_Header(lk_IndexFile, r_Marker::Undefined);
				lr_Header.Read();
				lk_IndexFile.seek(li_MarkerPosition);
				if (lr_Header.me_Type == r_Marker::DnaScaffold)
				{
					r_DnaScaffoldHeader lr_DnaScaffoldHeader(lk_IndexFile);
					lr_DnaScaffoldHeader.Read();
					r_Scaffold lr_Scaffold(QString(), lk_IndexFile.pos());
					lr_Scaffold.mui_NucleotideCount = lr_DnaScaffoldHeader.mui_NucleotideCount;
					mui_TotalNucleotideCount += lr_DnaScaffoldHeader.mui_NucleotideCount;
					if (lui_ScaffoldIndex == 0)
						lr_Scaffold.mui_NucleotideOffset = 0;
					else
						lr_Scaffold.mui_NucleotideOffset = mk_Scaffolds.last().mui_NucleotideOffset + mk_Scaffolds.last().mui_NucleotideCount;
					mk_Scaffolds.append(lr_Scaffold);
				}
				lk_IndexFile.seek(lr_Header.mi_NextEntry);
			}
		}
		else if (lr_Header.me_Type == r_Marker::ScaffoldLabels)
		{
			r_ScaffoldLabelsHeader lr_ScaffoldLabelsHeader(lk_IndexFile);
			lr_ScaffoldLabelsHeader.Read();
			unsigned int lui_MaxLabelLength = 1024;
			RefPtr<char> lr_pLabel(new char[lui_MaxLabelLength]);
			for (int i = 0; i < mk_Scaffolds.size(); ++i)
			{
				quint16 lui_Length;
				lk_IndexFile.read((char*)&lui_Length, sizeof(lui_Length));
				lk_IndexFile.read(lr_pLabel.get_Pointer(), lui_Length);
				lr_pLabel.get_Pointer()[lui_Length] = 0;
				mk_Scaffolds[i].ms_Id = QString(lr_pLabel.get_Pointer());
				mui_MaxScaffoldIdLength = std::max<unsigned int>(mui_MaxScaffoldIdLength, lui_Length);
			}
		}
		else if (lr_Header.me_Type == r_Marker::AminoAcids)
		{
			r_AminoAcidsHeader lr_AminoAcidsHeader(lk_IndexFile);
			lr_AminoAcidsHeader.Read();
			for (unsigned int lui_ScaffoldIndex = 0; lui_ScaffoldIndex < (unsigned int)mk_Scaffolds.size(); ++lui_ScaffoldIndex)
			{
				r_AminoAcidsScaffoldHeader lr_AminoAcidsScaffoldHeader(lk_IndexFile);
				lr_AminoAcidsScaffoldHeader.Read();
				for (unsigned int lui_FrameIndex = 0; lui_FrameIndex < 6; ++lui_FrameIndex)
				{
					r_AminoAcidsReadingFrameHeader lr_AminoAcidsReadingFrameHeader(lk_IndexFile);
					lr_AminoAcidsReadingFrameHeader.Read();
					mk_Scaffolds[lui_ScaffoldIndex].mr_ReadingFrames[lui_FrameIndex].mui_AminoAcidCount = lr_AminoAcidsReadingFrameHeader.mui_AminoAcidCount;
					mk_Scaffolds[lui_ScaffoldIndex].mr_ReadingFrames[lui_FrameIndex].mi_ReadingFrameFilePosition = lk_IndexFile.pos();
					mk_Scaffolds[lui_ScaffoldIndex].mui_TotalAminoAcidCount += lr_AminoAcidsReadingFrameHeader.mui_AminoAcidCount;
					lk_IndexFile.seek(lr_AminoAcidsReadingFrameHeader.mi_NextEntry);
				}
				lk_IndexFile.seek(lr_AminoAcidsScaffoldHeader.mi_NextEntry);
			}
		}
		else if (lr_Header.me_Type == r_Marker::PeptideSpanList)
		{
			r_PeptideSpanListHeader lr_PeptideSpanListHeader(lk_IndexFile);
			lr_PeptideSpanListHeader.Read();
			for (unsigned int lui_ReadingFrame = 0; lui_ReadingFrame < 6; ++lui_ReadingFrame)
			{
				r_PeptideSpanListReadingFrameHeader lr_PeptideSpanListReadingFrameHeader(lk_IndexFile);
				lr_PeptideSpanListReadingFrameHeader.Read();
				mui_PeptideSpanTablesSize_[lui_ReadingFrame] = lr_PeptideSpanListReadingFrameHeader.mui_PeptideCount;
				mr_PeptideSpanTables_[lui_ReadingFrame] = RefPtr<r_PeptideSpan>(new r_PeptideSpan[lr_PeptideSpanListReadingFrameHeader.mui_PeptideCount]);
				lk_IndexFile.read((char*)mr_PeptideSpanTables_[lui_ReadingFrame].get_Pointer(), sizeof(r_PeptideSpan) * lr_PeptideSpanListReadingFrameHeader.mui_PeptideCount);
				lk_IndexFile.seek(lr_PeptideSpanListReadingFrameHeader.mi_NextEntry);
			}
		}
		else if (lr_Header.me_Type == r_Marker::Index)
		{
			r_IndexHeader lr_IndexHeader(lk_IndexFile);
			lr_IndexHeader.Read();
			mui_TotalIndexedTrimerCount = lr_IndexHeader.mui_TotalIndexedTrimerCount;

			r_IndexOffsetsHeader lr_IndexOffsetsHeader(lk_IndexFile);
			lr_IndexOffsetsHeader.Read();
			lk_IndexFile.read((char*)mui_pTrimerOffsetTable.get_Pointer(), sizeof(unsigned int) * 8000);

			// derive trimer count table from trimer offset table
			for (unsigned int i = 0; i < 7999; ++i)
				mui_pTrimerCountTable.get_Pointer()[i] = mui_pTrimerOffsetTable.get_Pointer()[i + 1] - mui_pTrimerOffsetTable.get_Pointer()[i];
			mui_pTrimerCountTable.get_Pointer()[7999] = mui_TotalIndexedTrimerCount - mui_pTrimerOffsetTable.get_Pointer()[7999];

			// fetch four index table offsets
			r_IndexLeftTagMassesHeader lr_IndexLeftTagMassesHeader(lk_IndexFile);
			lr_IndexLeftTagMassesHeader.Read();
			if (lr_IndexLeftTagMassesHeader.me_Type != r_Marker::IndexLeftTagMasses)
			{
				printf("Error reading left tag mass table.\n");
				exit(1);
			}
			mi_LeftTagMassesFilePosition = lk_IndexFile.pos();
			lk_IndexFile.seek(lr_IndexLeftTagMassesHeader.mi_NextEntry);

			r_IndexLeftTagGnoHeader lr_IndexLeftTagGnoHeader(lk_IndexFile);
			lr_IndexLeftTagGnoHeader.Read();
			if (lr_IndexLeftTagGnoHeader.me_Type != r_Marker::IndexLeftTagGno)
			{
				printf("Error reading left tag GNO table.\n");
				exit(1);
			}

			mi_LeftTagGnoFilePosition = lk_IndexFile.pos();
			lk_IndexFile.seek(lr_IndexLeftTagGnoHeader.mi_NextEntry);

			r_IndexRightTagMassesHeader lr_IndexRightTagMassesHeader(lk_IndexFile);
			lr_IndexRightTagMassesHeader.Read();
			if (lr_IndexRightTagMassesHeader.me_Type != r_Marker::IndexRightTagMasses)
			{
				printf("Error reading right tag mass table.\n");
				exit(1);
			}
			mi_RightTagMassesFilePosition = lk_IndexFile.pos();
			lk_IndexFile.seek(lr_IndexRightTagMassesHeader.mi_NextEntry);

			r_IndexRightTagGnoHeader lr_IndexRightTagGnoHeader(lk_IndexFile);
			lr_IndexRightTagGnoHeader.Read();
			if (lr_IndexRightTagGnoHeader.me_Type != r_Marker::IndexRightTagGno)
			{
				printf("Error reading right tag GNO table.\n");
				exit(1);
			}
			mi_RightTagGnoFilePosition = lk_IndexFile.pos();
		}
		lk_IndexFile.seek(lr_Header.mi_NextEntry);
	}
}


k_GpfIndexFileInfo::k_GpfIndexFileInfo(k_GpfBase& ak_GpfBase)
	: mui_TotalNucleotideCount(0)
	, mui_TotalIndexedTrimerCount(0)
	, mi_LeftTagMassesFilePosition(0)
	, mi_LeftTagGnoFilePosition(0)
	, mi_RightTagMassesFilePosition(0)
	, mi_RightTagGnoFilePosition(0)
	, mui_pTrimerOffsetTable(RefPtr<unsigned int>(new unsigned int[8000]))
	, mui_pTrimerCountTable(RefPtr<unsigned int>(new unsigned int[8000]))
	, mui_MaxScaffoldIdLength(0)
	, mk_GpfBase(ak_GpfBase)
	, ms_Filename("n/a")
	, ms_Key("n/a")
{
}


k_GpfIndexFileInfo::~k_GpfIndexFileInfo()
{
}


QString k_GpfIndexFileInfo::get_Filename() const
{
	return ms_Filename;
}


QString k_GpfIndexFileInfo::get_Key() const
{
	return ms_Key;
}


QString k_GpfIndexFileInfo::get_Title() const
{
	return ms_Title;
}


int k_GpfIndexFileInfo::get_ScaffoldIndexForGno(unsigned int aui_Gno) const
{
	// mask out reverse flag in global nucleotide offset
	unsigned int lui_Offset = aui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag;

	if (lui_Offset >= mui_TotalNucleotideCount)
	{
		printf("Attention: invalid gno (%08x).\n", lui_Offset);
		return 0;
	}

	// TODO: change to log(n) time instead of n time
	int li_ScaffoldIndex = 0;
	while (li_ScaffoldIndex < mk_Scaffolds.size() && lui_Offset >= mk_Scaffolds[li_ScaffoldIndex].mui_NucleotideOffset + mk_Scaffolds[li_ScaffoldIndex].mui_NucleotideCount)
		++li_ScaffoldIndex;

	return li_ScaffoldIndex;
}


void k_GpfIndexFileInfo::FindPeptide(unsigned int& aui_PeptideStartGno, unsigned int& aui_PeptideLength, unsigned int aui_Gno, unsigned int aui_ReadingFrame) const
{
	unsigned int lui_First = 0;
	unsigned int lui_Last = mui_PeptideSpanTablesSize_[aui_ReadingFrame] - 1;
	r_PeptideSpan* lr_PeptideSpanList_ = mr_PeptideSpanTables_[aui_ReadingFrame].get_Pointer();

	if (aui_ReadingFrame < 3)
	{
		while (lui_Last > lui_First)
		{
			unsigned int lui_Mid = lui_First + (lui_Last - lui_First) / 2;
			if (aui_Gno < lr_PeptideSpanList_[lui_Mid].mui_Gno)
				// cut away right half
				lui_Last = lui_Mid - 1;
			else
			{
				if (aui_Gno >= lr_PeptideSpanList_[lui_Mid].mui_Gno && 
					aui_Gno < lr_PeptideSpanList_[lui_Mid].mui_Gno + lr_PeptideSpanList_[lui_Mid].muw_Length * 3)
				{
					// we have found the peptide span to which the GNO belongs
					lui_First = lui_Mid;
					break;
				}
				else
					// cut away left half
					lui_First = lui_Mid + 1;
			}
		}
	}
	else
	{
		while (lui_Last > lui_First)
		{
			unsigned int lui_Mid = lui_First + (lui_Last - lui_First) / 2;
			if (aui_Gno > lr_PeptideSpanList_[lui_Mid].mui_Gno)
				// cut away left half
				lui_First = lui_Mid + 1;
			else
			{
				if (aui_Gno <= lr_PeptideSpanList_[lui_Mid].mui_Gno && 
					aui_Gno > lr_PeptideSpanList_[lui_Mid].mui_Gno - lr_PeptideSpanList_[lui_Mid].muw_Length * 3)
				{
					// we have found the peptide span to which the GNO belongs
					lui_First = lui_Mid;
					break;
				}
				else
					// cut away right half
					lui_Last = lui_Mid - 1;
			}
		}
	}

	aui_PeptideStartGno = lr_PeptideSpanList_[lui_First].mui_Gno;
	aui_PeptideLength = lr_PeptideSpanList_[lui_First].muw_Length;
}


QList<QPair<unsigned int, unsigned int> > 
k_GpfIndexFileInfo::CreateSpansOfEqualScaffoldAndDirection(unsigned int* aui_Entries_, unsigned int aui_Size)
{
	// scan all entries and group by scaffold, then sort each scaffold according to reading direction and, again, offset
	unsigned int lui_ScaffoldStart = 0;
	unsigned int lui_CurrentScaffoldIndex = this->get_ScaffoldIndexForGno(*aui_Entries_);
	unsigned int lui_ForwardEntryCount = 0;
	QList<QPair<unsigned int, unsigned int> > lk_ScaffoldList;
	for (unsigned int i = 0; i < aui_Size; ++i)
	{
		if ((aui_Entries_[i] & gui_GlobalNucleotideOffsetBackwardFlag) == 0)
			++lui_ForwardEntryCount;

		if ((i == aui_Size - 1) || ((aui_Entries_[i + 1] & ~gui_GlobalNucleotideOffsetBackwardFlag) >= 
			mk_Scaffolds[lui_CurrentScaffoldIndex].mui_NucleotideOffset + mk_Scaffolds[lui_CurrentScaffoldIndex].mui_NucleotideCount))
		{
			// This is the last semi tag entry for the current scaffold,
			// sort the entries in this scaffold according to reading direction and GNO.
			SortUnsignedIntDirectionGno(aui_Entries_, lui_ScaffoldStart, i);

			unsigned int lui_BackwardEntryCount = i - lui_ScaffoldStart + 1 - lui_ForwardEntryCount;

			if (lui_ForwardEntryCount > 0)
				lk_ScaffoldList.append(QPair<unsigned int, unsigned int>(lui_CurrentScaffoldIndex * 2, lui_ScaffoldStart));

			if (lui_BackwardEntryCount > 0)
				lk_ScaffoldList.append(QPair<unsigned int, unsigned int>(
				lui_CurrentScaffoldIndex * 2 + 1, lui_ScaffoldStart + lui_ForwardEntryCount));

			if (i < aui_Size - 1)
			{
				// adjust values for next span
				lui_ForwardEntryCount = 0;
				lui_ScaffoldStart = i + 1;
				while ((aui_Entries_[i + 1] & ~gui_GlobalNucleotideOffsetBackwardFlag) >= 
					mk_Scaffolds[lui_CurrentScaffoldIndex].mui_NucleotideOffset + mk_Scaffolds[lui_CurrentScaffoldIndex].mui_NucleotideCount)
					++lui_CurrentScaffoldIndex;
			}
		}
	}

	return lk_ScaffoldList;
}


QString k_GpfIndexFileInfo::DecodeGno(unsigned int aui_Gno)
{
	QString ls_Result;
	unsigned int lui_ScaffoldIndex = get_ScaffoldIndexForGno(aui_Gno);
	r_Scaffold& lr_Scaffold = mk_Scaffolds[lui_ScaffoldIndex];
	ls_Result += QString("scaffold: '%1'\n").arg(lr_Scaffold.ms_Id);
	ls_Result += QString("position: %1\n").arg((aui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag) - lr_Scaffold.mui_NucleotideOffset);
	ls_Result += QString("frame: %1\n").arg(lr_Scaffold.getReadingFrameForGlobalNucleotideOffset(aui_Gno));
	return ls_Result;
}


QString k_GpfIndexFileInfo::Browse(unsigned int aui_StartPosition, unsigned int aui_Length)
{
	if (aui_StartPosition >= mui_TotalNucleotideCount)
	{
		// error
		return QString();
	}

	// cut length if it exceeds the DNA
	if (aui_StartPosition + aui_Length > mui_TotalNucleotideCount)
		aui_Length = mui_TotalNucleotideCount - aui_StartPosition;

	int li_ScaffoldIndex = get_ScaffoldIndexForGno(aui_StartPosition);
	r_Scaffold& lr_Scaffold = mk_Scaffolds[li_ScaffoldIndex];
	RefPtr<char> lc_pDna(new char[aui_Length]);

	qint64 li_Position = lr_Scaffold.mi_DnaFilePosition + aui_StartPosition - lr_Scaffold.mui_NucleotideOffset;

	QFile lk_File(ms_Filename);
	lk_File.open(QIODevice::ReadOnly | QIODevice::Unbuffered);
	lk_File.seek(li_Position);
	lk_File.read(lc_pDna.get_Pointer(), aui_Length);
	lk_File.close();

	QString ls_Dna;

	for (int i = 0; i < (int)aui_Length; ++i)
		ls_Dna.append(QChar(mk_GpfBase.mc_NucleotideIntToChar_[(unsigned char)lc_pDna.get_Pointer()[i]]));

	QString ls_Result;
	ls_Result += QString("dna: '%1'\n").arg(ls_Dna);
	ls_Result += QString("scaffoldId: '%1'\n").arg(lr_Scaffold.ms_Id);
	ls_Result += QString("scaffoldStart: %1\n").arg(lr_Scaffold.mui_NucleotideOffset);
	ls_Result += QString("scaffoldLength: %1\n").arg(lr_Scaffold.mui_NucleotideCount);
	ls_Result += QString("scaffoldPosition: %1\n").arg(aui_StartPosition - lr_Scaffold.mui_NucleotideOffset);
	if (li_ScaffoldIndex > 0)
		ls_Result += QString("previousScaffold: '%1'\n").arg(mk_Scaffolds[li_ScaffoldIndex - 1].ms_Id);
	if (li_ScaffoldIndex < mk_Scaffolds.size() - 1)
		ls_Result += QString("nextScaffold: '%1'\n").arg(mk_Scaffolds[li_ScaffoldIndex + 1].ms_Id);

	return ls_Result;
}


QString k_GpfIndexFileInfo::get_AssemblyInfoAsYaml(QString as_Assembly)
{
	QStringList lk_Parts = as_Assembly.split(QChar(';'));
	QString ls_Organism = lk_Parts[0];
	QString ls_Contigs = lk_Parts[1];
	bool lb_BackwardAssembly = ls_Contigs.startsWith(QChar('-'));

	// strip direction indicator from contigs string
	ls_Contigs = ls_Contigs.mid(1);

	QString ls_AssemblyInfo = QString("{ organism: '%1', forward: %2")
		.arg(ls_Organism).arg(lb_BackwardAssembly? "false": "true");
		
	lk_Parts = ls_Contigs.split(",");

	if (lk_Parts.size() > 1)
	{
		// this hit assembly contains an intron, determine the intron length
		QStringList lk_Elements = lk_Parts[0].split(":");
		unsigned int lui_Pos0 = lk_Elements[0].toUInt();
		unsigned int lui_Length0 = lk_Elements[1].toUInt();
		lk_Elements = lk_Parts[1].split(":");
		unsigned int lui_Pos1 = lk_Elements[0].toUInt();
		unsigned int lui_Length1 = lk_Elements[1].toUInt();
		unsigned int lui_IntronLength = 0;
		if (!lb_BackwardAssembly)
			lui_IntronLength = lui_Pos1 - (lui_Pos0 + lui_Length0);
		else
			lui_IntronLength = (lui_Pos0 - lui_Length0) - lui_Pos1;
		ls_AssemblyInfo += QString(", intronLength: %1").arg(lui_IntronLength);
	}

	ls_AssemblyInfo += ", parts: [";
	bool lb_First = true;
	foreach (QString ls_Part, lk_Parts)
	{
		QStringList lk_Elements = ls_Part.split(":");
		unsigned int lui_Gno = lk_Elements[0].toUInt();
		unsigned int lui_Length = lk_Elements[1].toUInt();
		if (lb_BackwardAssembly)
			lui_Gno |= gui_GlobalNucleotideOffsetBackwardFlag;

		unsigned int lui_ScaffoldIndex = get_ScaffoldIndexForGno(lui_Gno);
		r_Scaffold& lr_Scaffold = mk_Scaffolds[lui_ScaffoldIndex];

		QString ls_PartInfo = QString("{scaffold: '%1', position: %2, length: %3, frame: %4}")
			.arg(lr_Scaffold.ms_Id)
			.arg((lui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag) - lr_Scaffold.mui_NucleotideOffset)
			.arg(lui_Length)
			.arg(lr_Scaffold.getReadingFrameForGlobalNucleotideOffset(lui_Gno));

		if (!lb_First)
			ls_AssemblyInfo += ", ";
		else
			lb_First = false;

		ls_AssemblyInfo += ls_PartInfo;
	}
	ls_AssemblyInfo += "] }";

	return ls_AssemblyInfo;
}


RefPtr<unsigned char> k_GpfIndexFileInfo::readDnaSpanAsUnsignedChar(unsigned int aui_Gno, unsigned int aui_Length, unsigned int &aui_LeftBorder, unsigned int &aui_RightBorder, unsigned int &aui_ResultLength, QFile* ak_File_)
{
	// This functions reads the nucleotides associated with the given GNO, and then some adjacent nucleotides.
	// The characters of the return string are in the range of ascii 0..3, and the string is already 
	// reversed and transposed if necessary, that is, if the assembly was a backwards assembly.

	RefPtr<QFile> lk_pFile;

	// open file if no file was passed in as an argument
	if (ak_File_ == NULL)
	{
		lk_pFile = RefPtr<QFile>(new QFile(ms_Filename));
		lk_pFile->open(QIODevice::ReadOnly | QIODevice::Unbuffered);
		ak_File_ = lk_pFile.get_Pointer();
	}

	bool lb_Forwards = (aui_Gno & gui_GlobalNucleotideOffsetBackwardFlag) == 0;

	// determine which scaffold were are talking about
	r_Scaffold& lr_Scaffold = mk_Scaffolds[get_ScaffoldIndexForGno(aui_Gno)];

	if (aui_LeftBorder > 0)
	{
		// fix left border if it exceeds the scaffold
		if (lb_Forwards)
		{
			// check towards start of scaffold
			qint64 li_RemainingLeftNucleotides = (qint64)aui_Gno - lr_Scaffold.mui_NucleotideOffset;
			li_RemainingLeftNucleotides = std::max<qint64>(li_RemainingLeftNucleotides, 0);
			aui_LeftBorder = (unsigned int)std::min<qint64>(li_RemainingLeftNucleotides, (qint64)aui_LeftBorder);
			aui_Gno -= aui_LeftBorder;
		}
		else
		{
			// check towards end of scaffold
			qint64 li_RemainingLeftNucleotides = (qint64)(lr_Scaffold.mui_NucleotideOffset + lr_Scaffold.mui_NucleotideCount) - ((aui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag) + 1);
			li_RemainingLeftNucleotides = std::max<qint64>(li_RemainingLeftNucleotides, 0);
			aui_LeftBorder = (unsigned int)std::min<qint64>(li_RemainingLeftNucleotides, (qint64)aui_LeftBorder);
			aui_Gno += aui_LeftBorder;
		}
		aui_Length += aui_LeftBorder;
	}


	if (aui_RightBorder > 0)
	{
		// fix right border if it exceeds the scaffold
		if (lb_Forwards)
		{
			// check towards end of scaffold
			qint64 li_RemainingRightNucleotides = (qint64)(lr_Scaffold.mui_NucleotideOffset + lr_Scaffold.mui_NucleotideCount) - (aui_Gno + aui_Length);
			li_RemainingRightNucleotides = std::max<qint64>(li_RemainingRightNucleotides, 0);
			aui_RightBorder = (unsigned int)std::min<qint64>(li_RemainingRightNucleotides, (qint64)aui_RightBorder);
		}
		else
		{
			// check towards start of scaffold
			qint64 li_RemainingRightNucleotides = (qint64)((aui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag) + 1 - aui_Length) - lr_Scaffold.mui_NucleotideOffset;
			li_RemainingRightNucleotides = std::max<qint64>(li_RemainingRightNucleotides, 0);
			aui_RightBorder = (unsigned int)std::min<qint64>(li_RemainingRightNucleotides, (qint64)aui_RightBorder);
		}
		aui_Length += aui_RightBorder;
	}

	RefPtr<unsigned char> luc_pDna(new unsigned char[aui_Length]);

	qint64 li_Position = lr_Scaffold.mi_DnaFilePosition + (aui_Gno & ~gui_GlobalNucleotideOffsetBackwardFlag) - lr_Scaffold.mui_NucleotideOffset;
	if (!lb_Forwards)
		li_Position -= aui_Length - 1;

	ak_File_->seek(li_Position);
	ak_File_->read((char*)luc_pDna.get_Pointer(), aui_Length);

	if (!lb_Forwards)
	{
		// transpose the DNA in place
		int li_Center = aui_Length / 2;
		for (int i = 0; i < li_Center; ++i)
		{
			unsigned char luc_Left = luc_pDna.get_Pointer()[i];
			unsigned char luc_Right = luc_pDna.get_Pointer()[aui_Length - i - 1];
			luc_pDna.get_Pointer()[i] = luc_Right ^ 3;
			luc_pDna.get_Pointer()[aui_Length - i - 1] = luc_Left ^ 3;
		}
		if (aui_Length % 2 != 0)
			luc_pDna.get_Pointer()[li_Center] ^= 3;
	}

	aui_ResultLength = aui_Length;
	return luc_pDna;
}


QString k_GpfIndexFileInfo::readDnaSpan(unsigned int aui_Gno, unsigned int aui_Length, unsigned int &aui_LeftBorder, unsigned int &aui_RightBorder, QFile* ak_File_)
{
	unsigned int lui_Length;
	RefPtr<unsigned char> luc_pDna = this->readDnaSpanAsUnsignedChar(aui_Gno, aui_Length, aui_LeftBorder, aui_RightBorder, lui_Length, ak_File_);

	QString ls_Dna;
	for (int i = 0; i < (int)lui_Length; ++i)
		ls_Dna.append(QChar((char)luc_pDna.get_Pointer()[i]));

	return ls_Dna;
}


void k_GpfIndexFileInfo::readAssembly(QString as_Assembly, QString& as_Peptide, QString& as_AminoAcidsLeft, QString& as_AminoAcidsRight)
{
	// read the peptide described by the given assembly, and read a few surrounding amino acids
	as_Peptide = "";
	as_AminoAcidsLeft = "[some stupid amino acids to the left]";
	as_AminoAcidsRight = "[some stupid amino acids to the right]";
	QFile lk_File(ms_Filename);
	lk_File.open(QIODevice::ReadOnly | QIODevice::Unbuffered);

	as_Assembly = as_Assembly.split(";")[1];
	bool lb_BackwardAssembly = as_Assembly.startsWith("-");
	as_Assembly = as_Assembly.mid(1);

	RefPtr<char> lc_pTempString(new char[1024]);
	memset(lc_pTempString.get_Pointer(), 0, 1024);
	QString ls_TempString;

	QStringList lk_AssemblyParts = as_Assembly.split(",");
	unsigned int lui_LeftBorder = gui_MaxSurroundingAminoAcidCount * 3;
	unsigned int lui_RightBorder = gui_MaxSurroundingAminoAcidCount * 3;

	QString ls_Dna;
	for (int i = 0; i < lk_AssemblyParts.size(); ++i)
	{
		QString ls_Part = lk_AssemblyParts[i];
		QStringList lk_Parts = ls_Part.split(":");
		unsigned int lui_Gno = lk_Parts[0].toUInt();
		unsigned int lui_Length = lk_Parts[1].toUInt();
		if (lb_BackwardAssembly)
			lui_Gno |= gui_GlobalNucleotideOffsetBackwardFlag;

		unsigned int lui_LocalLeftBorder = (i == 0)? lui_LeftBorder: 0;
		unsigned int lui_LocalRightBorder = (i == lk_AssemblyParts.size() - 1)? lui_RightBorder: 0;

		ls_Dna += this->readDnaSpan(lui_Gno, lui_Length, lui_LocalLeftBorder, lui_LocalRightBorder, &lk_File);

		if (i == 0)
			lui_LeftBorder = lui_LocalLeftBorder;
		if (i == lk_AssemblyParts.size() - 1)
			lui_RightBorder = lui_LocalRightBorder;
	}

	// remove nucleotides from the front until left border is divisble by three
	while ((lui_LeftBorder % 3) != 0)
	{
		--lui_LeftBorder;
		ls_Dna = ls_Dna.mid(1);
	}

	// adjust borders for amino acids
	lui_LeftBorder /= 3;
	lui_RightBorder /= 3;

	// translate assembled DNA into amino acids
	for (int i = 0; i < ls_Dna.length() / 3; ++i)
	{
		unsigned int lui_Triplet = 
			((ls_Dna.at(i * 3 + 0).toAscii() & 7) << 6) |
			((ls_Dna.at(i * 3 + 1).toAscii() & 7) << 3) |
			(ls_Dna.at(i * 3 + 2).toAscii() & 7);
		r_AminoAcid::Enumeration le_AminoAcid = (r_AminoAcid::Enumeration)mk_GpfBase.mc_TripletToAminoAcidForward_[lui_Triplet];
		as_Peptide.append(mk_GpfBase.mc_AminoAcidToChar_[(unsigned char)le_AminoAcid]);
	}

	as_AminoAcidsLeft = as_Peptide.left(lui_LeftBorder);
	as_AminoAcidsRight = as_Peptide.right(lui_RightBorder);
	as_Peptide = as_Peptide.mid(lui_LeftBorder, as_Peptide.length() - lui_RightBorder - lui_LeftBorder);

	lk_File.close();
}


QString k_GpfIndexFileInfo::readAssemblyAsYaml(QString as_Assembly)
{
	QString ls_Peptide, ls_AminoAcidsLeft, ls_AminoAcidsRight;
	this->readAssembly(as_Assembly, ls_Peptide, ls_AminoAcidsLeft, ls_AminoAcidsRight);
	return QString("peptide: %1\nleft: %2\nright: %3").arg(ls_Peptide).arg(ls_AminoAcidsLeft).arg(ls_AminoAcidsRight);
}
