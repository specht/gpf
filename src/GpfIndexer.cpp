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

#include "GpfIndexer.h"
#include "GpfBase.h"
#include "RefPtr.h"
#include "StopWatch.h"


k_GpfIndexer::k_GpfIndexer(QString as_DnaPath, QString as_DnaIndexPath, QString as_Title)
	: ms_DnaPath(as_DnaPath)
	, ms_DnaIndexPath(as_DnaIndexPath)
	, ms_Title(as_Title)
	, mi_TotalNucleotideCount(0)
	, mi_OffsetBits(0)
	, mui_GnoBackwardsBit(0)
	, mi_MassBits(27)
	, mi_MassPrecision(10000)
	, mi_TagSize(5)
	, mi_DnaBufferLength(0)
	, mi_IndexBufferMaxLength(8 * 1024 * 1024)
	, mk_TempFileA_(NULL)
	, mk_TempFileB_(NULL)
{
}


k_GpfIndexer::~k_GpfIndexer()
{
}


void k_GpfIndexer::compileIndex()
{
	k_StopWatch lk_StopWatch("Compiling the DNA index took %1.\n");
	
	QFile lk_OutFile(ms_DnaIndexPath);
	lk_OutFile.open(QIODevice::WriteOnly);
	
	// parse FASTA DNA file
	parseDna(ms_DnaPath);
	
	mi_OffsetBits = 1;
	while ((((qint64)1) << mi_OffsetBits) < mi_TotalNucleotideCount)
		++mi_OffsetBits;
	
	// add one extra bit for recording the reading direction
	++mi_OffsetBits;
	mui_GnoBackwardsBit = (quint64)1 << (mi_OffsetBits - 1);
	
	// determine tag bits
	mi_TagCount =  1;
	for (int i = 0; i < mi_TagSize; ++i)
		mi_TagCount *= 20;
	mi_TagBits = 1;
	while ((((qint32)1) << mi_TagBits) < mi_TagCount)
		++mi_TagBits;
	
	mui_pTagDirectionCount = RefPtr<quint32>(new quint32[mi_TagCount * 2]);
	memset(mui_pTagDirectionCount.get_Pointer(), 0, mi_TagCount * 2 * 4);
	
	muc_pIndexBuffer = RefPtr<quint8>(new quint8[mi_IndexBufferMaxLength]);
	mi_IndexBufferOffset = 0;
	mi_IndexBufferBitOffset = 0;
	
	mi_PreIndexEntryBits = mi_TagBits + 1 + mi_MassBits + mi_OffsetBits;
	mi_PreIndexEntryBytes = mi_PreIndexEntryBits / 8 + 2;
	
	muc_pOnePreIndexEntry = RefPtr<quint8>(new quint8[mi_PreIndexEntryBytes]);
	
	writeIdentifierChunk(&lk_OutFile);
	writeInfoChunk(&lk_OutFile);
	writeDnaChunk(&lk_OutFile);
	
	mk_pTempFileA = RefPtr<QFile>(new QTemporaryFile(QFileInfo(ms_DnaIndexPath).absolutePath() + "/gpfindex.temp"));
	mk_pTempFileA->open(QIODevice::ReadWrite | QIODevice::Truncate);
	mk_pTempFileB = RefPtr<QFile>(new QTemporaryFile(QFileInfo(ms_DnaIndexPath).absolutePath() + "/gpfindex.temp"));
	mk_pTempFileB->open(QIODevice::ReadWrite | QIODevice::Truncate | QIODevice::Unbuffered);

	mk_TempFileA_ = mk_pTempFileA.get_Pointer();
	mk_TempFileB_ = mk_pTempFileB.get_Pointer();
	
	writeIndexChunk(&lk_OutFile);
	
	lk_OutFile.close();
}


void k_GpfIndexer::parseDna(QString as_DnaPath)
{
	printf("Parsing input file...");
	mk_ScaffoldLabels.clear();
	mk_ScaffoldLength.clear();
	
	QFile lk_DnaFile(as_DnaPath);
	lk_DnaFile.open(QIODevice::ReadOnly);
	QTextStream lk_Stream(&lk_DnaFile);
	
	muc_pDnaBuffer = RefPtr<quint8>(new quint8[lk_DnaFile.size() * 3 / 8 + 1]);
	qint64 li_DnaBufferOffset = 0;
	mi_DnaBufferLength = 0;
	
	quint16 lui_DnaBuffer = 0;
	// :CAREFUL: li_DnaBufferLength might be confused with mi_DnaBufferLength
	// li_DnaBufferLength is about our 16 bit mini buffer, mi_DnaBufferLength
	// is the whole DNA!
	int li_DnaBufferLength = 0;
	QString ls_CurrentLabel = "";
	mi_TotalNucleotideCount = 0;

	while (!lk_Stream.atEnd())
	{
		QString ls_Line = lk_Stream.readLine();
		ls_Line = ls_Line.trimmed();
		if (ls_Line[0] == QChar('>'))
		{
			ls_Line.remove(0, 1);
			QString ls_Label = ls_Line.trimmed();
			ls_CurrentLabel = ls_Label;
			mk_ScaffoldLabels.push_back(ls_Label);
			mk_ScaffoldLength[ls_Label] = 0;
			mk_ScaffoldStart[ls_Label] = mi_TotalNucleotideCount;
		}
		else
		{
			mk_ScaffoldLength[ls_CurrentLabel] += ls_Line.length();
			mi_TotalNucleotideCount += ls_Line.length();
			for (int i = 0; i < ls_Line.length(); ++i)
			{
				lui_DnaBuffer |= (gk_GpfBase.mk_DnaCharToNumber_[(unsigned char)(ls_Line.at(i).toAscii())] << li_DnaBufferLength);
				li_DnaBufferLength += 3;
				if (li_DnaBufferLength >= 8)
				{
					quint8 lui_DnaBufferBit = (quint8)(lui_DnaBuffer & 255);
					muc_pDnaBuffer.get_Pointer()[li_DnaBufferOffset] = lui_DnaBufferBit;
					++li_DnaBufferOffset;
					++mi_DnaBufferLength;
					li_DnaBufferLength -= 8;
					lui_DnaBuffer >>= 8;
				}
			}
		}
	}
	// write remaining nucleotides
	if (li_DnaBufferLength > 0)
	{
		quint8 lui_DnaBufferBit = (quint8)(lui_DnaBuffer & 255);
		muc_pDnaBuffer.get_Pointer()[li_DnaBufferOffset] = lui_DnaBufferBit;
		++li_DnaBufferOffset;
		++mi_DnaBufferLength;
	}
	printf(" done.\n");
}


qint64 k_GpfIndexer::writeChunkHeader(QFile* ak_OutFile_, r_DnaIndexChunkType::Enumeration ae_Type)
{
	quint32 lui_Type = (quint32)ae_Type;
	ak_OutFile_->write((char*)&lui_Type, 4);
	qint64 li_Pos = ak_OutFile_->pos();
	qint64 li_Size = 0;
	ak_OutFile_->write((char*)&li_Size, 8);
	return li_Pos;
}


void k_GpfIndexer::writeChunkSize(QFile* ak_OutFile_, qint64 ai_Position, qint64 ai_Size)
{
	qint64 li_OldPos = ak_OutFile_->pos();
	ak_OutFile_->seek(ai_Position);
	ak_OutFile_->write((char*)&ai_Size, 8);
	ak_OutFile_->seek(li_OldPos);
}


void k_GpfIndexer::writeIdentifierChunk(QFile* ak_OutFile_)
{
	ak_OutFile_->write(QByteArray("gpfindex"));
	qint16 li_VersionMajor = 3;
	qint16 li_VersionMinor = 0;
	ak_OutFile_->write((char*)&li_VersionMajor, 2);
	ak_OutFile_->write((char*)&li_VersionMinor, 2);
}


void k_GpfIndexer::writeInfoChunk(QFile* ak_OutFile_)
{
	qint64 li_SizeLocation = writeChunkHeader(ak_OutFile_, r_DnaIndexChunkType::Info);
	qint64 li_ChunkSize = ak_OutFile_->pos();
	
	// write genome title
	qint32 li_TitleLength = ms_Title.length();
	ak_OutFile_->write((char*)&li_TitleLength, 4);
	ak_OutFile_->write((char*)ms_Title.toStdString().c_str(), li_TitleLength);
	
	// write constants
	ak_OutFile_->write((char*)&mi_OffsetBits, 4);
	ak_OutFile_->write((char*)&mi_MassBits, 4);
	ak_OutFile_->write((char*)&mi_MassPrecision, 4);
	ak_OutFile_->write((char*)&mi_TagSize, 4);
	
	// write scaffold count, lengths and labels
	qint32 li_ScaffoldCount = mk_ScaffoldLabels.size();
	ak_OutFile_->write((char*)&li_ScaffoldCount, 4);
	foreach (QString ls_Label, mk_ScaffoldLabels)
	{
		qint64 li_ScaffoldSize = mk_ScaffoldLength[ls_Label];
		ak_OutFile_->write((char*)&li_ScaffoldSize, 8);
		qint32 li_ScaffoldLabelLength = ls_Label.length();
		ak_OutFile_->write((char*)&li_ScaffoldLabelLength, 4);
		ak_OutFile_->write((char*)ls_Label.toStdString().c_str(), li_ScaffoldLabelLength);
	}
	
	// finish info chunk
	li_ChunkSize = ak_OutFile_->pos() - li_ChunkSize;
	writeChunkSize(ak_OutFile_, li_SizeLocation, li_ChunkSize);
}


void k_GpfIndexer::writeDnaChunk(QFile* ak_OutFile_)
{
	qint64 li_SizeLocation = writeChunkHeader(ak_OutFile_, r_DnaIndexChunkType::Dna);
	ak_OutFile_->write((char*)muc_pDnaBuffer.get_Pointer(), mi_DnaBufferLength);
	writeChunkSize(ak_OutFile_, li_SizeLocation, mi_DnaBufferLength);
}


void k_GpfIndexer::writeIndexChunk(QFile* ak_OutFile_)
{
	qint64 li_SizeLocation = writeChunkHeader(ak_OutFile_, r_DnaIndexChunkType::Index);
	qint64 li_ChunkSize = ak_OutFile_->pos();
	
// calculate mass according to mass accuracy
	for (int i = 0; i < 256; ++i)
	{
		mi_AminoAcidMasses_[i] = 0;
		if (gk_GpfBase.mb_IsAminoAcid_[i])
			mi_AminoAcidMasses_[i] = (qint64)(gk_GpfBase.md_AminoAcidMasses_[i] * mi_MassPrecision);
	}
	qint64 li_MaxMass = ((qint64)1 << mi_MassBits) - 1;
	
	printf("Translating DNA, marking cleavage sites...");
	mi_PreIndexEntryCount = 0;
	// for each scaffold: translate all six reading frames
	foreach (QString ls_Label, mk_ScaffoldLabels)
	{
		RefPtr<char> lc_pOpenReadingFrame(new char[mk_ScaffoldLength[ls_Label] / 3 + 1]);
		qint64 li_ScaffoldStart = mk_ScaffoldStart[ls_Label];
		qint64 li_ScaffoldSize = mk_ScaffoldLength[ls_Label];
		for (int li_Direction = 0; li_Direction < 2; ++li_Direction)
		{
			for (int li_Frame = 0; li_Frame < 3; ++li_Frame)
			{
				// translate DNA into open reading frame, collect cleavage sites
				QSet<qint64> lk_CleavageSitesSet;
				lk_CleavageSitesSet << 0;
				qint64 li_OrfLength = 0;
				qint64 li_NucleotideStart = li_ScaffoldStart + li_Frame;
				qint64 li_NucleotideEnd = li_ScaffoldStart + li_ScaffoldSize - 2;
				if (li_Direction == 1)
				{
					li_NucleotideStart = li_ScaffoldStart + li_ScaffoldSize - 1 - li_Frame - 2;
					li_NucleotideEnd = li_ScaffoldStart - 1;
				}
				qint64 li_NucleotideStep = li_Direction == 0 ? 3 : -3;
				for (qint64 i = li_NucleotideStart; 
					 li_Direction == 0 ? 
						(i < li_NucleotideEnd) :
						(i > li_NucleotideEnd);
					 i += li_NucleotideStep)
				{
					quint16 lui_Triplet = readNucleotideTriplet(i);
					char lc_AminoAcid = li_Direction == 0 ? 
						gk_GpfBase.mk_DnaTripletToAminoAcid_[lui_Triplet] :
						gk_GpfBase.mk_DnaTripletToAminoAcidReverse_[lui_Triplet];
						
					lc_pOpenReadingFrame.get_Pointer()[li_OrfLength] = lc_AminoAcid;
					++li_OrfLength;
					if (lc_AminoAcid == '$' || lc_AminoAcid == 'R' || lc_AminoAcid == 'K')
						lk_CleavageSitesSet << li_OrfLength;
				}
				lk_CleavageSitesSet << li_OrfLength;
				QList<qint64> lk_CleavageSites = lk_CleavageSitesSet.toList();
				qSort(lk_CleavageSites);
				
				// now we have the ORF and all cleavage sites!
				for (int i = 0; i < lk_CleavageSites.size() - 1; ++i)
				{
					// here is an amino acid span (a tryptic peptide)
					qint64 li_SpanStart = lk_CleavageSites[i];
					qint64 li_SpanEnd = lk_CleavageSites[i + 1] - 1;
					
					// cut off trailing STOP if any
					if (lc_pOpenReadingFrame.get_Pointer()[li_SpanEnd] == '$')
						--li_SpanEnd;
					
					// check whether span is long enough
					if (li_SpanEnd - li_SpanStart + 1 < mi_TagSize)
						continue;
					
					// check whether unknown amino acids are in there
					bool lb_UnknownAminoAcids = false;
					for (qint64 k = li_SpanStart; !lb_UnknownAminoAcids && (k <= li_SpanEnd); ++k)
						if (!gk_GpfBase.mb_IsAminoAcid_[(int)lc_pOpenReadingFrame.get_Pointer()[k]])
							lb_UnknownAminoAcids = true;
					if (lb_UnknownAminoAcids)
						continue;
						
/*					for (qint64 k = li_SpanStart; (k <= li_SpanEnd); ++k)
						printf("%c", lc_pOpenReadingFrame.get_Pointer()[k]);
					printf("\n");*/
					
					qint64 li_LeftHalfMass = 0;
					for (qint64 k = li_SpanStart; ((k + mi_TagSize - 1) <= li_SpanEnd) && (li_LeftHalfMass <= li_MaxMass); ++k)
					{
						int li_Tag = gk_GpfBase.aminoAcidPolymerCode(lc_pOpenReadingFrame.get_Pointer() + k, mi_TagSize);
						quint64 lui_Gno = li_Direction == 0 ? 
							li_SpanStart * 3 + li_ScaffoldStart + li_Frame :
							(li_ScaffoldStart + li_ScaffoldSize - 1 - li_Frame) - li_SpanStart * 3;
						if (li_Direction == 1)
							lui_Gno |= mui_GnoBackwardsBit;
						recordHalfMassSequenceTag(li_Tag, false, li_LeftHalfMass, lui_Gno);
						li_LeftHalfMass += mi_AminoAcidMasses_[(int)lc_pOpenReadingFrame.get_Pointer()[k]];
					}

					qint64 li_RightHalfMass = 0;
					for (qint64 k = li_SpanEnd - mi_TagSize + 1; (k >= li_SpanStart) && (li_RightHalfMass <= li_MaxMass); --k)
					{
						int li_Tag = gk_GpfBase.aminoAcidPolymerCode(lc_pOpenReadingFrame.get_Pointer() + k, mi_TagSize);
						quint64 lui_Gno = li_Direction == 0 ?
							li_SpanEnd * 3 + 2 + li_ScaffoldStart + li_Frame :
							(li_ScaffoldStart + li_ScaffoldSize - 1 - li_Frame) - li_SpanEnd * 3 - 2;
						if (li_Direction == 1)
							lui_Gno |= mui_GnoBackwardsBit;
						recordHalfMassSequenceTag(li_Tag, true, li_RightHalfMass, lui_Gno);
						li_RightHalfMass += mi_AminoAcidMasses_[(int)lc_pOpenReadingFrame.get_Pointer()[k + mi_TagSize - 1]];
					}
				}
			}
		}
	}
	
	flushIndexBuffer();
	printf(" done.\n");
	
	// we now have the raw unsorted HMST in temp file A
	// now we already know how many entries we have for each tag/direction pair
	// rearrange the entries, so that all entries are already sorted by tag/direction
	mk_TempFileB_->resize(mi_PreIndexEntryBits * mi_PreIndexEntryCount / 8 + 1);
	RefPtr<quint64> lk_pTagDirectionOffset(new quint64[mi_TagCount * 2]);
	RefPtr<quint32> lk_pTagDirectionWritten(new quint32[mi_TagCount * 2]);
	memset(lk_pTagDirectionWritten.get_Pointer(), 0, mi_TagCount * 2 * 4);
	quint64 li_Offset = 0;
	for (int i = 0; i < mi_TagCount * 2; ++i)
	{
		lk_pTagDirectionOffset.get_Pointer()[i] = li_Offset;
		li_Offset += mui_pTagDirectionCount.get_Pointer()[i];
	}
	mk_TempFileA_->reset();
	mui_CurrentPreIndexByte = 0;
	mi_CurrentPreIndexByteBitsLeft = 0;
	for (qint64 i = 0; i < mi_PreIndexEntryCount; ++i)
	{
		if (i % 1000 == 0)
			printf("\rRearranging pre-index... %d%%", (int)(i * 100 / mi_PreIndexEntryCount));
		// read one pre index entry and write it to the correct position
		int li_Tag = readBitsFromTempFileA(mi_TagBits);
		int li_Direction = readBitsFromTempFileA(1);
		qint64 li_HalfMass = readBitsFromTempFileA(mi_MassBits);
		quint64 lui_Gno = readBitsFromTempFileA(mi_OffsetBits);
		qint64 li_OutPosition = 
			lk_pTagDirectionOffset.get_Pointer()[li_Tag * 2 + li_Direction] +
			lk_pTagDirectionWritten.get_Pointer()[li_Tag * 2 + li_Direction];
		writePreIndexEntryToFileB(li_OutPosition * mi_PreIndexEntryBits, li_Tag, li_Direction, li_HalfMass, lui_Gno);
		++lk_pTagDirectionWritten.get_Pointer()[li_Tag * 2 + li_Direction];
// 		printf("W %d %d %d %d\n", (int)li_Tag, (int)li_Direction, (int)li_HalfMass, (int)lui_Gno);
	}
	printf("\rRearranging pre-index... done.\n");
	
	(dynamic_cast<QTemporaryFile*>(mk_TempFileB_))->setAutoRemove(false);
	QString ls_TempFileBPath = mk_TempFileB_->fileName();
	mk_TempFileB_->close();
	mk_pTempFileB = RefPtr<QFile>(NULL);
	mk_pTempFileA = RefPtr<QFile>(new QFile(ls_TempFileBPath));
	mk_pTempFileA->open(QIODevice::ReadOnly);
	mk_TempFileB_ = NULL;
	mk_TempFileA_ = mk_pTempFileA.get_Pointer();

	int li_MaxEntryCount = 0;
	for (qint64 li_TagAndDir = 0; li_TagAndDir < mi_TagCount * 2; ++li_TagAndDir)
	{
		if ((int)mui_pTagDirectionCount.get_Pointer()[li_TagAndDir] > li_MaxEntryCount)
			li_MaxEntryCount = mui_pTagDirectionCount.get_Pointer()[li_TagAndDir];
	}
	
	RefPtr<quint8> luc_pEntryList(new quint8[li_MaxEntryCount * mi_PreIndexEntryBits / 8 + 1]);
	RefPtr<quint32> lui_pEntryListIndex(new quint32[li_MaxEntryCount]);
	if (!luc_pEntryList || !lui_pEntryListIndex)
	{
		printf("Error: Unable to allocate %d bytes for sorting.\n", li_MaxEntryCount * mi_PreIndexEntryBits / 8 + 1 + li_MaxEntryCount * 4);
		exit(1);
	}
	
	mk_TempFileA_->reset();
	mui_CurrentPreIndexByte = 0;
	mi_CurrentPreIndexByteBitsLeft = 0;
	
	for (qint64 li_TagAndDir = 0; li_TagAndDir < mi_TagCount * 2; ++li_TagAndDir)
	{
		if (li_TagAndDir % 1000 == 0)
			printf("\rCreating index... %d%%", (int)(li_TagAndDir * 100 / (mi_TagCount * 2)));
		int li_EntryListOffset = 0;
		int li_EntryListIndexOffset = 0;
		// read all entries from this tag/dir group into buffer
		for (quint32 li_Entry = 0; li_Entry < mui_pTagDirectionCount.get_Pointer()[li_TagAndDir]; ++li_Entry)
		{
			int li_Tag = readBitsFromTempFileA(mi_TagBits);
			int li_Direction = readBitsFromTempFileA(1);
			qint64 li_HalfMass = readBitsFromTempFileA(mi_MassBits);
			quint64 lui_Gno = readBitsFromTempFileA(mi_OffsetBits);
			if (li_TagAndDir != li_Tag * 2 + li_Direction)
			{
				printf("Error: expecting %d %d, got %d %d.\n", 
						(int)(li_TagAndDir / 2), (int)(li_TagAndDir % 2),
						(int)(li_Tag), (int)li_Direction
						);
				exit(1);
			}
			overwriteBitsInBuffer(luc_pEntryList.get_Pointer(), li_EntryListOffset, li_Tag, mi_TagBits);
			li_EntryListOffset += mi_TagBits;
			overwriteBitsInBuffer(luc_pEntryList.get_Pointer(), li_EntryListOffset, li_Direction, 1);
			li_EntryListOffset += 1;
			overwriteBitsInBuffer(luc_pEntryList.get_Pointer(), li_EntryListOffset, li_HalfMass, mi_MassBits);
			li_EntryListOffset += mi_MassBits;
			overwriteBitsInBuffer(luc_pEntryList.get_Pointer(), li_EntryListOffset, lui_Gno, mi_OffsetBits);
			li_EntryListOffset += mi_OffsetBits;
			lui_pEntryListIndex.get_Pointer()[li_EntryListIndexOffset] = li_EntryListIndexOffset;
			++li_EntryListIndexOffset;
		}
		// now sort the entries in RAM
	}
	printf("\rCreating index... done.\n");
	mk_pTempFileA->remove();
	
	// finish index chunk
	li_ChunkSize = ak_OutFile_->pos() - li_ChunkSize;
	writeChunkSize(ak_OutFile_, li_SizeLocation, li_ChunkSize);
}


inline quint16 k_GpfIndexer::readNucleotideTriplet(quint64 aui_Gno)
{
	// 9 bits are always contained in at most 2 bytes!
	quint64 lui_Offset = (aui_Gno * 3) / 8;
	quint16 lui_Bit;
	memcpy((char*)(&lui_Bit), muc_pDnaBuffer.get_Pointer() + lui_Offset, 1);
	memcpy((char*)(&lui_Bit) + 1, muc_pDnaBuffer.get_Pointer() + lui_Offset + 1, 1);
	lui_Bit >>= ((aui_Gno * 3) & 7);
	lui_Bit &= 511;
	return lui_Bit;
}


void k_GpfIndexer::recordHalfMassSequenceTag(int ai_Tag, bool ab_RightTag, qint64 ai_HalfMass, qint64 ai_Gno)
{
	// flush the HMST to disk, add it to the appropriate tag/direction hash
	int li_Index = ai_Tag * 2 + (ab_RightTag ? 1 : 0);
	++mui_pTagDirectionCount.get_Pointer()[li_Index];
	++mi_PreIndexEntryCount;
	writeBitsToIndexBuffer(ai_Tag, mi_TagBits);
	writeBitsToIndexBuffer(ab_RightTag, 1);
	writeBitsToIndexBuffer(ai_HalfMass, mi_MassBits);
	writeBitsToIndexBuffer(ai_Gno, mi_OffsetBits);
	//printf("W %d %d %d %d\n", (int)ai_Tag, (int)ab_RightTag, (int)ai_HalfMass, (int)ai_Gno);
}


void k_GpfIndexer::writeBitsToIndexBuffer(quint64 aui_Value, int ai_Size)
{
	while (ai_Size > 0)
	{
		int li_CopyBits = (8 - mi_IndexBufferBitOffset);
		if (li_CopyBits > ai_Size)
			li_CopyBits = ai_Size;
		quint8 lui_Byte = aui_Value & ((1 << li_CopyBits) - 1);
		quint8 lui_NullMask = (((quint32)1) << li_CopyBits) - 1;
		lui_NullMask <<= mi_IndexBufferBitOffset;
		lui_NullMask ^= 0xff;
		lui_Byte <<= mi_IndexBufferBitOffset;
		aui_Value >>= li_CopyBits;
		ai_Size -= li_CopyBits;
		muc_pIndexBuffer.get_Pointer()[mi_IndexBufferOffset] &= lui_NullMask;
		muc_pIndexBuffer.get_Pointer()[mi_IndexBufferOffset] |= lui_Byte;
		mi_IndexBufferBitOffset += li_CopyBits;
		if (mi_IndexBufferBitOffset > 7)
		{
			mi_IndexBufferBitOffset -= 8;
			++mi_IndexBufferOffset;
			if (mi_IndexBufferOffset >= mi_IndexBufferMaxLength)
				flushIndexBuffer();
		}
	}
}


void k_GpfIndexer::flushIndexBuffer()
{
	int li_Size = mi_IndexBufferOffset;
	if (mi_IndexBufferBitOffset > 0)
		++li_Size;
	mk_TempFileA_->write((char*)muc_pIndexBuffer.get_Pointer(), li_Size);
	mi_IndexBufferOffset = 0;
	mi_IndexBufferBitOffset = 0;
}


quint64 k_GpfIndexer::readBitsFromTempFileA(int ai_Size)
{
	quint64 lui_Result = 0;
	int li_BitsCopied = 0;
	while (ai_Size > 0)
	{
		if (mi_CurrentPreIndexByteBitsLeft == 0)
		{
			// read next byte from file
			mk_TempFileA_->read((char*)&mui_CurrentPreIndexByte, 1);
			mi_CurrentPreIndexByteBitsLeft = 8;
		}
		int li_CopyBits = mi_CurrentPreIndexByteBitsLeft;
		if (li_CopyBits > ai_Size)
			li_CopyBits = ai_Size;
		// extract bits
		quint64 lui_Byte = mui_CurrentPreIndexByte >> (8 - mi_CurrentPreIndexByteBitsLeft);
		lui_Byte &= ((quint64)1 << li_CopyBits) - 1;
		lui_Result |= lui_Byte << li_BitsCopied;
		li_BitsCopied += li_CopyBits;
		mi_CurrentPreIndexByteBitsLeft -= li_CopyBits;
		ai_Size -= li_CopyBits;
	}
	return lui_Result;
}


void k_GpfIndexer::writePreIndexEntryToFileB(qint64 ai_OutPosition, int ai_Tag, int ai_Direction, qint64 ai_HalfMass, quint64 aui_Gno)
{
	// write the pre index entry to the specified position (in bits)
	// read a small chunk to muc_pOnePreIndexEntry
	int li_FilePosition = ai_OutPosition / 8;
	int li_BufferPosition = ai_OutPosition - li_FilePosition * 8;
	mk_TempFileB_->seek(li_FilePosition);
	mk_TempFileB_->read((char*)muc_pOnePreIndexEntry.get_Pointer(), mi_PreIndexEntryBytes);
	// inject entry into these few bytes
	overwriteBitsInBuffer(muc_pOnePreIndexEntry.get_Pointer(), li_BufferPosition, ai_Tag, mi_TagBits);
	li_BufferPosition += mi_TagBits;
	overwriteBitsInBuffer(muc_pOnePreIndexEntry.get_Pointer(), li_BufferPosition, ai_Direction, 1);
	li_BufferPosition += 1;
	overwriteBitsInBuffer(muc_pOnePreIndexEntry.get_Pointer(), li_BufferPosition, ai_HalfMass, mi_MassBits);
	li_BufferPosition += mi_MassBits;
	overwriteBitsInBuffer(muc_pOnePreIndexEntry.get_Pointer(), li_BufferPosition, aui_Gno, mi_OffsetBits);
	li_BufferPosition += mi_OffsetBits;
	// write back the modified bytes
	mk_TempFileB_->seek(li_FilePosition);
	mk_TempFileB_->write((char*)muc_pOnePreIndexEntry.get_Pointer(), mi_PreIndexEntryBytes);
}


void k_GpfIndexer::overwriteBitsInBuffer(quint8* auc_Buffer_, int ai_Offset, quint64 aui_Value, int ai_Size)
{
	while (ai_Size > 0)
	{
		int li_ByteOffset = ai_Offset / 8;
		int li_BitOffset = ai_Offset % 8;
		int li_CopyBits = (8 - li_BitOffset);
		if (li_CopyBits > ai_Size)
			li_CopyBits = ai_Size;
		quint8 lui_CopyMask = ((((quint32)1) << li_CopyBits) - 1);
		quint8 lui_NullMask = ~(lui_CopyMask << li_BitOffset);
		auc_Buffer_[li_ByteOffset] &= lui_NullMask;
		quint8 lui_CopyByte = (aui_Value & lui_CopyMask) << li_BitOffset;
		auc_Buffer_[li_ByteOffset] |= lui_CopyByte;
		ai_Offset += li_CopyBits;
		aui_Value >>= li_CopyBits;
		ai_Size -= li_CopyBits;
	}
}
