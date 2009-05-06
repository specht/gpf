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
#include "BitWriter.h"
#include "GpfBase.h"
#include "HmstIterator.h"
#include "RefPtr.h"
#include "StopWatch.h"


// :UGLY: This is really ugly indeed. But it has been done in order to be able
// to use Qt's qsort() function with a simple less-than function pointer which
// can not be a class member function pointer. Therefore, we just store the current
// indexer in this global variable in order to do the comparison in a global
// less-than function. Otherwise, we would have to implement our own sorting function
// or, as an alternative, implement a new class containing the int and the pointer to
// the appropriate indexer, which means much more space. Because probably, only one
// genome file is indexed at a time, we can probably do this trick.
// Actually, it's not a pointer to the current indexer, but to the current list 
// containing all masses from the current tag/direction set while sorting out there
// final index.
quint8* guc_MassesBuffer_;
qint32 gi_MassBits;


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
	, mi_IndexBufferMaxLength(512 * 1024 * 1024)
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
	
	mi_HmstBits = mi_OffsetBits + mi_MassBits;
	
	printf("Allocating %s for tag/direction count (32 bits) list.\n", bytesToStr(mi_TagCount * 2 * 4).toStdString().c_str());
	mui_pTagDirectionCount = RefPtr<quint32>(new quint32[mi_TagCount * 2]);
	memset(mui_pTagDirectionCount.get_Pointer(), 0, mi_TagCount * 2 * 4);
	
	printf("Allocating %s for index buffer.\n", bytesToStr(mi_IndexBufferMaxLength).toStdString().c_str());
	muc_pIndexBuffer = RefPtr<quint8>(new quint8[mi_IndexBufferMaxLength]);
	memset(muc_pIndexBuffer.get_Pointer(), 0, mi_IndexBufferMaxLength);
	mi_IndexBufferOffset = 0;
	mi_IndexBufferBitOffset = 0;
	
	writeIdentifierChunk(&lk_OutFile);
	writeInfoChunk(&lk_OutFile);
	writeDnaChunk(&lk_OutFile);
	
	writeIndexChunk(&lk_OutFile);
	
	lk_OutFile.close();
}


void k_GpfIndexer::parseDna(QString as_DnaPath)
{
	//printf("Parsing input file...");
	mk_ScaffoldLabels.clear();
	mk_ScaffoldLength.clear();
	
	QFile lk_DnaFile(as_DnaPath);
	lk_DnaFile.open(QIODevice::ReadOnly);
	QTextStream lk_Stream(&lk_DnaFile);
	
	printf("Allocating %s for DNA.\n", bytesToStr(lk_DnaFile.size() * 3 / 8 + 1).toStdString().c_str());
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
			mk_ScaffoldLength.push_back(0);
			mk_ScaffoldStart.push_back(mi_TotalNucleotideCount);
		}
		else
		{
			mk_ScaffoldLength.last() += ls_Line.length();
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
	//printf(" done.\n");
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
	for (int i = 0; i < li_ScaffoldCount; ++i)
	{
		QString ls_Label = mk_ScaffoldLabels[i];
		qint64 li_ScaffoldSize = mk_ScaffoldLength[i];
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
	
//	mk_IndexBufferBitWriterFile_ = mk_TempFileA_;
	mi_IndexBufferOffset = 0;
	mi_IndexBufferBitOffset = 0;


	qint64 li_MaxLen = 0;
	foreach (qint64 li_Length, mk_ScaffoldLength)
		if (li_Length > li_MaxLen)
			li_MaxLen = li_Length;
		
	RefPtr<char> lc_pOpenReadingFrame(new char[li_MaxLen / 3 + 1]);
	//printf("Translating DNA, marking cleavage sites...");
	
	
	printf("HMST size: %d bits (%d for mass, %d for GNO)\n", (int)(mi_HmstBits), (int)mi_MassBits, (int)mi_OffsetBits);
	k_HmstIterator lk_HmstIterator(*this);
	r_Hmst lr_Hmst;
	qint64 li_TotalHmstCount = 0;
	RefPtr<qint64> li_pTagDirectionCount(new qint64[mi_TagCount * 2]);
	memset(li_pTagDirectionCount.get_Pointer(), 0, sizeof(qint64) * mi_TagCount * 2);
	qint64 li_BiggestBucketSize = 0;
	while (lk_HmstIterator.next(&lr_Hmst))
	{
		//printf("%d %d %d\n", (unsigned int)lr_Hmst.mui_TagDirectionIndex, (int)lr_Hmst.mi_HalfMass, (unsigned int)lr_Hmst.mui_Gno);
		// count
		++li_TotalHmstCount;
		if (li_TotalHmstCount % 10000 == 0)
			printf("\r%d %d", (int)(li_TotalHmstCount >> 32), (int)(li_TotalHmstCount & 0xffffffff));
		++li_pTagDirectionCount.get_Pointer()[lr_Hmst.mui_TagDirectionIndex];
		if (li_pTagDirectionCount.get_Pointer()[lr_Hmst.mui_TagDirectionIndex] > li_BiggestBucketSize)
			li_BiggestBucketSize = li_pTagDirectionCount.get_Pointer()[lr_Hmst.mui_TagDirectionIndex];
	}
	printf("\nfound %d %d HMST\n", (int)(li_TotalHmstCount >> 32), (int)(li_TotalHmstCount & 0xffffffff));
	
	qint64 li_MaxHmstPerIteration = ((qint64)mi_IndexBufferMaxLength * 8) / mi_HmstBits;
	
	if (li_BiggestBucketSize > li_MaxHmstPerIteration)
	{
		printf("Error: Unfortunately, we cannot continue here because the indexing buffer is to small.\n");
		printf("You need to specify at least %s for the index buffer.\n", bytesToStr(li_BiggestBucketSize * (mi_HmstBits) / 8 + 1).toStdString().c_str());
		exit(1);
	}
		
	// determine necessary iterations
	QList<QPair<int, int> > lk_IterationRanges;
	int li_TagDirection = 0;
	qint64 li_HmstForCurrentIteration = li_pTagDirectionCount.get_Pointer()[0];
	lk_IterationRanges << QPair<int, int>(0, 0);
	while (li_TagDirection < mi_TagCount * 2)
	{
		if (li_HmstForCurrentIteration + li_pTagDirectionCount.get_Pointer()[li_TagDirection] <= li_MaxHmstPerIteration)
		{
			lk_IterationRanges.last().second = li_TagDirection;
			li_HmstForCurrentIteration += li_pTagDirectionCount.get_Pointer()[li_TagDirection];
		}
		else
		{
			lk_IterationRanges << QPair<int, int>(li_TagDirection, li_TagDirection);
			li_HmstForCurrentIteration = li_pTagDirectionCount.get_Pointer()[li_TagDirection];
		}
		++li_TagDirection;
	}
	
	RefPtr<qint64> li_pTagDirectionOffset(new qint64[mi_TagCount * 2]);
	
	RefPtr<k_BitWriter> lk_pBitWriter(new k_BitWriter(ak_OutFile_));
	
	// perform iterations!
	printf("Required iterations: %d.\n", lk_IterationRanges.size());
	typedef QPair<int, int> tk_IntPair;
	foreach (tk_IntPair lk_Pair, lk_IterationRanges)
	{
		unsigned int lui_FirstTagDirectionIndex = lk_Pair.first;
		unsigned int lui_LastTagDirectionIndex = lk_Pair.second;
		
		// determine starting offsets for each tag/dir
		qint64 li_Offset = 0;
		for (unsigned int i = lui_FirstTagDirectionIndex; i <= lui_LastTagDirectionIndex; ++i)
		{
			li_pTagDirectionOffset.get_Pointer()[i] = li_Offset;
			li_Offset += li_pTagDirectionCount.get_Pointer()[i];
		}
		
		printf("range: %d - %d\n", lui_FirstTagDirectionIndex, lui_LastTagDirectionIndex);
		
		lk_HmstIterator.reset();
		while (lk_HmstIterator.next(&lr_Hmst))
		{
			if (lr_Hmst.mui_TagDirectionIndex >= lui_FirstTagDirectionIndex && lr_Hmst.mui_TagDirectionIndex <= lui_LastTagDirectionIndex)
			{
				//printf("HMST: %d, %08x\n", (int)lr_Hmst.mi_HalfMass, lr_Hmst.mui_Gno);
				// puts HMST into right place
				overwriteBitsInBuffer(muc_pIndexBuffer.get_Pointer(), li_pTagDirectionOffset.get_Pointer()[lr_Hmst.mui_TagDirectionIndex] * mi_HmstBits, lr_Hmst.mi_HalfMass, mi_MassBits);
				overwriteBitsInBuffer(muc_pIndexBuffer.get_Pointer(), li_pTagDirectionOffset.get_Pointer()[lr_Hmst.mui_TagDirectionIndex] * mi_HmstBits + mi_MassBits, lr_Hmst.mui_Gno, mi_OffsetBits);
/*				qint64 li_ControlMass = readBitsFromBuffer(muc_pIndexBuffer.get_Pointer(), li_pTagDirectionOffset.get_Pointer()[lr_Hmst.mui_TagDirectionIndex] * mi_HmstBits, mi_MassBits);
				qint64 li_ControlGno = readBitsFromBuffer(muc_pIndexBuffer.get_Pointer(), li_pTagDirectionOffset.get_Pointer()[lr_Hmst.mui_TagDirectionIndex] * mi_HmstBits + mi_MassBits, mi_OffsetBits);
				printf("%d/%d - %d/%d\n", (int)lr_Hmst.mi_HalfMass, (int)lr_Hmst.mui_Gno, 
						(int)li_ControlMass, (int)li_ControlGno);*/
				// advance tag/dir offset
				++li_pTagDirectionOffset.get_Pointer()[lr_Hmst.mui_TagDirectionIndex];
			}
		}
		
		// sort each tag/dir list
		li_Offset = 0;
		for (unsigned int i = lui_FirstTagDirectionIndex; i <= lui_LastTagDirectionIndex; ++i)
		{
			if (li_pTagDirectionCount.get_Pointer()[i] > 0)
			{
				// sort items from li_FirstItem to li_LastItem according to half mass
				qint64 li_FirstItem = li_Offset;
				qint64 li_LastItem = li_Offset + li_pTagDirectionCount.get_Pointer()[i] - 1;
				//printf("%d %d-%d\n", i, (int)li_FirstItem, (int)li_LastItem);
				quickSortHmst(muc_pIndexBuffer.get_Pointer(), li_FirstItem, li_LastItem, mi_MassBits, mi_OffsetBits);
				li_Offset += li_pTagDirectionCount.get_Pointer()[i];
			}
		}
		
		// flush range to index file
		li_Offset = 0;
		for (unsigned int i = lui_FirstTagDirectionIndex; i <= lui_LastTagDirectionIndex; ++i)
		{
			if (li_pTagDirectionCount.get_Pointer()[i] > 0)
			{
				// sort items from li_FirstItem to li_LastItem according to half mass
				qint64 li_FirstItem = li_Offset;
				qint64 li_LastItem = li_Offset + li_pTagDirectionCount.get_Pointer()[i] - 1;
				for (qint64 k = li_FirstItem; k <= li_LastItem; ++k)
					lk_pBitWriter->writeBits(readBitsFromBuffer(muc_pIndexBuffer.get_Pointer(), k * mi_HmstBits, mi_MassBits), mi_MassBits);
				for (qint64 k = li_FirstItem; k <= li_LastItem; ++k)
					lk_pBitWriter->writeBits(readBitsFromBuffer(muc_pIndexBuffer.get_Pointer(), k * mi_HmstBits + mi_MassBits, mi_OffsetBits), mi_OffsetBits);
				li_Offset += li_pTagDirectionCount.get_Pointer()[i];
			}
		}
	}
	
	lk_pBitWriter->flush();
	
	//printf(" done.\n");

	/*
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
	
	RefPtr<quint8> luc_pEntryMassList(new quint8[li_MaxEntryCount * mi_MassBits / 8 + 1]);
	RefPtr<quint8> luc_pEntryGnoList(new quint8[li_MaxEntryCount * mi_OffsetBits / 8 + 1]);
	QList<int> lk_EntryListIndex;
	if (!luc_pEntryMassList || !luc_pEntryGnoList)
	{
		printf("Error: Unable to allocate a few bytes for sorting.\n");
		exit(1);
	}
	guc_MassesBuffer_ = luc_pEntryMassList.get_Pointer();
	gi_MassBits = mi_MassBits;
	
	mk_TempFileA_->reset();
	mui_CurrentPreIndexByte = 0;
	mi_CurrentPreIndexByteBitsLeft = 0;
	
	// determine how many bits are required for storing the entry count
	quint32 lui_EntryCountBits = 1;
	while (((int)1 << lui_EntryCountBits) < li_MaxEntryCount + 1)
		++lui_EntryCountBits;

	// write entry count bit count
	ak_OutFile_->write((char*)&lui_EntryCountBits, 4);
	
	mk_IndexBufferBitWriterFile_ = ak_OutFile_;
	mi_IndexBufferOffset = 0;
	mi_IndexBufferBitOffset = 0;
	// write entry counts
	for (qint64 li_TagAndDir = 0; li_TagAndDir < mi_TagCount * 2; ++li_TagAndDir)
		writeBitsToIndexBuffer(mui_pTagDirectionCount.get_Pointer()[li_TagAndDir], lui_EntryCountBits);
	
	for (qint64 li_TagAndDir = 0; li_TagAndDir < mi_TagCount * 2; ++li_TagAndDir)
	{
		if (li_TagAndDir % 1000 == 0)
			printf("\rCreating index... %d%%", (int)(li_TagAndDir * 100 / (mi_TagCount * 2)));
		lk_EntryListIndex.clear();
		int li_EntryMassListOffset = 0;
		int li_EntryGnoListOffset = 0;
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
			overwriteBitsInBuffer(luc_pEntryMassList.get_Pointer(), li_EntryMassListOffset, li_HalfMass, mi_MassBits);
			li_EntryMassListOffset += mi_MassBits;
			overwriteBitsInBuffer(luc_pEntryGnoList.get_Pointer(), li_EntryGnoListOffset, lui_Gno, mi_OffsetBits);
			li_EntryGnoListOffset += mi_OffsetBits;
			lk_EntryListIndex.append(lk_EntryListIndex.size());
		}
		// now sort the entries in RAM
		qSort(lk_EntryListIndex.begin(), lk_EntryListIndex.end(), lessThanPreIndexEntry);
		for (int i = 0; i < lk_EntryListIndex.size(); ++i)
		{
			qint64 li_Mass = readBitsFromBuffer(luc_pEntryMassList.get_Pointer(), lk_EntryListIndex[i] * mi_MassBits, mi_MassBits);
			writeBitsToIndexBuffer(li_Mass, mi_MassBits);
		}
		for (int i = 0; i < lk_EntryListIndex.size(); ++i)
		{
			quint64 lui_Gno = readBitsFromBuffer(luc_pEntryGnoList.get_Pointer(), lk_EntryListIndex[i] * mi_OffsetBits, mi_OffsetBits);
			writeBitsToIndexBuffer(lui_Gno, mi_OffsetBits);
		}
	}
	flushIndexBuffer();
	printf("\rCreating index... done.\n");
	mk_pTempFileA->remove();
	*/
	
	// finish index chunk
	li_ChunkSize = ak_OutFile_->pos() - li_ChunkSize;
	writeChunkSize(ak_OutFile_, li_SizeLocation, li_ChunkSize);
}


quint16 k_GpfIndexer::readNucleotideTriplet(quint64 aui_Gno)
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


quint64 k_GpfIndexer::readBitsFromTempFileA(int ai_Size)
{
	quint64 lui_Result = 0;
	/*
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
	*/
	return lui_Result;
}


/*
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
*/


void overwriteBitsInBuffer(quint8* auc_Buffer_, qint64 ai_BitOffset, quint64 aui_Value, int ai_Size)
{
	while (ai_Size > 0)
	{
		int li_ByteOffset = ai_BitOffset / 8;
		int li_BitOffset = ai_BitOffset % 8;
		int li_CopyBits = (8 - li_BitOffset);
		if (li_CopyBits > ai_Size)
			li_CopyBits = ai_Size;
		quint8 lui_CopyMask = ((((quint32)1) << li_CopyBits) - 1);
		quint8 lui_NullMask = ~(lui_CopyMask << li_BitOffset);
		auc_Buffer_[li_ByteOffset] &= lui_NullMask;
		quint8 lui_CopyByte = (aui_Value & lui_CopyMask) << li_BitOffset;
		auc_Buffer_[li_ByteOffset] |= lui_CopyByte;
		ai_BitOffset += li_CopyBits;
		aui_Value >>= li_CopyBits;
		ai_Size -= li_CopyBits;
	}
}


quint64 readBitsFromBuffer(quint8* auc_Buffer_, qint64 ai_BitOffset, int ai_Size)
{
	quint64 lui_Result = 0;
	int li_BitsCopied = 0;
	while (ai_Size > 0)
	{
		int li_ByteOffset = ai_BitOffset / 8;
		int li_BitOffset = ai_BitOffset % 8;
		int li_CopyBits = (8 - li_BitOffset);
		if (li_CopyBits > ai_Size)
			li_CopyBits = ai_Size;
		quint8 lui_CopyMask = ((((quint32)1) << li_CopyBits) - 1);
		quint8 lui_CopyByte = (auc_Buffer_[li_ByteOffset] >> li_BitOffset) & lui_CopyMask;
		lui_Result |= (((quint64)lui_CopyByte) << li_BitsCopied);
		ai_BitOffset += li_CopyBits;
		ai_Size -= li_CopyBits;
		li_BitsCopied += li_CopyBits;
	}
	return lui_Result;
}


void quickSortHmst(quint8* auc_Buffer_, int ai_First, int ai_Last, int ai_MassBits, int ai_OffsetBits)
{
	if (ai_First < ai_Last)
	{
		int li_Pivot = quickSortHmstDivide(auc_Buffer_, ai_First, ai_Last, ai_MassBits, ai_OffsetBits);
		//printf("pivot is %d!\n", li_Pivot);
		quickSortHmst(auc_Buffer_, ai_First, li_Pivot - 1, ai_MassBits, ai_OffsetBits);
		quickSortHmst(auc_Buffer_, li_Pivot + 1, ai_Last, ai_MassBits, ai_OffsetBits);
	}
}


int quickSortHmstDivide(quint8* auc_Buffer_, int ai_First, int ai_Last, int ai_MassBits, int ai_OffsetBits)
{
	int i = ai_First;
	int j = ai_Last - 1;
	qint64 li_PivotMass = readBitsFromBuffer(auc_Buffer_, ai_Last * (ai_MassBits + ai_OffsetBits), ai_MassBits);
	//printf("pivot mass is %d!\n", (int)li_PivotMass);
	do 
	{
		while ((i < ai_Last) && (qint64)readBitsFromBuffer(auc_Buffer_, i * (ai_MassBits + ai_OffsetBits), ai_MassBits) <= li_PivotMass)
			++i;
		while ((j > ai_First) && (qint64)readBitsFromBuffer(auc_Buffer_, j * (ai_MassBits + ai_OffsetBits), ai_MassBits) >= li_PivotMass)
			--j;
		if (i < j)
		{
			qint64 li_MassA = readBitsFromBuffer(auc_Buffer_, j * (ai_MassBits + ai_OffsetBits), ai_MassBits);
			qint64 li_MassB = readBitsFromBuffer(auc_Buffer_, i * (ai_MassBits + ai_OffsetBits), ai_MassBits);
			overwriteBitsInBuffer(auc_Buffer_, i * (ai_MassBits + ai_OffsetBits), li_MassA, ai_MassBits);
			overwriteBitsInBuffer(auc_Buffer_, j * (ai_MassBits + ai_OffsetBits), li_MassB, ai_MassBits);
		}
	}
	while (i < j);
	
	if (readBitsFromBuffer(auc_Buffer_, i * (ai_MassBits + ai_OffsetBits), ai_MassBits) > li_PivotMass)
	{
		qint64 li_MassA = readBitsFromBuffer(auc_Buffer_, ai_Last * (ai_MassBits + ai_OffsetBits), ai_MassBits);
		qint64 li_MassB = readBitsFromBuffer(auc_Buffer_, i * (ai_MassBits + ai_OffsetBits), ai_MassBits);
		overwriteBitsInBuffer(auc_Buffer_, i * (ai_MassBits + ai_OffsetBits), li_MassA, ai_MassBits);
		overwriteBitsInBuffer(auc_Buffer_, ai_Last * (ai_MassBits + ai_OffsetBits), li_MassB, ai_MassBits);
	}
	return i;
}

/*
bool lessThanPreIndexEntry(const int ai_First, const int ai_Second)
{
	qint64 li_FirstMass = readBitsFromBuffer(guc_MassesBuffer_, ai_First * gi_MassBits, gi_MassBits);
	qint64 li_SecondMass = readBitsFromBuffer(guc_MassesBuffer_, ai_Second * gi_MassBits, gi_MassBits);
	return li_FirstMass < li_SecondMass;
}
*/


QString k_GpfIndexer::bytesToStr(qint64 ai_Size)
{
	if (ai_Size < (qint64)1024)
		return QString("%1 bytes").arg(ai_Size);
	else if (ai_Size < (qint64)1024 * 1024)
		return QString("%1 KB").arg(ai_Size / 1024);
	else if (ai_Size < (qint64)1024 * 1024 * 1024)
		return QString("%1 MB").arg(ai_Size / 1024 / 1024);
	else if (ai_Size < (qint64)1024 * 1024 * 1024 * 1024)
		return QString("%1 GB").arg(ai_Size / 1024 / 1024 / 1024);
	else 
		return QString("%1 TB").arg(ai_Size / 1024 / 1024 / 1024 / 1024);
		
}
