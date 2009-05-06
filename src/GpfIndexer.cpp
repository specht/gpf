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
	, mi_IndexBufferMaxLength(768 * 1024 * 1024)
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
	// allocate 9 extra bytes so that we're always safe if we should read several bytes at once
	muc_pIndexBuffer = RefPtr<quint8>(new quint8[mi_IndexBufferMaxLength + 9]);
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
	{
		k_StopWatch lk_SubStopWatch("One iteration took %1.\n");
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
	}
	
	qint64 li_MaxHmstPerIteration = ((qint64)mi_IndexBufferMaxLength * 8) / mi_HmstBits;
	
	RefPtr<quint32> lui_pIndicesToSort(new quint32[li_BiggestBucketSize]);
	
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
		mui_IndexBufferBitReader_ = (quint32*)(muc_pIndexBuffer.get_Pointer());
		mui_IndexBufferBitReaderOffset = 32;

		for (unsigned int i = lui_FirstTagDirectionIndex; i <= lui_LastTagDirectionIndex; ++i)
		{
			if (li_pTagDirectionCount.get_Pointer()[i] > 0)
			{
				typedef QPair<quint32, quint64> tk_MassGnoPair;
				QMap<quint32, tk_MassGnoPair> lk_Map;
				
				for (int k = 0; k < li_pTagDirectionCount.get_Pointer()[i]; ++k)
				{
					//quint32 lui_Mass = readBitsFromBuffer(muc_pIndexBuffer.get_Pointer(), (li_Offset + k) * mi_HmstBits, mi_MassBits);
					//quint64 lui_Gno = readBitsFromBuffer(muc_pIndexBuffer.get_Pointer(), (li_Offset + k) * mi_HmstBits + mi_MassBits, mi_OffsetBits);
					quint32 lui_Mass = readBitsFromIndexBuffer(mi_MassBits);
					quint64 lui_Gno = readBitsFromIndexBuffer(mi_OffsetBits);
					lk_Map.insertMulti(lui_Mass, tk_MassGnoPair(lui_Mass, lui_Gno));
				}
				
				QMapIterator<quint32, tk_MassGnoPair> lk_Iterator(lk_Map);
				while (lk_Iterator.hasNext())
				{
					lk_Iterator.next();
					lk_pBitWriter->writeBits(lk_Iterator.value().first, mi_MassBits);
				}

				lk_Iterator.toFront();
				while (lk_Iterator.hasNext())
				{
					lk_Iterator.next();
					lk_pBitWriter->writeBits(lk_Iterator.value().second, mi_OffsetBits);
				}
				
/*				// init indices
				for (unsigned int k = 0; k < li_pTagDirectionCount.get_Pointer()[i]; ++k)
					lui_pIndicesToSort.get_Pointer()[k] = k;
				
				// printf("sorting %d from %d to %d.\n", i, (int)li_FirstItem, (int)li_LastItem);
				// printf("%d %d-%d\n", i, (int)li_FirstItem, (int)li_LastItem);
				quickSortHmst(lui_pIndicesToSort.get_Pointer(), 
							   0, li_pTagDirectionCount.get_Pointer()[i] - 1, 
							   muc_pIndexBuffer.get_Pointer(), li_Offset, 
							   mi_MassBits, mi_OffsetBits);
				
				// flush to output file
				for (int k = 0; k < li_pTagDirectionCount.get_Pointer()[i]; ++k)
					lk_pBitWriter->writeBits(readBitsFromBuffer(muc_pIndexBuffer.get_Pointer(), (li_Offset + lui_pIndicesToSort.get_Pointer()[k]) * mi_HmstBits, mi_MassBits), mi_MassBits);
				for (int k = 0; k < li_pTagDirectionCount.get_Pointer()[i]; ++k)
					lk_pBitWriter->writeBits(readBitsFromBuffer(muc_pIndexBuffer.get_Pointer(), (li_Offset + lui_pIndicesToSort.get_Pointer()[k]) * mi_HmstBits + mi_MassBits, mi_OffsetBits), mi_OffsetBits);*/
				
				li_Offset += li_pTagDirectionCount.get_Pointer()[i];
			}
		}
		//printf("\n");
	}
	
	lk_pBitWriter->flush();
	
	// finish index chunk
	li_ChunkSize = ak_OutFile_->pos() - li_ChunkSize;
	writeChunkSize(ak_OutFile_, li_SizeLocation, li_ChunkSize);
}


quint16 k_GpfIndexer::readNucleotideTriplet(quint64 aui_Gno)
{
	// 9 bits are always contained in at most 2 bytes!
	quint64 lui_Offset = (aui_Gno * 3) / 8;
	quint16 lui_Bit;
	memcpy(&lui_Bit, muc_pDnaBuffer.get_Pointer() + lui_Offset, 2);
	lui_Bit >>= ((aui_Gno * 3) & 7);
	lui_Bit &= 511;
	return lui_Bit;
}


#define READ_BITS 32
#define READ_TYPE quint32


void overwriteBitsInBuffer(quint8* auc_Buffer_, qint64 ai_BitOffset, quint64 aui_Value, int ai_Size)
{
	while (ai_Size > 0)
	{
		int li_ByteOffset = ai_BitOffset / READ_BITS;
		int li_BitOffset = ai_BitOffset % READ_BITS;
		int li_CopyBits = (READ_BITS - li_BitOffset);
		if (li_CopyBits > ai_Size)
			li_CopyBits = ai_Size;
		READ_TYPE lui_CopyMask = ((((quint64)1) << li_CopyBits) - 1);
		READ_TYPE lui_NullMask = ~(lui_CopyMask << li_BitOffset);
		((READ_TYPE*)auc_Buffer_)[li_ByteOffset] &= lui_NullMask;
		READ_TYPE lui_CopyByte = (aui_Value & lui_CopyMask) << li_BitOffset;
		((READ_TYPE*)auc_Buffer_)[li_ByteOffset] |= lui_CopyByte;
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
		int li_ByteOffset = ai_BitOffset / READ_BITS;
		int li_BitOffset = ai_BitOffset % READ_BITS;
		int li_CopyBits = READ_BITS - li_BitOffset;
		if (li_CopyBits > ai_Size)
			li_CopyBits = ai_Size;
		READ_TYPE lui_CopyMask = ((((quint64)1) << li_CopyBits) - 1);
		//READ_TYPE lui_CopyByte = (auc_Buffer_[li_ByteOffset] >> li_BitOffset) & lui_CopyMask;
		READ_TYPE lui_CopyByte = (*((READ_TYPE*)auc_Buffer_ + li_ByteOffset) >> li_BitOffset) & lui_CopyMask;
		lui_Result |= (((quint64)lui_CopyByte) << li_BitsCopied);
		ai_BitOffset += li_CopyBits;
		ai_Size -= li_CopyBits;
		li_BitsCopied += li_CopyBits;
	}
	return lui_Result;
}


void quickSortHmst(quint32* aui_Indices_, int ai_First, int ai_Last, quint8* auc_Buffer_, qint64 ai_Offset, int ai_MassBits, int ai_OffsetBits)
{
	if (ai_First < ai_Last)
	{
		int li_Pivot = quickSortHmstDivide(aui_Indices_, ai_First, ai_Last, auc_Buffer_, ai_Offset, ai_MassBits, ai_OffsetBits);
		//printf("pivot is %d!\n", li_Pivot);
		quickSortHmst(aui_Indices_, ai_First, li_Pivot - 1, auc_Buffer_, ai_Offset, ai_MassBits, ai_OffsetBits);
		quickSortHmst(aui_Indices_, li_Pivot + 1, ai_Last, auc_Buffer_, ai_Offset, ai_MassBits, ai_OffsetBits);
	}
}


int quickSortHmstDivide(quint32* aui_Indices_, int ai_First, int ai_Last, quint8* auc_Buffer_, qint64 ai_Offset, int ai_MassBits, int ai_OffsetBits)
{
	int i = ai_First;
	int j = ai_Last - 1;
	qint64 li_PivotMass = readBitsFromBuffer(auc_Buffer_, (ai_Offset + aui_Indices_[ai_Last]) * (ai_MassBits + ai_OffsetBits), ai_MassBits);
	do 
	{
		while ((i < ai_Last) && (qint64)readBitsFromBuffer(auc_Buffer_, (ai_Offset + aui_Indices_[i]) * (ai_MassBits + ai_OffsetBits), ai_MassBits) <= li_PivotMass)
			++i;
		while ((j > ai_First) && (qint64)readBitsFromBuffer(auc_Buffer_, (ai_Offset + aui_Indices_[j]) * (ai_MassBits + ai_OffsetBits), ai_MassBits) >= li_PivotMass)
			--j;
		if (i < j)
		{
			quint32 lui_Temp = aui_Indices_[i];
			aui_Indices_[i] = aui_Indices_[j];
			aui_Indices_[j] = lui_Temp;
		}
	}
	while (i < j);
	
	if (readBitsFromBuffer(auc_Buffer_, (ai_Offset + aui_Indices_[i]) * (ai_MassBits + ai_OffsetBits), ai_MassBits) > li_PivotMass)
	{
		quint32 lui_Temp = aui_Indices_[ai_Last];
		aui_Indices_[ai_Last] = aui_Indices_[i];
		aui_Indices_[i] = lui_Temp;
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


quint64 k_GpfIndexer::readBitsFromIndexBuffer(quint8 ai_Bits)
{
	quint64 lui_Result = 0;
	while (ai_Bits > 0)
	{
		quint8 lui_CopyBits = 32 - mui_IndexBufferBitReaderOffset;
		if (lui_CopyBits > ai_Bits)
			lui_CopyBits = ai_Bits;
		quint32 lui_Mask = ((quint64)1 << lui_CopyBits) - 1;
		lui_Result <<= lui_CopyBits;
		lui_Result |= (mui_IndexBufferBitReaderCurrentByte >> mui_IndexBufferBitReaderOffset) & lui_Mask;
		ai_Bits -= lui_CopyBits;
		mui_IndexBufferBitReaderOffset += lui_CopyBits;
		if (mui_IndexBufferBitReaderOffset >= 32)
		{
			mui_IndexBufferBitReaderOffset -= 32;
			mui_IndexBufferBitReaderCurrentByte = *mui_IndexBufferBitReader_;
			++mui_IndexBufferBitReader_;
		}
	}
	return lui_Result;
}
