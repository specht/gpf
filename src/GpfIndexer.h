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

#pragma once

#include <QtCore>
#include "GpfBase.h"
#include "RefPtr.h"


class k_GpfIndexer
{
	friend class k_HmstIterator;
	
public:
	k_GpfIndexer(QString as_DnaPath, QString as_DnaIndexPath, QString as_Title);
	virtual ~k_GpfIndexer();
	
	virtual void compileIndex();
	
protected:
	virtual void parseDna(QString as_DnaPath);
	
	virtual qint64 writeChunkHeader(QFile* ak_OutFile_, r_DnaIndexChunkType::Enumeration ae_Type);
	virtual void writeChunkSize(QFile* ak_OutFile_, qint64 ai_Position, qint64 ai_Size);
	
	virtual void writeIdentifierChunk(QFile* ak_OutFile_);
	virtual void writeInfoChunk(QFile* ak_OutFile_);
	virtual void writeDnaChunk(QFile* ak_OutFile_);
	virtual void writeIndexChunk(QFile* ak_OutFile_);
	
	quint16 readNucleotideTriplet(quint64 aui_Gno);
	QString bytesToStr(qint64 ai_Size);
	quint64 readBitsFromIndexBuffer(quint8 ai_Bits);
	
	QString ms_DnaPath;
	QString ms_DnaIndexPath;
	QString ms_Title;
	
	QStringList mk_ScaffoldLabels;
	QList<qint64> mk_ScaffoldLength;
	QList<qint64> mk_ScaffoldStart;
	qint64 mi_TotalNucleotideCount;
	
	qint64 mi_OffsetBits;
	quint64 mui_GnoBackwardsBit;
	qint32 mi_MassBits;
	qint32 mi_MassPrecision;
	qint32 mi_TagSize;
	qint32 mi_TagCount;
	qint32 mi_HmstBits;
	
	qint64 mi_AminoAcidMasses_[256];
	
	// this buffer contains the nucleotides as 3-bit-entities
	RefPtr<quint8> muc_pDnaBuffer;
	qint64 mi_DnaBufferLength;
	
	RefPtr<quint8> muc_pIndexBuffer;
	qint32 mi_IndexBufferMaxLength;
	qint32 mi_IndexBufferOffset;
	qint32 mi_IndexBufferBitOffset;
	
	RefPtr<quint32> mui_pTagDirectionCount;
	
	quint32* mui_IndexBufferBitReader_;
	quint32 mui_IndexBufferBitReaderCurrentByte;
	quint8 mui_IndexBufferBitReaderOffset;
};

void overwriteBitsInBuffer(quint8* auc_Buffer_, qint64 ai_BitOffset, quint64 aui_Value, int ai_Size);
quint64 readBitsFromBuffer(quint8* auc_Buffer_, qint64 ai_BitOffset, int ai_Size);
void quickSortHmst(quint32* aui_Indices_, int ai_First, int ai_Last, quint8* auc_Buffer_, qint64 ai_Offset, int ai_MassBits, int ai_OffsetBits);
int quickSortHmstDivide(quint32* aui_Indices_, int ai_First, int ai_Last, quint8* auc_Buffer_, qint64 ai_Offset, int ai_MassBits, int ai_OffsetBits);
//bool lessThanPreIndexEntry(const int ai_First, const int ai_Second);
