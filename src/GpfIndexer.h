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
	
	inline quint16 readNucleotideTriplet(quint64 aui_Gno);
	void recordHalfMassSequenceTag(int ai_Tag, bool ab_RightTag, qint64 ai_HalfMass, qint64 ai_Gno);
	void writeBitsToIndexBuffer(quint64 aui_Value, int ai_Size);
	void flushIndexBuffer();
	quint64 readBitsFromTempFileA(int ai_Size);
	void writePreIndexEntryToFileB(qint64 ai_OutPosition, int ai_Tag, int ai_Direction, qint64 ai_HalfMass, quint64 aui_Gno);
	void overwriteBitsInBuffer(quint8* auc_Buffer_, int ai_Offset, quint64 aui_Value, int ai_Size);
	
	QString ms_DnaPath;
	QString ms_DnaIndexPath;
	QString ms_Title;
	
	QStringList mk_ScaffoldLabels;
	QHash<QString, qint64> mk_ScaffoldLength;
	QHash<QString, qint64> mk_ScaffoldStart;
	qint64 mi_TotalNucleotideCount;
	
	qint64 mi_OffsetBits;
	quint64 mui_GnoBackwardsBit;
	qint32 mi_MassBits;
	qint32 mi_MassPrecision;
	qint32 mi_TagSize;
	qint32 mi_TagBits;
	qint32 mi_PreIndexEntryBits;
	qint32 mi_PreIndexEntryBytes;
	qint32 mi_TagCount;
	
	qint64 mi_AminoAcidMasses_[256];
	
	// this buffer contains the nucleotides as 3-bit-entities
	RefPtr<quint8> muc_pDnaBuffer;
	qint64 mi_DnaBufferLength;
	
	RefPtr<quint8> muc_pIndexBuffer;
	qint32 mi_IndexBufferMaxLength;
	qint32 mi_IndexBufferOffset;
	qint32 mi_IndexBufferBitOffset;
	quint8 mui_CurrentPreIndexByte;
	quint32 mi_CurrentPreIndexByteBitsLeft;
	RefPtr<QFile> mk_pTempFileA;
	RefPtr<QFile> mk_pTempFileB;
	QFile* mk_TempFileA_;
	QFile* mk_TempFileB_;
	
	RefPtr<quint32> mui_pTagDirectionCount;
	qint64 mi_PreIndexEntryCount;
	RefPtr<quint8> muc_pOnePreIndexEntry;
};
