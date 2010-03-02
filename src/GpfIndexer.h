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
    k_GpfIndexer(QString as_DnaPath, QString as_DnaIndexPath, QString as_Title,
                 qint32 ai_TagSize, QString as_Enzyme, qint64 ai_IndexBufferAllocSize,
                 qint32 ai_MassPrecision, qint32 ai_MassBits, qint32 ai_GeneticCode);
    virtual ~k_GpfIndexer();
    
    virtual void compileIndex();
    
protected:
    virtual void parseDna(QString as_DnaPath);
    
    virtual qint64 writeChunkHeader(QFile* ak_OutFile_, r_DnaIndexChunkType::Enumeration ae_Type);
    virtual void writeChunkSize(QFile* ak_OutFile_, qint64 ai_Position, qint64 ai_Size);
    
    virtual void writeIdentifierChunk(QFile* ak_OutFile_);
    virtual void writeInfoChunk(QFile* ak_OutFile_);
    virtual void writeGeneticCodeChunk(QFile* ak_OutFile_);
    virtual void writeEnzymeChunk(QFile* ak_OutFile_);
    virtual void writeDnaChunk(QFile* ak_OutFile_);
    virtual void writeIndexChunk(QFile* ak_OutFile_);
    
    QString bytesToStr(qint64 ai_Size);
    
    QString ms_DnaPath;
    QString ms_DnaIndexPath;
    QString ms_Title;
    QString ms_Enzyme;
    
    QStringList mk_ScaffoldLabels;
    QList<qint64> mk_ScaffoldLength;
    QList<qint64> mk_ScaffoldStart;
    qint64 mi_TotalNucleotideCount;
    
    qint64 mi_OffsetBits;
    quint64 mui_GnoBackwardsBit;
    qint32 mi_MassBits;
    qint32 mi_MassPrecision;
    qint32 mi_TagSize;
    qint32 mi_GeneticCode;
    qint32 mi_TagCount;
    qint32 mi_HmstBits;
    qint64 mi_MaxMass;
    
    qint64 mi_AminoAcidMasses_[256];
    
    // this buffer contains the nucleotides as 3-bit-entities
    RefPtr<quint8> muc_pDnaBuffer;
    size_t mi_DnaBufferLength;
    
    RefPtr<quint8> muc_pIndexBuffer;
    size_t mi_IndexBufferMaxLength;
    size_t mi_IndexBufferOffset;
    size_t mi_IndexBufferBitOffset;
    
    RefPtr<quint32> mui_pTagDirectionCount;
    bool mb_CleaveBefore_[256];
    bool mb_CleaveAfter_[256];
};

