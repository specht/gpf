/*
Copyright (c) 2007-2010 Michael Specht

This file is part of GPF.

GPF is free software: you can redistribute it and/or modify it under 
the terms of the GNU Lesser General Public License as published by the 
Free Software Foundation, either version 3 of the License, or (at your 
option) any later version.

GPF is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GPF.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <QtCore>


class k_GpfIndexFile
{
public:
    k_GpfIndexFile(const QString& as_Path);
    virtual ~k_GpfIndexFile();
    
    bool isGood();
    quint64 readIndexBits(qint64 ai_Offset, qint32 ai_Size);
    qint64 peptideMass(const QString& as_Peptide);
    
    QFile mk_File;

    qint16 mi_VersionMajor;
    qint16 mi_VersionMinor;
    
    QString ms_Title;
    QString ms_ShortId;
    qint32 mi_OffsetBits;
    qint32 mi_MassBits;
    qint32 mi_MassPrecision;
    qint32 mi_TagSize;
    qint32 mi_HmstBits;
    qint64 mi_GnoBackwardsBit;
    qint32 mi_GeneticCode;

    QString ms_Enzyme;
    
    qint32 mi_TagCount;
    qint64 mi_MaxMass;

    QStringList mk_ScaffoldLabels;
    QList<qint64> mk_ScaffoldLength;
    QList<qint64> mk_ScaffoldStart;
    qint64 mi_TotalNucleotideCount;
    
    QSharedPointer<quint8> muc_pDnaBuffer;
    
    QList<qint64> mk_HmstOffset;
    QList<qint64> mk_HmstCount;
    qint64 mi_BiggestHmstCount;
    qint64 mi_TotalHmstCount;
    
    qint64 mi_IndexFilePosition;

    qint64 mi_AminoAcidMasses_[256];
    qint64 mi_MinAminoAcidMass;
    qint64 mi_WaterMass;
    qint64 mi_MassDecimalDigits;
    
protected:
    void parseGpfIndexFile(const QString& as_Path);

    bool mb_IsGood;
};
