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

#define WATER_MASS 18.01057


struct r_DnaIndexChunkType
{
    enum Enumeration
    {
        Info = 100,
        Dna,
        Index,
        GeneticCode,
        Enzyme
    };
};


class k_GpfBase
{
public:
    k_GpfBase();
    virtual ~k_GpfBase();
    
    quint16 readNucleotideTriplet(quint8* auc_Buffer_, quint64 aui_Gno);
    int aminoAcidPolymerCode(const char* ac_Buffer_, int ai_Length);
    QString aminoAcidSequenceForCode(int ai_Code, int ai_Length);
    QString nucleotideSequenceForCode(int ai_Code, int ai_Length);
    QString bytesToStr(qint64 ai_Size);
    
    quint16 mk_DnaCharToNumber_[256];
    
    QHash<int, QString> mk_TranslationTableTitle;
    QHash<int, QSharedPointer<char> > mk_TranslationTables;
    QHash<int, QSharedPointer<char> > mk_TranslationTablesReverse;
    
/*  char mk_DnaTripletToAminoAcid_[512];
    // takes the same triplet as the array above, but reverses and transposes
    char mk_DnaTripletToAminoAcidReverse_[512];*/
    
    double md_AminoAcidMasses_[256];
    bool mb_IsAminoAcid_[256];
    int mi_AminoAcidToNumber_[256];
    char mc_NumberToAminoAcid_[19];
};


extern k_GpfBase gk_GpfBase;


#define READ_BITS 32
#define READ_TYPE quint32


inline void overwriteBitsInBuffer(quint8* auc_Buffer_, qint64 ai_BitOffset, 
                                  quint64 aui_Value, int ai_Size)
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


inline quint64 readBitsFromBuffer(quint8* auc_Buffer_, qint64 ai_BitOffset, int ai_Size)
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


inline qint32 maskNucleotides(qint32 ai_Nucleotides, int ai_Length)
{
    return ai_Nucleotides & ((1 << (ai_Length * 3)) - 1);
}


inline qint32 reverseNucleotides(qint32 ai_Nucleotides, int ai_Length)
{
    // take all triplets and reverse their order
    qint32 li_Result = 0;
    for (int i = 0; i < ai_Length; ++i)
    {
        li_Result <<= 3;
        li_Result |= ai_Nucleotides & 7;
        ai_Nucleotides >>= 3;
    }
    return li_Result;
}


inline qint32 invertNucleotides(qint32 ai_Nucleotides, int ai_Length)
{
    qint32 li_Mask = 0;
    for (int i = 0; i < ai_Length; ++i)
    {
        li_Mask <<= 3;
        li_Mask |= 3;
    }
    return ai_Nucleotides ^ li_Mask;
}


inline qint32 transposeNucleotides(qint32 ai_Nucleotides, int ai_Length)
{
    return invertNucleotides(reverseNucleotides(ai_Nucleotides, ai_Length), ai_Length);
}


inline qint32 concatNucleotides(qint32 ai_FirstNucleotides, int ai_FirstLength,
                                qint32 ai_AppendNucleotides, int ai_AppendLength)
{
    qint32 li_Result = maskNucleotides(ai_FirstNucleotides, ai_FirstLength);
    ai_AppendNucleotides = maskNucleotides(ai_AppendNucleotides, ai_AppendLength);
    return li_Result | (ai_AppendNucleotides << (ai_FirstLength) * 3);
}
