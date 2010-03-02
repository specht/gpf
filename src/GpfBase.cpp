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

#include "GpfBase.h"


k_GpfBase gk_GpfBase;


k_GpfBase::k_GpfBase()
{
    Q_INIT_RESOURCE(libgpf);
    
    for (int i = 0; i < 256; ++i)
        mk_DnaCharToNumber_[i] = 7;
    mk_DnaCharToNumber_[(unsigned char)'A'] = 0;
    mk_DnaCharToNumber_[(unsigned char)'C'] = 1;
    mk_DnaCharToNumber_[(unsigned char)'G'] = 2;
    mk_DnaCharToNumber_[(unsigned char)'T'] = 3;
    mk_DnaCharToNumber_[(unsigned char)'a'] = 0;
    mk_DnaCharToNumber_[(unsigned char)'c'] = 1;
    mk_DnaCharToNumber_[(unsigned char)'g'] = 2;
    mk_DnaCharToNumber_[(unsigned char)'t'] = 3;
    
    {
        QFile lk_File(":ext/proteomics-knowledge-base/genetic-codes.txt");
        if (!lk_File.open(QIODevice::ReadOnly))
        {
            printf("Internal error: Unable to open translation table.\n");
            exit(1);
        }
        QTextStream lk_Stream(&lk_File);
        while (!lk_Stream.atEnd())
        {
            QString ls_Line = lk_Stream.readLine().trimmed();
            // skip comments or empty lines
            if (ls_Line.startsWith("#") || ls_Line.isEmpty())
                continue;
            
            // expect "1. The Standard Code (transl_table=1)"
            QString ls_Title = ls_Line;
            ls_Title = ls_Title.replace(QRegExp("\\d+\\.\\s+"), "");
            
            QRegExp lk_RegExp("\\(transl_table=(\\d+)\\)");
            if (lk_RegExp.indexIn(ls_Title) == -1)
            {
                printf("Internal error: Unable to extract transl_table id.\n");
                exit(1);
            }
            int li_Id = lk_RegExp.cap(1).toInt();
            ls_Title = ls_Title.replace(lk_RegExp, "").trimmed();
            
            int li_LinesRead = 0;
            // read five non-empty lines now
            QHash<QString, QString> lk_Lines;
            while (!lk_Stream.atEnd())
            {
                ls_Line = lk_Stream.readLine().trimmed();
                // skip comments or empty lines
                if (ls_Line.startsWith("#") || ls_Line.isEmpty())
                    continue;
                ++li_LinesRead;
                QStringList lk_Line = ls_Line.split("=");
                QString ls_Key = lk_Line[0].trimmed();
                QString ls_Value = lk_Line[1].trimmed();
                if (ls_Value.length() != 64)
                {
                    printf("Internal error: line length is not 64 chars in translation table.\n");
                    printf("The offending line is %s.\n", ls_Value.toStdString().c_str());
                    exit(1);
                }
                lk_Lines[ls_Key.toLower()] = ls_Value;
                if (li_LinesRead >= 5)
                    break;
            }
            
            mk_TranslationTableTitle[li_Id] = ls_Title;
            mk_TranslationTables[li_Id] = RefPtr<char>(new char[512]);
            memset(mk_TranslationTables[li_Id].get_Pointer(), 'X', 512);
            
            QHash<int, QSet<char> > lk_Ambiguities;
            
            for (int i = 0; i < 64; ++i)
            {
                int a = mk_DnaCharToNumber_[(unsigned char)lk_Lines["base1"].at(i).toAscii()];
                int b = mk_DnaCharToNumber_[(unsigned char)lk_Lines["base2"].at(i).toAscii()];
                int c = mk_DnaCharToNumber_[(unsigned char)lk_Lines["base3"].at(i).toAscii()];
                char lc_AminoAcid = lk_Lines["aas"].at(i).toAscii();
                mk_TranslationTables[li_Id].get_Pointer()[(a) | (b << 3) | (c << 6)] = lc_AminoAcid;
                // try all combinations of one or two unknown nucleotides
                for (int k = 1; k < 7; ++k)
                {
                    int sa = a;
                    int sb = b;
                    int sc = c;
                    if ((k & 1) != 0)
                        sa = 4;
                    if ((k & 2) != 0)
                        sb = 4;
                    if ((k & 4) != 0)
                        sc = 4;
                    int li_Code = (sa) | (sb << 3) | (sc << 6);
                    if (!lk_Ambiguities.contains(li_Code))
                        lk_Ambiguities[li_Code] = QSet<char>();
                    lk_Ambiguities[li_Code].insert(lc_AminoAcid);
                }
            }
            
            // :TODO: add incomplete triplets that still yield useful amino acids
            foreach (int li_Code, lk_Ambiguities.keys())
            {
                if (lk_Ambiguities[li_Code].size() == 1)
                {
                    // this incomplete triplet is unambiguous, so add it
                    char lc_AminoAcid = lk_Ambiguities[li_Code].toList().first();
                    int a = (li_Code >> 0) & 7;
                    int b = (li_Code >> 3) & 7;
                    int c = (li_Code >> 6) & 7;
                    int a0 = a;
                    int a1 = a;
                    int b0 = b;
                    int b1 = b;
                    int c0 = c;
                    int c1 = c;
                    if (a == 4)
                    {
                        a0 = 4;
                        a1 = 7;
                    }
                    if (b == 4)
                    {
                        b0 = 4;
                        b1 = 7;
                    }
                    if (c == 4)
                    {
                        c0 = 4;
                        c1 = 7;
                    }
                    for (int a = a0; a <= a1; ++a)
                        for (int b = b0; b <= b1; ++b)
                            for (int c = c0; c <= c1; ++c)
                                mk_TranslationTables[li_Id].get_Pointer()[(a) | (b << 3) | (c << 6)] = lc_AminoAcid;
                }
            }
            

            // add reverse table
            mk_TranslationTablesReverse[li_Id] = RefPtr<char>(new char[512]);
            memset(mk_TranslationTablesReverse[li_Id].get_Pointer(), 'X', 512);
            for (int i = 0; i < 512; ++i)
            {
                int a = i & 7;
                int b = (i >> 3) & 7;
                int c = (i >> 6) & 7;
                int li_Reverse = (a << 6) | (b << 3) | c;
                li_Reverse ^= 219;
                mk_TranslationTablesReverse[li_Id].get_Pointer()[i] = mk_TranslationTables[li_Id].get_Pointer()[li_Reverse];
            }
        }
        lk_File.close();
    }
    
    for (int i = 0; i < 256; ++i)
        mi_AminoAcidToNumber_[i] = -1;
    
    {
        QFile lk_File(":ext/proteomics-knowledge-base/amino-acids.csv");
        if (!lk_File.open(QIODevice::ReadOnly))
        {
            printf("Internal error: Unable to open amino acid table.\n");
            exit(1);
        }
        QTextStream lk_Stream(&lk_File);
        lk_Stream.readLine();
        for (int i = 0; i < 256; ++i)
        {
            md_AminoAcidMasses_[i] = 0.0;
            mb_IsAminoAcid_[i] = false;
        }
        while (!lk_Stream.atEnd())
        {
            QString ls_Line = lk_Stream.readLine().trimmed();
            QStringList lk_Line = ls_Line.split(",");
            QString ls_AminoAcid = lk_Line[3];
            QString ls_Mass = lk_Line[4];
            int li_Code = QVariant(lk_Line[0]).toInt();
            if (li_Code > 7)
                li_Code -= 1;
            if (li_Code <= 18)
            {
                mi_AminoAcidToNumber_[(int)ls_AminoAcid[0].toAscii()] = li_Code;
                mc_NumberToAminoAcid_[li_Code] = ls_AminoAcid[0].toAscii();
                mb_IsAminoAcid_[(int)ls_AminoAcid[0].toAscii()] = true;
                md_AminoAcidMasses_[(int)ls_AminoAcid[0].toAscii()] = QVariant(ls_Mass).toDouble();
            }
        }
        lk_File.close();
    }
}


k_GpfBase::~k_GpfBase()
{
}


quint16 k_GpfBase::readNucleotideTriplet(quint8* auc_Buffer_, quint64 aui_Gno)
{
    // 9 bits are always contained in at most 2 bytes!
    quint64 lui_Offset = (aui_Gno * 3) / 8;
    quint16 lui_Bit;
    memcpy(&lui_Bit, auc_Buffer_ + lui_Offset, 2);
    lui_Bit >>= ((aui_Gno * 3) & 7);
    lui_Bit &= 511;
    return lui_Bit;
}


int k_GpfBase::aminoAcidPolymerCode(const char* ac_Buffer_, int ai_Length)
{
    int li_Result = 0;
    for (int i = 0; i < ai_Length; ++i)
    {
        li_Result *= 19;
        int li_AminoAcidNumber = mi_AminoAcidToNumber_[(int)ac_Buffer_[i]];
        if (li_AminoAcidNumber < 0)
            return -1;
        li_Result += li_AminoAcidNumber;
    }
    return li_Result;
}


QString k_GpfBase::aminoAcidSequenceForCode(int ai_Code, int ai_Length)
{
    QString ls_Result;
    for (int i = 0; i < ai_Length; ++i)
    {
        ls_Result = QChar(mc_NumberToAminoAcid_[ai_Code % 19]) + ls_Result;
        ai_Code /= 19;
    }
    return ls_Result;
}


QString k_GpfBase::nucleotideSequenceForCode(int ai_Code, int ai_Length)
{
    QString ls_Result;
    for (int i = 0; i < ai_Length; ++i)
    {
        switch (ai_Code & 7)
        {
            case 0: ls_Result += "A"; break;
            case 1: ls_Result += "C"; break;
            case 2: ls_Result += "G"; break;
            case 3: ls_Result += "T"; break;
            default: ls_Result += "X"; break;
        };
        ai_Code >>= 3;
    }
    return ls_Result;
}
