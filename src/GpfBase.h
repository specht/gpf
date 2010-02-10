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

#define WATER_MASS 18.01057


struct r_DnaIndexChunkType
{
	enum Enumeration
	{
		Info = 100,
		Dna,
		Index
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
    qint64 reverseNucleotides(qint64 ai_Nucleotides, int ai_Length);
    qint64 invertNucleotides(qint64 ai_Nucleotides, int ai_Length);
    qint64 transposeNucleotides(qint64 ai_Nucleotides, int ai_Length);
	
	quint16 mk_DnaCharToNumber_[256];
	char mk_DnaTripletToAminoAcid_[512];
	// takes the same triplet as the array above, but reverses and transposes
	char mk_DnaTripletToAminoAcidReverse_[512];
	double md_AminoAcidMasses_[256];
	bool mb_IsAminoAcid_[256];
	int mi_AminoAcidToNumber_[256];
    char mc_NumberToAminoAcid_[19];
};


extern k_GpfBase gk_GpfBase;


void overwriteBitsInBuffer(quint8* auc_Buffer_, qint64 ai_BitOffset, quint64 aui_Value, int ai_Size);
quint64 readBitsFromBuffer(quint8* auc_Buffer_, qint64 ai_BitOffset, int ai_Size);
