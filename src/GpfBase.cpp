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
		QFile lk_File(":res/DnaToAminoAcid.csv");
		lk_File.open(QIODevice::ReadOnly);
		QTextStream lk_Stream(&lk_File);
		memset(mk_DnaTripletToAminoAcid_, 'X', 512);
		while (!lk_Stream.atEnd())
		{
			QString ls_Line = lk_Stream.readLine().trimmed();
			QStringList lk_Line = ls_Line.split(";");
			QString ls_Triplet = lk_Line[0];
			ls_Triplet.replace("U", "T");
			QString ls_AminoAcid = lk_Line.last();
			int a = mk_DnaCharToNumber_[(unsigned char)ls_Triplet.at(0).toAscii()];
			int b = mk_DnaCharToNumber_[(unsigned char)ls_Triplet.at(1).toAscii()];
			int c = mk_DnaCharToNumber_[(unsigned char)ls_Triplet.at(2).toAscii()];
			
			// the next lines make sure that when an unknown nucleotide pops up,
			// in some cases the correct amino acid is nevertheless returned,
			// even if the unknown is not encoded by 111 but any of 1xx.
			int a0 = a;
			int a1 = a;
			if (a > 3)
			{
				a0 = 4; a1 = 7;
			}
			int b0 = b;
			int b1 = b;
			if (b > 3)
			{
				b0 = 4; b1 = 7;
			}
			int c0 = c;
			int c1 = c;
			if (c > 3)
			{
				c0 = 4; c1 = 7;
			}
			for (int ia = a0; ia <= a1; ++ia)
				for (int ib = b0; ib <= b1; ++ib)
					for (int ic = c0; ic <= c1; ++ic)
						mk_DnaTripletToAminoAcid_[(ia) | (ib << 3) | (ic << 6)] = ls_AminoAcid.at(0).toAscii();
		}
		lk_File.close();
	}
	
	// build reverse triplet translation table
	for (int i = 0; i < 512; ++i)
	{
		int a = i & 7;
		int b = (i >> 3) & 7;
		int c = (i >> 6) & 7;
		int li_Reverse = (a << 6) | (b << 3) | c;
		li_Reverse ^= 219;
		mk_DnaTripletToAminoAcidReverse_[i] = mk_DnaTripletToAminoAcid_[li_Reverse];
	}
	
	for (int i = 0; i < 256; ++i)
		mi_AminoAcidToNumber_[i] = 20;
	
	{
		QFile lk_File(":res/AminoAcids.csv");
		lk_File.open(QIODevice::ReadOnly);
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
			QStringList lk_Line = ls_Line.split(";");
			QString ls_AminoAcid = lk_Line[3];
			QString ls_Mass = lk_Line[4];
			mb_IsAminoAcid_[(int)ls_AminoAcid[0].toAscii()] = true;
			md_AminoAcidMasses_[(int)ls_AminoAcid[0].toAscii()] = QVariant(ls_Mass).toDouble();
			mi_AminoAcidToNumber_[(int)ls_AminoAcid[0].toAscii()] = QVariant(lk_Line[0]).toInt();
		}
		lk_File.close();
	}
}


k_GpfBase::~k_GpfBase()
{
}


int k_GpfBase::aminoAcidPolymerCode(const char* ac_Buffer_, int ai_Length)
{
	int li_Result = 0;
	for (int i = 0; i < ai_Length; ++i)
	{
		li_Result *= 20;
		int li_AminoAcidNumber = mi_AminoAcidToNumber_[(int)ac_Buffer_[i]];
		if (li_AminoAcidNumber < 0 || li_AminoAcidNumber > 19)
			return -1;
		li_Result += li_AminoAcidNumber;
	}
	return li_Result;
}
