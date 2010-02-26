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

#include "GpfIndexFile.h"
#include <math.h>
#include "GpfBase.h"
#include <QtCore>


k_GpfIndexFile::k_GpfIndexFile(const QString& as_Path)
	: mb_IsGood(false)
{
	parseGpfIndexFile(as_Path);

    mi_MinAminoAcidMass = 0;
	// calculate mass according to mass accuracy
	for (int i = 0; i < 256; ++i)
	{
		mi_AminoAcidMasses_[i] = 0;
		if (gk_GpfBase.mb_IsAminoAcid_[i])
        {
            qint64 li_Mass = (qint64)(gk_GpfBase.md_AminoAcidMasses_[i] * mi_MassPrecision);
			mi_AminoAcidMasses_[i] = li_Mass;
            if (mi_MinAminoAcidMass == 0)
                mi_MinAminoAcidMass = li_Mass;
            else
            {
                if (li_Mass < mi_MinAminoAcidMass)
                    mi_MinAminoAcidMass = li_Mass;
            }
        }
	}
	mi_WaterMass = (qint64)(WATER_MASS * mi_MassPrecision);
    mi_MassDecimalDigits = (int)ceil(log(mi_MassPrecision) / log(10.0));
}


k_GpfIndexFile::~k_GpfIndexFile()
{
}


void k_GpfIndexFile::parseGpfIndexFile(const QString& as_Path)
{
    // generate short id from filename, remove colons and semicolons
    ms_ShortId = QFileInfo(as_Path).baseName();
    ms_ShortId.replace(":", "");
    ms_ShortId.replace(";", "");
    
	QFile lk_File(as_Path);
	if (!lk_File.open(QIODevice::ReadOnly))
	{
		printf("Error: Unable to open file.\n");
		return;
	}

	// check identifier
	QByteArray lk_Identifier = lk_File.read(8);
	if (lk_Identifier != "gpfindex")
	{
		printf("Error: gpfindex identifier missing.\n");
		return;
	}
	
	// check version
	memcpy(&mi_VersionMajor, lk_File.read(2).constData(), 2);
	memcpy(&mi_VersionMinor, lk_File.read(2).constData(), 2);
	if (mi_VersionMajor != 3 || mi_VersionMinor != 0)
	{
		printf("Error: Only version 3.0 is supported, this is %d.%d.\n", mi_VersionMajor, mi_VersionMinor);
		return;
	}
	
	while (!lk_File.atEnd())
	{
		quint32 lui_ChunkType;
		qint64 li_ChunkSize;
		memcpy(&lui_ChunkType, lk_File.read(4).constData(), 4);
		memcpy(&li_ChunkSize, lk_File.read(8).constData(), 8);
		if (lui_ChunkType == r_DnaIndexChunkType::Info)
		{
//             printf("Reading info chunk...\n");
			// read genome title
			qint32 li_TitleLength;
			memcpy(&li_TitleLength, lk_File.read(4).constData(), 4);
			QByteArray lk_Title = lk_File.read(li_TitleLength);
			lk_Title.append('\0');
			ms_Title = lk_Title;
			
			// read constants
			memcpy(&mi_OffsetBits, lk_File.read(4).constData(), 4);
			memcpy(&mi_MassBits, lk_File.read(4).constData(), 4);
			memcpy(&mi_MassPrecision, lk_File.read(4).constData(), 4);
			memcpy(&mi_TagSize, lk_File.read(4).constData(), 4);
			mi_HmstBits = mi_OffsetBits + mi_MassBits;
			mi_MaxMass = ((qint64)1 << mi_MassBits) - 1;
			mi_GnoBackwardsBit = (qint64)1 << (mi_OffsetBits - 1);
			
			// determine tag bits
			mi_TagCount =  1;
			for (int i = 0; i < mi_TagSize; ++i)
				mi_TagCount *= 19;
			
			// read scaffold sizes and labels
			qint32 li_ScaffoldCount;
			memcpy(&li_ScaffoldCount, lk_File.read(4).constData(), 4);
			mi_TotalNucleotideCount = 0;
			for (qint32 i = 0; i < li_ScaffoldCount; ++i)
			{
				qint64 li_ScaffoldSize;
				memcpy(&li_ScaffoldSize, lk_File.read(8).constData(), 8);
				qint32 li_ScaffoldLabelLength;
				memcpy(&li_ScaffoldLabelLength, lk_File.read(4).constData(), 4);
				QByteArray lk_Label = lk_File.read(li_ScaffoldLabelLength);
				lk_Label.append('\0');
				mk_ScaffoldLabels.append(QString(lk_Label));
				mk_ScaffoldStart.append(mi_TotalNucleotideCount);
				mk_ScaffoldLength.append(li_ScaffoldSize);
				mi_TotalNucleotideCount += li_ScaffoldSize;
			}
		}
		else if (lui_ChunkType == r_DnaIndexChunkType::Dna)
		{
//             printf("Reading DNA chunk...\n");
			muc_pDnaBuffer = RefPtr<quint8>(new quint8[li_ChunkSize]);
			// read DNA in 64M chunks
			qint64 li_Remaining = li_ChunkSize;
			qint64 li_Offset = 0;
			while (li_Remaining > 0)
			{
				qint64 li_Step = 64 * 1024 * 1024;
				if (li_Step > li_Remaining)
					li_Step = li_Remaining;
				memcpy(muc_pDnaBuffer.get_Pointer() + li_Offset, lk_File.read(li_Step).constData(), li_Step);
				li_Offset += li_Step;
				li_Remaining -= li_Step;
			}
		}
		else if (lui_ChunkType == r_DnaIndexChunkType::Index)
		{
//             printf("Reading index chunk...\n");
			qint64 li_OldPos = lk_File.pos();
			qint32 li_HmstCountBits;
			memcpy(&li_HmstCountBits, lk_File.read(4).constData(), 4);
			qint64 li_HmstCountEncodedBufferSize = mi_TagCount * 2 * li_HmstCountBits / 8 + 1;
			RefPtr<quint8> luc_pHmstCountEncoded(new quint8[li_HmstCountEncodedBufferSize]);
			memcpy(luc_pHmstCountEncoded.get_Pointer(), lk_File.read(li_HmstCountEncodedBufferSize).constData(), li_HmstCountEncodedBufferSize);
			mi_TotalHmstCount = 0;
			mi_BiggestHmstCount = 0;
			for (qint64 i = 0; i < mi_TagCount * 2; ++i)
			{
				qint64 li_HmstCount = readBitsFromBuffer(luc_pHmstCountEncoded.get_Pointer(), i * li_HmstCountBits, li_HmstCountBits);
				mk_HmstOffset.append(mi_TotalHmstCount);
				mk_HmstCount.append(li_HmstCount);
/*                if (li_HmstCount > 1)
                    printf("HMST COUNT %d\n", (unsigned int)li_HmstCount);*/
				mi_TotalHmstCount += li_HmstCount;
				if (li_HmstCount > mi_BiggestHmstCount)
					mi_BiggestHmstCount = li_HmstCount;
			}
			mi_IndexFilePosition = lk_File.pos();
			lk_File.seek(li_OldPos + li_ChunkSize);
		}
		else
		{
			printf("skipping type %d (%d bytes)\n", lui_ChunkType, (quint32)li_ChunkSize);
			lk_File.seek(lk_File.pos() + li_ChunkSize);
		}
	}
	
	lk_File.close();
	
	mk_File.setFileName(as_Path);
	// re-open the file unbuffered so that reading 4 bytes really means 
	// reading 4 bytes and not a few KiB more for good measure, because
	// we're in a bit of a hurry
	mk_File.open(QIODevice::ReadOnly | QIODevice::Unbuffered);
	
	mb_IsGood = true;
}


bool k_GpfIndexFile::isGood()
{
	return mb_IsGood;
}

#define READ_BITS 32
#define READ_TYPE quint32

quint64 k_GpfIndexFile::readIndexBits(qint64 ai_BitOffset, qint32 ai_Size)
{
	quint64 lui_Result = 0;
	int li_BitsCopied = 0;
	while (ai_Size > 0)
	{
		qint64 li_ByteOffset = ai_BitOffset / READ_BITS;
		int li_BitOffset = ai_BitOffset % READ_BITS;
		int li_CopyBits = READ_BITS - li_BitOffset;
		if (li_CopyBits > ai_Size)
			li_CopyBits = ai_Size;
		READ_TYPE lui_CopyMask = ((((quint64)1) << li_CopyBits) - 1);
		//READ_TYPE lui_CopyByte = (auc_Buffer_[li_ByteOffset] >> li_BitOffset) & lui_CopyMask;
		READ_TYPE lui_CopyByte;
		mk_File.seek(mi_IndexFilePosition + li_ByteOffset * 4);
		memcpy(&lui_CopyByte, mk_File.read(4).constData(), 4);
		lui_CopyByte >>= li_BitOffset;
		lui_CopyByte &= lui_CopyMask;
		lui_Result |= (((quint64)lui_CopyByte) << li_BitsCopied);
		ai_BitOffset += li_CopyBits;
		ai_Size -= li_CopyBits;
		li_BitsCopied += li_CopyBits;
	}
	return lui_Result;
}
