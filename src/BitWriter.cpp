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

#include "BitWriter.h"


k_BitWriter::k_BitWriter(QIODevice* ak_Device_)
	: mui_BufferSize(8 * 1024 * 1024)
	, mi_BufferOffset(0)
	, mi_BufferBitOffset(0)
	, muc_pBuffer(new quint8[mui_BufferSize])
	, mk_Device_(ak_Device_)
{
	memset(muc_pBuffer.get_Pointer(), 0, mui_BufferSize);
}


k_BitWriter::~k_BitWriter()
{
}


void k_BitWriter::writeBits(quint64 aui_Value, int ai_Bits)
{
	while (ai_Bits > 0)
	{
		int li_CopyBits = (8 - mi_BufferBitOffset);
		if (li_CopyBits > ai_Bits)
			li_CopyBits = ai_Bits;
		quint8 lui_Byte = aui_Value & ((1 << li_CopyBits) - 1);
		quint8 lui_NullMask = (((quint32)1) << li_CopyBits) - 1;
		lui_NullMask <<= mi_BufferBitOffset;
		lui_NullMask ^= 0xff;
		lui_Byte <<= mi_BufferBitOffset;
		aui_Value >>= li_CopyBits;
		ai_Bits -= li_CopyBits;
		muc_pBuffer.get_Pointer()[mi_BufferOffset] &= lui_NullMask;
		muc_pBuffer.get_Pointer()[mi_BufferOffset] |= lui_Byte;
		mi_BufferBitOffset += li_CopyBits;
		if (mi_BufferBitOffset > 7)
		{
			mi_BufferBitOffset -= 8;
			++mi_BufferOffset;
			if (mi_BufferOffset >= mui_BufferSize)
				flush();
		}
	}
}


void k_BitWriter::flush()
{
	size_t li_Size = mi_BufferOffset;
	if (mi_BufferBitOffset > 0)
		++li_Size;
	printf("writing %d bytes... to %p\n", li_Size, mk_Device_);
	mk_Device_->write((char*)muc_pBuffer.get_Pointer(), li_Size);
	printf("done.\n");
	mi_BufferOffset = 0;
	mi_BufferBitOffset = 0;
}
