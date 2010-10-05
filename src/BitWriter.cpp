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

#include "BitWriter.h"


k_BitWriter::k_BitWriter(QIODevice* ak_Device_)
    : mui_BufferSize(8 * 1024 * 1024)
    , mi_BufferOffset(0)
    , mi_BufferBitOffset(0)
    , muc_pBuffer(new quint8[mui_BufferSize])
    , mk_Device_(ak_Device_)
{
    memset(muc_pBuffer.data(), 0, mui_BufferSize);
}


k_BitWriter::~k_BitWriter()
{
}

#define READ_BITS 32
#define READ_TYPE quint32

void k_BitWriter::writeBits(quint64 aui_Value, int ai_Bits)
{
    while (ai_Bits > 0)
    {
        int li_CopyBits = READ_BITS - mi_BufferBitOffset;
        if (li_CopyBits > ai_Bits)
            li_CopyBits = ai_Bits;
        READ_TYPE lui_Byte = aui_Value & ((1 << li_CopyBits) - 1);
        READ_TYPE lui_NullMask = (((quint64)1) << li_CopyBits) - 1;
        lui_NullMask <<= mi_BufferBitOffset;
        lui_NullMask ^= (((quint64)1) << READ_BITS) - 1;
        lui_Byte <<= mi_BufferBitOffset;
        aui_Value >>= li_CopyBits;
        ai_Bits -= li_CopyBits;
        ((READ_TYPE*)(muc_pBuffer.data()))[mi_BufferOffset] &= lui_NullMask;
        ((READ_TYPE*)(muc_pBuffer.data()))[mi_BufferOffset] |= lui_Byte;
        mi_BufferBitOffset += li_CopyBits;
        if (mi_BufferBitOffset >= READ_BITS)
        {
            mi_BufferBitOffset -= READ_BITS;
            ++mi_BufferOffset;
            if (mi_BufferOffset * READ_BITS / 8 >= mui_BufferSize)
                flush();
        }
    }
}


void k_BitWriter::flush()
{
    size_t li_Size = mi_BufferOffset * READ_BITS / 8;
    while (mi_BufferBitOffset > 0)
    {
        ++li_Size;
        if (mi_BufferBitOffset >= 8)
            mi_BufferBitOffset -= 8;
        else
            mi_BufferBitOffset = 0;
    }
    mk_Device_->write((char*)muc_pBuffer.data(), li_Size);
    mi_BufferOffset = 0;
    mi_BufferBitOffset = 0;
}
