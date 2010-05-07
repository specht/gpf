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


class k_BitWriter
{
public:
    k_BitWriter(QIODevice* ak_Device_);
    virtual ~k_BitWriter();
    
    void writeBits(quint64 aui_Value, int ai_Bits);
    void flush();
    
protected:
    size_t mui_BufferSize;
    size_t mi_BufferOffset;
    size_t mi_BufferBitOffset;
    QSharedPointer<quint8> muc_pBuffer;
    QIODevice* mk_Device_;
};
