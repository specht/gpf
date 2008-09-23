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

#ifdef READ_FROM_RAM

#define READ_UINT32(aui_Target, ai_Offset) \
	aui_Target = *(unsigned int*)(muc_IndexFileBuffer_ + (ai_Offset));

#define READ_UINT32_CAST_INT32(ai_Target, ai_Offset) \
	ai_Target = (int)(*(unsigned int*)(muc_IndexFileBuffer_ + (ai_Offset)));

#else

#define READ_UINT32(aui_Target, ai_Offset) \
	{ \
		mk_IndexFile.seek(ai_Offset); \
		mk_IndexFile.read((char*)&aui_Target, 4); \
	}

#define READ_UINT32_CAST_INT32(ai_Target, ai_Offset) \
	{ \
		unsigned int lui_Value; \
		mk_IndexFile.seek(ai_Offset); \
		mk_IndexFile.read((char*)&lui_Value, 4); \
		ai_Target = (int)lui_Value; \
	}

#endif
