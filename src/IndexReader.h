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
