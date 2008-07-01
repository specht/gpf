#include "GpfConfig.h"


#define SortFunction(as_Name, T) \
	void as_Name(T* ar_Entries_, int ai_First, int ai_Last) \
	{ \
		int li_PivotIndex, li_Mid, i, j; \
		T lr_Temp; \
		T lr_PivotValue; \
		if (ai_Last > ai_First) \
		{ \
			li_Mid = ai_First + ((ai_Last - ai_First) >> 1); \
 \
			if (lt(ar_Entries_[ai_First], ar_Entries_[li_Mid])) \
			{ \
				if (lt(ar_Entries_[li_Mid], ar_Entries_[ai_Last])) \
					li_PivotIndex = li_Mid; \
				else \
				{ \
					if (lt(ar_Entries_[ai_First], ar_Entries_[ai_Last])) \
						li_PivotIndex = ai_Last; \
					else \
						li_PivotIndex = ai_First; \
				} \
			} \
			else \
			{ \
				if (lt(ar_Entries_[ai_Last], ar_Entries_[li_Mid])) \
					li_PivotIndex = li_Mid; \
				else \
				{ \
					if (lt(ar_Entries_[ai_First], ar_Entries_[ai_Last])) \
						li_PivotIndex = ai_First; \
					else \
						li_PivotIndex = ai_Last; \
				} \
			} \
 \
			lr_Temp = ar_Entries_[li_PivotIndex]; \
			ar_Entries_[li_PivotIndex] = ar_Entries_[ai_Last]; \
			ar_Entries_[ai_Last] = lr_Temp; \
 \
			lr_PivotValue = lr_Temp; \
 \
			i = ai_First - 1; \
			j = ai_Last; \
 \
			do \
			{ \
				do ++i; while (lt(ar_Entries_[i], lr_PivotValue)); \
				do --j; while (ai_First < j && lt(lr_PivotValue, ar_Entries_[j])); \
 \
				if (i < j) \
				{ \
					lr_Temp = ar_Entries_[i]; \
					ar_Entries_[i] = ar_Entries_[j]; \
					ar_Entries_[j] = lr_Temp; \
				} \
			} while (i < j); \
 \
			lr_Temp = ar_Entries_[i]; \
			ar_Entries_[i] = ar_Entries_[ai_Last]; \
			ar_Entries_[ai_Last] = lr_Temp; \
			as_Name(ar_Entries_, ai_First, i - 1); \
			as_Name(ar_Entries_, i + 1, ai_Last); \
		} \
	} \


#define lt(first, second) \
	(first < second)
SortFunction(SortPlainUnsignedIntGno, unsigned int);
#undef lt


#define lt(first, second) \
	((first & GNO_BACKWARD_INVERSE) < (second & GNO_BACKWARD_INVERSE))
SortFunction(SortUnsignedIntGno, unsigned int);
#undef lt


#define lt(first, second) \
	(((first & GNO_BACKWARD) && (!(second & GNO_BACKWARD)))? 0: \
	((!(first & GNO_BACKWARD)) && (second & GNO_BACKWARD))? 1: \
	(!(first & GNO_BACKWARD))? \
		((first & GNO_BACKWARD_INVERSE) < (second & GNO_BACKWARD_INVERSE)): \
		((second & GNO_BACKWARD_INVERSE) < (first & GNO_BACKWARD_INVERSE)) \
	)
SortFunction(SortUnsignedIntDirectionGno, unsigned int);
#undef lt

