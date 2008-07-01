#include "ValueEstimator.h"


k_ValueEstimator::k_ValueEstimator()
	: md_Sum(0)
	, mi_Count(0)
{
}


k_ValueEstimator::~k_ValueEstimator()
{
}


void k_ValueEstimator::AddValue(double ad_Value, int ai_Count)
{
	QMutexLocker lk_Locker(&mk_Mutex);
	md_Sum += ad_Value;
	mi_Count += ai_Count;

	// shrink estimate if necessary
	if (mi_Count >= mi_ShrinkAt)
	{
		md_Sum *= (double)mi_ShrinkTo / (double)mi_Count;
		mi_Count = mi_ShrinkTo;
	}
}


bool k_ValueEstimator::get_Estimate(double* ad_Result_)
{
	if (ad_Result_ == NULL)
		return false;

	QMutexLocker lk_Locker(&mk_Mutex);

	if (mi_Count < mi_StartReportingCount)
		return false;

	*ad_Result_ = md_Sum / mi_Count;
	return true;
}
