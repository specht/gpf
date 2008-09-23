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
