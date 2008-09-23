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


class k_ValueEstimator
{
public:
	k_ValueEstimator();
	~k_ValueEstimator();

	void AddValue(double ad_Value, int ai_Count = 1);
	bool get_Estimate(double* ad_Result_);
private:
	static const int mi_StartReportingCount = 10;
	static const int mi_ShrinkTo = 100;
	static const int mi_ShrinkAt = 1000;

	QMutex mk_Mutex;
	double md_Sum;
	int mi_Count;
};
