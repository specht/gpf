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
