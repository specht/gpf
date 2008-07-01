#pragma once

#include <QtCore>
#include <QMutex>
#include "GpfQuery.h"
#include "RefPtr.h"
#include "StopWatch.h"
#include "QueryBatch.h"


class k_GpfDaemon;
static const double gd_LowPriorityQueryPause = 5.0;


class k_QueryQueue
{
public:
	k_QueryQueue(k_GpfDaemon& ak_GpfDaemon);
	virtual ~k_QueryQueue();

	void Enqueue(RefPtr<k_GpfQuery> ak_pQuery);
	void AddBatch(QString as_Ticket);
	RefPtr<k_GpfQuery> get_NextQuery(RefPtr<k_QueryBatch>* ak_pQueryBatchPicked_);
	unsigned int get_Length();
	unsigned int get_BatchListLength();
	QString QueryBatchStatus(QString ls_Ticket);

	QMutex mk_BatchListMutex;

private:
	QList<RefPtr<k_GpfQuery> > mk_LowPriorityQueries;
	QList<RefPtr<k_QueryBatch> > mk_QueryBatches;
	QHash<QString, RefPtr<k_QueryBatch> > mk_QueryBatchByTicket;

	QMutex mk_QueueMutex;
	k_StopWatch mk_LowPriorityStopWatch;
	k_GpfDaemon& mk_GpfDaemon;

	// for batch processing remaining time estimation...
	k_StopWatch mk_BatchPerformanceStopWatch;
	bool mb_BatchPerformanceWatching;
	unsigned int mi_BatchPerformanceCount;
};
