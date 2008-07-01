#include "QueryQueue.h"
#include "GpfDaemon.h"


k_QueryQueue::k_QueryQueue(k_GpfDaemon& ak_GpfDaemon)
	: mk_GpfDaemon(ak_GpfDaemon)
	, mb_BatchPerformanceWatching(false)
	, mi_BatchPerformanceCount(-1)
{
}


k_QueryQueue::~k_QueryQueue()
{
}


void k_QueryQueue::Enqueue(RefPtr<k_GpfQuery> ak_pQuery)
{
	QMutexLocker lk_Locker(&mk_QueueMutex);
	mk_LowPriorityQueries.append(ak_pQuery);
}


void k_QueryQueue::AddBatch(QString as_Ticket)
{
	QMutexLocker lk_Locker(&mk_BatchListMutex);
	RefPtr<k_QueryBatch> lk_pQueryBatch(new k_QueryBatch(mk_GpfDaemon, as_Ticket));
	mk_QueryBatches.append(lk_pQueryBatch);
	mk_QueryBatchByTicket[as_Ticket] = lk_pQueryBatch;
}


unsigned int k_QueryQueue::get_Length()
{
	QMutexLocker lk_Locker(&mk_QueueMutex);
	return mk_LowPriorityQueries.size();
}


unsigned int k_QueryQueue::get_BatchListLength()
{
	QMutexLocker lk_Locker(&mk_BatchListMutex);
	return mk_QueryBatches.size();
}


QString k_QueryQueue::QueryBatchStatus(QString ls_Ticket)
{
	QString ls_Status("response: error\n");

	QMutexLocker lk_Locker(&mk_BatchListMutex);
	if (mk_QueryBatchByTicket.contains(ls_Ticket))
	{
		RefPtr<k_QueryBatch> lk_pQueryBatch = mk_QueryBatchByTicket[ls_Ticket];
		switch (lk_pQueryBatch->get_State())
		{
		case r_QueryBatchState::Waiting:
			ls_Status = QString("response: ok\nstate: waiting\n");
			ls_Status += QString("totalQueries: %1\n")
				.arg(lk_pQueryBatch->get_ExpectedQueryCount());
			break;
		case r_QueryBatchState::Executing:
			{
				ls_Status = QString("response: ok\nstate: executing\n");
				ls_Status += QString("totalQueries: %1\nprocessedQueries: %2\n")
					.arg(lk_pQueryBatch->get_ExpectedQueryCount())
					.arg(lk_pQueryBatch->get_ProcessedQueryCount());

				/*
				double ld_EstimatedSpeed = 0;
				if (mk_GpfDaemon.get_SpeedEstimator().get_Estimate(&ld_EstimatedSpeed))
				{
					double ld_EstimatedRemainingTime = ld_EstimatedSpeed * lk_pQueryBatch->get_RemainingQueryCount();
					QString ls_EstimatedRemainingTime = k_StopWatch::getTimeAsString(ld_EstimatedRemainingTime);
					ls_Status += QString("estimatedDuration: { seconds: %1, verbose: %2 }\n")
						.arg(ld_EstimatedRemainingTime)
						.arg(ls_EstimatedRemainingTime);
				}
				*/

				break;
			}
		case r_QueryBatchState::Finished:
			ls_Status = QString("response: ok\nstate: finished\n");
			break;
		}
	}
	else
	{
		// the query batch cannot be found but maybe a result file is waiting to be downloaded
		if (QFile::exists(mk_GpfDaemon.get_DownloadPath() + ls_Ticket + ".yaml"))
			ls_Status = QString("response: ok\nstate: finished\n");
	}

	return ls_Status;
}


RefPtr<k_GpfQuery> k_QueryQueue::get_NextQuery(RefPtr<k_QueryBatch>* ak_pQueryBatchPicked_)
{
	RefPtr<k_GpfQuery> lk_pQuery(NULL);
	*ak_pQueryBatchPicked_ = RefPtr<k_QueryBatch>(NULL);

	QMutexLocker lk_Locker(&mk_QueueMutex);
	QMutexLocker lk_BatchListLocker(&mk_BatchListMutex);

	// clean up high-priority batch list
	while (!mk_QueryBatches.empty() && mk_QueryBatches.first()->get_RemainingQueryCount() == 0)
	{
		mk_QueryBatchByTicket.remove(mk_QueryBatches.first()->get_Ticket());
		mk_QueryBatches.removeFirst();
	}

	if (mk_LowPriorityQueries.empty())
	{
		// no low priority queries, pick the next high priority query, if available
		if (!mk_QueryBatches.empty())
		{
			// TODO!!! ATTENTION ATTENTION THERE IS A POTENTIAL GLITCH HERE
			// What if the next query is null? Eh?
			lk_pQuery = mk_QueryBatches.first()->get_NextQuery();
			*ak_pQueryBatchPicked_ = mk_QueryBatches.first();
		}
	}
	else
	{
		// there are low priority queries waiting, do time check
		bool lb_PickLowPriorityQuery = true;
		if (!mk_QueryBatches.empty())
		{
			// there is also at least one high priority query waiting
			if (mk_LowPriorityStopWatch.get_Time() < gd_LowPriorityQueryPause)
				lb_PickLowPriorityQuery = false;
		}

		if (lb_PickLowPriorityQuery)
		{
			lk_pQuery = mk_LowPriorityQueries.takeFirst();
			mk_LowPriorityStopWatch.reset();
		}
		else
		{
			lk_pQuery = mk_QueryBatches.first()->get_NextQuery();
			*ak_pQueryBatchPicked_ = mk_QueryBatches.first();
		}
	}

	/*
	// watch batch processing performance
	if (mb_BatchPerformanceWatching)
	{
		++mi_BatchPerformanceCount;
		if (mk_QueryBatches.empty())
		{
			mb_BatchPerformanceWatching = false;
			mk_GpfDaemon.get_SpeedEstimator().AddValue(mk_BatchPerformanceStopWatch.get_Time(), mi_BatchPerformanceCount);
			mi_BatchPerformanceCount = -1;
			mk_BatchPerformanceStopWatch.reset();
		}
		else
		{
			if (mi_BatchPerformanceCount > 40)
			{
				mk_GpfDaemon.get_SpeedEstimator().AddValue(mk_BatchPerformanceStopWatch.get_Time(), mi_BatchPerformanceCount);
				mi_BatchPerformanceCount = -1;
				mk_BatchPerformanceStopWatch.reset();
			}
		}
	}
	else
	{
		if (!mk_QueryBatches.empty())
		{
			mb_BatchPerformanceWatching = true;
			mi_BatchPerformanceCount = -1;
			mk_BatchPerformanceStopWatch.reset();
		}
	}
	*/

	return lk_pQuery;
}
