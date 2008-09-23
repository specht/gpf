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

#include "GpfDaemon.h"
#include "GpfDaemonThread.h"
#include "GpfQuery.h"
#include <time.h>


k_GpfDaemon::k_GpfDaemon(QString as_HostName, quint16 auw_Port, QStringList ak_IndexFiles)
	: QTcpServer(NULL)
	, ms_HostName(as_HostName)
	, muw_Port(auw_Port)
	, mui_MaxWorkerThreadCount(2)
	, mk_GpfBase(ak_IndexFiles)
	, mk_QueryQueue(*this)
	, mui_TicketCounter(100)
{
	if (!listen(QHostAddress::Any, muw_Port))
	{
		printf("Error: Unable to listen on port %d.\n", muw_Port);
		exit(1);
	}
	srand(time(NULL));
}


k_GpfDaemon::~k_GpfDaemon()
{
}


quint16 k_GpfDaemon::get_Port() const
{
	return muw_Port;
}


QString k_GpfDaemon::get_HostName() const
{
	return ms_HostName;
}


k_QueryQueue& k_GpfDaemon::get_QueryQueue()
{
	return mk_QueryQueue;
}


k_GpfBase& k_GpfDaemon::get_GpfBase()
{
	return mk_GpfBase;
}


k_ValueEstimator& k_GpfDaemon::get_SpeedEstimator()
{
	return mk_SpeedEstimator;
}


QString k_GpfDaemon::get_UniqueTicket()
{
	QMutexLocker lk_Locker(&mk_UniqueTicketMutex);
	++mui_TicketCounter;
	return QString("x%1%2").arg(mui_TicketCounter, 0, 36).arg((int)(mk_StopWatch.get_Time() * 1000.0), 0, 36);
}


void k_GpfDaemon::DelegateQuery()
{
	QMutexLocker lk_FunctionLocker(&mk_DelegateQueryMutex);

	// delegate a query

	// only delegate a query if there is a worker thread available
	{
		QMutexLocker lk_Locker(&mk_WorkerThreadsMutex);
		if (mk_WorkerThreads.size() >= (int)mui_MaxWorkerThreadCount)
			return;
	}

	RefPtr<k_QueryBatch> lk_pQueryBatch;
	RefPtr<k_GpfQuery> lk_pQuery = mk_QueryQueue.get_NextQuery(&lk_pQueryBatch);

	if (lk_pQuery != NULL)
	{
		if (lk_pQueryBatch.get_Pointer() != NULL)
		{
			QMutexLocker lk_Locker(&mk_BatchQueriesMutex);
			mk_BatchQueries[lk_pQuery.get_Pointer()] = lk_pQueryBatch;
		}

		RefPtr<k_GpfWorkerThread> lk_pWorkerThread(new k_GpfWorkerThread(lk_pQuery));
		{
			QMutexLocker lk_Locker(&mk_WorkerThreadsMutex);
			mk_WorkerThreads[lk_pWorkerThread.get_Pointer()] = lk_pWorkerThread;
		}
		connect(lk_pWorkerThread.get_Pointer(), SIGNAL(finished()), this, SLOT(notifyFinished()));
		lk_pWorkerThread->start();
	}
}


void k_GpfDaemon::notifyFinished()
{
	// if this was a batch query, send result to appropriate query batch object
	/*
	RefPtr<k_GpfQuery> lk_pQuery = static_cast<k_GpfWorkerThread*>(sender())->get_Query();
	{
		QMutexLocker lk_Locker(&mk_BatchQueriesMutex);
		if (mk_BatchQueries.contains(lk_pQuery.get_Pointer()))
		{
			QString ls_Result = lk_pQuery->get_Result();
			// indent result
			ls_Result = "    " + ls_Result.replace("\n", "\n    ");
			QString ls_Key = QString("  '%1':").arg(lk_pQuery->get_Label());
			ls_Result = ls_Key + "\n" + ls_Result + "\n";
			k_GpfQuery* lk_GpfQuery_ = lk_pQuery.get_Pointer();
			RefPtr<k_QueryBatch> lk_pQueryBatch = mk_BatchQueries[lk_GpfQuery_];
			lk_pQueryBatch->AddResult(ls_Result);
			mk_BatchQueries.remove(lk_GpfQuery_);
		}
	}
	*/

	// remove worker thread from mk_WorkerThreads
	{
		QMutexLocker lk_Locker(&mk_WorkerThreadsMutex);
		mk_WorkerThreads.remove(sender());
	}

	// start next query processing, if necessary
	DelegateQuery();
}


void k_GpfDaemon::incomingConnection(int ai_SocketDescriptor)
{
	k_GpfDaemonThread* lk_Thread_ = new k_GpfDaemonThread(ai_SocketDescriptor, this, *this, mk_GpfBase);
	connect(lk_Thread_, SIGNAL(finished()), lk_Thread_, SLOT(deleteLater()));
	lk_Thread_->start();
}
