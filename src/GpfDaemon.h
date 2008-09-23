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
#include <QtNetwork>
#include "GpfBase.h"
#include "RefPtr.h"
#include "QueryQueue.h"
#include "GpfWorkerThread.h"
#include "StopWatch.h"
#include "ValueEstimator.h"


class k_GpfDaemon: public QTcpServer
{
	Q_OBJECT

public:
	k_GpfDaemon(QString as_HostName, quint16 auw_Port, QStringList ak_IndexFiles);
	virtual ~k_GpfDaemon();

	QString get_HostName() const;
	quint16 get_Port() const;
	k_QueryQueue& get_QueryQueue();
	k_GpfBase& get_GpfBase();
	k_ValueEstimator& get_SpeedEstimator();

	QString get_UniqueTicket();

	void DelegateQuery();

public slots:
	void notifyFinished();

protected:
	void incomingConnection(int ai_SocketDescriptor);

private:
	QMutex mk_DelegateQueryMutex;
	QMutex mk_WorkerThreadsMutex;
	QMutex mk_BatchQueriesMutex;
	QMutex mk_UniqueTicketMutex;

	QString ms_HostName;
	quint16 muw_Port;
	unsigned int mui_MaxWorkerThreadCount;
	k_GpfBase mk_GpfBase;

	k_QueryQueue mk_QueryQueue;
	QMap<QObject*, RefPtr<k_GpfWorkerThread> > mk_WorkerThreads;
	QMap<k_GpfQuery*, RefPtr<k_QueryBatch> > mk_BatchQueries;
	k_StopWatch mk_StopWatch;
	unsigned int mui_TicketCounter;
	k_ValueEstimator mk_SpeedEstimator;
};
