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
#include "GpfDaemon.h"

class k_GpfDaemonThread: public QThread
{
	Q_OBJECT

public:
	k_GpfDaemonThread(int ai_SocketDescriptor, QObject* ak_Parent_, k_GpfDaemon& ak_GpfDaemon, k_GpfBase& ak_GpfBase);
	virtual ~k_GpfDaemonThread();

protected:
	void run();

private slots:
	void EnqueueQuery();
	void finishedQuery();

private:
	int mi_SocketDescriptor;
	k_GpfDaemon& mk_GpfDaemon;
	k_GpfBase& mk_GpfBase;
	QString ms_ResponseContent;
	QString ms_ResponseContentType;
	RefPtr<k_GpfQuery> mk_pGpfQuery;
};
