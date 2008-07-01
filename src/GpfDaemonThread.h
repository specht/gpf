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
