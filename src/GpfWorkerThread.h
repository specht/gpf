#pragma once

#include <QtCore>
#include "RefPtr.h"


class k_GpfQuery;


class k_GpfWorkerThread: public QThread
{
	Q_OBJECT

public:
	k_GpfWorkerThread(RefPtr<k_GpfQuery> ak_pQuery);
	virtual ~k_GpfWorkerThread();
	RefPtr<k_GpfQuery> get_Query();

protected:
	virtual void run();

private:
	RefPtr<k_GpfQuery> mk_pQuery;
};
