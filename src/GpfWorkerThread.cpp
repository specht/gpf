#include "GpfWorkerThread.h"
#include "GpfQuery.h"


k_GpfWorkerThread::k_GpfWorkerThread(RefPtr<k_GpfQuery> ak_pQuery)
	: QThread(NULL)
	, mk_pQuery(ak_pQuery)
{
}


k_GpfWorkerThread::~k_GpfWorkerThread()
{
}


RefPtr<k_GpfQuery> k_GpfWorkerThread::get_Query()
{
	return mk_pQuery;
}


void k_GpfWorkerThread::run()
{
	mk_pQuery->Execute();
}
