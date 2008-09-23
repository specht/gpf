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
