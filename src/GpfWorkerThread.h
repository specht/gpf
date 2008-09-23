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
