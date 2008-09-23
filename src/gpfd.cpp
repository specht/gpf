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

#include <QCoreApplication>
#include <stdio.h>
#include "GpfDaemon.h"


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
	if (ai_ArgumentCount < 4)
	{
		printf("Usage: gpfd [own IP or hostname] [port] [GPF index files]\n");
		exit(1);
	}

	QString ls_HostName = ac_Arguments__[1];
	QString ls_Port = ac_Arguments__[2];
	QStringList lk_IndexFiles;
	for (int i = 3; i < ai_ArgumentCount; ++i)
		lk_IndexFiles += ac_Arguments__[i];

	QCoreApplication lk_Application(ai_ArgumentCount, ac_Arguments__);
	k_GpfDaemon lk_Daemon(ls_HostName, ls_Port.toUShort(), lk_IndexFiles);

	return lk_Application.exec();
}
