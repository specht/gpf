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
