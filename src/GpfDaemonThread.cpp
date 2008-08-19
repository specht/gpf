#include "GpfDaemonThread.h"
#include "GpfQuery.h"
#include "RefPtr.h"


k_GpfDaemonThread::k_GpfDaemonThread(int ai_SocketDescriptor, QObject* ak_Parent_, k_GpfDaemon& ak_GpfDaemon, k_GpfBase& ak_GpfBase)
	: QThread(ak_Parent_)
	, mi_SocketDescriptor(ai_SocketDescriptor)
	, mk_GpfDaemon(ak_GpfDaemon)
	, mk_GpfBase(ak_GpfBase)
{
}


k_GpfDaemonThread::~k_GpfDaemonThread()
{
}


void k_GpfDaemonThread::run()
{
	QTcpSocket mk_Socket;
	mk_Socket.setSocketDescriptor(mi_SocketDescriptor);

	// read HTTP request from socket
	QString ls_Request, ls_Headers, ls_Method, ls_Uri, ls_Content, ls_ContentType;
	int li_ContentLength = 0;
	int li_ContentReceived = 0;

	bool lb_HeadersComplete = false;
	bool lb_RequestComplete = false;
	bool lb_WriteRequestContentToDisk = false;
	bool lb_RemoveRequestContentFile = true;
	QFile lk_RequestContentFile;
	QTextStream lk_RequestContentStream(&lk_RequestContentFile);
	QString ls_Ticket;

	QString ls_Separator("\r\n\r\n");
	QString ls_ContentLengthKey("Content-Length:");
	QString ls_ContentTypeKey("Content-Type:");

	while (!lb_RequestComplete)
	{
		mk_Socket.waitForReadyRead();
		QString ls_Packet = QString(mk_Socket.readAll());
		if (!lb_HeadersComplete)
		{
			ls_Request += ls_Packet;
			// headers are not yet complete
			if (ls_Request.contains(ls_Separator))
			{
				// headers are complete now, process them and decide how to handle the remaining request data, if any
				lb_HeadersComplete = true;

				int li_SeparatorPosition = ls_Request.indexOf(ls_Separator);
				ls_Headers = ls_Request.left(li_SeparatorPosition);

				QStringList lk_Tokens = QString(ls_Headers).split("\r\n");
				QStringList lk_FirstLine = lk_Tokens[0].split(QChar(' '));
				ls_Method = lk_FirstLine[0];
				ls_Uri = lk_FirstLine[1];

				ls_Content = ls_Request.right(ls_Request.length() - li_SeparatorPosition - ls_Separator.length());
				if (ls_Headers.contains(ls_ContentLengthKey))
				{
					// there is content in the request, determine its length
					QString ls_ContentLength = ls_Headers.right(ls_Headers.length() - ls_Headers.indexOf(ls_ContentLengthKey) - ls_ContentLengthKey.length());
					if (ls_ContentLength.contains("\r\n"))
						ls_ContentLength = ls_ContentLength.left(ls_ContentLength.indexOf("\r\n"));
					li_ContentLength = ls_ContentLength.trimmed().toInt();

					/*
					if (ls_Method == "POST")
					{
						ls_Ticket = mk_GpfDaemon.get_UniqueTicket();
						lb_WriteRequestContentToDisk = true;
						lk_RequestContentFile.setFileName(mk_GpfDaemon.get_TemporaryPath() + ls_Ticket + ".request");
						lk_RequestContentFile.open(QIODevice::ReadWrite);
						lk_RequestContentStream << ls_Uri << "\r\n";

						// act as if the received packet was only made up of content,
						// so that the content is written to the content file further below
					}
					*/
					ls_Packet = ls_Content;
				}
				else
				{
					// there is no content, only headers
					li_ContentLength = 0;
					lb_RequestComplete = true;
				}
			}
		}

		if (lb_HeadersComplete)
		{
			if (li_ContentReceived < li_ContentLength)
			{
				if (lb_WriteRequestContentToDisk)
					lk_RequestContentStream << ls_Packet;
				else
					ls_Content += ls_Packet;

				li_ContentReceived += ls_Packet.length();

				if (li_ContentReceived >= li_ContentLength)
					lb_RequestComplete = true;
			}
			else
				lb_RequestComplete = true;
		}
	}

	if (lb_WriteRequestContentToDisk)
	{
		lk_RequestContentStream.flush();
		lk_RequestContentFile.close();
	}

	// fetch content type
	if (ls_Headers.contains(ls_ContentTypeKey))
	{
		ls_ContentType = ls_Headers.mid(ls_Headers.indexOf(ls_ContentTypeKey) + ls_ContentTypeKey.length());
		if (ls_ContentType.contains("\r\n"))
			ls_ContentType = ls_ContentType.mid(0, ls_ContentType.indexOf("\r\n"));
		ls_ContentType = ls_ContentType.trimmed();
	}

	if (ls_ContentType.startsWith("multipart/form-data"))
	{
		QString ls_BoundaryKey("boundary=");
		QString ls_Boundary = ls_ContentType.mid(ls_ContentType.indexOf(ls_BoundaryKey) + ls_BoundaryKey.length());
	}

	// fetch request variables
	QMap<QString, QString> lk_RequestVars;
	QString ls_Temp = ls_Uri.mid(2);
	QStringList ls_Variables = ls_Temp.split(QChar('&'));
	foreach (QString ls_Variable, ls_Variables)
	{
		QStringList ls_KeyAndValue = ls_Variable.split(QChar('='));
		if (ls_KeyAndValue.size() > 1)
			lk_RequestVars[ls_KeyAndValue[0]] = ls_KeyAndValue[1];
		else
			lk_RequestVars[ls_KeyAndValue[0]] = "";
	}

	QString ls_AdditionalResponseHeaders;
	QString ls_ResponseCode("200 OK");

	printf("-------------\nREQUEST START (from: %s)\n-------------\n%s\n-------------\nREQUEST END\n-------------\n", 
		mk_Socket.peerAddress().toString().toStdString().c_str(), 
		ls_Request.toStdString().c_str());
	// see if the query is specific to a certain genome, fetch that if possible
	RefPtr<k_GpfIndexFileInfo> lk_pIndexFileInfo = mk_GpfBase.get_DefaultIndexFileInfo();
	if (lk_RequestVars.contains("genome"))
	{
		printf("looking for genome: %s.\n", lk_RequestVars["genome"].toStdString().c_str());
		lk_pIndexFileInfo = mk_GpfBase.get_IndexFileInfo(lk_RequestVars["genome"]);
		if (lk_pIndexFileInfo.get_Pointer() == NULL)
			lk_pIndexFileInfo = mk_GpfBase.get_DefaultIndexFileInfo();
	}

	// we now have the headers and the content of the HTTP request, let's do some work!

	/*
	if (ls_Method == "POST")
	{
		if (lk_RequestVars.contains("queryBatch"))
		{
			lb_RemoveRequestContentFile = false;
			// process a query batch
			// check priority ticket

			// enqueue queries
			mk_GpfDaemon.get_QueryQueue().AddBatch(ls_Ticket);
			mk_GpfDaemon.DelegateQuery();

			ms_ResponseContentType = "text/yaml";

			ms_ResponseContent = QString("ticket: %1\n").arg(ls_Ticket);
			ms_ResponseContent += QString("status: http://%1:%2/?queryBatchStatus=%3\n")
				.arg(mk_GpfDaemon.get_HostName())
				.arg(mk_GpfDaemon.get_Port())
				.arg(ls_Ticket);
			ms_ResponseContent += QString("result: ftp://%1%2%3.yaml\n")
				.arg(mk_GpfDaemon.get_HostName())
				.arg(mk_GpfDaemon.get_UriDownloadPath())
				.arg(ls_Ticket);
		}
	}
	*/

	if (lk_RequestVars.contains("getParameters"))
	{
		// returns GPF parameters
		ms_ResponseContentType = "text/yaml";
		ms_ResponseContent = mk_GpfBase.get_GpfParameters();
	}
	/*
	else if (lk_RequestVars.contains("queryBatchStatus"))
	{
		QString ls_Ticket = lk_RequestVars["queryBatchStatus"];
		ms_ResponseContentType = "text/yaml";
		ms_ResponseContent = mk_GpfDaemon.get_QueryQueue().QueryBatchStatus(ls_Ticket);
	}
	else if (lk_RequestVars.contains("fetchBatchResult"))
	{
		QString ls_Ticket = lk_RequestVars["fetchBatchResult"];
		QString ls_Location = QString("http://%1/results/%2.yaml")
			.arg(mk_GpfDaemon.get_HostName()).arg(ls_Ticket);

		// return HTTP header redirect and delegate downloading to a patchy server!
		ls_ResponseCode = "302 FOUND";
		ls_AdditionalResponseHeaders = QString("Location: %1\r\n").arg(ls_Location);
		ms_ResponseContentType = "text/html";
		ms_ResponseContent = QString("Please download the result file at: <a href='%1'>%1</a>.").arg(ls_Location);
	}
	*/
	else if (lk_RequestVars.contains("query"))
	{
		// create a query and enqueue it
		double ld_PrecursorMass = 0.0;
		if (lk_RequestVars.contains("precursorMass"))
		{
			bool lb_Ok;
			ld_PrecursorMass = lk_RequestVars["precursorMass"].toDouble(&lb_Ok);
			if (!lb_Ok)
				ld_PrecursorMass = 0.0;
		}

		mk_pGpfQuery = RefPtr<k_GpfQuery>(new k_GpfQuery(mk_GpfBase, lk_RequestVars["query"], ld_PrecursorMass));
		mk_pGpfQuery->SetParameters(lk_RequestVars);

		// Use a timer to enqueue the query as soon as the event loop has been entered
		// because if we enqueue the query before entering the event loop, we might
		// miss the finished signal.
		QTimer lk_Timer;
		connect(&lk_Timer, SIGNAL(timeout()), this, SLOT(EnqueueQuery()));
		connect(mk_pGpfQuery.get_Pointer(), SIGNAL(finished()), this, SLOT(finishedQuery()));
		lk_Timer.setSingleShot(true);
		lk_Timer.start();

		// enter event loop and wait for query execution to finish
		this->exec();

		// query execution has finished, return results
		ms_ResponseContentType = "text/yaml";

		ms_ResponseContent = "results:\n";
		ms_ResponseContent += "  " + mk_pGpfQuery->get_Result().replace("\n", "\n  ");
		ms_ResponseContent += "\n";
		ms_ResponseContent += mk_pGpfQuery->get_Info();
	}
	/*
	else if (lk_RequestVars.contains("gno"))
	{
		// decode GNO
		ms_ResponseContentType = "text/yaml";
		bool lb_Ok;
		ms_ResponseContent = mk_GpfBase.get_IndexFileInfo()->DecodeGno(lk_RequestVars["gno"].toUInt(&lb_Ok, 16));
	}
	*/
	else if (lk_RequestVars.contains("browse"))
	{
		// browse DNA
		ms_ResponseContentType = "text/yaml";
		bool lb_Error = false;
		bool lb_Ok = false;
		unsigned int lui_Position = lk_RequestVars["position"].toUInt(&lb_Ok);
		if (!lb_Ok)
			lb_Error = true;
		unsigned int lui_Length = lk_RequestVars["length"].toUInt(&lb_Ok);
		if (!lb_Ok)
			lb_Error = true;
		if (lb_Error)
			ms_ResponseContent = "status: error\nmessage: You did not specify a proper position and length.";
		else
			ms_ResponseContent = lk_pIndexFileInfo->Browse(lui_Position, lui_Length);
	}
	else if (lk_RequestVars.contains("readAssembly"))
	{
		// read an assembly
		ms_ResponseContentType = "text/yaml";
		QString ls_Assembly = lk_RequestVars["readAssembly"];

		lk_pIndexFileInfo = mk_GpfBase.get_IndexFileInfoFromAssembly(ls_Assembly);
		ms_ResponseContent = lk_pIndexFileInfo->readAssemblyAsYaml(ls_Assembly);
	}

	QString ls_Response;
	ls_Response += QString("HTTP/1.0 %1\r\n").arg(ls_ResponseCode);
	ls_Response += ls_AdditionalResponseHeaders;
	ls_Response += "Content-Type: " + ms_ResponseContentType + "\r\n";
	ls_Response += QString("Content-Length: %1\r\n\r\n").arg(ms_ResponseContent.length());
	ls_Response += ms_ResponseContent;

	// return HTTP response

	QByteArray lk_ResponseByteArray;
	lk_ResponseByteArray.append(ls_Response);
	mk_Socket.write(lk_ResponseByteArray);
	mk_Socket.disconnectFromHost();
	mk_Socket.waitForDisconnected();

	if (lb_WriteRequestContentToDisk && lb_RemoveRequestContentFile)
		QFile::remove(lk_RequestContentFile.fileName());
}


void k_GpfDaemonThread::EnqueueQuery()
{
	mk_GpfDaemon.get_QueryQueue().Enqueue(mk_pGpfQuery);
	mk_GpfDaemon.DelegateQuery();
}


void k_GpfDaemonThread::finishedQuery()
{
	this->quit();
}
