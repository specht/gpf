#pragma once

#include <QtCore>
#include "RefPtr.h"
#include "StopWatch.h"


class k_GpfDaemon;
class k_GpfQuery;


struct r_QueryBatchState
{
	enum Enumeration
	{
		Waiting,
		Executing,
		Finished
	};
};


class k_QueryBatch
{
public:
	k_QueryBatch(k_GpfDaemon& ak_GpfDaemon, QString as_Ticket);
	virtual ~k_QueryBatch();

	RefPtr<k_GpfQuery> get_NextQuery();
	void AddResult(QString as_Result);

	QString get_Ticket() const;
	int get_ExpectedQueryCount() const;
	int get_ProcessedQueryCount() const;
	int get_RemainingQueryCount() const;

	r_QueryBatchState::Enumeration get_State() const;
	void set_State(r_QueryBatchState::Enumeration ae_State);


private:
	k_GpfDaemon& mk_GpfDaemon;
	QString ms_Ticket;
	QFile mk_RequestFile;
	QString ms_ResultFilename;
	QFile mk_ResultFile;
	QTextStream mk_RequestFileStream;
	int mi_ExpectedQueryCount;
	int mi_ProcessedQueryCount;
	QMap<QString, QString> mk_RequestVars;
	k_StopWatch mk_StopWatch;
	r_QueryBatchState::Enumeration me_State;
};
