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

#include "QueryBatch.h"
#include "GpfDaemon.h"
#include "GpfQuery.h"


k_QueryBatch::k_QueryBatch(k_GpfDaemon& ak_GpfDaemon, QString as_Ticket)
	: mk_GpfDaemon(ak_GpfDaemon)
	, ms_Ticket(as_Ticket)
	, mk_RequestFile(ak_GpfDaemon.get_TemporaryPath() + as_Ticket + ".request")
	, mk_RequestFileStream(&mk_RequestFile)
	, mi_ExpectedQueryCount(0)
	, mi_ProcessedQueryCount(0)
	, mk_StopWatch()
	, me_State(r_QueryBatchState::Waiting)
{
	mk_RequestFile.open(QIODevice::ReadOnly);

	// construct GPF options map from URI
	QString ls_Uri = mk_RequestFileStream.readLine();
	QString ls_Temp = ls_Uri.mid(2);
	QStringList ls_Variables = ls_Temp.split(QChar('&'));
	foreach (QString ls_Variable, ls_Variables)
	{
		QStringList ls_KeyAndValue = ls_Variable.split(QChar('='));
		if (ls_KeyAndValue.size() > 1)
			mk_RequestVars[ls_KeyAndValue[0]] = ls_KeyAndValue[1];
		else
			mk_RequestVars[ls_KeyAndValue[0]] = "";
	}

	// count lines == number of queries
	while (!mk_RequestFileStream.atEnd())
	{
		mk_RequestFileStream.readLine();
		++mi_ExpectedQueryCount;
	}

	// rewind request file and skip URI line
	mk_RequestFileStream.seek(0);
	mk_RequestFileStream.readLine();

	ms_ResultFilename = ms_Ticket + mk_GpfDaemon.get_UniqueTicket() + ".yaml";
	mk_ResultFile.setFileName(mk_GpfDaemon.get_DownloadPath() + ms_ResultFilename);
	mk_ResultFile.open(QIODevice::WriteOnly);

	// write query batch info, including parameters
	RefPtr<k_GpfQuery> lk_pQuery(new k_GpfQuery(mk_GpfDaemon.get_GpfBase(), "WLQYSEVIHAR", 0.0));
	lk_pQuery->SetParameters(mk_RequestVars);
	QString ls_Parameters = lk_pQuery->get_Parameters() + "\nresults:\n\n";
	mk_ResultFile.write(ls_Parameters.toAscii(), ls_Parameters.length());
}


k_QueryBatch::~k_QueryBatch()
{
	mk_RequestFile.close();
	mk_ResultFile.close();

	// remove request file
	QFile::remove(mk_GpfDaemon.get_TemporaryPath() + ms_Ticket + ".request");

	// rename result file
	mk_ResultFile.rename(mk_GpfDaemon.get_DownloadPath() + ms_Ticket + ".yaml");
}


RefPtr<k_GpfQuery> k_QueryBatch::get_NextQuery()
{
	RefPtr<k_GpfQuery> lk_pQuery(NULL);

	if (mk_RequestFileStream.atEnd())
	{
		mi_ProcessedQueryCount = mi_ExpectedQueryCount;
		return RefPtr<k_GpfQuery>(NULL);
	}

	// read next line from request file and construct query object
	QString ls_Line = mk_RequestFileStream.readLine();

	++mi_ProcessedQueryCount;
	if (me_State == r_QueryBatchState::Waiting)
		me_State = r_QueryBatchState::Executing;

	QStringList ls_Variables = ls_Line.split(QChar(';'));
	QMap<QString, QString> lk_Vars;

	foreach (QString ls_Variable, ls_Variables)
	{
		QStringList ls_KeyAndValue = ls_Variable.split(QChar('='));
		if (ls_KeyAndValue.size() > 1)
			lk_Vars[ls_KeyAndValue[0]] = ls_KeyAndValue[1];
		else
			lk_Vars[ls_KeyAndValue[0]] = "";
	}

	if (lk_Vars.contains("peptide"))
	{
		double ld_PrecursorMass = 0.0;
		QString ls_Label;
		if (lk_Vars.contains("precursorMass"))
		{
			bool lb_Ok;
			ld_PrecursorMass = lk_Vars["precursorMass"].toDouble(&lb_Ok);
			if (!lb_Ok)
				ld_PrecursorMass = 0.0;
		}

		lk_pQuery = RefPtr<k_GpfQuery>(new k_GpfQuery(mk_GpfDaemon.get_GpfBase(), lk_Vars["peptide"], ld_PrecursorMass, ls_Line));
		lk_pQuery->SetParameters(mk_RequestVars);
	}
	
	return lk_pQuery;
}


void k_QueryBatch::AddResult(QString as_Result)
{
	QByteArray lk_Result = as_Result.toAscii();
	mk_ResultFile.write(lk_Result, lk_Result.length());
	if (me_State == r_QueryBatchState::Executing && mi_ProcessedQueryCount == mi_ExpectedQueryCount)
		me_State = r_QueryBatchState::Finished;
}


QString k_QueryBatch::get_Ticket() const
{
	return ms_Ticket;
}


int k_QueryBatch::get_ExpectedQueryCount() const
{
	return mi_ExpectedQueryCount;
}


int k_QueryBatch::get_ProcessedQueryCount() const
{
	return mi_ProcessedQueryCount;
}


int k_QueryBatch::get_RemainingQueryCount() const
{
	return mi_ExpectedQueryCount - mi_ProcessedQueryCount;
}


r_QueryBatchState::Enumeration k_QueryBatch::get_State() const
{
	return me_State;
}


void k_QueryBatch::set_State(r_QueryBatchState::Enumeration ae_State)
{
	me_State = ae_State;
}
