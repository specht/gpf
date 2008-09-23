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

#include <QString>
#include <QReadWriteLock>
#include "GpfBase.h"
#include "GpfIndexFileInfo.h"


class k_Hit;

typedef QHash<QString, RefPtr<k_Hit> > tk_ResultList;
typedef QPair<unsigned int, unsigned int> tk_AssemblyPart;
typedef QList<tk_AssemblyPart> tk_Assembly;
typedef QPair<unsigned int, unsigned int> tk_PeptideSpan;
typedef QSet<tk_PeptideSpan> tk_PeptideLocations;


class k_GpfQuery: public QObject
{
	Q_OBJECT

public:
	k_GpfQuery(k_GpfBase& ak_GpfBase, QString as_Query, double ad_PrecursorMass = 0.0, QString as_Label = QString());
	virtual ~k_GpfQuery();

	void Execute();
	QString get_Query() const;
	QString get_CollapsedQuery() const;
	double get_PrecursorMass() const;
	QString get_Label() const;
	QString get_Result() const;
	QString get_Parameters();
	QString get_Info();
	bool get_IsFinished() const;

	k_GpfBase& get_GpfBase();
	QFile& get_IndexFile();
	k_GpfIndexFileInfo* get_GpfIndexFileInfo();

	void SetParameters(QMap<QString, QString> ak_Parameters);
	unsigned int mui_QueryMass;
	int mi_MassTolerance;
	unsigned int mui_ResultMassMinimum;
	unsigned int mui_ResultMassMaximum;

	QMap<r_GpfParameterName::Enumeration, k_GpfParameter*> mk_Parameters;

signals:
	void finished();

private:
	tk_ResultList QueryImmediateHit();
	//unsigned int QueryIntronHit();
	tk_ResultList QueryCsIntronHit();

	tk_PeptideLocations FindPentamers(bool ab_Left);
	tk_ResultList AssembleHalfHits(tk_PeptideLocations ak_Locations, bool ab_Left);

	void CropList(qint64 ai_FileOffset, unsigned int aui_Start, unsigned int aui_Count, int ai_MinMass, int ai_MaxMass, unsigned int& aui_CropStart, unsigned int& aui_CropCount);

	k_GpfBase& mk_GpfBase;
	RefPtr<k_GpfIndexFileInfo> mk_pIndexFileInfo;
	QString ms_Query;
	QString ms_CollapsedQuery;
	QString ms_Result;
	bool mb_IsFinished;
	mutable QReadWriteLock mk_IsFinishedLock;
	QFile mk_IndexFile;
	QString ms_QueryProcessingTime;
	double md_PrecursorMass;
	QString ms_Label;
};
