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
#include "GpfQuery.h"
#include "RefPtr.h"


struct r_ChainMarker
{
	enum Enumeration 
	{
		Left = 1,
		Right = 2
	};
};

typedef QHash<QString, QString> tk_StringHash;


class k_Hit
{
public:
	k_Hit(k_GpfQuery& ak_Query, bool ab_Forward, QList<QPair<unsigned int, unsigned int> > ak_Assembly);
	k_Hit(const k_Hit& ak_Other);
	virtual ~k_Hit();

	virtual QString get_Peptide() const;
	virtual QString get_Assembly() const;
	virtual unsigned int get_Mass() const;
	virtual void AddInformation(QString as_Key, QString as_Value);

	virtual void Finish();

	bool get_HasFullScore() const;
	bool get_IsDiscarded() const;
	void MarkDiscarded();

	virtual QString description();
    virtual tk_StringHash descriptionHash();

protected:
	virtual void CalculateScore();

	k_GpfQuery& mk_Query;
	bool mb_Forward;
	QList<QPair<unsigned int, unsigned int> > mk_Assembly;
	unsigned mui_QueryLength;

	QList<QPair<QString, QString> > mk_ResultItems;

	QString ms_Peptide;
	QString ms_Left;
	QString ms_Right;
	QString ms_Assembly;

	// derived values that are calculated on Finish().
	int mi_Score;
	QList<unsigned int> mk_PartScores;
	unsigned int mui_HitMass;
	bool mb_IsDiscarded;

private:
	void Initialize();
};
