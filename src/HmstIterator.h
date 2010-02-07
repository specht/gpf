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
#include "GpfBase.h"
#include "GpfIndexer.h"
#include "RefPtr.h"

struct r_Hmst
{
	unsigned int mui_TagDirectionIndex;
	qint64 mi_HalfMass;
	quint64 mui_Gno;
};


struct r_HmstIteratorLevel
{
	enum Enumeration
	{
		ScaffoldIndex = 0,
		OrfDirection,
		Frame,
		SpanIndex,
		MassDirection,
		TagOffset,
		Size
	};
};


class k_HmstIterator
{
public:
	k_HmstIterator(k_GpfIndexer& ak_GpfIndexer);
	virtual ~k_HmstIterator();
	
	void reset();
	bool next(r_Hmst* ar_Hmst_);
	
protected:
	bool advance(r_HmstIteratorLevel::Enumeration ae_Level);
	void updateOrfAndCleavageSites();
	void updateTagOffset();
	bool goodState();
	
	k_GpfIndexer& mk_GpfIndexer;
	
	qint64 mk_First_[(int)r_HmstIteratorLevel::Size];
	qint64 mk_Last_[(int)r_HmstIteratorLevel::Size];
	qint64 mk_Value_[(int)r_HmstIteratorLevel::Size];
	
	RefPtr<char> mc_pOrf;
    
    typedef QPair<qint64, qint64> tk_IntPair;
    QList<tk_IntPair> mk_Spans;
    
	bool mb_AtEnd;
	int mi_CurrentTag;
	qint64 mi_CurrentHalfMass;
	quint64 mui_CurrentGno;
	qint64 mi_CurrentAminoAcidSpanLength;
};
