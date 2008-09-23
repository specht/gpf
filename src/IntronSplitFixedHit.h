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

#include "Hit.h"


class k_IntronSplitFixedHit: public k_Hit
{
public:
	k_IntronSplitFixedHit();
	k_IntronSplitFixedHit(k_GpfQuery* ak_Query_, r_Scaffold ar_Scaffold, QString as_Query,
						  QString as_PeptideLeft, QString as_GapFix, QString as_PeptideRight,
						  unsigned int aui_GnoLeft, unsigned int aui_GnoRight,
						  unsigned int aui_NucleotideCountLeft, unsigned int aui_NucleotideCountRight,
						  unsigned int aui_IntronLength,
						  QString as_LeftSurroundings, QString as_RightSurroundings,
						  QString as_IntronEnds);
	virtual ~k_IntronSplitFixedHit();

	virtual QString description();

protected:
	virtual void CalculateScore();

	QString ms_PeptideLeft, ms_GapFix, ms_PeptideRight;
	unsigned int mui_GnoLeft, mui_GnoRight;
	unsigned int mui_NucleotideCountLeft, mui_NucleotideCountRight;
	unsigned int mui_IntronLength;
	QString ms_LeftSurroundings, ms_RightSurroundings;
	QString ms_IntronEnds;
	QString ms_Assembly;
};
