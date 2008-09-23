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

#include <QList>
#include "IntronSplitFixedHit.h"
#include "Hit.h"


class k_IntronSplitHit: public k_Hit
{
public:
	k_IntronSplitHit();
	k_IntronSplitHit(k_GpfQuery* ak_Query_, r_Scaffold ar_Scaffold, QString as_Query,
					 unsigned int aui_PeptideStartGno, unsigned int aui_LeftHalfLength, 
					 unsigned int aui_PeptideEndGno, unsigned int aui_RightHalfLength,
					 unsigned int aui_QueryLengthLeft, unsigned int aui_QueryLengthRight);
	virtual ~k_IntronSplitHit();

	void Update(unsigned int aui_LeftHalfLength, unsigned int aui_RightHalfLength,
				unsigned int aui_QueryLengthLeft, unsigned int aui_QueryLengthRight);

	virtual void Finish();
	const QList<k_IntronSplitFixedHit>& get_FixedHitList() const;

	virtual QString description();

protected:
	void FixGap();
	void ReadUnfixedPeptide();


	// how much of the query has been found (both 3 if only the outermost half sequence 
	// tags were found, higher if more half sequence tags have been found)
	unsigned int mui_QueryLengthLeft;
	unsigned int mui_QueryLengthRight;

	// the first nucleotide of the first amino acid
	unsigned int mui_PeptideStartGno;

	// left part length, in amino acids
	// This does not mean that so many matching amino acids have been found on the genome,
	// this just means that the gap starts at mui_PeptideStartGno + mui_LeftHalfLength * 3
	unsigned int mui_LeftHalfLength;

	// the last nucleotide of the last amino acid
	unsigned int mui_PeptideEndGno;

	// right part length, in amino acids
	unsigned int mui_RightHalfLength;

	// the part of the query that is missing in this result
	QString ms_QueryGap;

	// mass of the query gap
	unsigned int mui_QueryGapMass;

	// left and right results
	QString ms_PeptideLeft, ms_PeptideRight;
	QString ms_LeftSurroundings, ms_RightSurroundings;

	// number of nucleotides between left and right half sequence tag hit
	int mi_UnfixedIntronLength;

	bool mb_ReadUnfixedPeptide;

	QList<k_IntronSplitFixedHit> mk_FixedHitList;
};
