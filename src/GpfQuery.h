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
#include "RefPtr.h"


class k_GpfIndexFile;


class k_GpfQuery
{
public:
	k_GpfQuery(k_GpfIndexFile& ak_GpfIndexFile, QIODevice* ak_Output_ = NULL);
	virtual ~k_GpfQuery();
	
	void execute(const QString& as_Peptide);
	
protected:
    int reverseSpliceSequence(int ai_Sequence, int ai_Length);
    
	k_GpfIndexFile& mk_GpfIndexFile;
    QIODevice* mk_Output_;
	qint64 mi_Mass, mi_MinMass, mi_MaxMass;
	double md_MassAccuracy;
    // min intron length does not include slice donor/acceptor CS!
	int mi_MinIntronLength;
	int mi_MaxIntronLength;
	int mi_MaxNucleotideSpanLength;
	QString ms_IntronSpliceSites;
    bool mb_SimilaritySearch;
    bool mb_ImmediateHitsSufficient;
	
	typedef QPair<int, int> tk_IntPair;
    typedef QHash<tk_IntPair, QList<tk_IntPair> > tk_IntPairListHash;

	// GT/2: [AG/2]
	// GC/2: [AG/2]
	tk_IntPairListHash mk_IntronNTerm;
	
	// AG/2: [GT/2, GC/2]
	tk_IntPairListHash mk_IntronCTerm;
	
	int mi_IntronNTermMaxLength;
	int mi_IntronCTermMaxLength;
    RefPtr<QFile> mk_pStdOutFile;
};
