/*
Copyright (c) 2007-2010 Michael Specht

This file is part of GPF.

GPF is free software: you can redistribute it and/or modify it under 
the terms of the GNU Lesser General Public License as published by the 
Free Software Foundation, either version 3 of the License, or (at your 
option) any later version.

GPF is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GPF.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <QtCore>

class k_GpfIndexFile;

struct r_SearchIntronSplitAlignments
{
    enum Enumeration
    {
        Yes,
        No,
        Conditional
    };
};

struct r_IntronSearchType
{
    enum Enumeration
    {
        Exhaustive,
        Quick
    };
};

struct r_Query
{
    QString ms_Peptide;
    QString ms_Id;
    qint64 mi_PrecursorMass;
    qint64 mi_NTerminalMass;
    qint64 mi_CTerminalMass;

    r_Query(QString as_Peptide, qint64 ai_PrecursorMass)
        : ms_Peptide(as_Peptide)
        , ms_Id(QString())
        , mi_PrecursorMass(ai_PrecursorMass)
        , mi_NTerminalMass(0)
        , mi_CTerminalMass(0)
    {}

    r_Query(QString as_Peptide, qint64 ai_PrecursorMass, qint64 ai_NTerminalMass,
                 qint64 ai_CTerminalMass, QString as_Id)
        : ms_Peptide(as_Peptide)
        , ms_Id(as_Id)
        , mi_PrecursorMass(ai_PrecursorMass)
        , mi_NTerminalMass(ai_NTerminalMass)
        , mi_CTerminalMass(ai_CTerminalMass)
    {}
};

// mass direction: false == left, true == right
typedef QPair<qint64, bool> tk_GnoMassDirection;
typedef QMap<tk_GnoMassDirection, qint64> tk_GnoMap;
typedef QSet<QString> tk_StringSet;
typedef QPair<QString, qint64> tk_StringIntPair;
typedef QPair<int, int> tk_IntPair;
typedef QPair<qint64, qint64> tk_QInt64Pair;
typedef QHash<tk_IntPair, QList<tk_IntPair> > tk_IntPairListHash;
typedef QList<tk_IntPair> tk_IntPairList;


class k_GpfQuery
{
public:
    k_GpfQuery(k_GpfIndexFile& ak_GpfIndexFile, double ad_MassAccuracy,
               bool ab_SimilaritySearch, bool ab_DistinguishIL,
               bool ab_SearchImmediate, r_SearchIntronSplitAlignments::Enumeration 
               ae_SearchIntronSplitAlignments, r_IntronSearchType::Enumeration
               ae_IntronSearchType, int ai_MaxIntronLength, 
               QString as_IntronSpliceSites,
               int ai_PrintFlankingResidues, bool ab_Quiet,
               QIODevice* ak_CsvOutput_, QIODevice* ak_PeptidesOutput_);
    virtual ~k_GpfQuery();
    
    void execute(const QString& as_Peptide, qint64 ai_PrecursorMass, qint64 ai_NTerminalMass, qint64 ai_CTerminalMass, QString as_Id);
    void execute(const QList<r_Query> ak_Peptides);
    
protected:
    void findAlignments(const tk_GnoMap& ak_GnoMap,
                        bool ab_SearchImmediate,
                        bool ab_SearchIntronSplit,
                        tk_StringSet& ak_FoundAssemblies,
                        tk_StringSet& ak_AcceptedAssemblies);
    int reverseSpliceSequence(int ai_Sequence, int ai_Length);
    void readFlankingAminoAcids(qint64 ai_Start, qint64 ai_Stop,
                                QString& as_Left, QString& as_Right,
                                int ai_BStep1, qint64 ai_ScaffoldStart,
                                qint64 ai_ScaffoldEnd,
                                char* ac_TripletToAminoAcid_);
    
    k_GpfIndexFile& mk_GpfIndexFile;
    double md_MassAccuracy;
    bool mb_SimilaritySearch;
    bool mb_DistinguishIL;
    bool mb_SearchImmediateAlignments;
    r_SearchIntronSplitAlignments::Enumeration me_SearchIntronSplitAlignments;
    r_IntronSearchType::Enumeration me_IntronSearchType;
    int mi_MinIntronLength;
    int mi_MaxIntronLength;
    int mi_MinExonLength;
    int mi_FlankingSequenceLength;
    QString ms_IntronSpliceSites;
    bool mb_Quiet;
    QIODevice* mk_CsvOutput_;
    QTextStream mk_CsvOutStream;
    QIODevice* mk_PeptidesOutput_;
    QTextStream mk_PeptidesOutStream;
    QSet<QString> mk_ResultingPeptides;
    
    // GT/2: [AG/2]
    // GC/2: [AG/2]
    tk_IntPairListHash mk_IntronNTerm;
    // the sorted key list, sorted by descending splice site length
    tk_IntPairList mk_IntronNTermKeys;
    
    // AG/2: [GT/2, GC/2]
    tk_IntPairListHash mk_IntronCTerm;
    tk_IntPairList mk_IntronCTermKeys;
    
    // TG/2: [GA/2]
    // CG/2: [GA/2]
    tk_IntPairListHash mk_IntronNTermReverse;
    tk_IntPairList mk_IntronNTermReverseKeys;
    
    // GA/2: [TG/2, CG/2]
    tk_IntPairListHash mk_IntronCTermReverse;
    tk_IntPairList mk_IntronCTermReverseKeys;
    
    int mi_IntronNTermMaxLength;
    int mi_IntronCTermMaxLength;
    
    // per-query class variables
    int mi_AlignmentMinMass;
    int mi_AlignmentMaxMass;
    QString ms_QueryPeptide;
    QString ms_QueryPeptideIL;
    QString ms_QueryPeptideId;
};


int sortByDecreasingLength(const tk_IntPair& ak_First, const tk_IntPair& ak_Second);
