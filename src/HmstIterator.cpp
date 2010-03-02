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


#include "HmstIterator.h"
#include "GpfBase.h"


k_HmstIterator::k_HmstIterator(k_GpfIndexer& ak_GpfIndexer)
    : mk_GpfIndexer(ak_GpfIndexer)
{
    mk_First_[r_HmstIteratorLevel::ScaffoldIndex] = 0;
    mk_Last_[r_HmstIteratorLevel::ScaffoldIndex] = mk_GpfIndexer.mk_ScaffoldLabels.size() - 1;
    mk_First_[r_HmstIteratorLevel::OrfDirection] = 0;
    mk_Last_[r_HmstIteratorLevel::OrfDirection] = 1;
    mk_First_[r_HmstIteratorLevel::Frame] = 0;
    mk_Last_[r_HmstIteratorLevel::Frame] = 2;
    // CleavageSiteIndex is dependent
    mk_First_[r_HmstIteratorLevel::SpanIndex] = 0;
    mk_Last_[r_HmstIteratorLevel::SpanIndex] = -1;
    mk_First_[r_HmstIteratorLevel::MassDirection] = 0;
    mk_Last_[r_HmstIteratorLevel::MassDirection] = 1;
    // last TagOffset is dependent, start is always 0
    mk_First_[r_HmstIteratorLevel::TagOffset] = 0;
    mk_Last_[r_HmstIteratorLevel::TagOffset] = -1;
    
    // determine longest scaffold
    qint64 li_MaxLength = 0;
    foreach (qint64 li_Length, mk_GpfIndexer.mk_ScaffoldLength)
        if (li_Length > li_MaxLength)
            li_MaxLength = li_Length;
        
    mc_pOrf = RefPtr<char>(new char[li_MaxLength / 3 + 1]);
    
    this->reset();
}


k_HmstIterator::~k_HmstIterator()
{
}


void k_HmstIterator::reset()
{
    // SpanIndex is dependent
    mk_First_[r_HmstIteratorLevel::SpanIndex] = 0;
    mk_Last_[r_HmstIteratorLevel::SpanIndex] = -1;

    // TagOffset is dependent
    mk_First_[r_HmstIteratorLevel::TagOffset] = 0;
    mk_Last_[r_HmstIteratorLevel::TagOffset] = -1;
    
    for (int i = 0; i < r_HmstIteratorLevel::Size; ++i)
        mk_Value_[i] = mk_First_[i];

    updateOrfAndCleavageSites();
    updateTagOffset();
    
    mb_AtEnd = false;
}


bool k_HmstIterator::advance(r_HmstIteratorLevel::Enumeration ae_Level)
{
    if (mk_Value_[ae_Level] < mk_Last_[ae_Level])
    {
        ++mk_Value_[ae_Level];
        if ((int)ae_Level < r_HmstIteratorLevel::SpanIndex)
            updateOrfAndCleavageSites();
        if ((int)ae_Level < r_HmstIteratorLevel::TagOffset)
            updateTagOffset();
    }
    else
    {
        if ((int)ae_Level > 0)
        {
            mk_Value_[ae_Level] = mk_First_[ae_Level];
            bool lb_Result = advance((r_HmstIteratorLevel::Enumeration)(ae_Level - 1));
            return lb_Result;
        }
        else
            // we're at the end!
            return false;
    }
    return true;
}


void k_HmstIterator::updateOrfAndCleavageSites()
{
    char* lk_DnaTripletToAminoAcidTable_ = mk_Value_[r_HmstIteratorLevel::OrfDirection] == 0 ? 
        gk_GpfBase.mk_TranslationTables[mk_GpfIndexer.mi_GeneticCode].get_Pointer() :
        gk_GpfBase.mk_TranslationTablesReverse[mk_GpfIndexer.mi_GeneticCode].get_Pointer();
        
    QSet<qint64> lk_CleavageSites;
    lk_CleavageSites << 0;
    qint64 li_ScaffoldStart = mk_GpfIndexer.mk_ScaffoldStart[mk_Value_[r_HmstIteratorLevel::ScaffoldIndex]];
    qint64 li_ScaffoldSize = mk_GpfIndexer.mk_ScaffoldLength[mk_Value_[r_HmstIteratorLevel::ScaffoldIndex]];
    qint64 li_NucleotideStart = li_ScaffoldStart + mk_Value_[r_HmstIteratorLevel::Frame];
    qint64 li_NucleotideEnd = li_ScaffoldStart + li_ScaffoldSize - 3;
    if (mk_Value_[r_HmstIteratorLevel::OrfDirection] == 1)
    {
        li_NucleotideStart = li_ScaffoldStart + li_ScaffoldSize - mk_Value_[r_HmstIteratorLevel::Frame] - 3;
        li_NucleotideEnd = li_ScaffoldStart;
    }
    qint64 li_NucleotideStep = mk_Value_[r_HmstIteratorLevel::OrfDirection] == 0 ? 3 : -3;
    qint64 li_OrfLength = 0;
    for (qint64 i = li_NucleotideStart; 
            mk_Value_[r_HmstIteratorLevel::OrfDirection] == 0 ? 
            (i <= li_NucleotideEnd) :
            (i >= li_NucleotideEnd);
            i += li_NucleotideStep)
    {
        quint16 lui_Triplet = gk_GpfBase.readNucleotideTriplet(mk_GpfIndexer.muc_pDnaBuffer.get_Pointer(), i);
        char lc_AminoAcid = lk_DnaTripletToAminoAcidTable_[lui_Triplet];
        
        mc_pOrf.get_Pointer()[li_OrfLength] = lc_AminoAcid;
        ++li_OrfLength;
        // always cleave after a STOP codon
        if (lc_AminoAcid == '*')
            lk_CleavageSites << li_OrfLength;
        if (mk_GpfIndexer.mb_CleaveAfter_[(int)lc_AminoAcid])
            lk_CleavageSites << li_OrfLength;
        if (mk_GpfIndexer.mb_CleaveBefore_[(int)lc_AminoAcid])
            lk_CleavageSites << (li_OrfLength - 1);
    }
    lk_CleavageSites << li_OrfLength;
    
    QList<qint64> lk_CleavageSitesSorted = lk_CleavageSites.toList();
    qSort(lk_CleavageSitesSorted.begin(), lk_CleavageSitesSorted.end());
    
    mk_Spans.clear();
    
    for (int i = 0; i < lk_CleavageSitesSorted.size() - 1; ++i)
    {
        qint64 li_Start = lk_CleavageSitesSorted[i];
        qint64 li_End = lk_CleavageSitesSorted[i + 1] - 1;
        while (li_Start < li_OrfLength && mc_pOrf.get_Pointer()[li_Start] == '*')
            ++li_Start;
        while (li_End > 0 && mc_pOrf.get_Pointer()[li_End] == '*')
            --li_End;
        mk_Spans << tk_IntPair(li_Start, li_End);
    }
    
    mk_First_[r_HmstIteratorLevel::SpanIndex] = 0;
    mk_Last_[r_HmstIteratorLevel::SpanIndex] = mk_Spans.size() - 1;
    mk_Value_[r_HmstIteratorLevel::SpanIndex] = 0;
}


void k_HmstIterator::updateTagOffset()
{
    /*
    mk_First_[r_HmstIteratorLevel::TagOffset] = mk_CleavageSites[mk_Value_[r_HmstIteratorLevel::CleavageSiteIndex]];
    mk_Last_[r_HmstIteratorLevel::TagOffset] = mk_CleavageSites[mk_Value_[r_HmstIteratorLevel::CleavageSiteIndex] + 1] - mk_GpfIndexer.mi_TagSize;
    */
    mk_First_[r_HmstIteratorLevel::TagOffset] = 0;
    mk_Last_[r_HmstIteratorLevel::TagOffset] = mk_Spans[mk_Value_[r_HmstIteratorLevel::SpanIndex]].second - mk_Spans[mk_Value_[r_HmstIteratorLevel::SpanIndex]].first;
    
    mk_Value_[r_HmstIteratorLevel::TagOffset] = mk_First_[r_HmstIteratorLevel::TagOffset];
    mi_CurrentHalfMass = 0;
    mi_CurrentAminoAcidSpanLength = mk_GpfIndexer.mi_TagSize;
}


bool k_HmstIterator::goodState()
{
    for (int i = 0; i < r_HmstIteratorLevel::Size; ++i)
    {
        if ((mk_Value_[i] < mk_First_[i]) || (mk_Value_[i] > mk_Last_[i]))
            return false;
    }
    //  li_SpanEnd - (k - li_SpanStart) - mi_TagSize + 1
    qint32 li_SpanOffset;
    if (mk_Value_[r_HmstIteratorLevel::MassDirection] == 0)
        li_SpanOffset = (int)(mk_Spans[mk_Value_[r_HmstIteratorLevel::SpanIndex]].first + mk_Value_[r_HmstIteratorLevel::TagOffset]);
    else
        li_SpanOffset = (int)(mk_Spans[mk_Value_[r_HmstIteratorLevel::SpanIndex]].second - mk_GpfIndexer.mi_TagSize - (mk_Value_[r_HmstIteratorLevel::TagOffset] - 1));
    mi_CurrentTag = gk_GpfBase.aminoAcidPolymerCode(mc_pOrf.get_Pointer() + li_SpanOffset, mk_GpfIndexer.mi_TagSize);
    
    if (mi_CurrentTag < 0)
    {
        mi_CurrentHalfMass = mk_GpfIndexer.mi_MaxMass + 1;
        return false;
    }
    qint64 li_DeltaMass = 0;
    if (mk_Value_[r_HmstIteratorLevel::TagOffset] != mk_First_[r_HmstIteratorLevel::TagOffset])
    {
        int li_CurrentAminoAcidIndex = li_SpanOffset + ((mk_Value_[r_HmstIteratorLevel::MassDirection] == 0) ? -1 : mk_GpfIndexer.mi_TagSize);
        char lc_CurrentAminoAcid = mc_pOrf.get_Pointer()[li_CurrentAminoAcidIndex];
        li_DeltaMass = mk_GpfIndexer.mi_AminoAcidMasses_[(int)lc_CurrentAminoAcid];
    }
    // increase half mass if not greater than maximum allowed mass
    if (mi_CurrentHalfMass <= mk_GpfIndexer.mi_MaxMass)
        mi_CurrentHalfMass += li_DeltaMass;
    
    if (mi_CurrentHalfMass > mk_GpfIndexer.mi_MaxMass)
        return false;
    
    qint64 li_SpanAnchor = mk_Value_[r_HmstIteratorLevel::MassDirection] == 0 ? 
        mk_Spans[(int)mk_Value_[r_HmstIteratorLevel::SpanIndex]].first :
        mk_Spans[(int)mk_Value_[r_HmstIteratorLevel::SpanIndex]].second + 1;
    
    if (mk_Value_[r_HmstIteratorLevel::OrfDirection] == 0)
    {
        mui_CurrentGno = 
            mk_GpfIndexer.mk_ScaffoldStart[(int)mk_Value_[r_HmstIteratorLevel::ScaffoldIndex]] + 
            li_SpanAnchor * 3 +
            mk_Value_[r_HmstIteratorLevel::Frame] - mk_Value_[r_HmstIteratorLevel::MassDirection];
    }
    else
    {
        mui_CurrentGno = 
            mk_GpfIndexer.mk_ScaffoldStart[(int)mk_Value_[r_HmstIteratorLevel::ScaffoldIndex]] +
            mk_GpfIndexer.mk_ScaffoldLength[(int)mk_Value_[r_HmstIteratorLevel::ScaffoldIndex]] -
            li_SpanAnchor * 3 -
            mk_Value_[r_HmstIteratorLevel::Frame] - 1 + mk_Value_[r_HmstIteratorLevel::MassDirection];
        mui_CurrentGno |= mk_GpfIndexer.mui_GnoBackwardsBit;
    }
    return mi_CurrentTag >= 0;
}


bool k_HmstIterator::next(r_Hmst* ar_Hmst_)
{
    if (mb_AtEnd)
        return false;
    while (!goodState())
    {
        bool lb_Result = advance((r_HmstIteratorLevel::Enumeration)(r_HmstIteratorLevel::Size - 1));
        if (!lb_Result)
        {
            mb_AtEnd = true;
            return false;
        }
    }
    // now we are in a good state, return the current HMST
    ar_Hmst_->mui_TagDirectionIndex = mi_CurrentTag * 2 + mk_Value_[r_HmstIteratorLevel::MassDirection];
    ar_Hmst_->mi_HalfMass = mi_CurrentHalfMass;
    ar_Hmst_->mui_Gno = mui_CurrentGno;
    // advance at the end
    bool lb_Result = advance((r_HmstIteratorLevel::Enumeration)(r_HmstIteratorLevel::Size - 1));
    if (!lb_Result)
        mb_AtEnd = true;
    return true;
}
