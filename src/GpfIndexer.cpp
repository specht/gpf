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

#include "math.h"
#include "GpfIndexer.h"
#include "BitWriter.h"
#include "GpfBase.h"
#include "HmstIterator.h"
#include "StopWatch.h"


// :UGLY: This is really ugly indeed. But it has been done in order to be able
// to use Qt's qsort() function with a simple less-than function pointer which
// can not be a class member function pointer. Therefore, we just store the current
// indexer in this global variable in order to do the comparison in a global
// less-than function. Otherwise, we would have to implement our own sorting function
// or, as an alternative, implement a new class containing the int and the pointer to
// the appropriate indexer, which means much more space. Because probably, only one
// genome file is indexed at a time, we can probably do this trick.
// Actually, it's not a pointer to the current indexer, but to the current list 
// containing all masses from the current tag/direction set while sorting out there
// final index.
quint8* guc_MassesBuffer_;
qint32 gi_MassBits;


k_GpfIndexer::k_GpfIndexer(QString as_DnaPath, QString as_DnaIndexPath, 
                           QString as_Title, qint32 ai_TagSize, 
                           QString as_Enzyme, qint64 ai_IndexBufferAllocSize,
                           qint32 ai_MassPrecision, qint32 ai_MassBits,
                           qint32 ai_GeneticCode)
    : ms_DnaPath(as_DnaPath)
    , ms_DnaIndexPath(as_DnaIndexPath)
    , ms_Title(as_Title)
    , ms_Enzyme(as_Enzyme)
    , mi_TotalNucleotideCount(0)
    , mi_OffsetBits(0)
    , mui_GnoBackwardsBit(0)
    , mi_MassBits(ai_MassBits)
    , mi_MassPrecision(ai_MassPrecision)
    , mi_TagSize(ai_TagSize)
    , mi_GeneticCode(ai_GeneticCode)
    , mi_DnaBufferLength(0)
    , mi_IndexBufferMaxLength(ai_IndexBufferAllocSize)
{
    // parse enzyme
    QStringList lk_EnzymeList = as_Enzyme.split("|");
    if (lk_EnzymeList.size() < 2)
    {
        printf("Error: Pipe character (|) missing in specified enzyme: %s.\n",
               as_Enzyme.toStdString().c_str());
        exit(1);
    }
    if (lk_EnzymeList.size() > 2)
    {
        printf("Error: More than one pipe character (|) in specified enzyme: %s.\n",
               as_Enzyme.toStdString().c_str());
        exit(1);
    }
    for (int i = 0; i < 256; ++i)
    {
        mb_CleaveBefore_[i] = false;
        mb_CleaveAfter_[i] = false;
    }
    for (int i = 0; i < lk_EnzymeList[0].length(); ++i)
    {
        char lc_AminoAcid = (lk_EnzymeList[0].at(i).toAscii()) & 0xff;
        mb_CleaveAfter_[(int)lc_AminoAcid] = true;
    }
    for (int i = 0; i < lk_EnzymeList[1].length(); ++i)
    {
        char lc_AminoAcid = (lk_EnzymeList[1].at(i).toAscii()) & 0xff;
        mb_CleaveBefore_[(int)lc_AminoAcid] = true;
    }
}


k_GpfIndexer::~k_GpfIndexer()
{
}


void k_GpfIndexer::compileIndex()
{
    k_StopWatch lk_StopWatch("Compiling the DNA index took %1.\n");
    
    printf("Genome title is '%s', indexing with a tag size of %d.\n", 
           ms_Title.toStdString().c_str(),
           mi_TagSize);
    printf("Using %s for triplet translation.\n", gk_GpfBase.mk_TranslationTableTitle[mi_GeneticCode].toStdString().c_str());
    printf("Using a mass precision of %d, this corresponds to %1.2f decimal places.\n",
           mi_MassPrecision, log((double)mi_MassPrecision) / log(10.0));
    printf("Using %d bits for mass entries, the highest mass possible is %1.2f kDa.\n",
           mi_MassBits, (double)(((qint64)1 << mi_MassBits) - 1) / mi_MassPrecision / 1000.0);
    
    QFile lk_OutFile(ms_DnaIndexPath);
    lk_OutFile.open(QIODevice::WriteOnly);
    
    // parse FASTA DNA file
    parseDna(ms_DnaPath);
    
    mi_OffsetBits = 1;
    while ((((qint64)1) << mi_OffsetBits) < mi_TotalNucleotideCount)
        ++mi_OffsetBits;
    
    // add one extra bit for recording the reading direction
    ++mi_OffsetBits;
    mui_GnoBackwardsBit = (quint64)1 << (mi_OffsetBits - 1);
    
    // determine tag bits
    mi_TagCount =  1;
    for (int i = 0; i < mi_TagSize; ++i)
        mi_TagCount *= 19;
    
    mi_HmstBits = mi_OffsetBits + mi_MassBits;
    mi_MaxMass = ((qint64)1 << mi_MassBits) - 1;
    
    printf("Allocating %s for tag/direction count (32 bits) list.\n", gk_GpfBase.bytesToStr(mi_TagCount * 2 * 4).toStdString().c_str());
    mui_pTagDirectionCount = QSharedPointer<quint32>(new quint32[mi_TagCount * 2]);
    memset(mui_pTagDirectionCount.data(), 0, mi_TagCount * 2 * 4);

    // file size * 2 (all six reading frames) * 2 (left and right HMST)
    qint64 li_IndexBufferMaxRequiredSize = (QFileInfo(ms_DnaPath).size() * 2 * 2 * mi_HmstBits / 8) + 1;
    if (li_IndexBufferMaxRequiredSize < mi_IndexBufferMaxLength)
        mi_IndexBufferMaxLength = li_IndexBufferMaxRequiredSize;
    printf("Allocating %s for index buffer (this appears to be a %d bit system).\n", 
           gk_GpfBase.bytesToStr(mi_IndexBufferMaxLength).toStdString().c_str(),
           (int)(sizeof(size_t) * 8));
    // allocate 9 extra bytes so that we're always safe if we should read several bytes at once
    muc_pIndexBuffer = QSharedPointer<quint8>(new quint8[mi_IndexBufferMaxLength + 9]);
    memset(muc_pIndexBuffer.data(), 0, mi_IndexBufferMaxLength);
    mi_IndexBufferOffset = 0;
    mi_IndexBufferBitOffset = 0;
    
    writeIdentifierChunk(&lk_OutFile);
    writeInfoChunk(&lk_OutFile);
    
    // only write genetic code chunk if not standard code
    if (mi_GeneticCode != 1)
        writeGeneticCodeChunk(&lk_OutFile);
    
    // only write enzyme chunk if not 'RK|'
    if (ms_Enzyme != "RK|")
        writeEnzymeChunk(&lk_OutFile);
    
    writeDnaChunk(&lk_OutFile);
    
    writeIndexChunk(&lk_OutFile);
    
    lk_OutFile.close();
}


void k_GpfIndexer::parseDna(QString as_DnaPath)
{
    //printf("Parsing input file...");
    mk_ScaffoldLabels.clear();
    mk_ScaffoldLength.clear();
    
    QFile lk_DnaFile(as_DnaPath);
    lk_DnaFile.open(QIODevice::ReadOnly);
    QTextStream lk_Stream(&lk_DnaFile);
    
    printf("Allocating %s for DNA.\n", gk_GpfBase.bytesToStr(lk_DnaFile.size() * 3 / 8 + 1).toStdString().c_str());
    muc_pDnaBuffer = QSharedPointer<quint8>(new quint8[lk_DnaFile.size() * 3 / 8 + 1]);
    qint64 li_DnaBufferOffset = 0;
    mi_DnaBufferLength = 0;
    
    quint16 lui_DnaBuffer = 0;
    // :CAREFUL: li_DnaBufferLength might be confused with mi_DnaBufferLength
    // li_DnaBufferLength is about our 16 bit mini buffer, mi_DnaBufferLength
    // is the whole DNA!
    int li_DnaBufferLength = 0;
    QString ls_CurrentLabel = "";
    mi_TotalNucleotideCount = 0;
    
    int li_Percent = -1;
    
    while (!lk_Stream.atEnd())
    {
        int li_NowPercent = lk_DnaFile.pos() * 100 / lk_DnaFile.size();
        if (li_NowPercent != li_Percent)
        {
            printf("\rReading DNA... %d%%", li_NowPercent);
            li_Percent = li_NowPercent;
        }
        QString ls_Line = lk_Stream.readLine();
        ls_Line = ls_Line.trimmed();
        if (ls_Line[0] == QChar('>'))
        {
            ls_Line.remove(0, 1);
            QString ls_Label = ls_Line.trimmed();
            fflush(stdout);

            ls_CurrentLabel = ls_Label;
            mk_ScaffoldLabels.push_back(ls_Label);
            mk_ScaffoldLength.push_back(0);
            mk_ScaffoldStart.push_back(mi_TotalNucleotideCount);
        }
        else
        {
            mk_ScaffoldLength.last() += ls_Line.length();
            mi_TotalNucleotideCount += ls_Line.length();
            for (int i = 0; i < ls_Line.length(); ++i)
            {
                lui_DnaBuffer |= (gk_GpfBase.mk_DnaCharToNumber_[(unsigned char)(ls_Line.at(i).toAscii())] << li_DnaBufferLength);
                li_DnaBufferLength += 3;
                if (li_DnaBufferLength >= 8)
                {
                    quint8 lui_DnaBufferBit = (quint8)(lui_DnaBuffer & 255);
                    muc_pDnaBuffer.data()[li_DnaBufferOffset] = lui_DnaBufferBit;
                    ++li_DnaBufferOffset;
                    ++mi_DnaBufferLength;
                    li_DnaBufferLength -= 8;
                    lui_DnaBuffer >>= 8;
                }
            }
        }
    }
    // write remaining nucleotides
    if (li_DnaBufferLength > 0)
    {
        quint8 lui_DnaBufferBit = (quint8)(lui_DnaBuffer & 255);
        muc_pDnaBuffer.data()[li_DnaBufferOffset] = lui_DnaBufferBit;
        ++li_DnaBufferOffset;
        ++mi_DnaBufferLength;
    }
    printf(" done.\n");
}


qint64 k_GpfIndexer::writeChunkHeader(QFile* ak_OutFile_, r_DnaIndexChunkType::Enumeration ae_Type)
{
    quint32 lui_Type = (quint32)ae_Type;
    ak_OutFile_->write((char*)&lui_Type, 4);
    qint64 li_Pos = ak_OutFile_->pos();
    qint64 li_Size = 0;
    ak_OutFile_->write((char*)&li_Size, 8);
    return li_Pos;
}


void k_GpfIndexer::writeChunkSize(QFile* ak_OutFile_, qint64 ai_Position, qint64 ai_Size)
{
    qint64 li_OldPos = ak_OutFile_->pos();
    ak_OutFile_->seek(ai_Position);
    ak_OutFile_->write((char*)&ai_Size, 8);
    ak_OutFile_->seek(li_OldPos);
}


void k_GpfIndexer::writeIdentifierChunk(QFile* ak_OutFile_)
{
    ak_OutFile_->write(QByteArray("gpfindex"));
    qint16 li_VersionMajor = 3;
    qint16 li_VersionMinor = 0;
    ak_OutFile_->write((char*)&li_VersionMajor, 2);
    ak_OutFile_->write((char*)&li_VersionMinor, 2);
}


void k_GpfIndexer::writeInfoChunk(QFile* ak_OutFile_)
{
    qint64 li_SizeLocation = writeChunkHeader(ak_OutFile_, r_DnaIndexChunkType::Info);
    qint64 li_ChunkSize = ak_OutFile_->pos();
    
    // write genome title
    qint32 li_TitleLength = ms_Title.length();
    ak_OutFile_->write((char*)&li_TitleLength, 4);
    ak_OutFile_->write((char*)ms_Title.toStdString().c_str(), li_TitleLength);
    
    // write constants
    ak_OutFile_->write((char*)&mi_OffsetBits, 4);
    ak_OutFile_->write((char*)&mi_MassBits, 4);
    ak_OutFile_->write((char*)&mi_MassPrecision, 4);
    ak_OutFile_->write((char*)&mi_TagSize, 4);
    
    // write scaffold count, lengths and labels
    qint32 li_ScaffoldCount = mk_ScaffoldLabels.size();
    ak_OutFile_->write((char*)&li_ScaffoldCount, 4);
    for (int i = 0; i < li_ScaffoldCount; ++i)
    {
        QString ls_Label = mk_ScaffoldLabels[i];
        qint64 li_ScaffoldSize = mk_ScaffoldLength[i];
        ak_OutFile_->write((char*)&li_ScaffoldSize, 8);
        qint32 li_ScaffoldLabelLength = ls_Label.length();
        ak_OutFile_->write((char*)&li_ScaffoldLabelLength, 4);
        ak_OutFile_->write((char*)ls_Label.toStdString().c_str(), li_ScaffoldLabelLength);
    }
    
    // finish info chunk
    li_ChunkSize = ak_OutFile_->pos() - li_ChunkSize;
    writeChunkSize(ak_OutFile_, li_SizeLocation, li_ChunkSize);
}


void k_GpfIndexer::writeGeneticCodeChunk(QFile* ak_OutFile_)
{
    qint64 li_SizeLocation = writeChunkHeader(ak_OutFile_, r_DnaIndexChunkType::GeneticCode);
    qint64 li_ChunkSize = ak_OutFile_->pos();
    
    // write genetic code
    ak_OutFile_->write((char*)&mi_GeneticCode, 4);
    
    // finish genetic code chunk
    li_ChunkSize = ak_OutFile_->pos() - li_ChunkSize;
    writeChunkSize(ak_OutFile_, li_SizeLocation, li_ChunkSize);
}


void k_GpfIndexer::writeEnzymeChunk(QFile* ak_OutFile_)
{
    qint64 li_SizeLocation = writeChunkHeader(ak_OutFile_, r_DnaIndexChunkType::Enzyme);
    qint64 li_ChunkSize = ak_OutFile_->pos();
    
    // write enzyme
    qint32 li_EnzymeLength = ms_Enzyme.length();
    ak_OutFile_->write((char*)&li_EnzymeLength, 4);
    ak_OutFile_->write((char*)ms_Enzyme.toStdString().c_str(), li_EnzymeLength);
    
    // finish enzyme chunk
    li_ChunkSize = ak_OutFile_->pos() - li_ChunkSize;
    writeChunkSize(ak_OutFile_, li_SizeLocation, li_ChunkSize);
}


void k_GpfIndexer::writeDnaChunk(QFile* ak_OutFile_)
{
    qint64 li_SizeLocation = writeChunkHeader(ak_OutFile_, r_DnaIndexChunkType::Dna);
    ak_OutFile_->write((char*)muc_pDnaBuffer.data(), mi_DnaBufferLength);
    writeChunkSize(ak_OutFile_, li_SizeLocation, mi_DnaBufferLength);
}


void k_GpfIndexer::writeIndexChunk(QFile* ak_OutFile_)
{
    qint64 li_SizeLocation = writeChunkHeader(ak_OutFile_, r_DnaIndexChunkType::Index);
    qint64 li_ChunkSize = ak_OutFile_->pos();
    
// calculate mass according to mass accuracy
    for (int i = 0; i < 256; ++i)
    {
        mi_AminoAcidMasses_[i] = 0;
        if (gk_GpfBase.mb_IsAminoAcid_[i])
            mi_AminoAcidMasses_[i] = (qint64)(gk_GpfBase.md_AminoAcidMasses_[i] * mi_MassPrecision);
    }
    
    mi_IndexBufferOffset = 0;
    mi_IndexBufferBitOffset = 0;

    qint64 li_MaxLen = 0;
    foreach (qint64 li_Length, mk_ScaffoldLength)
        if (li_Length > li_MaxLen)
            li_MaxLen = li_Length;
        
    //printf("Translating DNA, marking cleavage sites...");
    
    
    printf("HMST size: %d bits (%d for mass, %d for GNO)\n", (int)(mi_HmstBits), (int)mi_MassBits, (int)mi_OffsetBits);
    k_HmstIterator lk_HmstIterator(*this);
    r_Hmst lr_Hmst;
    qint64 li_TotalHmstCount = 0;
    QSharedPointer<qint64> li_pTagDirectionCount(new qint64[mi_TagCount * 2]);
    memset(li_pTagDirectionCount.data(), 0, sizeof(qint64) * mi_TagCount * 2);
    qint64 li_BiggestBucketSize = 0;
    {
        k_StopWatch lk_SubStopWatch("One iteration took %1.\n");
        printf("Doing first iteration, counting HMST occurences...");
        fflush(stdout);
        while (lk_HmstIterator.next(&lr_Hmst))
        {
            //printf("%d %d %d\n", (unsigned int)lr_Hmst.mui_TagDirectionIndex, (int)lr_Hmst.mi_HalfMass, (unsigned int)lr_Hmst.mui_Gno);
            // count
            ++li_TotalHmstCount;
            ++li_pTagDirectionCount.data()[lr_Hmst.mui_TagDirectionIndex];
            if (li_pTagDirectionCount.data()[lr_Hmst.mui_TagDirectionIndex] > li_BiggestBucketSize)
                li_BiggestBucketSize = li_pTagDirectionCount.data()[lr_Hmst.mui_TagDirectionIndex];
        } 
        printf("done.\n");
    }
    
    qint64 li_MaxHmstPerIteration = ((qint64)mi_IndexBufferMaxLength * 8) / mi_HmstBits;
    
    QSharedPointer<quint32> lui_pIndicesToSort(new quint32[li_BiggestBucketSize]);
    
    if (li_BiggestBucketSize > li_MaxHmstPerIteration)
    {
        printf("Error: Unfortunately, we cannot continue here because the indexing buffer is to small.\n");
        printf("You need to specify at least %s for the index buffer.\n", gk_GpfBase.bytesToStr(li_BiggestBucketSize * (mi_HmstBits) / 8 + 1).toStdString().c_str());
        exit(1);
    }
    
    // write HMST count bits
    qint32 li_HmstCountBits = 1;
    while ((((qint64)1 << li_HmstCountBits) - 1) < li_BiggestBucketSize)
        ++li_HmstCountBits;
    
    printf("HMST count bits: %d.\n", li_HmstCountBits);
    ak_OutFile_->write((char*)&li_HmstCountBits, 4);
    // encode HMST counts
    QSharedPointer<quint8> lui_pHmstCounts(new quint8[mi_TagCount * 2 * li_HmstCountBits / 8 + 1]);
    for (qint64 i = 0; i < mi_TagCount * 2; ++i)
    {
/*        if (li_pTagDirectionCount.data()[i] > 0)
            printf("%x %d\n", (unsigned int)i, (unsigned int)li_pTagDirectionCount.data()[i]);*/
        overwriteBitsInBuffer(lui_pHmstCounts.data(), i * li_HmstCountBits, li_pTagDirectionCount.data()[i], li_HmstCountBits);
    }
    // write HMST counts
    ak_OutFile_->write((char*)lui_pHmstCounts.data(), mi_TagCount * 2 * li_HmstCountBits / 8 + 1);
    
    // determine necessary iterations
    QList<QPair<int, int> > lk_IterationRanges;
    int li_TagDirection = 0;
    qint64 li_HmstForCurrentIteration = li_pTagDirectionCount.data()[0];
    lk_IterationRanges << QPair<int, int>(0, 0);
    while (li_TagDirection < mi_TagCount * 2)
    {
        if (li_HmstForCurrentIteration + li_pTagDirectionCount.data()[li_TagDirection] <= li_MaxHmstPerIteration)
        {
            lk_IterationRanges.last().second = li_TagDirection;
            li_HmstForCurrentIteration += li_pTagDirectionCount.data()[li_TagDirection];
        }
        else
        {
            lk_IterationRanges << QPair<int, int>(li_TagDirection, li_TagDirection);
            li_HmstForCurrentIteration = li_pTagDirectionCount.data()[li_TagDirection];
        }
        ++li_TagDirection;
    }
    
    QSharedPointer<qint64> li_pTagDirectionOffset(new qint64[mi_TagCount * 2]);
    
    QSharedPointer<k_BitWriter> lk_pBitWriter;
    lk_pBitWriter = QSharedPointer<k_BitWriter>(new k_BitWriter(ak_OutFile_));
    
    // perform iterations!
    printf("Required iterations: %d.\n", lk_IterationRanges.size());
    fflush(stdout);
    typedef QPair<int, int> tk_IntPair;
    int li_IterationCount = 0;
    foreach (tk_IntPair lk_Pair, lk_IterationRanges)
    {
        ++li_IterationCount;
        printf("Performing iteration %d of %d.\n", li_IterationCount, lk_IterationRanges.size());
        fflush(stdout);
        unsigned int lui_FirstTagDirectionIndex = lk_Pair.first;
        unsigned int lui_LastTagDirectionIndex = lk_Pair.second;
        
        // determine starting offsets for each tag/dir
        qint64 li_Offset = 0;
        for (unsigned int i = lui_FirstTagDirectionIndex; i <= lui_LastTagDirectionIndex; ++i)
        {
            li_pTagDirectionOffset.data()[i] = li_Offset;
            li_Offset += li_pTagDirectionCount.data()[i];
        }
        
//      printf("range: %d - %d\n", lui_FirstTagDirectionIndex, lui_LastTagDirectionIndex);
        
        lk_HmstIterator.reset();
        while (lk_HmstIterator.next(&lr_Hmst))
        {
            if (lr_Hmst.mui_TagDirectionIndex >= lui_FirstTagDirectionIndex && lr_Hmst.mui_TagDirectionIndex <= lui_LastTagDirectionIndex)
            {
/*              printf("HMST: %10.4f, %8d (%s/%s)\n", 
                       (double)lr_Hmst.mi_HalfMass / mi_MassPrecision, 
                       (unsigned int)lr_Hmst.mui_Gno,
                       gk_GpfBase.aminoAcidSequenceForCode(lr_Hmst.mui_TagDirectionIndex / 2, mi_TagSize).toStdString().c_str(),
                       lr_Hmst.mui_TagDirectionIndex % 2 == 0 ? "L" : "R"
                );*/
                // puts HMST into right place
                overwriteBitsInBuffer(muc_pIndexBuffer.data(), li_pTagDirectionOffset.data()[lr_Hmst.mui_TagDirectionIndex] * mi_HmstBits, lr_Hmst.mi_HalfMass, mi_MassBits);
                overwriteBitsInBuffer(muc_pIndexBuffer.data(), li_pTagDirectionOffset.data()[lr_Hmst.mui_TagDirectionIndex] * mi_HmstBits + mi_MassBits, lr_Hmst.mui_Gno, mi_OffsetBits);
/*              qint64 li_ControlMass = readBitsFromBuffer(muc_pIndexBuffer.data(), li_pTagDirectionOffset.data()[lr_Hmst.mui_TagDirectionIndex] * mi_HmstBits, mi_MassBits);
                qint64 li_ControlGno = readBitsFromBuffer(muc_pIndexBuffer.data(), li_pTagDirectionOffset.data()[lr_Hmst.mui_TagDirectionIndex] * mi_HmstBits + mi_MassBits, mi_OffsetBits);
                printf("%d/%d - %d/%d\n", (int)lr_Hmst.mi_HalfMass, (int)lr_Hmst.mui_Gno, 
                        (int)li_ControlMass, (int)li_ControlGno);*/
                // advance tag/dir offset
                ++li_pTagDirectionOffset.data()[lr_Hmst.mui_TagDirectionIndex];
            }
        }
        
        // sort each tag/dir list
        li_Offset = 0;

        for (unsigned int i = lui_FirstTagDirectionIndex; i <= lui_LastTagDirectionIndex; ++i)
        {
            if (li_pTagDirectionCount.data()[i] > 0)
            {
/*                if (li_pTagDirectionCount.data()[i] > 1)
                    printf("HMST COUNT %d\n", (unsigned int)li_pTagDirectionCount.data()[i]);*/
                typedef QPair<quint32, quint64> tk_MassGnoPair;
                QMap<quint32, tk_MassGnoPair> lk_Map;
                
                for (int k = 0; k < li_pTagDirectionCount.data()[i]; ++k)
                {
                    quint32 lui_Mass = readBitsFromBuffer(muc_pIndexBuffer.data(), (li_Offset + k) * mi_HmstBits, mi_MassBits);
                    quint64 lui_Gno = readBitsFromBuffer(muc_pIndexBuffer.data(), (li_Offset + k) * mi_HmstBits + mi_MassBits, mi_OffsetBits);
                    
                    lk_Map.insertMulti(lui_Mass, tk_MassGnoPair(lui_Mass, lui_Gno));
                }
                
                QMapIterator<quint32, tk_MassGnoPair> lk_Iterator(lk_Map);
                while (lk_Iterator.hasNext())
                {
                    lk_Iterator.next();
/*                    if (i / 2 == gk_GpfBase.aminoAcidPolymerCode("QPAWQ", 5))
                        printf("QPAWQ %9.4f\n", (double)lk_Iterator.value().first / mi_MassPrecision);
                    if (i / 2 == gk_GpfBase.aminoAcidPolymerCode("QGEQG", 5))
                        printf("QGEQG %9.4f\n", (double)lk_Iterator.value().first / mi_MassPrecision);*/
                    lk_pBitWriter->writeBits(lk_Iterator.value().first, mi_MassBits);
//                     printf("%x M %d\n", i, (unsigned int)lk_Iterator.value().first);
                }
                lk_Iterator.toFront();
                while (lk_Iterator.hasNext())
                {
                    lk_Iterator.next();
                    lk_pBitWriter->writeBits(lk_Iterator.value().second, mi_OffsetBits);
//                     printf("%x O %d\n", i, (unsigned int)lk_Iterator.value().second);
                }
                
                li_Offset += li_pTagDirectionCount.data()[i];
            }
        }
    } 
    // 2278.1066 2282.2233 QPAWQ QGEQG
    
    lk_pBitWriter->flush();
    
    // finish index chunk
    li_ChunkSize = ak_OutFile_->pos() - li_ChunkSize;
    writeChunkSize(ak_OutFile_, li_SizeLocation, li_ChunkSize);
}
