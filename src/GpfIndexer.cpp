#include "GpfIndexer.h"
#include "GpfIndexHeader.h"
#include "StopWatch.h"
#include "RefPtr.h"


k_GpfIndexer::k_GpfIndexer(QString as_GenomeFilename, QString as_GpfIndexFilename, QString as_GenomeTitle)
	: ms_GenomeFilename(as_GenomeFilename)
	, ms_GenomeTitle(as_GenomeTitle)
	, mk_GpfBase()
	, mk_IndexFileInfo(mk_GpfBase)
	, ms_IndexFilename(as_GpfIndexFilename)
{
}


k_GpfIndexer::~k_GpfIndexer()
{
}


void k_GpfIndexer::CompileIndex()
{
	k_StopWatch lk_StopWatch("Compiling the index took %1.\n");

	printf("Compiling DNA index...\n");

	mk_IndexFile.setFileName(ms_IndexFilename);

	if (mk_IndexFile.exists())
	{
		printf("Removing existing index file...\n");
		mk_IndexFile.remove();
	}

	if (!mk_IndexFile.open(QFile::ReadWrite))
	{
		fprintf(stderr, "Error: Unable to open index file for r/w access.\n");
		exit(1);
	}

	WriteHeader();

	WriteGenomeTitle();

	// add all scaffolds (nucleotides) plus the scaffold labels
	RecodeDna(ms_GenomeFilename);
	WriteScaffoldLabels();

	// add all scaffolds, all reading frames (amino acids translated from nucleotides)
	TranslateDna();

	// create index
	CreateIndex();
}


void k_GpfIndexer::WriteHeader()
{
	mk_IndexFile.write(QByteArray("gpfdna"));
	quint16 luw_VersionNumberHigh = 0;
	quint16 luw_VersionNumberLow = 1;
	mk_IndexFile.write((const char*)(&luw_VersionNumberHigh), sizeof(luw_VersionNumberHigh));
	mk_IndexFile.write((const char*)(&luw_VersionNumberLow), sizeof(luw_VersionNumberLow));
}


void k_GpfIndexer::WriteGenomeTitle()
{
	r_GenomeTitleHeader lr_GenomeTitleHeader(mk_IndexFile);
	lr_GenomeTitleHeader.WritePreliminary();
	quint16 lui_Length = ms_GenomeTitle.length();
	mk_IndexFile.write((const char*)(&lui_Length), sizeof(lui_Length));
	mk_IndexFile.write(ms_GenomeTitle.toStdString().c_str(), lui_Length);
	lr_GenomeTitleHeader.WriteFinal();
}


void k_GpfIndexer::RecodeDna(QString as_GenomeFilename)
{
	QFile lk_GenomeFile(as_GenomeFilename);
	if (!lk_GenomeFile.open(QFile::ReadOnly))
		return;

	r_DnaHeader lr_DnaHeader(mk_IndexFile);
	lr_DnaHeader.WritePreliminary();

	RefPtr<r_DnaScaffoldHeader> lr_pDnaScaffoldHeader(new r_DnaScaffoldHeader(mk_IndexFile));
	int li_Percent = -1;

	forever
	{
		char lc_Buffer_[1024];
		qint64 li_Length = lk_GenomeFile.readLine(lc_Buffer_, sizeof(lc_Buffer_));

		if (li_Length == -1)
			break;

		QString ls_Line(lc_Buffer_);
		ls_Line = ls_Line.trimmed();
		if (ls_Line.length() > 0 && ls_Line.at(0) == '>')
		{
			if (!mk_IndexFileInfo.mk_Scaffolds.empty())
			{
				lr_pDnaScaffoldHeader->mui_NucleotideCount = mk_IndexFileInfo.mk_Scaffolds.last().mui_NucleotideCount;
				lr_pDnaScaffoldHeader->WriteFinal();
				if (mk_IndexFileInfo.mk_Scaffolds.size() > 1)
				{
					const r_Scaffold& lr_PreviousScaffold = mk_IndexFileInfo.mk_Scaffolds.at(mk_IndexFileInfo.mk_Scaffolds.size() - 2);
					mk_IndexFileInfo.mk_Scaffolds.last().mui_NucleotideOffset = lr_PreviousScaffold.mui_NucleotideOffset + lr_PreviousScaffold.mui_NucleotideCount;
				}
			}

			QString ls_Id = ls_Line.mid(1);
			mk_IndexFileInfo.mk_Scaffolds.append(r_Scaffold(ls_Id, 0));

			lr_pDnaScaffoldHeader = RefPtr<r_DnaScaffoldHeader>(new r_DnaScaffoldHeader(mk_IndexFile));
			lr_pDnaScaffoldHeader->WritePreliminary();
			mk_IndexFileInfo.mk_Scaffolds.last().mi_DnaFilePosition = mk_IndexFile.pos();
		}
		else
		{
			int li_LineLength = ls_Line.length();
			for (int i = 0; i < li_LineLength; ++i)
			{
				char lc_Char = ls_Line.at(i).toAscii();
				char lc_Code = mk_GpfBase.mc_NucleotideCharToInt_[(unsigned char)lc_Char];
				lc_Buffer_[i] = lc_Code;
			}
			mk_IndexFile.write(lc_Buffer_, li_LineLength);
			mk_IndexFileInfo.mk_Scaffolds.last().mui_NucleotideCount += li_LineLength;

			int li_NewPercent = (int)((double)lk_GenomeFile.pos() / lk_GenomeFile.size() * 100.0);
			if (li_NewPercent != li_Percent)
			{
				li_Percent = li_NewPercent;
				printf("\rRecoding DNA... %d%% done.", li_Percent);
			}
		}
	}

	// TODO: ugly code copy from a few lines above
	if (!mk_IndexFileInfo.mk_Scaffolds.empty())
	{
		lr_pDnaScaffoldHeader->mui_NucleotideCount = mk_IndexFileInfo.mk_Scaffolds.last().mui_NucleotideCount;
		lr_pDnaScaffoldHeader->WriteFinal();
		if (mk_IndexFileInfo.mk_Scaffolds.size() > 1)
		{
			const r_Scaffold& lr_PreviousScaffold = mk_IndexFileInfo.mk_Scaffolds.at(mk_IndexFileInfo.mk_Scaffolds.size() - 2);
			mk_IndexFileInfo.mk_Scaffolds.last().mui_NucleotideOffset = lr_PreviousScaffold.mui_NucleotideOffset + lr_PreviousScaffold.mui_NucleotideCount;
		}
	}

	lr_DnaHeader.mui_ScaffoldCount = mk_IndexFileInfo.mk_Scaffolds.size();
	lr_DnaHeader.WriteFinal();

	printf("\rRecoding DNA... finished.    \n");
}


void k_GpfIndexer::WriteScaffoldLabels()
{
	r_ScaffoldLabelsHeader lr_ScaffoldLabelsHeader(mk_IndexFile);
	lr_ScaffoldLabelsHeader.WritePreliminary();

	foreach (r_Scaffold lr_Scaffold, mk_IndexFileInfo.mk_Scaffolds)
	{
		quint16 lui_Length = lr_Scaffold.ms_Id.length();
		mk_IndexFile.write((const char*)(&lui_Length), sizeof(lui_Length));
		mk_IndexFile.write(lr_Scaffold.ms_Id.toStdString().c_str(), lui_Length);
	}

	lr_ScaffoldLabelsHeader.WriteFinal();
}


void k_GpfIndexer::TranslateDna()
{
	// At this point, we now know the total scaffold count and the number of nucleotides per scaffold.
	// From this information, we can calculate exactly how many amino acids there will be in each reading frame.

	// translate DNA into six reading frames

	unsigned int lui_TotalNucleotideCount = 0;
	foreach (r_Scaffold lr_Scaffold, mk_IndexFileInfo.mk_Scaffolds)
		lui_TotalNucleotideCount += lr_Scaffold.mui_NucleotideCount;
	unsigned int lui_ProcessedNucleotideCount = 0;

	r_AminoAcidsHeader lr_AminoAcidsHeader(mk_IndexFile);
	lr_AminoAcidsHeader.WritePreliminary();

	qint64 li_ExtraSize = 0;
	qint64 li_FirstScaffoldPosition = mk_IndexFile.pos();

	// precalculate the size of all six reading frames including six headers for all scaffolds

	QVector<r_Scaffold>::iterator lk_Iter = mk_IndexFileInfo.mk_Scaffolds.begin();
	for (; lk_Iter != mk_IndexFileInfo.mk_Scaffolds.end(); ++lk_Iter)
	{
		r_Scaffold& lr_Scaffold = *lk_Iter;
		int li_Count0, li_Count1, li_Count2;
		get_AminoAcidCount(lr_Scaffold.mui_NucleotideCount, li_Count0, li_Count1, li_Count2);
		lr_Scaffold.mr_ReadingFrames[0].mui_AminoAcidCount = li_Count0;
		lr_Scaffold.mr_ReadingFrames[1].mui_AminoAcidCount = li_Count1;
		lr_Scaffold.mr_ReadingFrames[2].mui_AminoAcidCount = li_Count2;
		lr_Scaffold.mr_ReadingFrames[3].mui_AminoAcidCount = li_Count0;
		lr_Scaffold.mr_ReadingFrames[4].mui_AminoAcidCount = li_Count1;
		lr_Scaffold.mr_ReadingFrames[5].mui_AminoAcidCount = li_Count2;

		lr_Scaffold.mui_TotalAminoAcidCount = 0;
		for (int i = 0; i < 6; ++i)
			lr_Scaffold.mui_TotalAminoAcidCount += lr_Scaffold.mr_ReadingFrames[i].mui_AminoAcidCount;

		lr_Scaffold.mi_AminoAcidScaffoldFilePosition = li_FirstScaffoldPosition + li_ExtraSize;

		lr_Scaffold.mr_ReadingFrames[0].mi_ReadingFrameFilePosition = lr_Scaffold.mi_AminoAcidScaffoldFilePosition + r_AminoAcidsScaffoldHeader::getSize();
		for (int i = 1; i < 6; ++i)
			lr_Scaffold.mr_ReadingFrames[i].mi_ReadingFrameFilePosition = lr_Scaffold.mr_ReadingFrames[i - 1].mi_ReadingFrameFilePosition + 
			lr_Scaffold.mr_ReadingFrames[i - 1].mui_AminoAcidCount + r_AminoAcidsReadingFrameHeader::getSize();

		li_ExtraSize += r_AminoAcidsScaffoldHeader::getSize() + r_AminoAcidsReadingFrameHeader::getSize() * 6 + lr_Scaffold.mui_TotalAminoAcidCount;
	}

	// preallocate disk space
	mk_IndexFile.resize(mk_IndexFile.size() + li_ExtraSize);
	mk_IndexFile.seek(li_FirstScaffoldPosition);

	QFile lk_DnaReader(mk_IndexFile.fileName());
	lk_DnaReader.open(QFile::ReadOnly);

	// allocate temporary buffers for DNA to amino acid translation
	RefPtr<char> lc_pNucleotideBatch(new char[gui_MaxNucleotideBatchSize]);
	unsigned int lui_AminoAcidBatchSize = gui_MaxNucleotideBatchSize / 3 + 1;
	RefPtr<char> lc_pAminoAcidBatch_[6];
	for (int i = 0; i < 6; ++i)
		lc_pAminoAcidBatch_[i] = RefPtr<char>(new char[lui_AminoAcidBatchSize]);

	// the locations of the amino acid scaffold headers and their reading 
	// frames are all determined by now and stored in the scaffolds as  
	// mi_AminoAcidScaffoldPosition and mi_ReadingFramePosition

	int li_Percent = -1;

	// translate whole DNA
	for (int li_ScaffoldIndex = 0; li_ScaffoldIndex < mk_IndexFileInfo.mk_Scaffolds.size(); ++li_ScaffoldIndex)
	{
		r_Scaffold& lr_Scaffold = mk_IndexFileInfo.mk_Scaffolds[li_ScaffoldIndex];

		// jump to start of amino acid scaffold, write header
		mk_IndexFile.seek(lr_Scaffold.mi_AminoAcidScaffoldFilePosition);
		r_AminoAcidsScaffoldHeader lr_AminoAcidsScaffoldHeader(mk_IndexFile);
		lr_AminoAcidsScaffoldHeader.WritePreliminary();

		RefPtr<r_AminoAcidsReadingFrameHeader> lr_pAminoAcidsReadingFrames_[6];
		for (int li_ReadingFrame = 0; li_ReadingFrame < 6; ++li_ReadingFrame)
		{
			lr_pAminoAcidsReadingFrames_[li_ReadingFrame] = RefPtr<r_AminoAcidsReadingFrameHeader>(new r_AminoAcidsReadingFrameHeader(mk_IndexFile));
			lr_pAminoAcidsReadingFrames_[li_ReadingFrame]->mui_FrameNumber = li_ReadingFrame;
			lr_pAminoAcidsReadingFrames_[li_ReadingFrame]->mui_AminoAcidCount = lr_Scaffold.mr_ReadingFrames[li_ReadingFrame].mui_AminoAcidCount;
			mk_IndexFile.seek(lr_Scaffold.mr_ReadingFrames[li_ReadingFrame].mi_ReadingFrameFilePosition);
			lr_pAminoAcidsReadingFrames_[li_ReadingFrame]->WritePreliminary();
		}

		int li_AminoAcidsWritten[6] = {0, 0, 0, 0, 0, 0};

		// scan whole scaffold and fill the prepared reading frames with amino acids
		lk_DnaReader.seek(lr_Scaffold.mi_DnaFilePosition);

		unsigned lui_NucleotidesLeft = lr_Scaffold.mui_NucleotideCount;
		unsigned int lui_CurrentForwardFrame = 0;
		unsigned int lui_CurrentBackwardFrame = lr_Scaffold.mui_NucleotideCount % 3;

		unsigned int lui_Triplet = 0;
		unsigned int lui_ScaffoldProcessedNucleotideCount = 0;

		while (lui_NucleotidesLeft > 0)
		{
			// translate nucleotide batch from current scaffold
			char* lc_AminoAcidBatchForwardPointer_[3] = 
			{
				lc_pAminoAcidBatch_[0].get_Pointer(),
				lc_pAminoAcidBatch_[1].get_Pointer(),
				lc_pAminoAcidBatch_[2].get_Pointer()
			};

			char* lc_AminoAcidBatchBackwardPointer_[3] = 
			{
				lc_pAminoAcidBatch_[3].get_Pointer() + lui_AminoAcidBatchSize - 1,
				lc_pAminoAcidBatch_[4].get_Pointer() + lui_AminoAcidBatchSize - 1,
				lc_pAminoAcidBatch_[5].get_Pointer() + lui_AminoAcidBatchSize - 1
			};

			unsigned int lui_NucleotideBatchSize = std::min<unsigned int>(lui_NucleotidesLeft, gui_MaxNucleotideBatchSize);
			lui_NucleotidesLeft -= lui_NucleotideBatchSize;

			// read next batch of DNA
			lk_DnaReader.read(lc_pNucleotideBatch.get_Pointer(), lui_NucleotideBatchSize);

			lui_ProcessedNucleotideCount += lui_NucleotideBatchSize;
			
			int li_NewPercent = (int)((double)lui_ProcessedNucleotideCount / lui_TotalNucleotideCount * 100.0);
			if (li_NewPercent != li_Percent)
			{
				li_Percent = li_NewPercent;
				printf("\rTranslating DNA... %d%% done.", li_Percent);
			}

			char* lc_NucleotidePointer_ = lc_pNucleotideBatch.get_Pointer();
			for (unsigned int i = 0; i < lui_NucleotideBatchSize; ++i)
			{
				lui_Triplet <<= 3;
				lui_Triplet |= (*(lc_NucleotidePointer_++) & 7);
				lui_Triplet &= 511;

				if (lui_ScaffoldProcessedNucleotideCount < 3)
					++lui_ScaffoldProcessedNucleotideCount;

				// don't translate if we have less than three nucleotides in our triplet
				if (lui_ScaffoldProcessedNucleotideCount < 3)
					continue;

				// we now have a triplet in lui_Triplet (n1 came first): 
				// ... n1 n1 n1 n2 n2 n2 n3 n3 n3
				// each nucleotide is represented by three bits.
				// if bit 2 is set, the nucleotide is unknown, if it is cleared, bit 0 and 1 denote the nucleotide (A, C, G or T)

				// determine amino acid depending on current triplet, forward direction
				char lc_ForwardAminoAcid = mk_GpfBase.mc_TripletToAminoAcidForward_[lui_Triplet];
				// place current forward translated amino acid in current reading frame buffer
				*((lc_AminoAcidBatchForwardPointer_[lui_CurrentForwardFrame])++) = lc_ForwardAminoAcid;
				// advance current forward reading frame (plus 1)
				lui_CurrentForwardFrame = (lui_CurrentForwardFrame + 1) % 3;

				// determine amino acid depending on current triplet, backward direction
				char lc_BackwardAminoAcid = mk_GpfBase.mc_TripletToAminoAcidBackward_[lui_Triplet];
				// place current backward translated amino acid in current reading frame buffer
				*((lc_AminoAcidBatchBackwardPointer_[lui_CurrentBackwardFrame])--) = lc_BackwardAminoAcid;
				// advance current backward reading frame (minus 1, that is, plus 2)
				lui_CurrentBackwardFrame = (lui_CurrentBackwardFrame + 2) % 3;
			}

			// write amino acid forward batches
			for (int i = 0; i < 3; ++i)
			{
				unsigned int lui_AminoAcidCount = lc_AminoAcidBatchForwardPointer_[i] - lc_pAminoAcidBatch_[i].get_Pointer();
				mk_IndexFile.seek(lr_Scaffold.mr_ReadingFrames[i].mi_ReadingFrameFilePosition + r_AminoAcidsReadingFrameHeader::getSize() + li_AminoAcidsWritten[i]);
				mk_IndexFile.write(lc_pAminoAcidBatch_[i].get_Pointer(), lui_AminoAcidCount);
				li_AminoAcidsWritten[i] += lui_AminoAcidCount;
			}

			// write amino acid backward batches
			for (int i = 3; i < 6; ++i)
			{
				unsigned int lui_AminoAcidCount = lc_pAminoAcidBatch_[i].get_Pointer() + lui_AminoAcidBatchSize - 1 - lc_AminoAcidBatchBackwardPointer_[i - 3];
				mk_IndexFile.seek(lr_Scaffold.mr_ReadingFrames[i].mi_ReadingFrameFilePosition + r_AminoAcidsReadingFrameHeader::getSize() + lr_Scaffold.mr_ReadingFrames[i].mui_AminoAcidCount - li_AminoAcidsWritten[i] - lui_AminoAcidCount);
				mk_IndexFile.write(lc_AminoAcidBatchBackwardPointer_[i - 3] + 1, lui_AminoAcidCount);
				li_AminoAcidsWritten[i] += lui_AminoAcidCount;
			}
		}

		// seek to end position so that the calls to WriteFinal() won't be messed up
		mk_IndexFile.seek(lr_Scaffold.mr_ReadingFrames[5].mi_ReadingFrameFilePosition + r_AminoAcidsReadingFrameHeader::getSize() + lr_Scaffold.mr_ReadingFrames[5].mui_AminoAcidCount);

		for (int li_ReadingFrame = 0; li_ReadingFrame < 5; ++li_ReadingFrame)
			lr_pAminoAcidsReadingFrames_[li_ReadingFrame]->WriteFinal(lr_Scaffold.mr_ReadingFrames[li_ReadingFrame + 1].mi_ReadingFrameFilePosition);

		//qint64 li_NextScaffoldPosition = (mk_Scaffolds[li_ScaffoldIndex + 1].mi_AminoAcidScaffoldFilePosition)
		lr_pAminoAcidsReadingFrames_[5]->WriteFinal();
		lr_AminoAcidsScaffoldHeader.WriteFinal();
	}
	lr_AminoAcidsHeader.WriteFinal();
	printf("\rTranslating DNA... finished.   \n");
}


void k_GpfIndexer::CreateIndex()
{
	// index 20 amino acids with 3 digits, this makes for 8000 different amino acid trimers
	int li_States = 20;
	int li_Digits = 3;

	int li_OffsetCount = 1;
	for (int i = 0; i < li_Digits; ++i)
		li_OffsetCount *= li_States;

	memset(mk_IndexFileInfo.mui_pTrimerOffsetTable.get_Pointer(), 0, sizeof(int) * li_OffsetCount);
	memset(mk_IndexFileInfo.mui_pTrimerCountTable.get_Pointer(), 0, sizeof(int) * li_OffsetCount);

	RefPtr<unsigned int> lui_pMassTable(new unsigned int[li_OffsetCount]);

	for (int i = 0; i < li_OffsetCount; ++i)
	{
		r_AminoAcid::Enumeration le_First = (r_AminoAcid::Enumeration)((i / 400) % 20);
		r_AminoAcid::Enumeration le_Second = (r_AminoAcid::Enumeration)((i / 20) % 20);
		r_AminoAcid::Enumeration le_Third = (r_AminoAcid::Enumeration)(i % 20);
		lui_pMassTable.get_Pointer()[i] = mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)le_First] + 
			mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)le_Second] + 
			mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)le_Third];
	}

	QFile lk_TrimerReader(mk_IndexFile.fileName());
	lk_TrimerReader.open(QFile::ReadOnly);

	QFile lk_AminoAcidsReader(mk_IndexFile.fileName());
	lk_AminoAcidsReader.open(QFile::ReadOnly);

	RefPtr<char> lc_pAminoAcidBuffer(new char[gui_MaxNucleotideBatchSize]);

	unsigned int lui_TotalPeptideCount = 0;
	unsigned int lui_TotalIndexedTrimerCount = 0;
	qint64 li_AveragePeptideLength = 0;
	unsigned int lui_MinPeptideLength = 0;
	unsigned int lui_MaxPeptideLength = 0;

	unsigned int lui_TotalProcessedAminoAcidCount = 0;

	unsigned int lui_TotalAminoAcidCount = 0;
	foreach (r_Scaffold lr_Scaffold, mk_IndexFileInfo.mk_Scaffolds)
		lui_TotalAminoAcidCount += lr_Scaffold.mui_TotalAminoAcidCount;

	RefPtr<r_PeptideSpan> lr_pPeptideSpans(new r_PeptideSpan[gui_SortBatchSize]);

	QTemporaryFile lk_TempFile(mk_IndexFile.fileName());
	lk_TempFile.open();

	RefPtr<r_IndexEntry> lr_pIndexEntries(new r_IndexEntry[gui_SortBatchSize]);
	unsigned int lui_IndexedTrimerBatchCount = 0;

	int li_Percent = -1;

	r_PeptideSpanListHeader lr_PeptideSpanListHeader(mk_IndexFile);
	lr_PeptideSpanListHeader.WritePreliminary();

	// iterate through all amino acids, count and index trimers on the fly, find cleavage sites
	for (unsigned int lui_ReadingFrame = 0; lui_ReadingFrame < 6; ++lui_ReadingFrame)
	{
		r_PeptideSpanListReadingFrameHeader lr_PeptideSpanListReadingFrameHeader(mk_IndexFile);
		lr_PeptideSpanListReadingFrameHeader.WritePreliminary();

		QTemporaryFile lk_PeptideSpansTempFile(mk_IndexFile.fileName());
		unsigned int lui_PeptideSpansBatchCount = 0;

		foreach (r_Scaffold lr_Scaffold, mk_IndexFileInfo.mk_Scaffolds)
		{
			r_AminoAcid::Enumeration le_PreviousAminoAcid = r_AminoAcid::Unknown;

			lk_TrimerReader.seek(lr_Scaffold.mr_ReadingFrames[lui_ReadingFrame].mi_ReadingFrameFilePosition + r_AminoAcidsReadingFrameHeader::getSize());
			unsigned int lui_Count = lr_Scaffold.mr_ReadingFrames[lui_ReadingFrame].mui_AminoAcidCount;

			unsigned int lui_CharTrimer = 0;
			unsigned int lui_ProcessedAminoAcidCount = 0;

			QString ls_Peptide = "";
			unsigned int lui_PeptideOffset = 0;
			unsigned int lui_PeptideUnknownAminoAcidCount = 0;

			while (lui_Count > 0)
			{
				unsigned int lui_BatchSize = std::min<unsigned int>(lui_Count, gui_MaxNucleotideBatchSize);

				lui_TotalProcessedAminoAcidCount += lui_BatchSize;
				int li_NewPercent = (int)((double)lui_TotalProcessedAminoAcidCount / lui_TotalAminoAcidCount * 100.0);
				if (li_NewPercent != li_Percent)
				{
					li_Percent = li_NewPercent;
					printf("\rIndexing peptides... %d%% done.", li_Percent);
				}

				lui_Count -= lui_BatchSize;
				lk_TrimerReader.read((char*)(lc_pAminoAcidBuffer.get_Pointer()), lui_BatchSize);

				char* lc_AminoAcidBufferPointer_ = lc_pAminoAcidBuffer.get_Pointer();

				for (unsigned int lui_BatchIndex = 0; lui_BatchIndex < lui_BatchSize; ++lui_BatchIndex)
				{
					++lui_ProcessedAminoAcidCount;
					r_AminoAcid::Enumeration le_AminoAcid = (r_AminoAcid::Enumeration)(*((unsigned char*)lc_AminoAcidBufferPointer_++));
					char lc_AminoAcid = mk_GpfBase.mc_AminoAcidToChar_[(unsigned char)le_AminoAcid];

					if (le_AminoAcid != r_AminoAcid::Stop)
						ls_Peptide += lc_AminoAcid;

					if (lc_AminoAcid == 'X')
						++lui_PeptideUnknownAminoAcidCount;

					lui_CharTrimer <<= 8;
					lui_CharTrimer |= (char)le_AminoAcid;
					lui_CharTrimer &= 0xffffff;

					if ((le_AminoAcid == r_AminoAcid::Stop) || (/*(le_PreviousAminoAcid != r_AminoAcid::Pro) && */((le_AminoAcid == r_AminoAcid::Lys) || (le_AminoAcid == r_AminoAcid::Arg))))
					{
						if (lui_PeptideUnknownAminoAcidCount <= gui_MaxUnknownAminoAcidCount && ls_Peptide.length() > (int)gui_MaxPeptideLength)
							printf("\nAttention: A peptide was ignored because of it is too long (%d amino acids).\n", ls_Peptide.length());

						// yay, we found a tryptic cleavage site
						if (ls_Peptide.length() >= 3 && ls_Peptide.length() <= (int)gui_MaxPeptideLength && lui_PeptideUnknownAminoAcidCount <= gui_MaxUnknownAminoAcidCount)
						{
							unsigned int lui_Gno = lr_Scaffold.getGlobalNucleotideOffset(lui_ReadingFrame, lui_PeptideOffset);

							// we have a peptide with a reasonable length and at most one unkown amino acid --> add to index!
							unsigned int lui_Length = ls_Peptide.length();

							// add this peptide span to peptide span list
							r_PeptideSpan& lr_PeptideSpan = lr_pPeptideSpans.get_Pointer()[lui_PeptideSpansBatchCount];

							lr_PeptideSpan.mui_Gno = lui_Gno;
							lr_PeptideSpan.muw_Length = lui_Length;

							++lui_PeptideSpansBatchCount;
							++lr_PeptideSpanListReadingFrameHeader.mui_PeptideCount;

							if (lui_PeptideSpansBatchCount == gui_SortBatchSize)
							{
								// sort array and flush to peptide span temp file
								Sort<r_PeptideSpan>(lr_pPeptideSpans.get_Pointer(), 0, lui_PeptideSpansBatchCount - 1, &lt_PeptideSpan_Gno);
								lk_PeptideSpansTempFile.write((const char*)lr_pPeptideSpans.get_Pointer(), sizeof(r_PeptideSpan) * lui_PeptideSpansBatchCount);
								lui_PeptideSpansBatchCount = 0;
							}

							unsigned int lui_EnumTrimer = ((unsigned int)mk_GpfBase.me_CharToAminoAcid_[(unsigned char)ls_Peptide.at(0).toAscii()]);
							lui_EnumTrimer <<= 8;
							lui_EnumTrimer |= ((unsigned int)mk_GpfBase.me_CharToAminoAcid_[(unsigned char)ls_Peptide.at(1).toAscii()]);

							unsigned int lui_LeftMass = 0;
							unsigned int lui_RightMass = 0;
							for (unsigned int i = 3; i < lui_Length; ++i)
								lui_RightMass += mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)mk_GpfBase.me_CharToAminoAcid_[(unsigned char)ls_Peptide.at(i).toAscii()]];

							// iterate over all trimers
							for (unsigned int i = 2; i < lui_Length; ++i)
							{
								lui_EnumTrimer <<= 8;
								lui_EnumTrimer |= ((unsigned int)mk_GpfBase.me_CharToAminoAcid_[(unsigned char)ls_Peptide.at(i).toAscii()]);
								lui_EnumTrimer &= 0xffffff;
								if (((lui_EnumTrimer & 0xff) < 20) && (((lui_EnumTrimer >> 8) & 0xff) < 20) && (((lui_EnumTrimer >> 16) & 0xff) < 20))
								{
									// We have a trimer that should be indexed, calculate the trimer id

									// Fix this trimer to accomodate for mass ambiguities that result from equal 
									// or close to equal masses (look up in me_AminoAcidMassSimilarityCollapse).
									//unsigned int lui_Trimer = ((lui_EnumTrimer >> 16) & 0xff) * 400 + ((lui_EnumTrimer >> 8) & 0xff) * 20 + (lui_EnumTrimer & 0xff);
									unsigned int lui_Trimer = mk_GpfBase.me_AminoAcidMassSimilarityCollapse[(lui_EnumTrimer >> 16) & 0xff] * 400 + 
										mk_GpfBase.me_AminoAcidMassSimilarityCollapse[(lui_EnumTrimer >> 8) & 0xff] * 20 + 
										mk_GpfBase.me_AminoAcidMassSimilarityCollapse[lui_EnumTrimer & 0xff];

									++mk_IndexFileInfo.mui_pTrimerCountTable.get_Pointer()[lui_Trimer];

									unsigned int lui_Gno = lr_Scaffold.getGlobalNucleotideOffset(lui_ReadingFrame, lui_PeptideOffset + i - 2);

									r_IndexEntry& lr_IndexEntry = lr_pIndexEntries.get_Pointer()[lui_IndexedTrimerBatchCount];

									lr_IndexEntry.muw_Trimer = (unsigned short)lui_Trimer;
									lr_IndexEntry.mui_Gno = lui_Gno;
									lr_IndexEntry.mui_LeftMass = lui_LeftMass;
									lr_IndexEntry.mui_RightMass = lui_RightMass;

									++lui_IndexedTrimerBatchCount;
									++lui_TotalIndexedTrimerCount;

									if (lui_IndexedTrimerBatchCount == gui_SortBatchSize)
									{
										// sort array and flush to temp file
										Sort<r_IndexEntry>(lr_pIndexEntries.get_Pointer(), 0, lui_IndexedTrimerBatchCount - 1, &lt_IndexEntry_Trimer_LeftMass_Gno);
										lk_TempFile.write((const char*)lr_pIndexEntries.get_Pointer(), sizeof(r_IndexEntry) * lui_IndexedTrimerBatchCount);
										lui_IndexedTrimerBatchCount = 0;
									}
								}
								if (i < lui_Length - 1)
								{
									lui_LeftMass += mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)(r_AminoAcid::Enumeration)((lui_EnumTrimer >> 16) & 0xff)];
									lui_RightMass -= mk_GpfBase.mui_AminoAcidWeight_[(unsigned char)mk_GpfBase.me_CharToAminoAcid_[(unsigned char)ls_Peptide.at(i + 1).toAscii()]];
								}
							}

							if (lui_TotalPeptideCount == 0)
								lui_MaxPeptideLength = lui_MinPeptideLength = lui_Length;
							else
							{
								if (lui_Length > lui_MaxPeptideLength)
									lui_MaxPeptideLength = lui_Length;
								if (lui_Length < lui_MinPeptideLength)
									lui_MinPeptideLength = lui_Length;
							}
							li_AveragePeptideLength += lui_Length;
							++lui_TotalPeptideCount;
						}
						ls_Peptide = "";
						lui_PeptideUnknownAminoAcidCount = 0;
						lui_PeptideOffset = lui_ProcessedAminoAcidCount;
					}

					le_PreviousAminoAcid = le_AminoAcid;
				}
			}
		}

		// finish sorting of peptide spans for current reading frame
		if (lui_PeptideSpansBatchCount > 0)
		{
			// sort array and flush to peptide span temp file
			Sort<r_PeptideSpan>(lr_pPeptideSpans.get_Pointer(), 0, lui_PeptideSpansBatchCount - 1, &lt_PeptideSpan_Gno);
			lk_PeptideSpansTempFile.write((const char*)lr_pPeptideSpans.get_Pointer(), sizeof(r_PeptideSpan) * lui_PeptideSpansBatchCount);
		}

		// merge sorted chunks in the peptide spans temp file into the index file
		QVector<unsigned int> lk_ChunkOffsets;
		unsigned int lui_Offset = 0;
		while (lui_Offset < lr_PeptideSpanListReadingFrameHeader.mui_PeptideCount)
		{
			lk_ChunkOffsets.append(lui_Offset);
			lui_Offset += gui_SortBatchSize;
		}

		QTemporaryFile lk_TempFile(mk_IndexFile.fileName());
		QFile* lk_ResultFile_ = this->SortEntriesOnDisk<r_PeptideSpan>(lk_PeptideSpansTempFile, lk_TempFile, lk_ChunkOffsets, lr_PeptideSpanListReadingFrameHeader.mui_PeptideCount, &lt_PeptideSpan_Gno, false);

		lk_ResultFile_->open(QIODevice::ReadOnly);
		lk_ResultFile_->reset();
		unsigned int lui_Count = lr_PeptideSpanListReadingFrameHeader.mui_PeptideCount;
		while (lui_Count > 0)
		{
			unsigned int lui_BatchSize = std::min<unsigned int>(gui_SortBatchSize, lui_Count);
			lui_Count -= lui_BatchSize;
			lk_ResultFile_->read((char*)lr_pPeptideSpans.get_Pointer(), sizeof(r_PeptideSpan) * lui_BatchSize);
			mk_IndexFile.write((const char*)lr_pPeptideSpans.get_Pointer(), sizeof(r_PeptideSpan) * lui_BatchSize);
		}

		lr_PeptideSpanListReadingFrameHeader.WriteFinal();
	}

	lr_PeptideSpanListHeader.WriteFinal();

	printf("\rIndexing peptides... finished.     \n");
	
	li_AveragePeptideLength /= lui_TotalPeptideCount;

	if (lui_IndexedTrimerBatchCount > 0)
	{
		// sort array and flush to temp file
		Sort<r_IndexEntry>(lr_pIndexEntries.get_Pointer(), 0, lui_IndexedTrimerBatchCount - 1, &lt_IndexEntry_Trimer_LeftMass_Gno);
		lk_TempFile.write((const char*)lr_pIndexEntries.get_Pointer(), sizeof(r_IndexEntry) * lui_IndexedTrimerBatchCount);
	}

	// derive trimer offset array from trimer count array
	mk_IndexFileInfo.mui_pTrimerOffsetTable.get_Pointer()[0] = 0;
	for (int i = 1; i < li_OffsetCount; ++i)
		mk_IndexFileInfo.mui_pTrimerOffsetTable.get_Pointer()[i] = 
			mk_IndexFileInfo.mui_pTrimerOffsetTable.get_Pointer()[i - 1] + mk_IndexFileInfo.mui_pTrimerCountTable.get_Pointer()[i - 1];

	// TODO: crash if too many entries for a certain trimer, we should be able to handle that case
	bool lb_SizesOk = true;
	for (int i = 0; i < 8000; ++i)
	{
		if (mk_IndexFileInfo.mui_pTrimerCountTable.get_Pointer()[i] > gui_SortBatchSize)
		{
			printf("Attention: more entries for trimer %d (%d) than currently supported... (to be implemented).\n", i, mk_IndexFileInfo.mui_pTrimerCountTable.get_Pointer()[i]);
			lb_SizesOk = false;
		}
	}
	if (!lb_SizesOk)
	{
		mk_IndexFile.remove();
		exit(1);
	}

	r_IndexHeader lr_IndexHeader(mk_IndexFile);
	lr_IndexHeader.WritePreliminary();
	lr_IndexHeader.mui_TotalIndexedTrimerCount = lui_TotalIndexedTrimerCount;

	r_IndexOffsetsHeader lr_IndexOffsetsHeader(mk_IndexFile);
	lr_IndexOffsetsHeader.WritePreliminary();
	mk_IndexFile.write((const char*)mk_IndexFileInfo.mui_pTrimerOffsetTable.get_Pointer(), sizeof(unsigned int) * li_OffsetCount);
	lr_IndexOffsetsHeader.WriteFinal();

	// The temporary file now contains chunks which are each correctly sorted,
	// now we must merge all sorted chunks into one huge sorted chunk.

	QVector<unsigned int> lk_ChunkOffsets;
	unsigned int lui_Offset = 0;
	while (lui_Offset < lui_TotalIndexedTrimerCount)
	{
		lk_ChunkOffsets.append(lui_Offset);
		lui_Offset += gui_SortBatchSize;
	}

	QTemporaryFile lk_TempFile2(lk_TempFile.fileName());
	QFile* lk_ResultFile_ = this->SortEntriesOnDisk<r_IndexEntry>(lk_TempFile, lk_TempFile2, lk_ChunkOffsets, lui_TotalIndexedTrimerCount, &lt_IndexEntry_Trimer_LeftMass_Gno);
	this->WriteMassTables(*lk_ResultFile_, lui_TotalIndexedTrimerCount);

	lr_IndexHeader.WriteFinal();

	mk_IndexFile.seek(mk_IndexFile.size());

	printf("\n");
	printf("Total number of peptides  : %d\n", lui_TotalPeptideCount);
	printf("Number of indexed trimers : %d\n", lui_TotalIndexedTrimerCount);
	printf("Minimum peptide length    : %d\n", lui_MinPeptideLength);
	printf("Maximum peptide length    : %d\n", lui_MaxPeptideLength);
	printf("Average peptide length    : %d\n", (int)li_AveragePeptideLength);
}


void k_GpfIndexer::get_AminoAcidCount(int ai_NucleotideCount, int& ai_Count0, int& ai_Count1, int& ai_Count2)
{
	int n = ai_NucleotideCount / 3;
	switch (ai_NucleotideCount % 3)
	{
		case 0:
			ai_Count0 = n;
			ai_Count1 = n - 1;
			ai_Count2 = n - 1;
			break;
		case 1:
			ai_Count0 = n;
			ai_Count1 = n;
			ai_Count2 = n - 1;
			break;
		case 2:
			ai_Count0 = n;
			ai_Count1 = n;
			ai_Count2 = n;
			break;
	}
}


template <typename T> QFile* k_GpfIndexer::SortEntriesOnDisk(QFile& ak_SourceFile, QFile& ak_TemporaryFile, QVector<unsigned int> ak_ChunkOffsets, unsigned int aui_EntryCount, bool (*af_LessThan_)(T*, T*), bool ab_Verbose)
{
	if (ak_ChunkOffsets.size() < 2)
		return &ak_SourceFile;

	QFile* lk_SourceFile_ = &ak_SourceFile;
	QFile* lk_DestinationFile_ = &ak_TemporaryFile;

	RefPtr<T> lr_pLeftChunk(new T[gui_SortBatchSize]);
	RefPtr<T> lr_pRightChunk(new T[gui_SortBatchSize]);
	RefPtr<T> lr_pMergedChunk(new T[gui_SortBatchSize]);

	unsigned int lui_MergedChunkOffset = 0;

	// precalculate number of sorting steps for progress notification
	unsigned int lui_TotalStepCount = 0;
	int li_TempCount = ak_ChunkOffsets.size();
	do
	{
		++lui_TotalStepCount;
		li_TempCount = ((li_TempCount - 1) / 2) + 1;
	} while (li_TempCount > 1);

	int li_StepsPerformed = 0;
	// perform ping-pong style merging of sorted chunks (sorted by trimer id, then left mass)
	do
	{
		++li_StepsPerformed;
		if (ab_Verbose)
			printf("\rSorting on disk... performing step %d of %d.", li_StepsPerformed, lui_TotalStepCount);

		// rewind destination file
		lk_SourceFile_->open(QIODevice::ReadOnly);
		lk_DestinationFile_->open(QIODevice::WriteOnly);
		lk_DestinationFile_->seek(0);

		// merge every pair of consequtive chunks
		int li_Index = 0;
		while (li_Index + 1 < ak_ChunkOffsets.size())
		{
			// we have a pair that should be merged
			qint64 li_LeftCount = ak_ChunkOffsets[li_Index + 1] - ak_ChunkOffsets[li_Index];
			qint64 li_RightCount;
			if (li_Index + 2 < ak_ChunkOffsets.size())
				li_RightCount = ak_ChunkOffsets[li_Index + 2] - ak_ChunkOffsets[li_Index + 1];
			else
				li_RightCount = aui_EntryCount - ak_ChunkOffsets[li_Index + 1];

			qint64 li_LeftOffset = 0;
			qint64 li_RightOffset = 0;
			unsigned int lui_LeftBufferSize = 0;
			unsigned int lui_LeftBufferOffset = 0;
			unsigned int lui_RightBufferSize = 0;
			unsigned int lui_RightBufferOffset = 0;

			// repeat while not at end in left and right chunk
			while (li_LeftOffset < li_LeftCount && li_RightOffset < li_RightCount)
			{
				// refill left and right buffers if necessary
				if (lui_LeftBufferSize == lui_LeftBufferOffset)
				{
					unsigned int lui_BatchSize = std::min<unsigned int>(gui_SortBatchSize, (unsigned int)(li_LeftCount - li_LeftOffset));
					lk_SourceFile_->seek((ak_ChunkOffsets[li_Index] + li_LeftOffset) * sizeof(T));
					lk_SourceFile_->read((char*)lr_pLeftChunk.get_Pointer(), lui_BatchSize * sizeof(T));
					lui_LeftBufferSize = lui_BatchSize;
					lui_LeftBufferOffset = 0;
				}
				if (lui_RightBufferSize == lui_RightBufferOffset)
				{
					unsigned int lui_BatchSize = std::min<unsigned int>(gui_SortBatchSize, (unsigned int)(li_RightCount - li_RightOffset));
					lk_SourceFile_->seek((ak_ChunkOffsets[li_Index + 1] + li_RightOffset) * sizeof(T));
					lk_SourceFile_->read((char*)lr_pRightChunk.get_Pointer(), lui_BatchSize * sizeof(T));
					lui_RightBufferSize = lui_BatchSize;
					lui_RightBufferOffset = 0;
				}

				if ((*af_LessThan_)(&lr_pLeftChunk.get_Pointer()[lui_LeftBufferOffset], &lr_pRightChunk.get_Pointer()[lui_RightBufferOffset]))
				{
					// choose left entry
					lr_pMergedChunk.get_Pointer()[lui_MergedChunkOffset++] = lr_pLeftChunk.get_Pointer()[lui_LeftBufferOffset++];
					if (lui_MergedChunkOffset == gui_SortBatchSize)
					{
						lk_DestinationFile_->write((const char*)lr_pMergedChunk.get_Pointer(), lui_MergedChunkOffset * sizeof(T));
						lui_MergedChunkOffset = 0;
					}
					++li_LeftOffset;
				}
				else
				{
					// choose right entry
					lr_pMergedChunk.get_Pointer()[lui_MergedChunkOffset++] = lr_pRightChunk.get_Pointer()[lui_RightBufferOffset++];
					if (lui_MergedChunkOffset == gui_SortBatchSize)
					{
						lk_DestinationFile_->write((const char*)lr_pMergedChunk.get_Pointer(), lui_MergedChunkOffset * sizeof(T));
						lui_MergedChunkOffset = 0;
					}
					++li_RightOffset;
				}
			}

			// flush remaining entries in merge buffer
			if (lui_MergedChunkOffset > 0)
			{
				lk_DestinationFile_->write((const char*)lr_pMergedChunk.get_Pointer(), lui_MergedChunkOffset * sizeof(T));
				lui_MergedChunkOffset = 0;
			}

			// copy remaining contents of left chunk, if any
			while (li_LeftOffset < li_LeftCount)
			{
				if (lui_LeftBufferSize == lui_LeftBufferOffset)
				{
					unsigned int lui_BatchSize = std::min<unsigned int>(gui_SortBatchSize, (unsigned int)(li_LeftCount - li_LeftOffset));
					lk_SourceFile_->seek((ak_ChunkOffsets[li_Index] + li_LeftOffset) * sizeof(T));
					lk_SourceFile_->read((char*)lr_pLeftChunk.get_Pointer(), lui_BatchSize * sizeof(T));
					lui_LeftBufferSize = lui_BatchSize;
					lui_LeftBufferOffset = 0;
				}

				lk_DestinationFile_->write((const char*)&lr_pLeftChunk.get_Pointer()[lui_LeftBufferOffset], sizeof(T) * (lui_LeftBufferSize - lui_LeftBufferOffset));

				li_LeftOffset += (lui_LeftBufferSize - lui_LeftBufferOffset);
				lui_LeftBufferOffset = lui_LeftBufferSize;
			}

			// copy remaining contents of right chunk, if any
			while (li_RightOffset < li_RightCount)
			{
				if (lui_RightBufferSize == lui_RightBufferOffset)
				{
					unsigned int lui_BatchSize = std::min<unsigned int>(gui_SortBatchSize, (unsigned int)(li_RightCount - li_RightOffset));
					lk_SourceFile_->seek((ak_ChunkOffsets[li_Index + 1] + li_RightOffset) * sizeof(T));
					lk_SourceFile_->read((char*)lr_pRightChunk.get_Pointer(), lui_BatchSize * sizeof(T));
					lui_RightBufferSize = lui_BatchSize;
					lui_RightBufferOffset = 0;
				}

				lk_DestinationFile_->write((const char*)&lr_pRightChunk.get_Pointer()[lui_RightBufferOffset], sizeof(T) * (lui_RightBufferSize - lui_RightBufferOffset));

				li_RightOffset += (lui_RightBufferSize - lui_RightBufferOffset);
				lui_RightBufferOffset = lui_RightBufferSize;
			}

			// advance to next pair of consecutive chunks
			li_Index += 2;
		}

		// copy remaining last chunk, if any
		if (li_Index < ak_ChunkOffsets.size())
		{
			lk_SourceFile_->seek(ak_ChunkOffsets[li_Index] * sizeof(r_IndexEntry));
			qint64 li_Count = aui_EntryCount - ak_ChunkOffsets[li_Index];
			while (li_Count > 0)
			{
				unsigned int lui_BatchSize = std::min<unsigned int>(gui_SortBatchSize, (unsigned int)li_Count);
				lk_SourceFile_->read((char*)lr_pLeftChunk.get_Pointer(), lui_BatchSize * sizeof(T));
				lk_DestinationFile_->write((const char*)lr_pLeftChunk.get_Pointer(), lui_BatchSize * sizeof(T));
				li_Count -= lui_BatchSize;
			}
		}

		// remove every second entry from chunk offsets vector
		li_Index = 1;
		while (li_Index < ak_ChunkOffsets.size())
			ak_ChunkOffsets.remove(li_Index++);

		// swap file source and destination file
		QFile* lk_Temp_ = lk_SourceFile_;
		lk_SourceFile_ = lk_DestinationFile_;
		lk_DestinationFile_ = lk_Temp_;
	} while (ak_ChunkOffsets.size() > 1);

	if (ab_Verbose)
		printf("\rSorting on disk... finished.                          \n");
	return lk_SourceFile_;
}


void k_GpfIndexer::WriteMassTables(QFile& ak_SourceFile, unsigned int aui_EntryCount)
{
	// Now we have all entries sorted by trimer, then left mass in the source file.
	// Next thing is to convert the sorted r_IndexEntry entries to r_IndexFileEntries,
	// that is, chuck out unimportant data

	ak_SourceFile.open(QIODevice::ReadOnly);
	RefPtr<r_IndexEntry> lr_pLeftChunk(new r_IndexEntry[gui_SortBatchSize]);

	printf("Writing left mass index... ");

	// copy left tag sorted entries to output file (masses and GNO)
	RefPtr<unsigned int> lui_pMasses(new unsigned int[gui_SortBatchSize]);
	RefPtr<unsigned int> lui_pGno(new unsigned int[gui_SortBatchSize]);

	r_IndexLeftTagMassesHeader lr_IndexLeftTagMassesHeader(mk_IndexFile);
	lr_IndexLeftTagMassesHeader.WritePreliminary();
	mk_IndexFileInfo.mi_LeftTagMassesFilePosition = mk_IndexFile.pos();
	mk_IndexFile.resize(mk_IndexFileInfo.mi_LeftTagMassesFilePosition + (qint64)aui_EntryCount * sizeof(unsigned int));
	mk_IndexFile.seek(mk_IndexFileInfo.mi_LeftTagMassesFilePosition + (qint64)aui_EntryCount * sizeof(unsigned int));
	r_IndexLeftTagGnoHeader lr_IndexLeftTagGnoHeader(mk_IndexFile);
	lr_IndexLeftTagGnoHeader.WritePreliminary();
	mk_IndexFileInfo.mi_LeftTagGnoFilePosition = mk_IndexFile.pos();

	ak_SourceFile.seek(0);

	int li_Percent = -1;
	unsigned int lui_Count = aui_EntryCount;
	unsigned int lui_ProcessedTrimerCount = 0;
	unsigned int lui_WrittenEntries = 0;
	while (lui_Count > 0)
	{
		unsigned int lui_BatchSize = std::min<unsigned int>(gui_SortBatchSize, (unsigned int)lui_Count);
		lui_Count -= lui_BatchSize;

		ak_SourceFile.read((char*)lr_pLeftChunk.get_Pointer(), lui_BatchSize * sizeof(r_IndexEntry));
		r_IndexEntry* lr_IndexEntry_ = lr_pLeftChunk.get_Pointer();
		r_IndexEntry* lr_IndexEntryEnd_ = lr_IndexEntry_ + lui_BatchSize;
		unsigned int* lui_Mass_ = lui_pMasses.get_Pointer();
		unsigned int* lui_Gno_ = lui_pGno.get_Pointer();

		while (lr_IndexEntry_ < lr_IndexEntryEnd_)
		{
			++lui_ProcessedTrimerCount;
			int li_NewPercent = (int)((double)lui_ProcessedTrimerCount / aui_EntryCount * 100.0);
			if (li_NewPercent != li_Percent)
			{
				printf("\rWriting left mass index... %d%% done.", li_NewPercent);
				li_Percent = li_NewPercent;
			}
			*(lui_Mass_++) = (lr_IndexEntry_)->mui_LeftMass;
			*(lui_Gno_++) = (lr_IndexEntry_)->mui_Gno;
			++lr_IndexEntry_;
		}
		mk_IndexFile.seek(mk_IndexFileInfo.mi_LeftTagMassesFilePosition + lui_WrittenEntries * sizeof(unsigned int));
		mk_IndexFile.write((const char*)lui_pMasses.get_Pointer(), lui_BatchSize * sizeof(unsigned int));
		mk_IndexFile.seek(mk_IndexFileInfo.mi_LeftTagGnoFilePosition + lui_WrittenEntries * sizeof(unsigned int));
		mk_IndexFile.write((const char*)lui_pGno.get_Pointer(), lui_BatchSize * sizeof(unsigned int));
		lui_WrittenEntries += lui_BatchSize;
	}

	lr_IndexLeftTagMassesHeader.WriteFinal(mk_IndexFileInfo.mi_LeftTagGnoFilePosition - r_IndexLeftTagGnoHeader::getSize());
	lr_IndexLeftTagGnoHeader.WriteFinal();

	printf("\rWriting left mass index... finished.      \n");

	// copy right tag sorted entries to output file (right mass and GNO)
	r_IndexRightTagMassesHeader lr_IndexRightTagMassesHeader(mk_IndexFile);
	lr_IndexRightTagMassesHeader.WritePreliminary();
	mk_IndexFileInfo.mi_RightTagMassesFilePosition = mk_IndexFile.pos();
	mk_IndexFile.resize(mk_IndexFileInfo.mi_RightTagMassesFilePosition + (qint64)aui_EntryCount * sizeof(unsigned int));
	mk_IndexFile.seek(mk_IndexFileInfo.mi_RightTagMassesFilePosition + (qint64)aui_EntryCount * sizeof(unsigned int));
	r_IndexRightTagGnoHeader lr_IndexRightTagGnoHeader(mk_IndexFile);
	lr_IndexRightTagGnoHeader.WritePreliminary();
	mk_IndexFileInfo.mi_RightTagGnoFilePosition = mk_IndexFile.pos();

	ak_SourceFile.seek(0);

	li_Percent = -1;
	lui_ProcessedTrimerCount = 0;
	for (int li_TrimerIndex = 0; li_TrimerIndex < 8000; ++li_TrimerIndex)
	{
		int li_NewPercent = (int)((double)lui_ProcessedTrimerCount / aui_EntryCount * 100.0);
		if (li_NewPercent != li_Percent)
		{
			printf("\rWriting right mass index... %d%% done.", li_NewPercent);
			li_Percent = li_NewPercent;
		}

		unsigned int lui_EntryCount = mk_IndexFileInfo.mui_pTrimerCountTable.get_Pointer()[li_TrimerIndex];

		if (lui_EntryCount > 0)
		{
			ak_SourceFile.read((char*)lr_pLeftChunk.get_Pointer(), lui_EntryCount * sizeof(r_IndexEntry));

			Sort<r_IndexEntry>(lr_pLeftChunk.get_Pointer(), 0, lui_EntryCount - 1, &lt_IndexEntry_Trimer_RightMass_Gno);

			r_IndexEntry* lr_IndexEntry_ = lr_pLeftChunk.get_Pointer();
			r_IndexEntry* lr_IndexEntryEnd_ = lr_IndexEntry_ + lui_EntryCount;
			unsigned int* lui_Mass_ = lui_pMasses.get_Pointer();
			unsigned int* lui_Gno_ = lui_pGno.get_Pointer();

			while (lr_IndexEntry_ < lr_IndexEntryEnd_)
			{
				*(lui_Mass_++) = lr_IndexEntry_->mui_RightMass;
				*(lui_Gno_++) = lr_IndexEntry_->mui_Gno;
				++lr_IndexEntry_;
			}
			mk_IndexFile.seek(mk_IndexFileInfo.mi_RightTagMassesFilePosition + lui_ProcessedTrimerCount * sizeof(unsigned int));
			mk_IndexFile.write((const char*)lui_pMasses.get_Pointer(), lui_EntryCount * sizeof(unsigned int));
			mk_IndexFile.seek(mk_IndexFileInfo.mi_RightTagGnoFilePosition + lui_ProcessedTrimerCount * sizeof(unsigned int));
			mk_IndexFile.write((const char*)lui_pGno.get_Pointer(), lui_EntryCount * sizeof(unsigned int));
		}
		lui_ProcessedTrimerCount += lui_EntryCount;
	}

	lr_IndexRightTagMassesHeader.WriteFinal(mk_IndexFileInfo.mi_RightTagGnoFilePosition - r_IndexRightTagGnoHeader::getSize());
	lr_IndexRightTagGnoHeader.WriteFinal();

	printf("\rWriting right mass index... finished.      \n");
}


template <typename T> void k_GpfIndexer::Sort(T* ar_Entries_, unsigned int aui_First, unsigned int aui_Last, bool (*af_LessThan_)(T*, T*), T* ar_TempEntries_)
{
	bool lb_AllocatedTempEntries = false;
	if (ar_TempEntries_ == NULL)
	{
		ar_TempEntries_ = new T[aui_Last - aui_First + 1];
		lb_AllocatedTempEntries = true;
	}

	if (aui_Last - aui_First > 1)
	{
		// divide...

		unsigned int lui_Mid = aui_First + (aui_Last - aui_First) / 2;
		Sort<T>(ar_Entries_, aui_First, lui_Mid, af_LessThan_, ar_TempEntries_);
		Sort<T>(ar_Entries_, lui_Mid + 1, aui_Last, af_LessThan_, ar_TempEntries_);

		// ... and conquer!

		// merge step, use ar_Temp_ as the tmeporary desitnation here
		// merge [first .. mid] with [mid + 1 .. last]
		T* lr_Left_ = &(ar_Entries_[aui_First]);
		T* lr_Right_ = &(ar_Entries_[lui_Mid + 1]);
		T* lr_LastLeft_ = &(ar_Entries_[lui_Mid + 1]);
		T* lr_LastRight_ = &(ar_Entries_[aui_Last + 1]);
		T* lr_Result_ = ar_TempEntries_;

		unsigned int lui_Count = aui_Last - aui_First + 1;
		unsigned int lui_Counter = lui_Count;
		// pick smaller element while not at end of left or right list
		while ((lr_Left_ < lr_LastLeft_) && (lr_Right_ < lr_LastRight_) && ((lui_Counter--) != 0))
		{
			if ((*af_LessThan_)(lr_Left_, lr_Right_))
				*(lr_Result_++) = *(lr_Left_++);
			else
				*(lr_Result_++) = *(lr_Right_++);
		}

		if (lui_Counter > 0)
		{
			// there are items remaining in either the left or the right 
			// list, append them to the end of the temporary array
			if (lr_Left_ < lr_LastLeft_)
				memcpy(lr_Result_, lr_Left_, sizeof(T) * lui_Counter);
			else
				memcpy(lr_Result_, lr_Right_, sizeof(T) * lui_Counter);
		}

		// overwrite unsorted entries with freshly sorted entries from the temp array
		memcpy(&ar_Entries_[aui_First], ar_TempEntries_, sizeof(T) * lui_Count);
	}
	else if (aui_First != aui_Last)
	{
		// our task here is to sort an array consisting of two entries
		if ((*af_LessThan_)(&ar_Entries_[aui_Last], &ar_Entries_[aui_First]))
		{
			// swap entries
			T lr_Temp = ar_Entries_[aui_First];
			ar_Entries_[aui_First] = ar_Entries_[aui_Last];
			ar_Entries_[aui_Last] = lr_Temp;
		}
	}

	if (lb_AllocatedTempEntries)
		delete [] ar_TempEntries_;
}


bool get_ValueLessThan(const QMap<unsigned int, unsigned int>::const_iterator& ak_Element1, 
					   const QMap<unsigned int, unsigned int>::const_iterator& ak_Element2)
{
	return ak_Element1.value() < ak_Element2.value();
}


