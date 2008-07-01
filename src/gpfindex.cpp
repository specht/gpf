#include "GpfIndexer.h"


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
	if (ai_ArgumentCount < 4)
	{
		printf("Usage: gpfindex [genome fasta file] [GPF index out file] [genome title]\n");
		exit(1);
	}

	// construct index filename
	QString ls_GenomeFilename = ac_Arguments__[1];
	QString ls_IndexFilename = ac_Arguments__[2];
	QString ls_GenomeTitle = ac_Arguments__[3];

	k_GpfIndexer lk_GpfIndexer(ls_GenomeFilename, ls_IndexFilename, ls_GenomeTitle);
	lk_GpfIndexer.CompileIndex();
}
