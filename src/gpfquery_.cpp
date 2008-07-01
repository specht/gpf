#include <stdio.h>
#include <QtCore>
#include "GpfBase.h"
#include "GpfQuery.h"


int main(int ai_ArgumentCount, char **ac_Arguments__) 
{
 	if (ai_ArgumentCount < 3)
	{
		fprintf(stderr, "Usage:   gpfquery [GPF index file] [query] [precursor mass (optional)]\n");
		exit(1);
	}

	QString ls_IndexFile = ac_Arguments__[1];
	QString ls_Peptide = ac_Arguments__[2];

	double ld_PrecursorMass = 0.0;
	if (ai_ArgumentCount > 3)
	{
		bool lb_Ok;
		ld_PrecursorMass = QString(ac_Arguments__[3]).toDouble(&lb_Ok);
		if (!lb_Ok)
			ld_PrecursorMass = 0.0;
	}

	k_GpfBase lk_GpfBase(QStringList() << ls_IndexFile);
	k_GpfQuery lk_Query(lk_GpfBase, ls_Peptide, ld_PrecursorMass);
	lk_Query.Execute();
	printf(lk_Query.get_Result().toStdString().c_str());
}
