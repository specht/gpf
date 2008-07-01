#pragma once

#include <QString>
#include "GpfQuery.h"
#include "RefPtr.h"


struct r_ChainMarker
{
	enum Enumeration 
	{
		Left = 1,
		Right = 2
	};
};


class k_Hit
{
public:
	k_Hit(k_GpfQuery& ak_Query, bool ab_Forward, QList<QPair<unsigned int, unsigned int> > ak_Assembly);
	k_Hit(const k_Hit& ak_Other);
	virtual ~k_Hit();

	virtual QString get_Peptide() const;
	virtual QString get_Assembly() const;
	virtual unsigned int get_Mass() const;
	virtual void AddInformation(QString as_Key, QString as_Value);

	virtual void Finish();

	bool get_HasFullScore() const;
	bool get_IsDiscarded() const;
	void MarkDiscarded();

	virtual QString description();

protected:

	virtual unsigned int CalculateScore() const;

	k_GpfQuery& mk_Query;
	bool mb_Forward;
	QList<QPair<unsigned int, unsigned int> > mk_Assembly;
	unsigned mui_QueryLength;

	QList<QPair<QString, QString> > mk_ResultItems;

	QString ms_Peptide;
	QString ms_Left;
	QString ms_Right;
	QString ms_Assembly;

	// derived values that are calculated on Finish().
	unsigned int mui_MaxChainLength;
	unsigned int mui_HitMass;
	bool mb_IsDiscarded;

private:
	void Initialize();
};
