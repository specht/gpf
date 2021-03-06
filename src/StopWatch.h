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

#ifdef _WIN32
#include <windows.h>
#endif


class k_StopWatch
{
public:
    k_StopWatch();
    k_StopWatch(QString as_Message, QTextStream* ak_OutputStream_ = NULL);
    virtual ~k_StopWatch();

    static QString getTimeAsString(double ad_Time);

    QString getTimeAsString();
    double get_Time();
    void reset();
    void setExitMessage(QString as_Message);

private:

    double get_AbsoluteTime();

#ifdef _WIN32
    LARGE_INTEGER ml_Frequency;
#endif

    double md_StartTime;
    QString ms_Message;
    QTextStream mk_StdOutputStream;
    QTextStream* mk_OutputStream_;
    bool mb_Print;
};
