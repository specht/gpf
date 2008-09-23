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


#define GNO_BACKWARD 0x80000000
#define GNO_BACKWARD_INVERSE 0x7fffffff

static const unsigned int gui_GlobalNucleotideOffsetBackwardFlag = 0x80000000;
static const double gd_WaterMass = 18.0;

static const unsigned int gui_MassPrecision = 10000;
static const unsigned int gui_MaxUnknownAminoAcidCount = 1;
static const unsigned int gui_MaxPeptideLength = 1000;

static const unsigned int gui_MaxNucleotideBatchSize = 64 * 1024; // 64 K
static const unsigned int gui_SortBatchSize = 8 * 1024 * 1024; // 8 M (multiplied by 14 == 112 MB for one table, and we need two!)
static const unsigned int gui_InvertTripletMask = 0xdb;

static const unsigned int gui_MaxSurroundingAminoAcidCount = 5;
