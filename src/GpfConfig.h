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
