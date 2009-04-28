include(../base.pro)

TARGET = gpfindex
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

HEADERS += \
	../../src/GpfIndexer.h \
	
SOURCES += \
	../../src/gpfindex.cpp \
	../../src/GpfIndexer.cpp \
