include(../base.pro)

TARGET = gpfcheck
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

HEADERS += \
	../../src/GpfIndexFile.h \
	
SOURCES += \
	../../src/gpfcheck.cpp \
	../../src/GpfIndexFile.cpp \
