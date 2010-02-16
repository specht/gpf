include(../base.pro)

TARGET = gpfquery
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

HEADERS += \
	../../src/GpfIndexFile.h \
	../../src/GpfQuery.h \
	
SOURCES += \
	../../src/gpfquery.cpp \
	../../src/GpfIndexFile.cpp \
	../../src/GpfQuery.cpp \
