include(../base.pro)

TARGET = gpfindex
CONFIG(debug, debug|release) {
    TARGET = $$join(TARGET,,,_debug)
}

HEADERS += \
    ../../src/GpfIndexer.h \
    ../../src/HmstIterator.h \
    
SOURCES += \
    ../../src/gpfindex.cpp \
    ../../src/GpfIndexer.cpp \
    ../../src/HmstIterator.cpp \
