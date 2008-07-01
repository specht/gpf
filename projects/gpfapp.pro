TEMPLATE = app
win32 {
    CONFIG += embed_manifest_exe
}

macx {
	CONFIG -= app_bundle
	CONFIG += x86 ppc
}

CONFIG += debug_and_release

OBJECTS_DIR = ../../obj/
MOC_DIR = ../../obj/
RCC_DIR = ../../obj/
DESTDIR = ../../

QT -= gui

CONFIG += console

RESOURCES += ../../src/libgpf.qrc

HEADERS += \
	../../src/GpfBase.h \
	../../src/GpfBrowser.h \
	../../src/GpfConfig.h \
	../../src/GpfIndexer.h \
	../../src/GpfIndexFileInfo.h \
	../../src/GpfIndexHeader.h \
	../../src/GpfQuery.h \
	../../src/Hit.h \
	../../src/IndexReader.h \
	../../src/RefPtr.h \
	../../src/Sorter.h \
	../../src/StopWatch.h \

SOURCES += \
	../../src/GpfBase.cpp \
	../../src/GpfBrowser.cpp \
	../../src/GpfIndexer.cpp \
	../../src/GpfIndexFileInfo.cpp \
	../../src/GpfQuery.cpp \
	../../src/Hit.cpp \
	../../src/Sorter.c \
	../../src/StopWatch.cpp \
