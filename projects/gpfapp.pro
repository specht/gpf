TEMPLATE = app

macx {
	CONFIG -= app_bundle
	CONFIG += x86 ppc
	QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.4
}

CONFIG += debug_and_release

CONFIG(debug, debug|release) {
	OBJECTS_DIR = ../../obj/debug/
	MOC_DIR = ../../obj/debug/
	RCC_DIR = ../../obj/debug/
}
else {
	OBJECTS_DIR = ../../obj/release/
	MOC_DIR = ../../obj/release/
	RCC_DIR = ../../obj/release/
}

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
