TEMPLATE = app
    
macx {
    CONFIG -= app_bundle
    CONFIG += x86 ppc
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.4
}

; linux-g++ {
;     CONFIG += static
;     QMAKE_LFLAGS += -L/usr/local/lib
;     LIBS += -ldl -lrt
;     message("STATIC LINKAGE")
; }

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
    ../../src/BitWriter.h \
    ../../src/GpfBase.h \
    ../../src/StopWatch.h \

SOURCES += \
    ../../src/BitWriter.cpp \
    ../../src/GpfBase.cpp \
    ../../src/StopWatch.cpp \
