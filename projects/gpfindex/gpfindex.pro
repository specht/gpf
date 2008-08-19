include(../gpfapp.pro)

TARGET = gpfindex
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

SOURCES += \
	../../src/gpfindex.cpp \
	