include(../gpfapp.pro)

TARGET = gpfbatch
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

SOURCES += \
	../../src/gpfbatch.cpp \
	