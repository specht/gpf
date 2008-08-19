include(../gpfapp.pro)

TARGET = gpfdump
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

SOURCES += \
	../../src/gpfdump.cpp \
	