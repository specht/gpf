include(../gpfapp.pro)

TARGET = gpfquery
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

SOURCES += \
	../../src/gpfquery_.cpp \
	