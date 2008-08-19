include(../gpfapp.pro)

TARGET = gpfd
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

QT += network

HEADERS += \
	../../src/GpfDaemon.h \
	../../src/GpfDaemonThread.h \
	../../src/GpfWorkerThread.h \
#	../../src/QueryBatch.h \
	../../src/QueryQueue.h \
	../../src/ValueEstimator.h \

SOURCES += \
	../../src/gpfd.cpp \
	../../src/GpfDaemon.cpp \
	../../src/GpfDaemonThread.cpp \
	../../src/GpfWorkerThread.cpp \
#	../../src/QueryBatch.cpp \
	../../src/QueryQueue.cpp \
	../../src/ValueEstimator.cpp \
