TARGET = g4blgui
COMPANY = muonsinc

SOURCES += *.cc  ../g4bl/Util.cc
HEADERS += *.h

linux-g++:QMAKE_LFLAGS_RPATH=
linux-g++:QMAKE_LFLAGS += "-Wl,-rpath,\'\$$ORIGIN\'/../lib"

INCLUDEPATH += . ../g4bl
TEMPLATE = app
QMAKE_CXXFLAGS_WARN_ON = -Wall -Wno-unused-parameter

