TARGET = g4bldata
COMPANY = muonsinc

SOURCES += *.cc  ../g4bl/Util.cc
HEADERS += *.h

# libarchive
INCLUDEPATH += $$(PWD)/libarchive-3.1.2/include
macx:LIBS += $$(PWD)/libarchive-3.1.2/lib/libarchive.dylib

linux-g++:QMAKE_LFLAGS_RPATH=
linux-g++:QMAKE_LFLAGS += "-Wl,-rpath,\'\$$ORIGIN\'/../lib"

INCLUDEPATH += .
TEMPLATE = app
QT += network
QMAKE_CXXFLAGS_WARN_ON = -Wall -Wno-unused-parameter

