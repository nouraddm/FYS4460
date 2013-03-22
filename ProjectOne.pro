TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt


SOURCES += main.cpp \
    atom.cpp \
    lib.cpp \
    generatequantities.cpp \
    cell.cpp \
    potentials.cpp \
    pressurecells.cpp

HEADERS += \
    atom.h \
    lib.h \
    generatequantities.h \
    cell.h \
    potentials.h \
    pressurecells.h

# MPI Settings
QMAKE_CXX = /opt/local/bin/mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = /opt/local/bin/mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
