TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt


SOURCES += main.cpp \
    atom.cpp \
    lib.cpp \
    generatequantities.cpp \
    cell.cpp \
    potentials.cpp

HEADERS += \
    atom.h \
    lib.h \
    generatequantities.h \
    cell.h \
    potentials.h
