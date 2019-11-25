TEMPLATE = lib
QMAKE_CXXFLAGS += -std=c++0x
QMAKE_RPATHDIR += /usr/autodesk/maya2018/lib
CONFIG += lib
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    implicitskin.cpp \
    hrbfgenerator.cpp \
    hrbfmanager.cpp

HEADERS += \
    implicitskin.h \
    hrbfgenerator.h \
    hrbfmanager.h

INCLUDEPATH += /usr/autodesk/maya2018/include/

LIBS += -L/usr/autodesk/maya2018/lib -lOpenMaya -lFoundation -lIMFbase -lOpenMayaAnim
