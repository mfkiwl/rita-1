## Process this file with automake to produce Makefile.in

bin_PROGRAMS = rita
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y
rita_LDADD = $(top_builddir)/muparserx/libmp.a

rita_SOURCES = rita.cpp \
               rita.h \
               ritaException.h \
               approximation.cpp \
               approximation.h \
               calc.cpp \
               calc.h \
               cmd.cpp \
               cmd.h \
               configure.cpp \
               configure.h \
               data.cpp \
               data.h \
               eigen.cpp \
               eigen.h \
               equa.cpp \
               equa.h \
               help.h \
               integration.cpp \
               integration.h \
               mesh.cpp \
               mesh.h \
               optim.cpp \
               optim.h \
               runAE.cpp \
               runPDE.cpp \
               runODE.cpp \
               solve.cpp \
               solve.h \
               stationary.cpp \
               stationary.h \
               transient.cpp \
               transient.h

clean-local:
	-rm -f .rita.log .rita.his
