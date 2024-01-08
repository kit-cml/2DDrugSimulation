BLDIR=/opt/prog/petsc
PETSC_DIR=$(BLDIR)/petsc-3.1-p5
PETSC_ARCH=arch-intel-opt
PETSC_LIB= $(PETSC_DIR)/lib
#include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/conf/variables
#PETSC_LIB = -Wl,-rpath,/opt/Lib/dev/petsc/lib -L/opt/Lib/dev/petsc/lib -lpetsc
#PETSC_LIB=-Wl,-rpath,/opt/Lib/dev/petsc/lib -L/opt/Lib/dev/petsc/lib -lpetsc -Wl,-rpath,/opt/Lib/dev/petsc/lib -L/opt/Lib/dev/petsc/lib -lf2clapack -lf2cblas -Wl,-rpath,/usr/local/lib -Wl,-rpath,/opt/Lib/dev/openmpi/lib -ldl /usr/local/lib64/libsundials_nvecserial.a /usr/local/lib64/libsundials_cvode.a 
LEX      = flex
YACC     = yacc
#CXXFLAGS = -fpermissive -std=c++11
CXXFLAGS = -fpermissive -DTISSUE -DVTK_COMBINED
#uncomment this if you want to use Tomek model
#CXXFLAGS += -DTOMEK_2019
LEXFLAGS =
YACCFLAGS= -d
INCPATH  = -I./ -I/ $(PETSC_DIR)/include 
LFLAGS   =  -pg
LIBS	=  $(PETSC_LIB) -lpetsc -lz -lm


AR       = ar cqs
RANLIB   =
TAR      = tar -cf
GZIP     = gzip -9f
COPY     = cp -f
COPY_FILE= $(COPY)
COPY_DIR = $(COPY) -r
INSTALL_FILE= $(COPY_FILE)
INSTALL_DIR = $(COPY_DIR)
DEL_FILE = rm -f
SYMLINK  = ln -sf
DEL_DIR  = rmdir
MOVE     = mv -f
CHK_DIR_EXISTS= test -d
MKDIR    = mkdir -p

####### Output directory

OBJECTS_DIR = ./

####### Files

HEADERS := fem.h \
	  heart.h \
	  mesh.h \
	fiber.h \
	bspm.h \
    ccalcexception.h \
    DrugSimulation/modules/glob_type.hpp \
    DrugSimulation/modules/cipa_t.hpp \
    DrugSimulation/modules/commons.hpp \
    DrugSimulation/libs/zip.h \
    libs/ode.hpp \
    $(wildcard DrugSimulation/cellmodels/*.hpp)


SOURCES := bench3d.cpp \
	  fem.cpp \
	  heart.cpp \
	  mesh.cpp \
	  fiber.cpp \
    bspm.cpp \
    ccalcexception.cpp \
    DrugSimulation/modules/glob_type.cpp \
    DrugSimulation/modules/cipa_t.cpp \
    DrugSimulation/modules/commons.cpp \
    DrugSimulation/libs/zip.c \
    libs/ode.cpp \
    $(wildcard DrugSimulation/cellmodels/*.cpp)

OBJECTS := $(SOURCES:%.cpp=%.o)

FORMS =
UICDECLS =
UICIMPLS =
SRCMOC   =
OBJMOC =
DESTDIR  = 
#TARGET   = ./bin/CEP6_A1656D 
TARGET   = ./bin/CEP6_2D_Drug
#TARGET   = ./bin/CEP6_3D
#TARGET   = ./bin/CEP4_bspm

first: all
####### Implicit rules

.SUFFIXES: .c .o .cpp .cc .cxx .C

.cpp.o:
	$(PETSC_COMPILE_SINGLE) $(CXXFLAGS) $(INCPATH) -o $@ $<

.cc.o:
	$(PETSC_COMPILE_SINGLE) $(CXXFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(PETSC_COMPILE_SINGLE) $(CXXFLAGS) $(INCPATH) -o $@ $<

.C.o:
	$(PETSC_COMPILE_SINGLE) $(CXXFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(PETSC_COMPILE_SINGLE) $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)
	cep2017_heart -d bin/ || mkdir -p ./bin
	$(info PCC="$(PCC)"  )
	$(info PETSC_LIB="$(PETSC_LIB)")
	$(PCC) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

dist:
	@mkdir -p .tmp/two && $(COPY_FILE) --parents $(SOURCES) $(HEADERS) $(FORMS) $(DIST) .tmp/two/ && ( cd `dirname .tmp/two` && $(TAR) two.tar two && $(GZIP) two.tar ) && $\(MOVE) `dirname .tmp/two`/two.tar.gz . && $(DEL_FILE) -r .tmp/two


clean:
	-$(DEL_FILE) -f *.o
	-$(DEL_FILE) -f CEPCell_NEO/*.o
	-$(DEL_FILE) -f CEPCell_NEO/**/*.o
	-$(DEL_FILE) *~ core *.core
	-$(DEL_FILE) -f *.out
	-$(DEL_FILE) $(TARGET)
