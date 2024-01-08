## Makefile for C++ project using Boost
#
# @author Cedric "levif" Le Dillau
# @modified Aroli Marcellinus 
#
# Some notes:
# - Using ':=' instead of '=' assign the value at Makefile parsing time,
#   others are evaluated at usage time. This discards
# - Use ':set list' in Vi/Vim to show tabs (Ctrl-v-i force tab insertion)
#

# List to '.PHONY' all fake targets, those that are neither files nor folders.
# "all" and "clean" are good candidates.
.PHONY: all, clean

# Define the final program name
PROGNAME := cell_CRN1998

# Pre-processor flags to be used for includes (-I) and defines (-D) 
CPPFLAGS := -I/opt/Lib/dev/petsc-new/include

# CFLAGS is used for C compilation options.
CFLAGS := -Wall -O0

# CXX to set the compiler
CXX := icpc

# CXXFLAGS is used for C++ compilation options.
#CXXFLAGS += -Wall -O0 -fpermissive -std=c++11
CXXFLAGS += -Wall -O2 -fno-alias -fpermissive -DSINGLE_CELL -DDEBUG_COMPLETE

# LDFLAGS is used for linker (-g enables debug symbols)
LDFLAGS  += -g -static /opt/Lib/dev/petsc-new/lib/libsundials_cvode.a /opt/Lib/dev/petsc-new/lib/libsundials_nvecserial.a

# Which Boost modules to use (all)
BOOST_MODULES := \
  filesystem \
  iostreams \
  regex \
  system

# Boost libraries' type (a suffix)
BOOST_MODULES_TYPE := -mt

# Define library names with their type
BOOST_MODULES_LIBS := $(addsuffix $(BOOT_MODULES_TYPE),$(BOOST_MODULES))

# Define the linker argument to use the Boost libraries.
BOOST_LDFLAGS := $(addprefix -lboost_,$(BOOST_MODULES_LIBS))

# Feed compiler/linker flags with Boost's
#CPPFLAGS += $(BOOST_CPPFLAGS)
#LDFLAGS += $(BOOST_LDFLAGS)

# List the project' sources to compile or let the Makefile recognize
# them for you using 'wildcard' function.
#
#SOURCES = simple_ls.cpp rawr.cpp converter.cpp
#SOURCES := $(wildcard *.cpp) \
	$(wildcard **/*.cpp)
SOURCES	= $(wildcard **/*.cpp) main.cpp

# List the project' headers or let the Makefile recognize
# them for you using 'wildcard' function.
#
#HEADERS = simple_ls.h 2dquicksort.h rawr.h
HEADERS	=	$(wildcard **/*.hpp)

# Construct the list of object files based on source files using
# simple extension substitution.
OBJECTS := $(SOURCES:%.cpp=%.o)

#
# Now declare the dependencies rules and targets
#
# Starting with 'all' make it  becomes the default target when none 
# is specified on 'make' command line.
all : $(PROGNAME)

# Declare that the final program depends on all objects and the Makfile
$(PROGNAME) : $(OBJECTS) Makefile
	$(CXX) -o bin/$@ $(OBJECTS) $(LDFLAGS)

# Now the choice of using implicit rules or not (my choice)...
#
# Choice 1: use implicit rules and then we only need to add some dependencies
#           to each object.
#
## Tells make that each object file depends on all headers and this Makefile.
#$(OBJECTS) : $(HEADERS) Makefile
#
# Choice 2: don't use implicit rules and specify our will
%.o: %.cpp $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(OUTPUT_OPTION) $<


# Simple clean-up target
# notes:
# - the '@' before 'echo' informs make to hide command invocation.
# - the '-' before 'rm' command to informs make to ignore errors.
clean :
	@echo "Clean."
	-rm -rf *.o **/*.o bin/$(PROGNAME)
