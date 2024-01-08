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
#PROGNAME := cell_sim_3
PROGNAME := cell_sim_qnet_last
#PROGNAME := cell_sim

# Pre-processor flags to be used for includes (-I) and defines (-D) 
CPPFLAGS := -I/opt/Lib/dev/sundials/include -DSINGLE_CELL

# CXX to set the compiler
#CXX := mpicxx
CXX := /opt/Lib/dev/openmpi/bin/mpicxx


# CXXFLAGS is used for C++ compilation options.
#CXXFLAGS += -Wall -O0 -fpermissive -std=c++11
#CXXFLAGS += -Wall -O2 -fno-alias -fpermissive 

# LDFLAGS is used for linker (-g enables debug symbols)
LDFLAGS  += -g -lmagic /opt/Lib/dev/sundials/lib/libsundials_cvode.a /opt/Lib/dev/sundials/lib/libsundials_nvecserial.a 

# List the project' sources to compile or let the Makefile recognize
# them for you using 'wildcard' function.
#
SOURCES	= $(wildcard **/*.cpp) $(wildcard **/*.c) main.cpp

# List the project' headers or let the Makefile recognize
# them for you using 'wildcard' function.
#
HEADERS	= $(wildcard **/*.hpp) $(wildcard **/*.h)

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
