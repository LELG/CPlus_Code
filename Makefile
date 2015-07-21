# Compiler and Linker
CXX       := g++

# Ensure we have g++ 4.7 or greater
# This assumes that the next version of gcc will NOT be 4.10, but 4.9x, or 5.x
# See Confluence for an explanation of how this command works
GCC_GTE47 := $(shell awk -v gcc_ver=`gcc -dumpversion | cut -f1-2 -d.` 'BEGIN {if (gcc_ver >= 4.7) {print "1"} else {print "0"}}')

ifeq ($(GCC_GTE47), 0)
    $(error g++ version too old, >= 4.7 required)
endif

# Target executable
TARGET    := tumourEvoSim

# Directories
BUILDDIR  := obj
TARGETDIR := bin

# Flags, libraries, and includes
CPPFLAGS  := -Wall -std=c++11
LIB       := -L/usr/local/lib -fopenmp -lgsl -lgslcblas -lm
INC       := -I/usr/local/include 

# Specify object files in OBJECTS; OBJ just adds the directory prefix
OBJECTS   := V1_Main_Tumour_Evoulution.o
OBJ       := $(patsubst %,$(BUILDDIR)/%,$(OBJECTS))

# -----------------------------------------------------------------------------
# TARGETS
# -----------------------------------------------------------------------------

# Default make
all: directories $(TARGET)

# Create directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

# Remake
remake: clean all

# Clean up objects and executables
clean:
	@rm -rf $(BUILDDIR)
	@rm -rf $(TARGETDIR)

# Link
$(TARGET): $(OBJ)
	$(CXX) $(LIB) -o $(TARGETDIR)/$(TARGET) $^

# Compile
$(BUILDDIR)/%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(INC) -c -o $@ $<
