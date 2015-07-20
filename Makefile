# Compiler and Linker
CXX       := g++

# Ensure we have g++-4.7 or greater
GCC_GTE47 := $(shell echo `g++ -dumpversion | cut -f1-2 -d.` \>= 4.7 | sed -e 's/\./*100+/g' | bc )

ifeq ($(GCC_GTE47), 0)
    $(error g++ version too old. >= 4.7 required)
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
