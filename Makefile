include make.inc

# CC := clang --analyze # and comment out the linker last line for sanity
GITREV = -D'GIT_REVISION="$(shell git log --pretty=format:'%h' -n 1)"'
SRCDIR := src
BUILDDIR := build

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

INC := -I ./ -I ./include -I ./include/libgeodecomp -I $(GEODECOMP_DIR)/include -I $(BOOST_INCLUDE_DIR)
CFLAGS += -DLIBGEODECOMP_WITH_MPI -Wfatal-errors -fopenmp -std=c++11 $(GITREV)
LDFLAGS := -fopenmp -L $(GEODECOMP_DIR)/lib -L $(BOOST_LIB_DIR) 
LIBS := -lgeodecomp -lboost_date_time 

TYPEMAP_TEST_OBJECTS := src/catchmentmodel/LSDCatchmentModel.o src/libgeodecomp/typemaps.o test/typemaptest.o

TARGET := bin/HAIL-CAESAR.mpi

ifdef $(MPI_DIR)
INC += -I $(MPI_DIR)/include 
LDFLAGS += -L $(MPI_DIR)/lib
LIBS += -lmpi
endif


$(TARGET): $(OBJECTS) 
	@echo -e " \n Linking... \n"
	@echo " $(CXX) $(LDFLAGS) $(INC) $(LIBS) $^ -o $(TARGET)"; $(CXX) $(LDFLAGS) $(INC) $(LIBS) $^ -o $(TARGET) 

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) 
	@mkdir -p bin
	@mkdir -p $(BUILDDIR)/topotools
	@mkdir -p $(BUILDDIR)/catchmentmodel
	@mkdir -p $(BUILDDIR)/libgeodecomp
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

typemaps: # only necessary to generate new MPI typemaps run if the model has been changed, not needed during normal build process
	@echo " Generating xml using doxygen..."; echo "doxygen ./include/libgeodecomp/doxygen.conf"; doxygen ./include/libgeodecomp/doxygen.conf
	@echo " Generating MPI typemaps...";
	@mkdir typemaps
	@$(GEODECOMP_SRC_DIR)/tools/typemapgenerator/generate.rb -S typemaps-doxygen-docs/xml typemaps
	@sed -i '' 's/CellType/int/g' typemaps/typemaps.cpp  # correct mistake in LibGeoDecomp typemap generation script for enums, see https://github.com/gentryx/libgeodecomp/issues/72
	@cp typemaps/*.h include/libgeodecomp/
	@cp typemaps/*.cpp src/libgeodecomp/
	@echo " Typemaps saved in src/libgeodecomp and include/libgeodecomp"
	@rm -rf typemaps
	@rm -rf typemaps-doxygen-docs

typemaptest: $(TYPEMAP_TEST_OBJECTS) $(OBJECTS)
	@echo -e " \n Linking... \n"
	@echo " $(CXX) $(LDFLAGS) $(INC) $(LIBS) $(TYPEMAP_SRC_OBJECTS) -o bin/typemaptest"; $(CXX) $(LDFLAGS) $(INC) $(LIBS) $(TYPEMAP_TEST_OBJECTS) -o bin/typemaptest

typemaptest.o : test/typemaptest.cpp include/catchmentmodel/LSDCatchmentModel.hpp include/libgeodecomp/typemaps.h
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o test/typemaptest.o test/typemaptest.cpp"; $(CXX) $(CFLAGS) $(INC) -c -o test/typemaptest.o test/typemaptest.cpp

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -rf $(BUILDDIR) $(TARGET) typemaps typemaps-doxygen-docs"; $(RM) -r $(BUILDDIR) $(TARGET) typemaps typemaps-doxygen-docs



.PHONY: clean
