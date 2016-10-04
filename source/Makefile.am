AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = averager
bin_SCRIPTS = averager.so

noinst_LIBRARIES = libaverage.a 
libaverage_a_SOURCES = fillarrays.f toblockdiag.f error_logging.f output.f statrecalc.f \
					 common_tools.f averaging.f initave.f readdata.f 

# all source files which go into the main target
averager_SOURCES = fillarrays.f  averaging.f  averager.f \
                      readdata.f  statrecalc.f  toblockdiag.f \
                      initave.f output.f error_logging.f covartonui.f common_tools.f

# Force recompilation if common.inc is modified:
*.o: ../include/common.inc ../include/data.inc  ../include/settings.inc

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = -I$(srcdir)/../include -fno-automatic -fPIC
averager.so:
# 
AM_CPPFLAGS = -I$(srcdir)/../include 


# libraries to be linked to the main target
LDADD =  @FLIBS@ -L../num_utils/ -lnum_utils  $(LAPACK_LIBS) $(BLAS_LIBS)
# \ -lgsl -lgslcblas -lstdc++

GFRTLLIB = $(shell gfortran  -print-file-name=libgfortran.a )
GFRTLIB  = -L $(subst /libgfortran.a, ,$(GFRTLLIB)) -lgfortran

if ENABLE_PYTHON
  clean-averager: 
	rm  averager.so
  all: all-am averager.so
  clean: clean-am clean-average
  averager.so: average_py.f90  libaverage.a ../num_utils/libnum_utils.a
	f2py --quiet -c ../num_utils/libnum_utils.a libaverage.a average_py.f90 ../num_utils/libnum_utils.a   $(LAPACK_LIBS) $(BLAS_LIBS) -I../include -m averager --fcompiler=$(F77)
endif

