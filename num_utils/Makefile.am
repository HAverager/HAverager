AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libnum_utils.a

libnum_utils_a_SOURCES = deqinv.F  deqn.F  dfact.F  dfeqn.F  dfinv.F  dinv.F  ucopy.F mydsyevd.f

AM_FFLAGS = -I$(srcdir)/../include -ffixed-line-length-132   -fPIC
#AM_FFLAGS = -I$(srcdir)/../include -fno-automatic -finit-local-zero -ffixed-line-length-132

#if ENABLE_OPENMP
#AM_FFLAGS+= -fopenmp
#endif

