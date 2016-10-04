ACLOCAL_AMFLAGS=-I m4

AUTOMAKE_OPTIONS = foreign
SUBDIRS =  num_utils source

# DOXYGEN SUPPORT
include $(top_srcdir)/aminclude.am

bin_SCRIPTS = average-config

# ensure the distribution of the doxygen configuration file
EXTRA_DIST = doxygen.cfg
