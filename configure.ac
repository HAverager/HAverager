# Process this file with autoconf to produce a configure script.
AC_INIT([averager], [1.0], [glazov@mail.desy.de], [averager])
m4_ifndef([AC_PACKAGE_URL], AC_DEFINE([PACKAGE_URL],["http://xfitter.org"]))

#Suppress verbose output when compiling (use make V=99 for verbose output)
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AM_INIT_AUTOMAKE
AC_PREFIX_DEFAULT(${PWD})

# Checks for programs.
AC_PROG_FC([gfortran])
AC_PROG_F77([gfortran])
AC_PROG_CPP
# AC_PROG_RANLIB
AC_PROG_LIBTOOL

# Fortran configures
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS
AC_F77_WRAPPERS
#AC_FC_FREEFORM


#Use libtoolize to build libraries
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

# lapack macros
m4_include([m4/ax_blas.m4])
m4_include([m4/ax_lapack.m4])




# Check for trapFPE flag
AC_ARG_ENABLE([trapFPE],
  [AC_HELP_STRING(--enable-trapFPE, [Stop of floating point errors (default=no)])])
if test x$enable_trapFPE == xyes; then
  echo "Add floating point trap"
  FFLAGS="$FFLAGS -ffpe-trap='invalid,overflow' "
fi

# Check for check bounds flag

AC_ARG_ENABLE([checkBounds],
  [AC_HELP_STRING(--enable-checkBounds, [add -fbounds-check flag for compilation (default=no)])])
if test x$enable_checkBounds == xyes; then
  echo "Add check bounds flag"
  FFLAGS="$FFLAGS -fbounds-check"
fi

AC_ARG_ENABLE([python],
  [AC_HELP_STRING(--enable-python, [add python support (default=no)])])
if test x$enable_python == xyes; then
        AC_MSG_CHECKING([F2PY using: f2py])
	f2py -v  >/dev/null 2>&1

	# Check if f2py3 is called using "f2py"
	if test "$?" != "0"; then
		AC_MSG_RESULT([no])
		HAVE_F2PY=no
   		AM_CONDITIONAL(ENABLE_PYTHON, [false]) 
	else
		HAVE_F2PY=yes
		F2PY=`which f2py`
		F2PYV=`$F2PY -v`
		AC_DEFINE([HAVE_F2PY], [1])
		AC_DEFINE_UNQUOTED([F2PY], ["${F2PY}"], [Path to F2PY])
		AC_SUBST([F2PY], ["${F2PY}"])
		AC_MSG_RESULT([yes v${F2PYV}])
		AM_CONDITIONAL(ENABLE_PYTHON, [true]) 
	fi
  	echo "Add python support"
else
   		AM_CONDITIONAL(ENABLE_PYTHON, [false]) 
fi



# check for lapack and blas
AX_BLAS
if test "$ax_blas_ok" = "no"; then
   AC_MSG_ERROR([Cannot find BLAS libraries])
fi

AX_LAPACK
if test "$ax_lapack_ok" = "no"; then
   AC_MSG_ERROR([Cannot find LAPACK libraries])
fi



# DOXYGEN SUPPORT

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN([$PACKAGE_NAME],[doxygen.cfg])



# Output
AC_CONFIG_FILES([source/Makefile
                 num_utils/Makefile
                 average-config
                 Makefile])

AC_OUTPUT
