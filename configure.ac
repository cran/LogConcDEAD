AC_INIT([LogConcDEAD], 1.0-9) dnl package name, version

dnl find what R uses for the C compiler and initialize flags
: ${R_HOME=`R RHOME`}
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=""

dnl check that it has the flags we're looking for
for option in -fno-strict-aliasing; do
	SAVE_CFLAGS="$CFLAGS"
	CFLAGS="$CFLAGS $option"
	AC_TRY_COMPILE([], [],
		has_option=yes,
		has_option=no,)
	if test $has_option = no; then
		CFLAGS="$SAVE_CFLAGS"
	fi
	AC_MSG_CHECKING([whether $CC understands $option])
	AC_MSG_RESULT($has_option)
	unset has_option
	unset SAVE_CFLAGS
done
unset option

dnl write out the Makevars file with the right CFLAGS
AC_CONFIG_FILES([src/Makevars])

dnl done
AC_OUTPUT
