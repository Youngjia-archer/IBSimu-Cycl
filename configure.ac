AC_INIT([libibsimu-1.0.6dev], [1.0.6dev], [taneli.kalvas@jyu.fi])
AC_CANONICAL_HOST
#AM_CONFIG_HEADER(src/config.h)
AC_CONFIG_HEADERS(src/config.h)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
AC_LANG(C++)

dnl AC_USE_SYSTEM_EXTENSIONS

AC_DEFUN(
	AC_GCC_SUPPORTS_RDYNAMIC,
	[
		AC_PATH_PROG([grep], [grep])
                echo -n "checking if the compiler supports -rdynamic... "
                if test "X${ac_compiler_gnu}X" = "XyesX"; then
                        if test "X`${CC} -dumpspecs | ${grep} rdynamic`X" != "XX"; then
                                echo "yes"
                                AC_SUBST([rdynamic_flag], ["-rdynamic"])
                        else
                                echo "no"
                                AC_SUBST([rdynamic_flag], [])
                        fi
                else
                        echo "no"
                        AC_SUBST([rdynamic_flag], [])
                fi
        ]
)
AC_GCC_SUPPORTS_RDYNAMIC([rdynamic_flag])

lt_enable_auto_import=""
case "$host_os" in
   mingw* | cegcc*)
      lt_enable_auto_import="-Wl,--enable-auto-import"
esac
AC_SUBST(lt_enable_auto_import)

dnl ---------------------------------
dnl Check for programs
dnl ---------------------------------

dnl m4_pattern_allow([AC_PROG_LIBTOOL])
dnl AC_PROG_LIBTOOL
LT_INIT
AC_PROG_CXX
dnl AC_PROG_RANLIB

dnl ---------------------------------
dnl OpenMP (IBSimu-Cycl: particle-level parallelism)
dnl ---------------------------------
AC_OPENMP

dnl ---------------------------------
dnl CUDA GPU backend (IBSimu-Cycl: optional, NVRTC + Driver API)
dnl
dnl   --without-cuda    禁用
dnl   --with-cuda=DIR   指定 CUDA 前缀（默认自动探测 /opt/cuda, /usr/local/cuda, /usr）
dnl
dnl   构建期只需要 CUDA 头文件与 libnvrtc/libcuda，**不需要 nvcc**：
dnl   内核源码以字符串内嵌，由 NVRTC 在运行期编译（见 src/cyclotron/cuda/）。
dnl ---------------------------------
AC_ARG_WITH([cuda],
    [AS_HELP_STRING([--with-cuda@<:@=DIR@:>@],
                    [CUDA toolkit prefix for the optional GPU backend (auto-detected)])],
    [], [with_cuda=auto])

have_cuda=no
CUDA_CPPFLAGS=""
CUDA_LDFLAGS=""
CUDA_LIBS=""

if test "x$with_cuda" != "xno"; then
    if test "x$with_cuda" = "xauto" || test "x$with_cuda" = "xyes"; then
        cuda_candidates="/opt/cuda /usr/local/cuda /usr"
    else
        cuda_candidates="$with_cuda"
    fi

    for cuda_p in $cuda_candidates; do
        for cuda_inc in "$cuda_p/include" "$cuda_p/targets/x86_64-linux/include" \
                        "$cuda_p/targets/aarch64-linux/include"; do
            if test -f "$cuda_inc/nvrtc.h" && test -f "$cuda_inc/cuda.h"; then
                for cuda_libd in "$cuda_p/lib64" "$cuda_p/lib" \
                                 "$cuda_p/targets/x86_64-linux/lib" \
                                 "$cuda_p/targets/aarch64-linux/lib"; do
                    if test -f "$cuda_libd/libnvrtc.so" || test -f "$cuda_libd/libnvrtc.a"; then
                        have_cuda=yes
                        CUDA_CPPFLAGS="-I$cuda_inc"
                        CUDA_LDFLAGS="-L$cuda_libd"
                        CUDA_LIBS="-lnvrtc -lcuda"
                        break
                    fi
                done
            fi
            if test "x$have_cuda" = "xyes"; then break; fi
        done
        if test "x$have_cuda" = "xyes"; then break; fi
    done
fi

if test "x$have_cuda" = "xyes"; then
    AC_DEFINE([IBSIMU_CYCL_HAVE_CUDA], [1],
              [Define to 1 if the optional CUDA GPU backend is enabled])
    AC_MSG_NOTICE([CUDA GPU backend: enabled ($CUDA_CPPFLAGS $CUDA_LDFLAGS $CUDA_LIBS)])
else
    if test "x$with_cuda" = "xno"; then
        AC_MSG_NOTICE([CUDA GPU backend: disabled by --without-cuda])
    else
        AC_MSG_NOTICE([CUDA GPU backend: disabled (toolkit not found; use --with-cuda=DIR)])
    fi
fi
AC_SUBST([CUDA_CPPFLAGS])
AC_SUBST([CUDA_LDFLAGS])
AC_SUBST([CUDA_LIBS])

dnl ---------------------------------
dnl pkgconfig
dnl ---------------------------------

PKG_PROG_PKG_CONFIG()
if test "x$PKG_CONFIG" = "x" ; then
   AC_MSG_ERROR([pkg-config >= $PKGCONFIG_REQUIRED required but not
found (http://pkgconfig.freedesktop.org/)])
fi

dnl -------------
dnl zlib

CHECK_ZLIB()

dnl ---------------------------------
dnl Check for selectable compilation
dnl ---------------------------------

AC_ARG_ENABLE([sigsegv_stack],
  AS_HELP_STRING([--disable-sigsegv_stack], [Disable stack printing on SEGV]))

AS_IF([test "x$enable_sigsegv_stack" != "xno"], [
  CPPFLAGS="$CPPFLAGS -DSIGSEGV_STACK"
])

dnl Gtk
dnl Check by default, can be disabled
dnl ---------------------------------
AC_ARG_WITH([gtk3],
  [AS_HELP_STRING([--with-gtk3],
                  [enable interactive gtk plotter @<:@default=check@:>@])],
                  [],
                  [with_gtk3=check])
AS_IF([test "x$with_gtk3" != "xno"],
      [PKG_CHECK_MODULES([GTK3], [gtk+-3.0], [], [echo "Disabling functionality"])
       if $PKG_CONFIG --exists gtk+-3.0 ; then
         USEGTK3="yes"
         AC_DEFINE([GTK3], [1], [Define to 1 if Gtk3 dependent code is to be compiled.])
         LIBS="$LIBS $GTK3_LIBS"
         CPPFLAGS="$CPPFLAGS $GTK3_CFLAGS"
         AC_SUBST(GTK3_COND,[gtk+-3.0])
       elif test "x$with_gtk3" != "xcheck" ; then
         AC_MSG_ERROR([--with-gtk3 was given, but pkgconfig check for gtk+-3.0 library failed])
       else
         AC_SUBST(GTK3_COND,"")
       fi]
)
AM_CONDITIONAL([GTK3], [test x$USEGTK3 = xyes])

dnl OpenGL/GtkGLExt
dnl Check by default, can be disabled
dnl ---------------------------------
AC_ARG_WITH([opengl],
  [AS_HELP_STRING([--with-opengl],
                  [enable opengl 3d geometry plotter @<:@default=check@:>@])],
                  [],
                  [with_opengl=check])
AS_IF([test "x$with_opengl" != "xno"],
      [PKG_CHECK_MODULES([GTKGLEXT], [gtkglext-3.0], [], [echo "Using software rendering"])
       if $PKG_CONFIG --exists gtkglext-3.0 ; then
         USEOPENGL="yes"
         AC_DEFINE([OPENGL], [1], [Define to 1 if GtkGLExt/OpenGL dependent code is to be compiled.])
         LIBS="$LIBS $GTKGLEXT_LIBS"
         CPPFLAGS="$CPPFLAGS $GTKGLEXT_CFLAGS"
         AC_SUBST(OPENGL_COND,[gtkglext-3.0])
       elif test "x$with_opengl" != "xcheck" ; then
         AC_MSG_ERROR([--with-opengl was given, but pkgconfig check for gtkglext-3.0 library failed])
       else
         AC_SUBST(OPENGL_COND,"")
       fi]
)
AM_CONDITIONAL([OPENGL], [test x$USEOPENGL = xyes])


CHECK_UMFPACK()
AM_CONDITIONAL([UMFPACK], [test "x$USEUMFPACK" = "xyes"])


dnl CSG
dnl Check by default, can be disabled
dnl ---------------------------------
AC_ARG_WITH([csg],
  [AS_HELP_STRING([--with-csg],
                  [enable support for csg solid @<:@default=check@:>@])],
                  [],
                  [with_csg=check])
AS_IF([test "x$with_csg" != "xno"],
      [PKG_CHECK_MODULES([CSG], [csg], [], [echo "Disabling functionality"])
       if $PKG_CONFIG --exists csg ; then
         USECSG="yes"
         AC_DEFINE([CSG], [1], [Define to 1 if CSG solid is to be compiled.])
         LIBS="$LIBS $CSG_LIBS"
         CPPFLAGS="$CPPFLAGS $CSG_CFLAGS"
         AC_SUBST(CSG_COND,[csg])
       elif test "x$with_csg" != "xcheck" ; then
         AC_MSG_ERROR([--with-csg was given, but pkgconfig check for csg library failed])
       else
         AC_SUBST(CSG_COND,"")
       fi]
)
AM_CONDITIONAL([CSG], [test x$USECSG = xyes])

dnl ---------------------------------
dnl Check for libraries
dnl ---------------------------------

dnl -------------
dnl libm

AC_CHECK_LIBM
CPPFLAGS="$CFLAGS $CPPFLAGS"
LIBS="$LIBS $LIBM"

dnl -------------
dnl gtk+-3.0

dnl AM_PATH_GTK_3_0([3.6.0])

LIBS="$LIBS $GTK_LIBS"
CPPFLAGS="$CPPFLAGS $GTK_CFLAGS"


dnl -------------
dnl libpng

if $PKG_CONFIG --exists libpng ; then
  PKG_CHECK_MODULES([LIBPNG], [libpng])
else
  AC_MSG_ERROR([required library libpng not found])
fi
LIBS="$LIBS $LIBPNG_LIBS"
CPPFLAGS="$CPPFLAGS $LIBPNG_CFLAGS"


dnl -------------
dnl fontconfig

if $PKG_CONFIG --exists fontconfig ; then
  PKG_CHECK_MODULES([FONTCONFIG], [fontconfig])
else
  AC_MSG_ERROR([required library fontconfig not found])
fi
LIBS="$LIBS $FONTCONFIG_LIBS"
CPPFLAGS="$CPPFLAGS $FONTCONFIG_CFLAGS"


dnl -------------
dnl Freetype2

if $PKG_CONFIG --exists freetype2 ; then
  PKG_CHECK_MODULES([FREETYPE], [freetype2])
else
  AC_MSG_ERROR([required library freetype2 not found])
fi
LIBS="$LIBS $FREETYPE_LIBS"
CPPFLAGS="$CPPFLAGS $FREETYPE_CFLAGS"


dnl -------------
dnl Cairo
dnl
dnl No accurate information about required cairo version. 
dnl Tested to work down to 1.2.4.
dnl Postscript support may need higher version

if $PKG_CONFIG --exists cairo ; then
  PKG_CHECK_MODULES([CAIRO], [cairo >= 1.2.4])
else
  AC_MSG_ERROR([required library cairo not found])
fi
LIBS="$LIBS $CAIRO_LIBS"
CPPFLAGS="$CPPFLAGS $CAIRO_CFLAGS"

dnl -------------
dnl POSIX threads 

ACX_PTHREAD()
LIBS="$LIBS $PTHREAD_LIBS $PTHREAD_CFLAGS"
CPPFLAGS="$CPPFLAGS $PTHREAD_CFLAGS"

dnl AC_CHECK_HEADERS([pthread.h])
dnl if $HAVE_PTHREAD_H ; then
dnl   LIBS="$LIBS -lpthread"
dnl else
dnl   AC_MSG_ERROR([required library pthread not found])
dnl fi

dnl -------------
dnl GSL

if $PKG_CONFIG --exists gsl ; then
  PKG_CHECK_MODULES([GSL], [gsl >= 1.12])
else
  AC_MSG_ERROR([required version of gsl not found])
fi
LIBS="$LIBS $GSL_LIBS"
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"

dnl AC_CHECK_LIB([gslcblas],[cblas_dgemm])
dnl AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

dnl -------------
dnl getcwd

AC_CHECK_HEADERS([unistd.h])
AC_CHECK_FUNCS([getcwd])

dnl -------------
dnl siginfo_t

AC_CHECK_TYPES([siginfo_t], [], [], [[#include <signal.h>]])
AC_CHECK_TYPES([clockid_t], [], [], [[#include <time.h>]])
AC_CHECK_TYPES([struct timespec], [], [], [[#include <time.h>]])

dnl -------------
dnl Librt

AC_CHECK_LIB(rt, clock_gettime)

dnl -------------
dnl Time

AC_CHECK_HEADERS([sys/time.h])
AC_HEADER_TIME
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_FUNCS([clock_gettime])

dnl -------------
dnl strerror_r

AC_FUNC_STRERROR_R

dnl ---------------------------------
dnl HDF5 二进制 IO 后端 (IBSimu-Cycl: optional)
dnl
 dnl   --without-hdf5    禁用
dnl   --with-hdf5=DIR   指定 HDF5 前缀
 dnl
 dnl   用于大规模数据（>1e8 点）的场图/粒子导出（ROADMAP R3-④）。
 dnl   未探测到时 vtkwriter 的零依赖 VTK XML 通路仍然可用。
dnl ---------------------------------
AC_ARG_WITH([hdf5],
    [AS_HELP_STRING([--with-hdf5@<:@=DIR@:>@],
                    [HDF5 prefix for the optional binary IO backend (auto-detected)])],
    [], [with_hdf5=auto])

have_hdf5=no
HDF5_CPPFLAGS=""
HDF5_LDFLAGS=""
HDF5_LIBS=""

if test "x$with_hdf5" != "xno"; then
    if test "x$with_hdf5" = "xauto" || test "x$with_hdf5" = "xyes"; then
        dnl 优先用 pkg-config：Debian/Ubuntu 把头文件放在
        dnl /usr/include/hdf5/serial 下，只有 pkg-config 知道这个路径。
        for hdf5_pc in hdf5 hdf5-serial; do
            if $PKG_CONFIG --exists "$hdf5_pc" 2>/dev/null; then
                HDF5_CPPFLAGS=`$PKG_CONFIG --cflags "$hdf5_pc"`
                HDF5_LIBS=`$PKG_CONFIG --libs "$hdf5_pc"`
                have_hdf5=yes
                break
            fi
        done
        if test "x$have_hdf5" = "xno"; then
            for hdf5_p in /usr /usr/local; do
                if test -f "$hdf5_p/include/hdf5.h"; then
                    for hdf5_libd in "$hdf5_p/lib" "$hdf5_p/lib64"; do
                        if test -f "$hdf5_libd/libhdf5.so" || test -f "$hdf5_libd/libhdf5.a"; then
                            have_hdf5=yes
                            HDF5_CPPFLAGS="-I$hdf5_p/include"
                            HDF5_LDFLAGS="-L$hdf5_libd"
                            HDF5_LIBS="-lhdf5"
                            break
                        fi
                    done
                fi
                if test "x$have_hdf5" = "xyes"; then break; fi
            done
        fi
    else
        dnl 显式前缀：必须**真的**找到头文件与库才置 have_hdf5=yes，
        dnl 否则 configure 会假装成功、到编译时才炸（也允许 Debian 的
        dnl include/hdf5/serial 布局）。
        hdf5_inc="$with_hdf5/include"
        if test -f "$with_hdf5/include/hdf5/serial/hdf5.h"; then
            hdf5_inc="$with_hdf5/include/hdf5/serial"
        fi
        hdf5_libd=""
        for hdf5_d in "$with_hdf5/lib" "$with_hdf5/lib64"; do
            if test -f "$hdf5_d/libhdf5.so" || test -f "$hdf5_d/libhdf5.a"; then
                hdf5_libd="$hdf5_d"
                break
            fi
        done
        if test -f "$hdf5_inc/hdf5.h" && test -n "$hdf5_libd"; then
            have_hdf5=yes
            HDF5_CPPFLAGS="-I$hdf5_inc"
            HDF5_LDFLAGS="-L$hdf5_libd"
            HDF5_LIBS="-lhdf5"
        fi
    fi
fi

if test "x$have_hdf5" = "xyes"; then
    AC_DEFINE([IBSIMU_CYCL_HAVE_HDF5], [1],
              [Define to 1 if the optional HDF5 binary IO backend is enabled])
    AC_MSG_NOTICE([HDF5 binary IO backend: enabled ($HDF5_CPPFLAGS $HDF5_LDFLAGS $HDF5_LIBS)])
else
    if test "x$with_hdf5" = "xno"; then
        AC_MSG_NOTICE([HDF5 binary IO backend: disabled by --without-hdf5])
    else
        AC_MSG_NOTICE([HDF5 binary IO backend: disabled (library not found; use --with-hdf5=DIR)])
    fi
fi
AC_SUBST([HDF5_CPPFLAGS])
AC_SUBST([HDF5_LDFLAGS])
AC_SUBST([HDF5_LIBS])

dnl ---------------------------------
dnl Output
dnl ---------------------------------

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile doc/Makefile doc/Doxyfile doc/Doxyfile_web ibsimu-1.0.6dev.pc])
AC_OUTPUT

echo ""
echo "Optional compilation"
echo "--------------------"
echo ""
AS_IF([test "x$USEGTK3" = "xyes"],
      [echo "Gtk plotter     = yes"],
      [echo "Gtk plotter     = no"])
AS_IF([test "x$USEOPENGL" = "xyes"],
      [echo "OpenGL renderer = yes"],
      [echo "OpenGL renderer = no"])
AS_IF([test "x$USEUMFPACK" = "xyes"],
      [echo "UMFPACK solver  = yes"],
      [echo "UMFPACK solver  = no"])
AS_IF([test "x$USECSG" = "xyes"],
      [echo "CSG solid       = yes"],
      [echo "CSG solid       = no"])
echo ""
