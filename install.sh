#! /bin/sh
#
## Customization section

# comment any of the following lines if you already have
# the required version (or a newer one) installed in the
# standard PATH/LD_LIBRARY_PATH
#
PYTHON=2.7.2 # http://www.python.org/ftp/python/2.7.2/Python-2.7.2.tar.bz2
CYTHON=0.15.1
SWIG=1.3.40

# LinBox requires GMP, Givaro and ATLAS
LINBOX=1.1.7
GMP=5.0.2
GIVARO=3.3.3 # later versions in the 3.3.x series give incorrect results
ATLAS=3.8.3 # 3.8.4 does not compile on Ubuntu 11.04 w/ GCC 4.5.2

# this is really optional
TCMALLOC=1.8.3 # http://google-perftools.googlecode.com/files/google-perftools-1.8.3.tar.gz


## No customization should be necessary further down here

PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG [DIR [CFLAGS ...]]

Download and install all required software for FatGHoL into
directory DIR.  If omitted, DIR defaults to '`pwd`/sw'.

Any additional arguments are placed directly on the 'configure'
command line of every package, so it can be used to set CFLAGS,
CXXFLAGS etc.

EOF
}


## helper functions
die () {
  rc="$1"
  shift
  (echo -n "$PROG: ERROR: ";
      if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
  exit $rc
}

have_command () {
  type "$1" >/dev/null 2>/dev/null
}

require_command () {
  if ! have_command "$1"; then
    die 1 "Could not find required command '$1' in system PATH. Aborting."
  fi
}

is_absolute_path () {
  expr "$1" : '/' >/dev/null 2>/dev/null
}

_ () {
    echo
    echo ==== "$@" ...;
}


## parse command-line 

case "$1" in
    -h|--help)
        usage
        exit 0
        ;;
    *)
esac


## main

require_command bzip2
require_command gzip
require_command make
require_command tar
require_command wget

set -e

# target directory
if [ -n "$1" ]; then
    sw="$1"
    shift
else
    sw='sw'
fi
if ! is_absolute_path "$sw"; then 
    sw="$(pwd)/${sw}"
fi
mkdir -p "${sw}"
mkdir -p "${sw}/src"
cd "${sw}"

# paths
PATH=${sw}/bin:$PATH
LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
PYTHONPATH=${sw}/lib/python

# misc
ncpus="$(grep -c '^processor' /proc/cpuinfo)"
concurrent_make="make -j $ncpus"

# the `cpu MHz` line in /proc/cpuinfo states the *current* CPU clock,
# whereas "bogomips" is approximately twice the maximal clock on all Intel/AMD CPUs,
# see: http://tldp.org/HOWTO/BogoMips/bogo-faq.html#AEN192
bogomips=`fgrep bogomips /proc/cpuinfo | head -1 | cut -d: -f2 | cut -d. -f1`
mhz=`expr $bogomips / 2`

case `uname -m` in
        i?86) bits=32 ;;
        x86_64) bits=64 ;;
        *) die 1 "Unknown architecture `uname -m`: is it 32-bit or 64-bit?" ;;
esac


# Python
if [ -n "${PYTHON}" ]; then
    _ Installing Python ${PYTHON}
    cd ${sw}/src/
    wget -N http://www.python.org/ftp/python/${PYTHON}/Python-${PYTHON}.tar.bz2
    set -x
    tar -xjf Python-${PYTHON}.tar.bz2
    cd Python-${PYTHON}
    ./configure --prefix=${sw} "$@"
    $concurrent_make
    make install
    set +x
fi # PYTHON


# SWIG
if [ -n "${SWIG}" ]; then
    _ Downloading SWIG ${SWIG}
    cd ${sw}/src/
    wget -N http://downloads.sourceforge.net/project/swig/swig/swig-${SWIG}/swig-${SWIG}.tar.gz
    set -x
    tar -xzf swig-${SWIG}.tar.gz
    cd swig-${SWIG}
    ./configure --prefix=${sw} \
        --with-python=${sw}/bin \
        --without-allegrocl \
        --without-chicken \
        --without-clisp \
        --without-csharp \
        --without-gcj \
        --without-guile \
        --without-java \
        --without-lua \
        --without-mzscheme \
        --without-ocaml \
        --without-octave \
        --without-perl5 \
        --without-php4 \
        --without-pike \
        --without-r \
        --without-ruby \
        --without-rxspencer \
        --without-tcl \
        "$@";
    $concurrent_make
    make install
    set +x
fi # SWIG


# Cython
if [ -n "$CYTHON" ]; then
    _ Installing Cython-${CYTHON}
    cd ${sw}/src/
    wget -N http://www.cython.org/release/Cython-${CYTHON}.tar.gz
    set -x
    tar -xzf Cython-${CYTHON}.tar.gz
    cd Cython-${CYTHON}
    python setup.py build
    python setup.py install --home=${sw}
    
    PYTHONPATH=$PYTHONPATH:`pwd`; export PYTHONPATH
    PATH=$PATH:`pwd`/bin; export PATH
    set +x
fi # CYTHON


# GMP
if [ -n "$GMP" ]; then
    _ Installing GMP ${GMP}
    cd ${sw}/src
    wget -N ftp://sunsite.cnlab-switch.ch/mirror/gnu/gmp/gmp-${GMP}.tar.bz2
    set -x
    tar -xjf gmp-${GMP}.tar.bz2
    cd gmp-${GMP}
    ./configure --prefix=${sw} --enable-cxx "$@"
    $concurrent_make
    make install
    set +x
fi # GMP


# GIVARO (cfr. http://groups.google.com/group/linbox-use/browse_thread/thread/82673844f6921271)
if [ -n "$GIVARO" ]; then
    _ Installing Givaro ${GIVARO}
    cd ${sw}/src
    case "$GIVARO" in
        [0-9]*) 
            # it's a version number, use fixed-format URL
            wget -N http://www-lmc.imag.fr/CASYS/LOGICIELS/givaro/Downloads/givaro-${GIVARO}.tar.gz
            ;;
        *) 
            # presume it's a URL: download it and extract version number
            wget --no-check-certificate -N "$GIVARO" 
            GIVARO=$(expr match "$GIVARO" '.*/givaro-\([0-9\.]\+\).tar.gz')
            ;;
    esac
    set -x
    tar -xzf givaro-${GIVARO}.tar.gz
    cd givaro-${GIVARO%%rc[0-9]}
    # work around bug in ./configure: the test for GMP cannot find it
    # unless it's in the LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
    ./configure  --prefix=${sw} --enable-shared ${GMP:+"--with-gmp=${sw}"} "$@"
    $concurrent_make
    make install
    set +x
fi # GIVARO


# ATLAS
if [ -n "$ATLAS" ]; then
    cd ${sw}/src
    wget "http://switch.dl.sourceforge.net/sourceforge/math-atlas/atlas${ATLAS}.tar.bz2" \
         -O atlas${ATLAS}.tar.bz2
    set -x
    tar -xjf atlas${ATLAS}.tar.bz2
    cd ATLAS
    mkdir -p BLDdir
    cd BLDdir
    ../configure -v 2 \
        -b ${bits} \
        -m ${mhz} \
        -D c -DPentiumCPS=${mhz} \
        -Si cputhrchk 0 \
        -Ss pmake "$concurrent_make" \
        --prefix=${sw}
    make build
    # for shared library support, the following flags
    # must be added to the `../configure` line above:
    #             -Fa alg -fPIC 
    # then the shared libraries can be compiled with this:
    #(cd lib; make cshared cptshared && cp -a *.so ${sw}/lib)
    make install
    set +x
fi # ATLAS


# LinBox
if [ -n "$LINBOX" ]; then
    _ Installing LinBox ${LINBOX}
    cd ${sw}/src
    wget -N http://linalg.org/linbox-${LINBOX}.tar.gz
    set -x
    tar -xzf linbox-${LINBOX}.tar.gz
    cd linbox-${LINBOX}
    ./configure --prefix=${sw} \
        ${ATLAS:+"--with-blas=${sw}"} \
        ${GMP:+"--with-gmp=${sw}"} \
        ${GIVARO:+"--with-givaro=${sw}"} \
        "$@";
    if [ "$LINBOX" = '1.1.7' ]; then
        # patch error in distribute sources
        sed -r -i -e 's/^namespace/using namespace std;\nnamespace/;' \
            linbox/algorithms/rational-reconstruction-base.h \
            linbox/algorithms/rational-reconstruction.h
    fi
    $concurrent_make
    make install
    set +x
fi # LINBOX


# Google perftools
if [ -n "$TCMALLOC" ]; then
    _ Installing Google PerfTools $TCMALLOC ...
    cd ${sw}/src
    set -x
    wget -N http://google-perftools.googlecode.com/files/google-perftools-${TCMALLOC}.tar.gz
    tar -xzf "google-perftools-${TCMALLOC}.tar.gz"
    cd google-perftools-${TCMALLOC}
    ./configure --prefix=${sw} \
        --enable-frame-pointers --disable-debugalloc \
        "$@";
    $concurrent_make
    make install
    set +x
fi


_ All done.
