#! /bin/sh

basedir="`dirname $0`"
if [ "$basedir" = '.' ]; then
  basedir="`pwd`"
fi

## set up Python module search path

PATH="$basedir"/sw/bin:$PATH; export PATH
if type hash >/dev/null; then
    # bash might have stored the 'python' location -- reset it
    hash -d python 2>/dev/null
fi

arch="`uname -m`"
v="`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f1,2`"
os="`uname -s | tr A-Z a-z`"

PYTHONPATH="$basedir":"$basedir"/build/lib.${os}-${arch}-${v}:$PYTHONPATH; export PYTHONPATH
LD_LIBRARY_PATH="$basedir"/build/lib.${os}-${arch}-${v}:"$basedir"/sw/lib/python${v}:"$basedir"/sw/lib:"$basedir":$LD_LIBRARY_PATH; export LD_LIBRARY_PATH

## try to use GNU time
if [ -x /usr/bin/time ]; then
    # usual location in Linux systems
    time_cmd=/usr/bin/time
else
    # whatever is in PATH
    time_cmd=time
fi
if ($time_cmd --format='%M' /bin/true) >/dev/null 2>&1; then
    :   # time is GNU time, nothing to do
else
    time_cmd=''
fi

## run the main program, with timing information
if [ -n "$time_cmd" ]; then
    exec $time_cmd --format='\nUSED: wctime=%e cputime.usr=%U cputime.sys=%S maxmem=%MkB cs.involuntary=%c cs.voluntary=%w majflt=%F minflt=%R' python "$basedir"/fatghol/mgn.py "$@"
else
    exec  python "$basedir"/fatghol/mgn.py "$@"
fi

