#! /bin/sh

basedir="`dirname $0`"

arch="`uname -m`"
v="`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f1,2`"
os="`uname -s | tr A-Z a-z`"

PYTHONPATH="$basedir":"$basedir"/build/lib.${os}-${arch}-${v}:$PYTHONPATH; export PYTHONPATH
LD_LIBRARY_PATH="$basedir"/build/lib.${os}-${arch}-${v}:"$basedir"/sw/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH

#exec python -O /usr/lib/python2.5/cProfile.py -o profile.bin "$basedir"/main.py "$@"
exec python "$basedir"/main.py "$@"
