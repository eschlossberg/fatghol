#! /bin/sh

basedir="`dirname $0`"

arch="`uname -m`"
v="`uname -r | cut -d. -f1-2`"
os="`uname -s`"

PYTHONPATH="$basedir":"$basedir"/build/lib.${os}-${arch}-${v}:$PYTHONPATH; export PYTHONPATH
LD_LIBRARY_PATH="$basedir"/build/lib.${os}-${arch}-${v}:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH

#exec python -O /usr/lib/python2.5/cProfile.py -o profile.bin "$basedir"/main.py "$@"
exec python "$basedir"/main.py "$@"
