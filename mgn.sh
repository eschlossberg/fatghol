#! /bin/sh

basedir="`dirname $0`"

PYTHONPATH="$basedir":"$basedir"/simplematrix:$PYTHONPATH; export PYTHONPATH
LD_LIBRARY_PATH="$basedir"/simplematrix/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH

exec python "$basedir"/main.py "$@"
