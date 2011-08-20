
# for debug:
#CXXFLAGS=-O2 -g

# for production runs:
CXXFLAGS=-O2 -DNDEBUG

default: python

python: simplematrix.py _simplematrix.so
simplematrix.py simplematrix_wrap.cxx: simplematrix.i simplematrix.hpp
	swig -python -c++ simplematrix.i
_simplematrix.so: simplematrix_wrap.o
	g++ -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions simplematrix_wrap.o -o _simplematrix.so -L./lib -lgivaro -lgmp -llinbox
simplematrix_wrap.o: simplematrix_wrap.cxx simplematrix.hpp
	g++ $(CXXFLAGS) -pthread -fno-strict-aliasing -fwrapv -Wall -fPIC -I/usr/include/python2.5 -I./include -c simplematrix_wrap.cxx -o simplematrix_wrap.o

linbox:
	cd linbox-1.1.6 && \
	./configure --prefix=$$(dirname `pwd`) && \
	make && \
	make install

install: _simplematrix.so
	install -m 0755 _simplematrix.so lib/

clean:
	rm -f *~ \
		_simplematrix.so \
		simplematrix.py \
		simplematrix.pyo \
		simplematrix.pyc \
		simplematrix_wrap.o \
		simplematrix_wrap.cxx
