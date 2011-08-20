
default: python

python: simplematrix.py _simplematrix.so
simplematrix.py simplematrix_wrap.cxx: simplematrix.i
	swig -python -c++ simplematrix.i
_simplematrix.so: simplematrix_wrap.o
	g++ -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions simplematrix_wrap.o -o _simplematrix.so -L./lib -lgivaro -lgmp -llinbox
simplematrix_wrap.o: simplematrix_wrap.cxx
	g++ -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -fPIC -I/usr/include/python2.5 -I./include -c simplematrix_wrap.cxx -o simplematrix_wrap.o

linbox:
	cd linbox-1.1.6 && \
	./configure --prefix=$$(dirname `pwd`) && \
	make && \
	make install

install: _simplematrix.so
	install -m 0755 _simplematrix.so lib/

clean:
	rm -f *~ \
		_simplematrix.so
		simplematrix_wrap.o \
		simplematrix_wrap.cxx
