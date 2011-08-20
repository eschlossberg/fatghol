// --- python interface ---

#include <boost/python.hpp>
using namespace boost::python;

#include "simplematrix.hpp"

BOOST_PYTHON_MODULE(simplematrix)
{
  class_<SimpleMatrix>("simplematrix", init<size_t, size_t>())
    .def("entry", &SimpleMatrix::entry)
    .def("rank", &SimpleMatrix::rank)
    ;
}


