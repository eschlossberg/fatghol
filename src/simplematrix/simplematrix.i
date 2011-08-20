/**
 * @file   simple.i
 *
 * Interface of the matrix class.
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */

%module simplematrix

%{
#define SWIG_FILE_WITH_INIT
#include "simplematrix.hpp"
%}

%include simplematrix.hpp
