#ifndef AMG_CORE_VECTOR_HPP_
#define AMG_CORE_VECTOR_HPP_

#include "types.h"

// Vector Class
//
// This class constructs a vector, supporting simple linear
// algebra operations.
//
// Attributes
// -------------
// values : std::vector<double>
//    stl vector of vector values
// size : index_t
//    Dimension of vector
//
// Methods
// -------
// set_const_value(data_t alpha)
//    Sets the vector to a constant value
// set_rand_values()
//    Sets each element of the vector to a random value
// axpy(Vector& y, data_t alpha)
//    Multiplies each element by a constant, alpha, and then
//    adds corresponding values from y
// scale(data_t alpha)
//    Multiplies entries of vector by a constant
// norm(index_t p)
//    Calculates the p-norm of the vector
// print()
//    Prints the nonzero values and positions
// data()
//    Returns the data values as a data_t*
//
namespace AMG
{
    class Vector
    {

    public:
        /**************************************************************
        *****   Vector Class Constructor
        **************************************************************
        ***** Initializes an empty vector of the given size
        *****
        ***** Parameters
        ***** -------------
        ***** len : index_t
        *****    Size of the vector
        **************************************************************/
        Vector(int len)
        {
            resize(len);
        }

        /**************************************************************
        *****   Vector Class Constructor
        **************************************************************
        ***** Initializes an empty vector without setting the size
        **************************************************************/
        Vector()
        {
            num_values = 0;
        }

        Vector(const Vector& v)
        {
            copy(v);
        }

        void resize(int len)
        {
            values.resize(len);
            num_values = len;
        }

        /**************************************************************
        *****   Vector Set Constant Value
        **************************************************************
        ***** Initializes the vector to a constant value
        *****
        ***** Parameters
        ***** -------------
        ***** alpha : data_t
        *****    Constant value to set each element of vector to
        **************************************************************/
        void set_const_value(data_t alpha);

        /**************************************************************
        *****   Vector Set Random Values
        **************************************************************
        ***** Initializes each element of the vector to a random
        ***** value
        **************************************************************/
        void set_rand_values(data_t min = 0.0, data_t max = 1.0);

        /**************************************************************
        *****   Vector AXPY
        **************************************************************
        ***** Multiplies the vector by a constant, alpha, and then
        ***** sums each element with corresponding entry of Y
        *****
        ***** Parameters
        ***** -------------
        ***** y : Vector&
        *****    Vector to be summed with
        ***** alpha : data_t
        *****    Constant value to multiply each element of vector by
        **************************************************************/
        void axpy(Vector& y, data_t alpha);

        /**************************************************************
        *****   Vector Copy
        **************************************************************
        ***** Copies each vector value of y into values
        *****
        ***** Parameters
        ***** -------------
        ***** y : Vector&
        *****    Vector to be copied
        **************************************************************/
        void copy(const Vector& y);

        /**************************************************************
        *****   Vector Scale
        **************************************************************
        ***** Multiplies each element of the vector by a constant value
        *****
        ***** Parameters
        ***** -------------
        ***** alpha : data_t
        *****    Constant value to set multiply element of vector by
        **************************************************************/
        void scale(data_t alpha);

        /**************************************************************
        *****   Vector Norm
        **************************************************************
        ***** Calculates the P norm of the vector (for a given P)
        *****
        ***** Parameters
        ***** -------------
        ***** p : index_t
        *****    Determines which p-norm to calculate
         **************************************************************/
        data_t norm(index_t p);

        /**************************************************************
        *****   Print Vector
        **************************************************************
        ***** Prints all nonzero elements in vector
        *****
        ***** Parameters
        ***** -------------
        ***** vec_name : const char* (optional)
        *****    Name to be printed.  Default prints Vec[%d] = %e.
        **************************************************************/
        void print(const char* vec_name = "Vec");

        /**************************************************************
        *****   Vector Element Access
        **************************************************************
        ***** Function overload for element access
        *****
        ***** Returns
        ***** ------------
        ***** data_t& element at position passed
        **************************************************************/
        data_t& operator[](const int index);

        /**************************************************************
        *****   Vector Data
        **************************************************************
        ***** Returns pointer to vector entries
        *****
        ***** Returns
        ***** -------------
        ***** data_t*
        *****    Pointer to values of vector
        **************************************************************/
        data_t* data()
        {
            return values.data();
        }

        index_t size()
        {
            return num_values;
        }

        data_t inner_product(Vector& x);

        std::vector<double> values;
        index_t num_values;
    }; // class Vector

}   // namespace AMG


#endif
