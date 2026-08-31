//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes
#include <pybind11/complex.h>

// Project includes
#include "includes/define_python.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"
#ifdef KRATOS_USE_EIGEN_BACKEND
#include "includes/kratos_eigen_interface.h"
#endif
#include "containers/array_1d.h"
#include "python/add_vector_to_python.h"

namespace Kratos::Python
{

namespace py = pybind11;

// The fixed-size array_1d types get value (copy) semantics for python slice
// reads: unlike the dynamic uBLAS vectors they are not required to model the
// uBLAS vector-expression concept (under the Eigen backend they do not), so
// no live vector_slice view can be built on them.
template<class T> struct IsArray1d : std::false_type {};
template<class T, std::size_t N> struct IsArray1d<array_1d<T, N>> : std::true_type {};

template< typename TVectorType >
py::class_< TVectorType > CreateVectorInterface(pybind11::module& m, std::string Name )
{
    py::class_< TVectorType, std::shared_ptr<TVectorType> > binder(m,Name.c_str(), py::buffer_protocol());
    binder.def(py::init<>());
    binder.def(py::init<const TVectorType&>());

    //binder.def(py::init<std::TVectorType& >())
    binder.def("Size", [](const TVectorType& self)
    {
        return self.size();
    } );
    binder.def("Resize", [](TVectorType& self, const typename TVectorType::size_type  new_size)
    {
        if(self.size() != new_size) self.resize(new_size, false);
    } );
    binder.def("__len__", [](const TVectorType& self)
    {
        return self.size();
    } );

    //operating on the object itself, +=, -=, *=, etc
    binder.def("__iadd__", [](TVectorType& self, const typename TVectorType::value_type scalar)
    {
        for(unsigned int i=0; i<self.size(); ++i) self[i]+=scalar;
        return self;
    }, py::is_operator());
    binder.def("__isub__", [](TVectorType& self, const typename TVectorType::value_type scalar)
    {
        for(unsigned int i=0; i<self.size(); ++i) self[i]-=scalar;
        return self;
    }, py::is_operator());
    binder.def("__imul__", [](TVectorType& self, const typename TVectorType::value_type scalar)
    {
        for(unsigned int i=0; i<self.size(); ++i) self[i]*=scalar;
        return self;
    }, py::is_operator());
    binder.def("__itruediv__", [](TVectorType& self, const typename TVectorType::value_type scalar)
    {
        for(unsigned int i=0; i<self.size(); ++i) self[i]/=scalar;
        return self;
    }, py::is_operator());

    binder.def("__iadd__", [](TVectorType& self, const TVectorType& other_vec)
    {
        noalias(self) += other_vec;
        return self;
    }, py::is_operator());
    binder.def("__isub__", [](TVectorType& self, const TVectorType& other_vec)
    {
        noalias(self) -= other_vec;
        return self;
    }, py::is_operator());

    //returning a different object
//         binder.def("__add__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]+=scalar; return vec1;}, py::is_operator());
//         binder.def("__sub__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]-=scalar; return vec1;}, py::is_operator());
    binder.def("__mul__", [](TVectorType vec1, const typename TVectorType::value_type scalar)
    {
        for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar;
        return vec1;
    }, py::is_operator());
    binder.def("__div__", [](TVectorType vec1, const typename TVectorType::value_type scalar)
    {
        for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar;
        return vec1;
    }, py::is_operator());
    binder.def("__truediv__", [](TVectorType vec1, const typename TVectorType::value_type scalar)
    {
        for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar;
        return vec1;
    }, py::is_operator());
//         binder.def("__radd__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]+=scalar; return vec1;}, py::is_operator());
//         binder.def("__rsub__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]-=scalar; return vec1;}, py::is_operator());
    binder.def("__rmul__", [](TVectorType vec1, const typename TVectorType::value_type scalar)
    {
        for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar;
        return vec1;
    }, py::is_operator());
    binder.def("__rdiv__", [](TVectorType vec1, const typename TVectorType::value_type scalar)
    {
        for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar;
    }, py::is_operator());
    binder.def("__add__", [](const TVectorType& vec1, const TVectorType& vec2)
    {
        TVectorType aux(vec1);
        aux += vec2;
        return aux;
    }, py::is_operator());
    binder.def("__sub__", [](const TVectorType& vec1, const TVectorType& vec2)
    {
        TVectorType aux(vec1);
        aux -= vec2;
        return aux;
    }, py::is_operator());

    //access operators
    binder.def("__setitem__", [](TVectorType& self, const unsigned int i, const typename TVectorType::value_type value)
    {
        self[i] = value;
    } );
    binder.def("__getitem__", [](const TVectorType& self, const unsigned int i)
    {
        return self[i];
    } );

    //access using slices
    binder.def("__setitem__", [](TVectorType &self, pybind11::slice this_slice, const TVectorType &value)
    {
        size_t start, stop, step, slicelength;
        if (!this_slice.compute(self.size(), &start, &stop, &step, &slicelength))
            throw pybind11::error_already_set();
        if (slicelength != value.size())
            throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
        for (size_t i = 0; i < slicelength; ++i)
        {
            self[start] = value[i];
            start += step;
        }
    });

    binder.def("__getitem__", [](TVectorType &self, pybind11::slice this_slice)
    {
        size_t start, stop, step, slicelength;
        if (!this_slice.compute(self.size(), &start, &stop, &step, &slicelength))
            throw pybind11::error_already_set();
        if constexpr (IsArray1d<TVectorType>::value) {
            // fixed-size arrays: return a dense copy of the slice
            DenseVector<typename TVectorType::value_type> sliced_copy(slicelength);
            for (size_t i = 0; i < slicelength; ++i) {
                sliced_copy[i] = self[start];
                start += step;
            }
            return sliced_copy;
        } else {
#ifdef KRATOS_USE_EIGEN_BACKEND
            // No uBLAS vector_slice over the Eigen-backed dynamic vector:
            // return a dense copy of the slice instead.
            TVectorType sliced_copy(slicelength);
            for (size_t i = 0; i < slicelength; ++i) {
                sliced_copy[i] = self[start];
                start += step;
            }
            return sliced_copy;
#else
            boost::numeric::ublas::slice ublas_slice(start, step, slicelength);
            boost::numeric::ublas::vector_slice<TVectorType> sliced_self(self, ublas_slice);
            return sliced_self;
#endif
        }
    });
    binder.def("fill", [](TVectorType& self, const typename TVectorType::value_type value)
    {
        noalias(self) = TVectorType(self.size(),value);
    });
    // unqualified so the calls resolve per vector type: to the uBLAS norms
    // for the dynamic vectors and to the backend overloads for array_1d
    binder.def("norm_1", [](TVectorType& self)
    {
        return norm_1(self);
    });
    binder.def("norm_2", [](TVectorType& self)
    {
        return norm_2(self);
    });
    binder.def("norm_inf", [](TVectorType& self)
    {
        return norm_inf(self);
    });
    binder.def("__iter__", [](TVectorType& self)
    {
        return py::make_iterator(self.begin(), self.end(), py::return_value_policy::reference_internal);
    }, py::keep_alive<0,1>() ) ;
    binder.def("__str__", PrintObject<TVectorType>);

    return std::move(binder);
}

template< std::size_t TSize >
void CreateArray1DInterface(pybind11::module& m, const std::string& Name )
{
    auto binder = CreateVectorInterface< array_1d<double,TSize> >(m,Name);
    binder.def(py::init( [](double value)
    {
        array_1d<double,TSize> tmp;
        for(std::size_t i=0; i < TSize; ++i)
            tmp[i] = value;
        return tmp;
    }));
    binder.def(py::init( [](const Vector& input)
    {
        KRATOS_ERROR_IF(input.size() != TSize)
                << "Attempting to initialize an array_1d<double," << TSize << "> from a Vector of size "
                << input.size() << ". Input should have size " << TSize <<"." << std::endl;

        array_1d<double,TSize> tmp(input);
        return tmp;
    }));
    binder.def(py::init<array_1d<double,TSize>>());
    binder.def(py::init( [](const py::list& input)
    {
        KRATOS_ERROR_IF(input.size() != TSize)
                << "Attempting to initialize an array_1d<double," << TSize << "> from a Python list of size "
                << input.size() << ". Input should have size " << TSize <<"." << std::endl;

        array_1d<double,TSize> tmp;
        for(std::size_t i=0; i<TSize; ++i)
        {
            tmp[i] = py::cast<double>(input[i]);
        }
        return tmp;
    }));
    binder.def(py::init( [](py::buffer b)
    {
        py::buffer_info info = b.request();
        KRATOS_ERROR_IF(info.shape[0] != TSize)
                << "Attempting to initialize an array_1d<double," << TSize << "> from a Python buffer of size "
                << info.shape[0] << ". Input should have size " << TSize <<"." << std::endl;
        KRATOS_ERROR_IF( info.ndim != 1 ) << "Buffer dimension of 1 is required, got: " << info.ndim << std::endl;
        array_1d<double,TSize> vec(info.shape[0]);

        for( int i=0; i<info.shape[0]; ++i )
        {
            vec[i]= static_cast<typename array_1d<double,TSize>::value_type *>(info.ptr)[i];
        }

        return vec;
    }));
    /*binder.def_buffer( [](array_1d<double,TSize>& self)-> py::buffer_info
    {
        return py::buffer_info(
            &self.data()[0], // pointer for the Eigen data(), unbounded_array for uBLAS
            sizeof(typename array_1d<double,TSize>::value_type),
            py::format_descriptor<typename array_1d<double,TSize>::value_type>::format(),
            1,
        {self.size()},
        {sizeof(typename array_1d<double,TSize>::value_type)}
        );
    });*/

    py::implicitly_convertible<py::buffer, array_1d<double,TSize>>();
    py::implicitly_convertible<py::list, array_1d<double,TSize>>();
    py::implicitly_convertible<Vector, array_1d<double,TSize>>();
}

void  AddVectorToPython(pybind11::module& m)
{
#ifndef KRATOS_USE_EIGEN_BACKEND
    // VectorSlice is a boost::numeric::ublas proxy over the Vector alias, so
    // it only exists when Vector is the uBLAS container (Eigen-backed Vectors
    // are sliced through numpy views on the python side instead).
    typedef boost::numeric::ublas::vector_slice<Vector> VectorSlice;
    py::class_< VectorSlice >(m, "VectorSlice")
    .def("Size", [](const VectorSlice& self)
    {
        return self.size();
    } )
    .def("__len__", [](const VectorSlice& self)
    {
        return self.size();
    } )
    .def("__iadd__", [](VectorSlice& self, const double scalar)
    {
        for(unsigned int i=0; i<self.size(); ++i) self[i]+=scalar;
        return self;
    }, py::is_operator())
    .def("__isub__", [](VectorSlice& self, const double scalar)
    {
        for(unsigned int i=0; i<self.size(); ++i) self[i]-=scalar;
        return self;
    }, py::is_operator())
    .def("__imul__", [](VectorSlice& self, const double scalar)
    {
        for(unsigned int i=0; i<self.size(); ++i) self[i]*=scalar;
        return self;
    }, py::is_operator())
    .def("__itruediv__", [](VectorSlice& self, const double scalar)
    {
        for(unsigned int i=0; i<self.size(); ++i) self[i]/=scalar;
        return self;
    }, py::is_operator())
    .def("__iadd__", [](VectorSlice& self, const VectorSlice& other_vec)
    {
        noalias(self) += other_vec;
        return self;
    }, py::is_operator())
    .def("__isub__", [](VectorSlice& self, const VectorSlice& other_vec)
    {
        noalias(self) -= other_vec;
        return self;
    }, py::is_operator())
    .def("__mul__", [](VectorSlice vec1, const double scalar)
    {
        for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar;
        return vec1;
    }, py::is_operator())
    .def("__div__", [](VectorSlice vec1, const double scalar)
    {
        for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar;
        return vec1;
    }, py::is_operator())
    .def("__rmul__", [](VectorSlice vec1, const double scalar)
    {
        for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar;
        return vec1;
    }, py::is_operator())
    .def("__rdiv__", [](VectorSlice vec1, const double scalar)
    {
        for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar;
    }, py::is_operator())
    .def("__add__", [](const VectorSlice& vec1, const VectorSlice& vec2)
    {
        Vector aux(vec1);
        aux += vec2;
        return aux;
    }, py::is_operator())
    .def("__sub__", [](const VectorSlice& vec1, const VectorSlice& vec2)
    {
        Vector aux(vec1);
        aux -= vec2;
        return aux;
    }, py::is_operator())
    .def("__setitem__", [](VectorSlice& self, const unsigned int i, const typename VectorSlice::value_type value)
    {
        self[i] = value;
    })
    .def("__getitem__", [](const VectorSlice& self, const unsigned int i)
    {
        return self[i];
    } )
    .def("__setitem__", [](VectorSlice &self, pybind11::slice this_slice, const VectorSlice &value)
    {
        size_t start, stop, step, slicelength;
        if (!this_slice.compute(self.size(), &start, &stop, &step, &slicelength))
            throw pybind11::error_already_set();
        if (slicelength != value.size())
            throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
        for (size_t i = 0; i < slicelength; ++i)
        {
            self[start] = value[i];
            start += step;
        }
    })
    .def("__setitem__", [](VectorSlice &self, pybind11::slice this_slice, const Vector &value)
    {
        size_t start, stop, step, slicelength;
        if (!this_slice.compute(self.size(), &start, &stop, &step, &slicelength))
            throw pybind11::error_already_set();
        if (slicelength != value.size())
            throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
        for (size_t i = 0; i < slicelength; ++i)
        {
            self[start] = value[i];
            start += step;
        }
    })
    .def("__iter__", [](VectorSlice& self)
    {
        return py::make_iterator(self.begin(), self.end(), py::return_value_policy::reference_internal);
    }, py::keep_alive<0,1>() )
    .def("__str__", PrintObject<VectorSlice>)
    ;
#endif // !KRATOS_USE_EIGEN_BACKEND

    auto vector_binder = CreateVectorInterface<Vector>(m, "Vector");
    vector_binder.def(py::init<typename Vector::size_type>());
    vector_binder.def(py::init<typename Vector::size_type, double>());
    vector_binder.def(py::init<Vector>());
    vector_binder.def(py::init( [](const array_1d<double,3>& input)
    {
        // element-wise so no vector-expression concept is required of array_1d
        Vector tmp(3);
        for(std::size_t i=0; i<3; ++i)
            tmp[i] = input[i];
        return tmp;
    }));
    vector_binder.def(py::init( [](const py::list& input)
    {
        Vector tmp(input.size());
        for(unsigned int i=0; i<tmp.size(); ++i)
            tmp[i] = py::cast<double>(input[i]);
        return tmp;
    }));
#ifndef KRATOS_USE_EIGEN_BACKEND
    vector_binder.def(py::init([](const VectorSlice& input) -> Vector
    {
        return input;
    }));
#endif
    vector_binder.def(py::init( [](py::buffer b)
    {
        py::buffer_info info = b.request();
        KRATOS_ERROR_IF( info.format != py::format_descriptor<typename Vector::value_type >::value ) << "Expected a double array\n";
        KRATOS_ERROR_IF( info.ndim != 1 ) << "Buffer dimension of 1 is required, got: " << info.ndim << std::endl;
        Vector vec(info.shape[0]);

        for( int i=0; i<info.shape[0]; ++i )
        {
            vec[i]= static_cast<typename Vector::value_type *>(info.ptr)[i];
        }

        return vec;
    }));
    vector_binder.def_buffer( [](Vector& self)-> py::buffer_info
    {
        return py::buffer_info(
            &self.data()[0], // pointer for the Eigen data(), unbounded_array for uBLAS
            sizeof(typename Vector::value_type),
            py::format_descriptor<typename Vector::value_type>::format(),
            1,
        {self.size()},
        {sizeof(typename Vector::value_type)}
        );
    });
    py::implicitly_convertible<py::buffer, Vector>();
    py::implicitly_convertible<py::list, Vector>();
    py::implicitly_convertible<array_1d<double,3>, Vector>();

#ifdef KRATOS_USE_EIGEN_BACKEND
    // Under the Eigen backend the dynamic Vector alias IS EigenVector<double>,
    // so the historical "EigenVector" python name (the strategies' system
    // vector type) is an alias of the "Vector" class rather than a second
    // registration of the same C++ type.
    m.attr("EigenVector") = m.attr("Vector");
#endif

    //***********************************************************************************
        //***********************************************************************************
    py::class_< DenseVector<bool>, std::shared_ptr<DenseVector<bool>> > bool_vector_binder(m,"BoolVector", py::buffer_protocol());
    bool_vector_binder.def(py::init<>());
    bool_vector_binder.def(py::init<const DenseVector<bool>&>());
    bool_vector_binder.def(py::init< DenseVector<bool>::size_type>());
    bool_vector_binder.def(py::init< DenseVector<bool>::size_type, bool>());
    bool_vector_binder.def(py::init( [](const py::list& input)
    {
        DenseVector<bool> tmp(input.size());
        for(unsigned int i=0; i<tmp.size(); ++i)
            tmp[i] = py::cast<bool>(input[i]);
        return tmp;
    }));
    bool_vector_binder.def("Size", [](const DenseVector<bool>& self)
    {
        return self.size();
    } );
    bool_vector_binder.def("Resize", [](DenseVector<bool>& self, const typename DenseVector<bool>::size_type  new_size)
    {
        if(self.size() != new_size) self.resize(new_size, false);
    } );
    bool_vector_binder.def("__len__", [](const DenseVector<bool>& self)
    {
        return self.size();
    } );
    bool_vector_binder.def("__setitem__", [](DenseVector<bool>& self,
                                            const unsigned int i,
                                            const bool value)
    {
        self[i] = value;
    } );
    bool_vector_binder.def("__getitem__", [](const DenseVector<bool>& self, const unsigned int i)
    {
        return self[i];
    } );
    bool_vector_binder.def("__iter__", [](DenseVector<bool>& self)
    {
        return py::make_iterator(self.begin(), self.end(), py::return_value_policy::reference_internal);
    }, py::keep_alive<0,1>() ) ;
    bool_vector_binder.def("__str__", PrintObject<DenseVector<bool>>);

    py::implicitly_convertible<py::list, DenseVector<bool>>();

    //***********************************************************************************
    auto int_vector_binder = CreateVectorInterface<DenseVector<int>>(m, "DenseVectorInt");
    int_vector_binder.def(py::init<typename DenseVector<int>::size_type>());
    int_vector_binder.def(py::init<typename DenseVector<int>::size_type, int>());
    int_vector_binder.def(py::init<DenseVector<int>>());
    int_vector_binder.def(py::init( [](const py::list& input)
    {
        DenseVector<int> tmp(input.size());
        for(unsigned int i=0; i<tmp.size(); ++i)
            tmp[i] = py::cast<int>(input[i]);
        return tmp;
    }));
    int_vector_binder.def(py::init( [](py::buffer b)
    {
        py::buffer_info info = b.request();
        KRATOS_ERROR_IF( info.format != py::format_descriptor<typename DenseVector<int>::value_type >::value ) << "Expected a double array\n";
        KRATOS_ERROR_IF( info.ndim != 1 ) << "Buffer dimension of 1 is required, got: " << info.ndim << std::endl;
        DenseVector<int> vec(info.shape[0]);

        for( int i=0; i<info.shape[0]; ++i )
        {
            vec[i]= static_cast<typename DenseVector<int>::value_type *>(info.ptr)[i];
        }

        return vec;
    }));
    int_vector_binder.def_buffer( [](DenseVector<int>& self)-> py::buffer_info
    {
        return py::buffer_info(
            &self.data()[0], // pointer for the Eigen data(), unbounded_array for uBLAS
            sizeof(typename DenseVector<int>::value_type),
            py::format_descriptor<typename DenseVector<int>::value_type>::format(),
            1,
        {self.size()},
        {sizeof(typename DenseVector<int>::value_type)}
        );
    });
    py::implicitly_convertible<py::list, DenseVector<int>>();

//***********************************************************************************
    auto unsigned_int_vector_binder = CreateVectorInterface<DenseVector<unsigned int>>(m, "DenseVectorUnsignedInt");
    unsigned_int_vector_binder.def(py::init<typename DenseVector<unsigned int>::size_type>());
    unsigned_int_vector_binder.def(py::init<typename DenseVector<unsigned int>::size_type, int>());
    unsigned_int_vector_binder.def(py::init<DenseVector<unsigned int>>());
    unsigned_int_vector_binder.def(py::init( [](const py::list& input)
    {
        DenseVector<unsigned int> tmp(input.size());
        for(unsigned int i=0; i<tmp.size(); ++i)
            tmp[i] = py::cast<int>(input[i]);
        return tmp;
    }));
    unsigned_int_vector_binder.def(py::init( [](py::buffer b)
    {
        py::buffer_info info = b.request();
        KRATOS_ERROR_IF( info.format != py::format_descriptor<typename DenseVector<unsigned int>::value_type >::value ) << "Expected a double array\n";
        KRATOS_ERROR_IF( info.ndim != 1 ) << "Buffer dimension of 1 is required, got: " << info.ndim << std::endl;
        DenseVector<unsigned int> vec(info.shape[0]);

        for( unsigned int i=0; i<info.shape[0]; ++i ) {
            vec[i]= static_cast<typename DenseVector<unsigned int>::value_type *>(info.ptr)[i];
        }

        return vec;
    }));
    unsigned_int_vector_binder.def_buffer( [](DenseVector<unsigned int>& self)-> py::buffer_info
    {
        return py::buffer_info(
            &self.data()[0], // pointer for the Eigen data(), unbounded_array for uBLAS
            sizeof(typename DenseVector<unsigned int>::value_type),
            py::format_descriptor<typename DenseVector<unsigned int>::value_type>::format(),
            1,
        {self.size()},
        {sizeof(typename DenseVector<unsigned int>::value_type)}
        );
    });
    py::implicitly_convertible<py::list, DenseVector<unsigned int>>();


    //***********************************************************************************
    auto cplx_vector_binder = CreateVectorInterface<ComplexVector>(m, "ComplexVector");
    cplx_vector_binder.def(py::init<typename ComplexVector::size_type>());
    cplx_vector_binder.def(py::init<typename ComplexVector::size_type, double>());
    cplx_vector_binder.def(py::init<typename ComplexVector::size_type, std::complex<double>>());
    cplx_vector_binder.def(py::init<Vector>());
    cplx_vector_binder.def(py::init<ComplexVector>());
    // vector_binder.def(py::init<array_1d<double,3>>());
    cplx_vector_binder.def(py::init( [](const py::list& input)
    {
        ComplexVector tmp(input.size());
        for(unsigned int i=0; i<tmp.size(); ++i)
            tmp[i] = py::cast<std::complex<double>>(input[i]);
        return tmp;
    }));
    py::implicitly_convertible<py::list, ComplexVector>();

    CreateArray1DInterface< 3 >(m,"Array3");
    CreateArray1DInterface< 4 >(m,"Array4");
    CreateArray1DInterface< 6 >(m,"Array6");
    CreateArray1DInterface< 9 >(m,"Array9");

}
}  // namespace Kratos::Python.
