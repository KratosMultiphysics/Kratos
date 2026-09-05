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

namespace
{

#ifdef KRATOS_USE_EIGEN_BACKEND
// Python-side counterpart of the uBLAS vector_slice proxy for the Eigen-backed
// Vector: a strided, writable, non-owning view into the parent's storage, so
// slice reads and writes go through to the parent (numpy-like behaviour). As
// with numpy, resizing the parent reallocates its storage and invalidates any
// outstanding view.
using PythonVectorSlice = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>>;
#else
using PythonVectorSlice = boost::numeric::ublas::vector_slice<Vector>;
#endif

/// Live (write-through) slice view over a uBLAS container.
template<class TVectorType>
auto MakeVectorSliceView(TVectorType& rVector, const std::size_t Start, const std::size_t Step, const std::size_t Length)
{
    return boost::numeric::ublas::vector_slice<TVectorType>(rVector, boost::numeric::ublas::slice(Start, Step, Length));
}

#ifdef KRATOS_USE_EIGEN_BACKEND
/// Live (write-through) slice view over the Eigen-backed Vector. A non-template
/// overload, so it is an exact match for Vector and the uBLAS proxy above (not
/// constructible over an Eigen container) is never selected for it.
inline PythonVectorSlice MakeVectorSliceView(Vector& rVector, const std::size_t Start, const std::size_t Step, const std::size_t Length)
{
    // data() is nullptr for an empty Vector; a zero-length map over it is valid.
    // Step comes from py::slice::compute(size_t): a negative python step is
    // wrapped there, and the cast to the signed Eigen::Index restores it.
    return PythonVectorSlice(rVector.data() + Start, static_cast<Eigen::Index>(Length), Eigen::InnerStride<>(static_cast<Eigen::Index>(Step)));
}
#endif

/// size() as an unsigned count for both slice proxies (Eigen's is signed).
template<class TSliceType>
std::size_t SliceSize(const TSliceType& rSlice)
{
    return static_cast<std::size_t>(rSlice.size());
}

} // namespace

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
            Vector sliced_copy(slicelength);
            for (size_t i = 0; i < slicelength; ++i) {
                sliced_copy[i] = self[start];
                start += step;
            }
            return sliced_copy;
        } else {
            // dynamic vectors: return a live view, so writes go through to self
            return MakeVectorSliceView(self, start, step, slicelength);
        }
    }, py::keep_alive<0,1>()); // the view must not outlive its parent
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
            GetContiguousDataPointer(self), // raw storage pointer for both backends, valid (nullptr) for an empty container
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

template<class TSliceType>
void CreateVectorSliceInterface(pybind11::module& m, const std::string& Name)
{
    using value_type = typename TSliceType::value_type;

    // Both slice proxies (uBLAS vector_slice, Eigen strided Map) are views: a
    // by-value copy still refers to the parent's storage, so the in-place
    // operators return self by value and keep the parent alive through the
    // returned object. All element-wise operations are written as index loops
    // so the same interface serves both proxies.
    py::class_< TSliceType >(m, Name.c_str())
    .def("Size", [](const TSliceType& self)
    {
        return SliceSize(self);
    } )
    .def("__len__", [](const TSliceType& self)
    {
        return SliceSize(self);
    } )
    .def("__iadd__", [](TSliceType& self, const value_type scalar)
    {
        for(std::size_t i=0; i<SliceSize(self); ++i) self[i]+=scalar;
        return self;
    }, py::is_operator(), py::keep_alive<0,1>())
    .def("__isub__", [](TSliceType& self, const value_type scalar)
    {
        for(std::size_t i=0; i<SliceSize(self); ++i) self[i]-=scalar;
        return self;
    }, py::is_operator(), py::keep_alive<0,1>())
    .def("__imul__", [](TSliceType& self, const value_type scalar)
    {
        for(std::size_t i=0; i<SliceSize(self); ++i) self[i]*=scalar;
        return self;
    }, py::is_operator(), py::keep_alive<0,1>())
    .def("__itruediv__", [](TSliceType& self, const value_type scalar)
    {
        for(std::size_t i=0; i<SliceSize(self); ++i) self[i]/=scalar;
        return self;
    }, py::is_operator(), py::keep_alive<0,1>())
    .def("__iadd__", [](TSliceType& self, const TSliceType& other_vec)
    {
        KRATOS_ERROR_IF(SliceSize(self) != SliceSize(other_vec)) << "Left and right hand side of slice addition have different sizes!" << std::endl;
        for(std::size_t i=0; i<SliceSize(self); ++i) self[i]+=other_vec[i];
        return self;
    }, py::is_operator(), py::keep_alive<0,1>())
    .def("__isub__", [](TSliceType& self, const TSliceType& other_vec)
    {
        KRATOS_ERROR_IF(SliceSize(self) != SliceSize(other_vec)) << "Left and right hand side of slice subtraction have different sizes!" << std::endl;
        for(std::size_t i=0; i<SliceSize(self); ++i) self[i]-=other_vec[i];
        return self;
    }, py::is_operator(), py::keep_alive<0,1>())
    .def("__mul__", [](TSliceType vec1, const value_type scalar)
    {
        for(std::size_t i=0; i<SliceSize(vec1); ++i) vec1[i]*=scalar;
        return vec1;
    }, py::is_operator(), py::keep_alive<0,1>())
    .def("__div__", [](TSliceType vec1, const value_type scalar)
    {
        for(std::size_t i=0; i<SliceSize(vec1); ++i) vec1[i]/=scalar;
        return vec1;
    }, py::is_operator(), py::keep_alive<0,1>())
    .def("__rmul__", [](TSliceType vec1, const value_type scalar)
    {
        for(std::size_t i=0; i<SliceSize(vec1); ++i) vec1[i]*=scalar;
        return vec1;
    }, py::is_operator(), py::keep_alive<0,1>())
    .def("__rdiv__", [](TSliceType vec1, const value_type scalar)
    {
        for(std::size_t i=0; i<SliceSize(vec1); ++i) vec1[i]/=scalar;
    }, py::is_operator())
    .def("__add__", [](const TSliceType& vec1, const TSliceType& vec2)
    {
        KRATOS_ERROR_IF(SliceSize(vec1) != SliceSize(vec2)) << "Left and right hand side of slice addition have different sizes!" << std::endl;
        Vector aux(SliceSize(vec1));
        for(std::size_t i=0; i<SliceSize(vec1); ++i) aux[i] = vec1[i] + vec2[i];
        return aux;
    }, py::is_operator())
    .def("__sub__", [](const TSliceType& vec1, const TSliceType& vec2)
    {
        KRATOS_ERROR_IF(SliceSize(vec1) != SliceSize(vec2)) << "Left and right hand side of slice subtraction have different sizes!" << std::endl;
        Vector aux(SliceSize(vec1));
        for(std::size_t i=0; i<SliceSize(vec1); ++i) aux[i] = vec1[i] - vec2[i];
        return aux;
    }, py::is_operator())
    .def("__setitem__", [](TSliceType& self, const unsigned int i, const value_type value)
    {
        self[i] = value;
    })
    .def("__getitem__", [](const TSliceType& self, const unsigned int i)
    {
        return self[i];
    } )
    .def("__setitem__", [](TSliceType &self, pybind11::slice this_slice, const TSliceType &value)
    {
        size_t start, stop, step, slicelength;
        if (!this_slice.compute(SliceSize(self), &start, &stop, &step, &slicelength))
            throw pybind11::error_already_set();
        if (slicelength != SliceSize(value))
            throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
        for (size_t i = 0; i < slicelength; ++i)
        {
            self[start] = value[i];
            start += step;
        }
    })
    .def("__setitem__", [](TSliceType &self, pybind11::slice this_slice, const Vector &value)
    {
        size_t start, stop, step, slicelength;
        if (!this_slice.compute(SliceSize(self), &start, &stop, &step, &slicelength))
            throw pybind11::error_already_set();
        if (slicelength != value.size())
            throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
        for (size_t i = 0; i < slicelength; ++i)
        {
            self[start] = value[i];
            start += step;
        }
    })
    .def("__iter__", [](TSliceType& self)
    {
        return py::make_iterator(self.begin(), self.end(), py::return_value_policy::reference_internal);
    }, py::keep_alive<0,1>() )
    .def("__str__", [](const TSliceType& self)
    {
        // uBLAS text format, so str() is the same under both backends
        std::stringstream ss;
        ss << "[" << SliceSize(self) << "](";
        for(std::size_t i=0; i<SliceSize(self); ++i)
        {
            if (i > 0) ss << ",";
            ss << self[i];
        }
        ss << ")";
        return ss.str();
    })
    ;
}

void  AddVectorToPython(pybind11::module& m)
{
    // Live slice view of the dynamic Vector (see PythonVectorSlice above);
    // registered before Vector so the Vector(VectorSlice) constructor and the
    // implicit conversion below can refer to it.
    CreateVectorSliceInterface<PythonVectorSlice>(m, "VectorSlice");

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
    vector_binder.def(py::init([](const PythonVectorSlice& input)
    {
        return Vector(input);
    }));
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
            GetContiguousDataPointer(self), // raw storage pointer for both backends, valid (nullptr) for an empty container
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
    py::implicitly_convertible<PythonVectorSlice, Vector>();

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
            GetContiguousDataPointer(self), // raw storage pointer for both backends, valid (nullptr) for an empty container
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
            GetContiguousDataPointer(self), // raw storage pointer for both backends, valid (nullptr) for an empty container
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
