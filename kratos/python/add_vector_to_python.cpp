//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//



// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"
#include "python/add_vector_to_python.h"

namespace Kratos
{

namespace Python
{

    namespace py = pybind11;

    template< typename TVectorType > py::class_< TVectorType > CreateVectorInterface(pybind11::module& m, std::string Name )
    {
        py::class_< TVectorType, std::shared_ptr<TVectorType> > binder(m,Name.c_str());
        binder.def(py::init<>());

        //binder.def(py::init<std::TVectorType& >())
        binder.def("Size", [](const TVectorType& self){return self.size();} );
        binder.def("Resize", [](TVectorType& self, const typename TVectorType::size_type  new_size){if(self.size() != new_size) self.resize(new_size, false);} );
        binder.def("__len__", [](const TVectorType& self){return self.size();} );

        //operating on the object itself, +=, -=, *=, etc
        binder.def("__iadd__", [](TVectorType& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]+=scalar; return self;}, py::is_operator());
        binder.def("__isub__", [](TVectorType& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]-=scalar; return self;}, py::is_operator());
        binder.def("__imul__", [](TVectorType& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]*=scalar; return self;}, py::is_operator());
        binder.def("__itruediv__", [](TVectorType& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]/=scalar; return self;}, py::is_operator());

        binder.def("__iadd__", [](TVectorType& self, const TVectorType& other_vec){noalias(self) += other_vec; return self;}, py::is_operator());
        binder.def("__isub__", [](TVectorType& self, const TVectorType& other_vec){noalias(self) -= other_vec; return self; }, py::is_operator());

        //returning a different object
//         binder.def("__add__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]+=scalar; return vec1;}, py::is_operator());
//         binder.def("__sub__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]-=scalar; return vec1;}, py::is_operator());
         binder.def("__mul__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar; return vec1;}, py::is_operator());
         binder.def("__div__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar; return vec1;}, py::is_operator());
         binder.def("__truediv__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar; return vec1;}, py::is_operator());
//         binder.def("__radd__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]+=scalar; return vec1;}, py::is_operator());
//         binder.def("__rsub__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]-=scalar; return vec1;}, py::is_operator());
         binder.def("__rmul__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar; return vec1;}, py::is_operator());
         binder.def("__rdiv__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar;}, py::is_operator());
        binder.def("__add__", [](const TVectorType& vec1, const TVectorType& vec2){Vector aux(vec1); aux += vec2; return aux;}, py::is_operator());
        binder.def("__sub__", [](const TVectorType& vec1, const TVectorType& vec2){Vector aux(vec1); aux -= vec2; return aux;}, py::is_operator());

        //access operators
        binder.def("__setitem__", [](TVectorType& self, const unsigned int i, const typename TVectorType::value_type value){self[i] = value;} );
        binder.def("__getitem__", [](const TVectorType& self, const unsigned int i){return self[i];} );

        //access using slices
        binder.def("__setitem__", [](TVectorType &self, pybind11::slice this_slice, const TVectorType &value) {
            size_t start, stop, step, slicelength;
            if (!this_slice.compute(self.size(), &start, &stop, &step, &slicelength))
                throw pybind11::error_already_set();
            if (slicelength != value.size())
                throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
            for (size_t i = 0; i < slicelength; ++i) {
                self[start] = value[i]; start += step;
            }
        });
    #ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        binder.def("__getitem__", [](TVectorType &self, pybind11::slice this_slice) -> AMatrix::SubVector<TVectorType> {
          size_t start, stop, step, slicelength;
          if (!this_slice.compute(self.size(), &start, &stop, &step, &slicelength))
            throw pybind11::error_already_set();
          KRATOS_ERROR_IF(step != 1) << "The AMatrix only supports continuous slices with step == 1" << std::endl;
          AMatrix::SubVector<TVectorType> sliced_self(self, start, slicelength);
          return sliced_self;
        });
        binder.def("fill", [](TVectorType& self, const typename TVectorType::value_type value) { self.fill(value); });
        binder.def("norm_1", [](TVectorType& self) { return norm_1(self); });
        binder.def("norm_2", [](TVectorType& self) { return norm_2(self); });
    #else
        binder.def("__getitem__", [](TVectorType &self, pybind11::slice this_slice) -> boost::numeric::ublas::vector_slice<TVectorType> {
          size_t start, stop, step, slicelength;
          if (!this_slice.compute(self.size(), &start, &stop, &step, &slicelength))
            throw pybind11::error_already_set();
          boost::numeric::ublas::slice ublas_slice(start, step, slicelength);
          boost::numeric::ublas::vector_slice<TVectorType> sliced_self(self, ublas_slice);
          return sliced_self;
        });
        binder.def("fill", [](TVectorType& self, const typename TVectorType::value_type value) { noalias(self) = TVectorType(self.size(),value); });
        binder.def("norm_1", [](TVectorType& self) { return boost::numeric::ublas::norm_1(self); });
        binder.def("norm_2", [](TVectorType& self) { return boost::numeric::ublas::norm_2(self); });
    #endif // KRATOS_USE_AMATRIX

        binder.def("__iter__", [](TVectorType& self){ return py::make_iterator(self.begin(), self.end(), py::return_value_policy::reference_internal); } , py::keep_alive<0,1>() ) ;
        binder.def("__str__", PrintObject<TVectorType>);

        return std::move(binder);
    }

    template< std::size_t TSize >
    void CreateArray1DInterface(pybind11::module& m, const std::string& Name )
    {
        auto binder = CreateVectorInterface< array_1d<double,TSize> >(m,Name);
        binder.def(py::init( [](double value){
            array_1d<double,TSize> tmp;
            for(std::size_t i=0; i < TSize; ++i)
                tmp[i] = value;
            return tmp;
        }));
        binder.def(py::init( [](const Vector& input){
            KRATOS_ERROR_IF(input.size() != TSize)
            << "Attempting to initialize an array_1d<double," << TSize << "> from a Vector of size "
            << input.size() << ". Input should have size " << TSize <<"." << std::endl;

            array_1d<double,TSize> tmp(input);
            return tmp;
        }));
        binder.def(py::init<array_1d<double,TSize>>());
        binder.def(py::init( [](const py::list& input){
            KRATOS_ERROR_IF(input.size() != TSize)
            << "Attempting to initialize an array_1d<double," << TSize << "> from a Python list of size "
            << input.size() << ". Input should have size " << TSize <<"." << std::endl;

            array_1d<double,TSize> tmp;
            for(std::size_t i=0; i<TSize; ++i) {
                tmp[i] = py::cast<double>(input[i]);
            }
            return tmp;
        }));


        py::implicitly_convertible<py::list, array_1d<double,TSize>>();
        py::implicitly_convertible<Vector, array_1d<double,TSize>>();
    }

    void  AddVectorToPython(pybind11::module& m)
    {
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        using VectorSlice = AMatrix::SubVector<Vector>;
#else
        typedef boost::numeric::ublas::vector_slice<Vector> VectorSlice;
#endif // KRATOS_USE_AMATRIX
        py::class_< VectorSlice >(m, "VectorSlice")
        .def("Size", [](const VectorSlice& self){return self.size();} )
        .def("__len__", [](const VectorSlice& self){return self.size();} )
        .def("__iadd__", [](VectorSlice& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]+=scalar; return self;}, py::is_operator())
        .def("__isub__", [](VectorSlice& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]-=scalar; return self;}, py::is_operator())
        .def("__imul__", [](VectorSlice& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]*=scalar; return self;}, py::is_operator())
        .def("__itruediv__", [](VectorSlice& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]/=scalar; return self;}, py::is_operator())
        .def("__iadd__", [](VectorSlice& self, const VectorSlice& other_vec){noalias(self) += other_vec; return self;}, py::is_operator())
        .def("__isub__", [](VectorSlice& self, const VectorSlice& other_vec){noalias(self) -= other_vec; return self; }, py::is_operator())
        .def("__mul__", [](VectorSlice vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar; return vec1;}, py::is_operator())
        .def("__div__", [](VectorSlice vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar; return vec1;}, py::is_operator())
        .def("__rmul__", [](VectorSlice vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar; return vec1;}, py::is_operator())
        .def("__rdiv__", [](VectorSlice vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar;}, py::is_operator())
        .def("__add__", [](const VectorSlice& vec1, const VectorSlice& vec2){Vector aux(vec1); aux += vec2; return aux;}, py::is_operator())
        .def("__sub__", [](const VectorSlice& vec1, const VectorSlice& vec2){Vector aux(vec1); aux -= vec2; return aux;}, py::is_operator())
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
            .def("__setitem__", [](VectorSlice& self, const unsigned int i, const typename VectorSlice::data_type value) {self[i] = value; })
#else
            .def("__setitem__", [](VectorSlice& self, const unsigned int i, const typename VectorSlice::value_type value) {self[i] = value; })
#endif // KRATOS_USE_AMATRIX
        .def("__getitem__", [](const VectorSlice& self, const unsigned int i){return self[i];} )
        .def("__setitem__", [](VectorSlice &self, pybind11::slice this_slice, const VectorSlice &value) {
            size_t start, stop, step, slicelength;
            if (!this_slice.compute(self.size(), &start, &stop, &step, &slicelength))
                throw pybind11::error_already_set();
            if (slicelength != value.size())
                throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
            for (size_t i = 0; i < slicelength; ++i) {
                self[start] = value[i]; start += step;
            }
        })
        .def("__setitem__", [](VectorSlice &self, pybind11::slice this_slice, const Vector &value) {
            size_t start, stop, step, slicelength;
            if (!this_slice.compute(self.size(), &start, &stop, &step, &slicelength))
                throw pybind11::error_already_set();
            if (slicelength != value.size())
                throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
            for (size_t i = 0; i < slicelength; ++i) {
                self[start] = value[i]; start += step;
            }
        })
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        .def("__iter__", [](VectorSlice& self){ return py::make_iterator(self.data(), self.data() + self.size(), py::return_value_policy::reference_internal); } , py::keep_alive<0,1>() )
#else
        .def("__iter__", [](VectorSlice& self){ return py::make_iterator(self.begin(), self.end(), py::return_value_policy::reference_internal); } , py::keep_alive<0,1>() )
#endif // ifdef KRATOS_USE_AMATRIX
        .def("__str__", PrintObject<VectorSlice>)
        ;

        auto vector_binder = CreateVectorInterface<Vector>(m, "Vector");
        vector_binder.def(py::init<typename Vector::size_type>());
        vector_binder.def(py::init<typename Vector::size_type, double>());
        vector_binder.def(py::init<Vector>());
        vector_binder.def(py::init<array_1d<double,3>>());
        vector_binder.def(py::init( [](const py::list& input){
                                Vector tmp(input.size());
                                for(unsigned int i=0; i<tmp.size(); ++i)
                                    tmp[i] = py::cast<double>(input[i]);
                                return tmp;
                                }));
        py::implicitly_convertible<py::list, Vector>();
        py::implicitly_convertible<array_1d<double,3>, Vector>();

        auto int_vector_binder = CreateVectorInterface<DenseVector<int>>(m, "DenseVectorInt");
        int_vector_binder.def(py::init<typename DenseVector<int>::size_type>());
        int_vector_binder.def(py::init<typename DenseVector<int>::size_type, int>());
        int_vector_binder.def(py::init<DenseVector<int>>());
        int_vector_binder.def(py::init( [](const py::list& input){
                                DenseVector<int> tmp(input.size());
                                for(unsigned int i=0; i<tmp.size(); ++i)
                                    tmp[i] = py::cast<int>(input[i]);
                                return tmp;
                                }));
        py::implicitly_convertible<py::list, DenseVector<int>>();

        CreateArray1DInterface< 3 >(m,"Array3");
        CreateArray1DInterface< 4 >(m,"Array4");
        CreateArray1DInterface< 6 >(m,"Array6");
        CreateArray1DInterface< 9 >(m,"Array9");

    }
}  // namespace Python.

} // Namespace Kratos

