//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//



// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"
#include "python/add_vector_to_python.h"

namespace Kratos
{

namespace Python
{

    using namespace pybind11;

    template< typename TVectorType > class_< TVectorType > CreateVectorInterface(pybind11::module& m, std::string Name )
        {
        class_< TVectorType > binder(m,Name.c_str());

        //binder.def(init<std::TVectorType& >())
        binder.def("Size", [](const TVectorType& self){return self.size();} );
        binder.def("Resize", [](TVectorType& self, unsigned int new_size){if(self.size() != new_size) self.resize(new_size, false);} );
        binder.def("__len__", [](const TVectorType& self){return self.size();} );
        
        //operating on the object itself, +=, -=, *=, etc
        binder.def("__iadd__", [](TVectorType& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]+=scalar; return self;}, is_operator());
        binder.def("__isub__", [](TVectorType& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]-=scalar; return self;}, is_operator());
        binder.def("__imul__", [](TVectorType& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]*=scalar; return self;}, is_operator());
        binder.def("__itruediv__", [](TVectorType& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]/=scalar; return self;}, is_operator());
        
        binder.def("__iadd__", [](TVectorType& self, const TVectorType& other_vec){noalias(self) += other_vec; return self;}, is_operator());
        binder.def("__isub__", [](TVectorType& self, const TVectorType& other_vec){noalias(self) -= other_vec; return self; }, is_operator());
       
        //returning a different object
//         binder.def("__add__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]+=scalar; return vec1;}, is_operator());
//         binder.def("__sub__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]-=scalar; return vec1;}, is_operator());
//         binder.def("__mul__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar; return vec1;}, is_operator());
//         binder.def("__div__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar; return vec1;}, is_operator());
//         binder.def("__radd__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]+=scalar; return vec1;}, is_operator());
//         binder.def("__rsub__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]-=scalar; return vec1;}, is_operator());
//         binder.def("__rmul__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]*=scalar; return vec1;}, is_operator());
//         binder.def("__rdiv__", [](TVectorType vec1, const double scalar){for(unsigned int i=0; i<vec1.size(); ++i) vec1[i]/=scalar;}, is_operator());
        binder.def("__add__", [](TVectorType& vec1, const TVectorType& vec2){Vector aux(vec1); aux += vec2; return aux;}, is_operator());
        binder.def("__sub__", [](TVectorType& vec1, const TVectorType& vec2){Vector aux(vec1); aux -= vec2; return aux;}, is_operator());
         
        //access operators
        binder.def("__setitem__", [](TVectorType& self, const unsigned int i, const double value){self[i] = value;} );
        binder.def("__getitem__", [](const TVectorType& self, const unsigned int i){return self[i];} );
        
        binder.def("__iter__", [](TVectorType& self){ return make_iterator(self.begin(), self.end(), return_value_policy::reference_internal); } , keep_alive<0,1>() ) ;
        binder.def("__repr__", [](const TVectorType& self)
                { std::stringstream out;
                out << ( self );
                  return out.str();
                }); 
        
        return binder;
        }

    void  AddVectorToPython(pybind11::module& m)
    {
        auto vector_binder = CreateVectorInterface<Vector>(m, "Vector");
        vector_binder.def(init<typename Vector::size_type>());
        vector_binder.def(init<typename Vector::size_type, double>());
        vector_binder.def(init<Vector>());
        vector_binder.def(init<array_1d<double,3>>());
        vector_binder.def(init( [](const list& input){
                                Vector* tmp = new Vector(input.size());
                                for(unsigned int i=0; i<tmp->size(); ++i)
                                    (*tmp)[i] = cast<double>(input[i]);
                                return tmp;
                                }));
        
        auto array3_binder = CreateVectorInterface< Kratos::array_1d<double,3> >(m, "Array3");
        array3_binder.def(init<>());
        array3_binder.def(init( [](double value){
                                array_1d<double,3>* tmp = new array_1d<double,3>();
                                for(unsigned int i=0; i<tmp->size(); ++i)
                                    (*tmp)[i] = value;
                                return tmp;
                                }));
    }
}  // namespace Python.

} // Namespace Kratos

