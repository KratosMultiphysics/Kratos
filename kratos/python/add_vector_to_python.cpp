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
#include "python/add_vector_to_python.h"
// #include "python/vector_python_interface.h"
// #include "python/vector_scalar_operator_python.h"
// #include "python/vector_vector_operator_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

// template<class TContainerType>
// struct UblasVectorModifierRenamed
// {
//     typedef typename TContainerType::size_type index_type;
//     static void Resize(TContainerType& ThisContainer, typename TContainerType::size_type NewSize)
//     {
//         ThisContainer.resize(NewSize, true);
//     }
//     static void MoveSlice(TContainerType& ThisContainer, index_type Index, index_type From, index_type To)
//     {
//         if(Index > From)
//         {
//             ThisContainer.resize(ThisContainer.size() + Index - From, true);
//             std::copy_backward(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index + To - From);
//         }
//         else
//         {
//             std::copy(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index);
//             ThisContainer.resize(ThisContainer.size() + Index - From, true);
//         }
//     }
// };


void  AddVectorToPython(pybind11::module& m)
{

    class_< vector<double> >(m,"Vector")
    .def(init<vector<double>::size_type>())
    .def(init<vector<double>::size_type, double>())
    .def("__add__", [](Vector& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]+=scalar;} )
    .def("__sub__", [](Vector& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]-=scalar;})
    .def("__mul__", [](Vector& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]*=scalar;})
    .def("__div__", [](Vector& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]/=scalar;})
    .def("__radd__", [](Vector& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]+=scalar;})
    .def("__rsub__", [](Vector& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]-=scalar;})
    .def("__rmul__", [](Vector& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]*=scalar;})
    .def("__rdiv__", [](Vector& self, const double scalar){for(unsigned int i=0; i<self.size(); ++i) self[i]/=scalar;})
    .def("__add__", [](Vector& self, const Vector& other_vec){noalias(self) += other_vec; } )
    .def("__sub__", [](Vector& self, const Vector& other_vec){noalias(self) -= other_vec; } )
    ;
    

}
}  // namespace Python.

} // Namespace Kratos

