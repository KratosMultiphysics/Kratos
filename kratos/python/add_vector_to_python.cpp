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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "python/add_vector_to_python.h"
#include "python/vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"
#include "python/vector_vector_operator_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

template<class TContainerType>
struct UblasVectorModifierRenamed
{
    typedef typename TContainerType::size_type index_type;
    static void Resize(TContainerType& ThisContainer, typename TContainerType::size_type NewSize)
    {
        ThisContainer.resize(NewSize, true);
    }
    static void MoveSlice(TContainerType& ThisContainer, index_type Index, index_type From, index_type To)
    {
        if(Index > From)
        {
            ThisContainer.resize(ThisContainer.size() + Index - From, true);
            std::copy_backward(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index + To - From);
        }
        else
        {
            std::copy(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index);
            ThisContainer.resize(ThisContainer.size() + Index - From, true);
        }
    }
};


void  AddVectorToPython()
{

    ReadonlyVectorPythonInterface<zero_vector<double> >::CreateInterface("ZeroVector")
    .def(init<zero_vector<double>::size_type>())
//       .def(VectorScalarOperatorPython<zero_vector<double>, double, vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, unit_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, scalar_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, mapped_vector<double>, mapped_vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, compressed_vector<double>, compressed_vector<double> >())
//       .def(VectorVectorOperatorPython<zero_vector<double>, coordinate_vector<double>, coordinate_vector<double> >())
    ;

    ReadonlyVectorPythonInterface<unit_vector<double> >::CreateInterface("UnitVector")
    .def(init<unit_vector<double>::size_type, vector<double>::size_type>())
//       .def(VectorScalarOperatorPython<unit_vector<double>, double, vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, zero_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, scalar_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, mapped_vector<double>, mapped_vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, compressed_vector<double>, compressed_vector<double> >())
//       .def(VectorVectorOperatorPython<unit_vector<double>, coordinate_vector<double>, coordinate_vector<double> >())
    ;

    ReadonlyVectorPythonInterface<scalar_vector<double> >::CreateInterface("ScalarVector")
    .def(init<scalar_vector<double>::size_type, scalar_vector<double>::value_type>())
//       .def(VectorScalarOperatorPython<scalar_vector<double>, double, vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, zero_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, unit_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, mapped_vector<double>, mapped_vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, compressed_vector<double>, compressed_vector<double> >())
//       .def(VectorVectorOperatorPython<scalar_vector<double>, coordinate_vector<double>, coordinate_vector<double> >())
    ;

    VectorPythonInterface<vector<double>, UblasVectorModifierRenamed<vector<double> > >::CreateInterface("Vector")
    .def(init<vector<double>::size_type>())
    .def(init<vector_expression<vector<double> > >())
    .def(VectorScalarOperatorPython<vector<double>, double, vector<double> >())
//       .def(VectorVectorOperatorPython<vector<double>, zero_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<vector<double>, unit_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<vector<double>, scalar_vector<double>, vector<double> >())
//       .def(VectorVectorOperatorPython<vector<double>, mapped_vector<double>, vector<double> >())
   ;
    
      VectorPythonInterface<vector<int>, UblasVectorModifierRenamed<vector<int> > >::CreateInterface("IntegerVector")
      .def(init<vector<int>::size_type>())
   ;

}
}  // namespace Python.

} // Namespace Kratos

