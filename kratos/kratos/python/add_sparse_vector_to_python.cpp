// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "python/add_sparse_vector_to_python.h"
#include "python/vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"
#include "python/vector_vector_operator_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

template<class TContainerType>
struct UblasSparseVectorModifier
{
    typedef typename TContainerType::size_type index_type;
    typedef typename TContainerType::value_type data_type;

    static void Resize(TContainerType& ThisContainer, typename TContainerType::size_type NewSize)
    {
        index_type size = std::min(ThisContainer.size(), NewSize);
        std::vector<std::pair<index_type, data_type> > temp;
        data_type value;

        for(index_type i = 0 ; i < size ; i++)
            if((value = ThisContainer(i)) != data_type())
                temp.push_back(std::pair<index_type, data_type>(i,value));

        ThisContainer.clear(); // There is no way to know which resize hold the data and which not. So better to make it certain! :-)
        ThisContainer.resize(NewSize, false);

        for(typename std::vector<std::pair<index_type, data_type> >::iterator j = temp.begin() ; j != temp.end() ; j++)
            ThisContainer.insert_element(j->first, j->second);

    }

    static void MoveSlice(TContainerType& ThisContainer, index_type Index, index_type From, index_type To)
    {
        index_type i;
        index_type size = std::min(Index, From);
        index_type new_size = ThisContainer.size() + Index - From;
        std::vector<std::pair<index_type, data_type> > temp;
        data_type value;

        for(i = 0 ; i < size ; i++)
            if((value = ThisContainer(i)) != data_type())
                temp.push_back(std::pair<index_type, data_type>(i,value));

        for(i = From ; i < To ; i++)
            if((value = ThisContainer(i)) != data_type())
                temp.push_back(std::pair<index_type, data_type>(i + Index - From,value));

        ThisContainer.clear(); // There is no way to know which resize hold the data and which not. So better to make it certain! :-)
        ThisContainer.resize(new_size, false);

        for(typename std::vector<std::pair<index_type, data_type> >::iterator j = temp.begin() ; j != temp.end() ; j++)
            ThisContainer.insert_element(j->first, j->second);
    }
};


void  AddSparseVectorToPython()
{

    VectorPythonInterface<mapped_vector<double>, UblasSparseVectorModifier<mapped_vector<double> > >::CreateInterface("SparseVector")
    .def(init<mapped_vector<double>::size_type>())
    .def("NonZeros", &mapped_vector<double>::nnz)
    .def(VectorScalarOperatorPython<mapped_vector<double>, double, mapped_vector<double> >())
    .def(VectorVectorOperatorPython<mapped_vector<double>, zero_vector<double>, mapped_vector<double> >())
    .def(VectorVectorOperatorPython<mapped_vector<double>, unit_vector<double>, mapped_vector<double> >())
    .def(VectorVectorOperatorPython<mapped_vector<double>, scalar_vector<double>, mapped_vector<double> >())
    .def(VectorVectorOperatorPython<mapped_vector<double>, vector<double>, vector<double> >())
    ;

    VectorPythonInterface<compressed_vector<double>, UblasSparseVectorModifier<compressed_vector<double> > >::CreateInterface("CompressedVector")
    .def(init<compressed_vector<double>::size_type>())
    .def("NonZeros", &compressed_vector<double>::nnz)
    .def(VectorScalarOperatorPython<compressed_vector<double>, double, compressed_vector<double> >())
    .def(VectorVectorOperatorPython<compressed_vector<double>, zero_vector<double>, compressed_vector<double> >())
    .def(VectorVectorOperatorPython<compressed_vector<double>, unit_vector<double>, compressed_vector<double> >())
    .def(VectorVectorOperatorPython<compressed_vector<double>, scalar_vector<double>, compressed_vector<double> >())
    .def(VectorVectorOperatorPython<compressed_vector<double>, vector<double>, vector<double> >())
    .def(VectorVectorOperatorPython<compressed_vector<double>, mapped_vector<double>, compressed_vector<double> >())
    .def(VectorVectorOperatorPython<compressed_vector<double>, coordinate_vector<double>, compressed_vector<double> >())
    ;

    VectorPythonInterface<coordinate_vector<double>, UblasSparseVectorModifier<coordinate_vector<double> > >::CreateInterface("CoordinateVector")
    .def(init<coordinate_vector<double>::size_type>())
    .def("NonZeros", &coordinate_vector<double>::nnz)
    .def(VectorScalarOperatorPython<coordinate_vector<double>, double, coordinate_vector<double> >())
    .def(VectorVectorOperatorPython<coordinate_vector<double>, zero_vector<double>, coordinate_vector<double> >())
    .def(VectorVectorOperatorPython<coordinate_vector<double>, unit_vector<double>, coordinate_vector<double> >())
    .def(VectorVectorOperatorPython<coordinate_vector<double>, scalar_vector<double>, coordinate_vector<double> >())
    .def(VectorVectorOperatorPython<coordinate_vector<double>, vector<double>, vector<double> >())
    ;




}

}  // namespace Python.

} // Namespace Kratos

