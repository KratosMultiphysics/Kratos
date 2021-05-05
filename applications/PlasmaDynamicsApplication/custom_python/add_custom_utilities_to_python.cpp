/*
==============================================================================
KratosTestApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: G.Casas$
//   Date:                $Date: 2011-6-27 16:07:33 $
//   Revision:            $Revision: 1.2 $
//
//

// System includes

// External includes

// Project includes

#include "add_custom_utilities_to_python.h"


namespace Kratos{

namespace Python{

// typedef ModelPart::NodesContainerType::iterator PointIterator;
// typedef std::vector<array_1d<double, 3 > > ComponentVectorType;
// typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;

class VariableChecker{
public:
    // VariableChecker(){}
    // virtual ~VariableChecker(){}

    // template<class TDataType>
    // bool ModelPartHasNodalVariableOrNot(ModelPart& r_model_part, const Variable<TDataType>& rThisVariable)
    // {
    //     // return (r_model_part.GetNodalSolutionStepVariablesList()).Has(rThisVariable);
    // }
};

// template<class TDataType>
// bool ModelPartHasNodalVariableOrNot(VariableChecker& rChecker, ModelPart& rModelPart, Variable<TDataType> const& rThisVariable)
// {
//     // return rChecker.ModelPartHasNodalVariableOrNot(rModelPart, rThisVariable);
// }

namespace py = pybind11;



void  AddCustomUtilitiesToPython(pybind11::module& m){

    }

}  // namespace Python.

} // Namespace Kratos
