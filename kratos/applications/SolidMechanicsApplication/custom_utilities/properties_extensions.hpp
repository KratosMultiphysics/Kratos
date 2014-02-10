/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-06-06 10:37:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(PROPERTIES_EXTENSIONS_H_INCLUDED)
#define PROPERTIES_EXTENSIONS_H_INCLUDED

#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES \
template<class TVariableType> \
static void AddToProperties(TVariableType const& rV, typename TVariableType::Type const& rValue, Properties::Pointer& p) \
{ \
	p->SetValue(rV, rValue); \
}

#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON(TClassName) \
.def("AddToProperties", &TClassName::AddToProperties< Variable< TClassName > >).staticmethod("AddToProperties")

#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON_AS_POINTER(TClassName) \
.def("AddToProperties", &TClassName::AddToProperties< Variable< TClassName::Pointer > >).staticmethod("AddToProperties")



#define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES \
template<class TVariableType> \
static typename TVariableType::Type GetFromProperties(TVariableType const& rV, Properties::Pointer& p) \
{ \
	return p->GetValue(rV); \
}

#define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON(TClassName) \
.def("GetFromProperties", &TClassName::GetFromProperties< Variable< TClassName > >).staticmethod("GetFromProperties")

#define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON_AS_POINTER(TClassName) \
.def("GetFromProperties", &TClassName::GetFromProperties< Variable< TClassName::Pointer > >).staticmethod("GetFromProperties")



#endif // PROPERTIES_EXTENSIONS_H_INCLUDED
