// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Massimo Petracca
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
