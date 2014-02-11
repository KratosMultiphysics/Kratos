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
