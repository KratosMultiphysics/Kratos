//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            December 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(PROCESS_INFO_EXTENSIONS_H_INCLUDED)
#define PROCESS_INFO_EXTENSIONS_H_INCLUDED


#define DECLARE_HAS_THIS_TYPE_PROCESS_INFO \
template<class TVariable> \
static bool HasProcessInfo(TVariable const& rV, ProcessInfo::Pointer& p) \
{ \
	return p->Has(rV); \
}

#define DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON(TClassName) \
.def("HasProcessInfo", &TClassName::HasProcessInfo< Variable< TClassName > >).staticmethod("HasProcessInfo")

#define DECLARE_HAS_THIS_TYPE_PROCESS_INFO_PYTHON_AS_POINTER(TClassName) \
.def("HasProcessInfo", &TClassName::HasProcessInfo< Variable< TClassName::Pointer > >).staticmethod("HasProcessInfo")



#define DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO \
template<class TVariable> \
static void AddToProcessInfo(TVariable const& rV, typename TVariable::Type const& rValue, ProcessInfo::Pointer& p) \
{ \
	p->SetValue(rV, rValue); \
}

#define DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON(TClassName) \
.def("AddToProcessInfo", &TClassName::AddToProcessInfo< Variable< TClassName > >).staticmethod("AddToProcessInfo")

#define DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO_PYTHON_AS_POINTER(TClassName) \
.def("AddToProcessInfo", &TClassName::AddToProcessInfo< Variable< TClassName::Pointer > >).staticmethod("AddToProcessInfo")



#define DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO \
template<class TVariable> \
static typename TVariable::Type GetFromProcessInfo(TVariable const& rV, ProcessInfo::Pointer& p) \
{ \
	return p->GetValue(rV); \
}

#define DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON(TClassName) \
.def("GetFromProcessInfo", &TClassName::GetFromProcessInfo< Variable< TClassName > >).staticmethod("GetFromProcessInfo")

#define DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO_PYTHON_AS_POINTER(TClassName) \
.def("GetFromProcessInfo", &TClassName::GetFromProcessInfo< Variable< TClassName::Pointer > >).staticmethod("GetFromProcessInfo")



#endif // PROCESS_INFO_EXTENSIONS_H_INCLUDED
