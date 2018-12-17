//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                 October 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(PROPERTIES_EXTENSIONS_H_INCLUDED)
#define PROPERTIES_EXTENSIONS_H_INCLUDED

#if !defined(DECLARE_ADD_THIS_TYPE_TO_PROPERTIES)
#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES                             \
  template<class TVariable>                                         \
  static void AddToProperties(TVariable const& rV, typename TVariable::Type const& rValue, Properties::Pointer& p) \
  {                                                                     \
    p->SetValue(rV, rValue);                                            \
  }

#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON(TClassName)          \
  .def_static("AddToProperties", &TClassName::AddToProperties< Variable< TClassName > >)

#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON_AS_POINTER(TClassName) \
  .def_static("AddToProperties", &TClassName::AddToProperties< Variable< TClassName::Pointer > >)
#endif

#if !defined(DECLARE_GET_THIS_TYPE_FROM_PROPERTIES)
#define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES                           \
  template<class TVariable>                                         \
  static typename TVariable::Type GetFromProperties(TVariable const& rV, Properties::Pointer& p) \
  {                                                                     \
    return p->GetValue(rV);                                             \
  }

#define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON(TClassName)        \
  .def_static("GetFromProperties", &TClassName::GetFromProperties< Variable< TClassName > >)

#define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON_AS_POINTER(TClassName) \
  .def_static("GetFromProperties", &TClassName::GetFromProperties< Variable< TClassName::Pointer > >)
#endif


#endif // PROPERTIES_EXTENSIONS_H_INCLUDED
