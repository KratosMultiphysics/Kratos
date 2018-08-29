//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(PROPERTIES_EXTENSIONS_H_INCLUDED)
#define PROPERTIES_EXTENSIONS_H_INCLUDED

#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES                             \
  template<class TVariableType>                                         \
  static void AddToProperties(TVariableType const& rV, typename TVariableType::Type const& rValue, Properties::Pointer& p) \
  {                                                                     \
    p->SetValue(rV, rValue);                                            \
  }

#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON(TClassName)          \
  .def_static("AddToProperties", &TClassName::AddToProperties< Variable< TClassName > >)

#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON_AS_POINTER(TClassName) \
  .def_static("AddToProperties", &TClassName::AddToProperties< Variable< TClassName::Pointer > >)



#define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES                           \
  template<class TVariableType>                                         \
  static typename TVariableType::Type GetFromProperties(TVariableType const& rV, Properties::Pointer& p) \
  {                                                                     \
    return p->GetValue(rV);                                             \
  }

#define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON(TClassName)        \
  .def_static("GetFromProperties", &TClassName::GetFromProperties< Variable< TClassName > >)

#define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON_AS_POINTER(TClassName) \
  .def_static("GetFromProperties", &TClassName::GetFromProperties< Variable< TClassName::Pointer > >)



#endif // PROPERTIES_EXTENSIONS_H_INCLUDED
