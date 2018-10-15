//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                 October 2018 $
//   Revision:            $Revision:                      0.0 $
//
//


#if !defined(KRATOS_TABLE_KEY_VARIABLES_H_INCLUDED )
#define  KRATOS_TABLE_KEY_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "custom_utilities/properties_extensions.hpp"

// #if !defined(PROPERTIES_EXTENSIONS_H_INCLUDED)
// #define PROPERTIES_EXTENSIONS_H_INCLUDED

// #if !defined(DECLARE_ADD_THIS_TYPE_TO_PROPERTIES)
// #define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES                             \
//   template<class TVariable>                                             \
//   static void AddToProperties(TVariable const& rV, typename TVariable::Type const& rValue, Properties::Pointer& p) \
//   {                                                                     \
//     p->SetValue(rV, rValue);                                            \
//   }

// #define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON(TClassName)          \
//   .def_static("AddToProperties", &TClassName::AddToProperties< Variable< TClassName > >)

// #define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES_PYTHON_AS_POINTER(TClassName) \
//   .def_static("AddToProperties", &TClassName::AddToProperties< Variable< TClassName::Pointer > >)
// #endif

// #if !defined(DECLARE_GET_THIS_TYPE_FROM_PROPERTIES)
// #define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES                           \
//   template<class TVariable>                                             \
//   static typename TVariable::Type GetFromProperties(TVariable const& rV, Properties::Pointer& p) \
//   {                                                                     \
//     return p->GetValue(rV);                                             \
//   }

// #define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON(TClassName)        \
//   .def_static("GetFromProperties", &TClassName::GetFromProperties< Variable< TClassName > >)

// #define DECLARE_GET_THIS_TYPE_FROM_PROPERTIES_PYTHON_AS_POINTER(TClassName) \
//   .def_static("GetFromProperties", &TClassName::GetFromProperties< Variable< TClassName::Pointer > >)
// #endif


// #endif // PROPERTIES_EXTENSIONS_H_INCLUDED


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// TableKeyVariables 
/**
   Configures the variables set in the tables of the properties
*/
template<class TArgumentType, class TResultType = TArgumentType>
class TableKeyVariables
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of TableKeyVariables
  KRATOS_CLASS_POINTER_DEFINITION(TableKeyVariables);

#ifdef  _WIN32 // work around for windows int64_t error
  typedef __int64 int64_t;
#endif

  typedef Variable<TArgumentType> XVariableType;

  typedef Variable<TResultType> YVariableType;
  
  typedef std::pair<XVariableType*, YVariableType*> TableVariablesType;
  
  typedef std::vector<TableVariablesType> TableVariablesContainerType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  TableKeyVariables() {}

  /// Copy constructor.
  TableKeyVariables(const TableKeyVariables& rOther)
      :mKeys(rOther.mKeys)
      ,mData(rOther.mData)
  {}

  /// Clone.
  TableKeyVariables::Pointer Clone() const
  {
    return Kratos::make_shared< TableKeyVariables >(*this);
  }
  
  /// Destructor.
  ~TableKeyVariables() {}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{
  
  void SetTable(const XVariableType& XVariable, const YVariableType& YVariable)
  {
    mKeys.push_back(Key(XVariable.Key(), YVariable.Key()));
    mData.push_back(TableVariablesType(&XVariable,&YVariable));
  }

  XVariableType& GetXVariable(std::size_t rTableKey)
  {
    return *(mData[rTableKey].first);
  }
  
  YVariableType& GetYVariable(std::size_t rTableKey)
  {
    return *(mData[rTableKey].second);
  }
  
  ///@}
  ///@name Access
  ///@{

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  std::string Info() const
  {
    return "TableKeyVariables";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const
  {
    rOStream << "TableKeyVariables";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const
  {
    rOStream << "This table keys contain " << mData.size() << " variables";
  }

  ///@}
  ///@name Friends
  ///@{

  ///@}

 protected:
  ///@name Protected static Member Variables
  ///@{

  ///@}
  ///@name Protected member Variables
  ///@{

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

  ///@}
  ///@name Protected  Access
  ///@{

  ///@}
  ///@name Protected Inquiry
  ///@{

  ///@}
  ///@name Protected LifeCycle
  ///@{

  ///@}

 private:
  ///@name Static Member Variables
  ///@{

  ///@}
  ///@name Member Variables
  ///@{

  std::vector<std::size_t> mKeys;

  TableVariablesContainerType mData;
   
  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  // must be the same method as in properties.h
  int64_t Key(std::size_t XKey, std::size_t YKey) const
  {
    int64_t result_key = XKey;
    result_key = result_key << 32;
    result_key |= YKey; // I know that the key is less than 2^32 so I don't need zeroing the upper part
    return result_key;
  }
  
  ///@}
  ///@name Serialization
  ///@{

  friend class Serializer;

  void save(Serializer& rSerializer) const
  {
    rSerializer.save("mKeys", mKeys);
    rSerializer.save("mData", mData);
  }

  void load(Serializer& rSerializer)
  {
    rSerializer.load("mKeys", mKeys);
    rSerializer.load("mData", mData);
  }

  ///@}
  ///@name Private  Access
  ///@{


  ///@}
  ///@name Private Inquiry
  ///@{


  ///@}
  ///@name Un accessible methods
  ///@{
  
  ///@}

 public:

  DECLARE_ADD_THIS_TYPE_TO_PROPERTIES
  DECLARE_GET_THIS_TYPE_FROM_PROPERTIES

}; // Class TableKeyVariables

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<class TArgumentType, class TResultType>
inline std::istream& operator >> (std::istream& rIStream,
                                  TableKeyVariables<TArgumentType, TResultType>& rThis);

/// output stream function
template<class TArgumentType, class TResultType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TableKeyVariables<TArgumentType, TResultType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}


}  // namespace Kratos.

#endif // KRATOS_TABLE_KEY_VARIABLES_H_INCLUDED  defined
