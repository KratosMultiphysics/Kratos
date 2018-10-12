//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                 October 2018 $
//   Revision:            $Revision:                      0.0 $
//
//


#if !defined(KRATOS_PROPERTIES_LAYOUT_H_INCLUDED )
#define  KRATOS_PROPERTIES_LAYOUT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/properties.h"


#if !defined(PROPERTIES_EXTENSIONS_H_INCLUDED)
#define PROPERTIES_EXTENSIONS_H_INCLUDED

#if !defined(DECLARE_ADD_THIS_TYPE_TO_PROPERTIES)
#define DECLARE_ADD_THIS_TYPE_TO_PROPERTIES                             \
  template<class TVariable>                                             \
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
  template<class TVariable>                                             \
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

/// PropertiesLayout 
/**
   Configures the propeties data and supplies it propertly according to problem configuration
*/
class PropertiesLayout
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of PropertiesLayout
  KRATOS_CLASS_POINTER_DEFINITION(PropertiesLayout);

#ifdef  _WIN32 // work around for windows int64_t error
  typedef __int64 int64_t;
#endif

  /// Type of container used for variables
  typedef DataValueContainer ContainerType;

  typedef Node<3> NodeType;

  typedef Geometry<NodeType> GeometryType;

  typedef NodeType::IndexType IndexType;

  /// Type of container used for tables
  typedef Table<double> TableType;

  typedef std::unordered_map<std::size_t, TableType> TablesContainerType; 

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  PropertiesLayout()
  {
    mpData = nullptr;
    mpTables = nullptr;
  }

  /// Copy constructor.
  PropertiesLayout(const PropertiesLayout& rOther)
      :mpData(rOther.mpData)
      ,mpTables(rOther.mpTables)
      ,mTableArguments(rOther.mTableArguments)
  {}

  /// Clone.
  PropertiesLayout::Pointer Clone() const
  {
    return Kratos::make_shared< PropertiesLayout >(*this);
  }
  
  /// Destructor.
  ~PropertiesLayout() {}

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  template<class TVariableType>
  typename TVariableType::Type& operator()(const TVariableType& rV)
  {
    return GetValue(rV);
  }

  template<class TVariableType>
  typename TVariableType::Type const& operator()(const TVariableType& rV) const
  {
    return GetValue(rV);
  }

  template<class TVariableType>
  typename TVariableType::Type& operator[](const TVariableType& rV)
  {
    return GetValue(rV);
  }

  template<class TVariableType>
  typename TVariableType::Type const& operator[](const TVariableType& rV) const
  {
    return GetValue(rV);
  }

  ///@}
  ///@name Operations
  ///@{

  void Configure(const Properties& rProperties, const GeometryType& rGeometry, const Vector& rShapeFunctions)
  {
    mpData = &(rProperties.Data());
    mpTables = &(rProperties.Tables());

    std::size_t max_size = 0;
    for(auto it = mpTables->begin(); it != mpTables->end(); ++it)
      if( it->first > max_size )
        max_size = it->first;

    mTableArguments.resize(max_size);

    for(auto it = mpTables->begin(); it != mpTables->end(); ++it)
    {
      double Variable = 0.0;
      for(std::size_t j=1; j<rShapeFunctions.size(); ++j)
        Variable = rShapeFunctions[j] * rGeometry[j].FastGetSolutionStepValue(it->second.GetYVariable());
      mTableArguments[it->first] = Variable;
    }
  }

  template<class TVariableType>
  typename TVariableType::Type& GetValue(const TVariableType& rV)
  {
    typename TablesContainerType::iterator i;
    
    if((i = std::find_if(mpTables->begin(), mpTables->end(), VariableCheck(rV.Key()))) != mpTables->end())
    {
      return *static_cast<const typename TVariableType::Type*>((i->second)[mTableArguments[(i->first)]]);
    }
    return mpData->GetValue(rV);
  }

  template<class TVariableType>
  typename TVariableType::Type const& GetValue(const TVariableType& rV) const
  {
    typename TablesContainerType::iterator i;
    
    if((i = std::find_if(mpTables->begin(), mpTables->end(), VariableCheck(rV.Key()))) != mpTables->end())
    {
      return *static_cast<const typename TVariableType::Type*>((i->second)[mTableArguments[(i->first)]]);
    }
    return mpData->GetValue(rV);
  }

  bool HasVariables() const
  {
    return !mpData->IsEmpty();
  }

  bool HasTables() const
  {
    return !mpTables->empty();
  }

  bool IsEmpty() const
  {
    return !( HasVariables() || HasTables() );
  }

  int64_t Key(std::size_t XKey, std::size_t YKey) const
  {
    int64_t result_key = XKey;
    result_key = result_key << 32;
    result_key |= YKey; // I know that the key is less than 2^32 so I don't need zeroing the upper part
    return result_key;
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
    return "PropertiesLayout";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const
  {
    rOStream << "PropertiesLayout";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const
  {
    if(mpData != nullptr)
      mpData->PrintData(rOStream);
    if(mpTables != nullptr)
      rOStream << "This properties contains " << mpTables->size() << " tables";
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

  const ContainerType* mpData;

  const TablesContainerType* mpTables;

  Vector mTableArguments; // all table arguments considered as "double" type

  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  class VariableCheck
  {
    std::size_t mI;
   public:
    VariableCheck(std::size_t I) : mI(I) {}
    bool operator()(const DataValueContainer::ValueType& I)
    {
      return I.second->GetYVariable()->Key() == mI;
    }
  };
  
  ///@}
  ///@name Serialization
  ///@{

  friend class Serializer;

  void save(Serializer& rSerializer) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject);
    rSerializer.save("mTableArguments", mTableArguments);
  }

  void load(Serializer& rSerializer)
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject);
    rSerializer.load("mTableArguments", mTableArguments);
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

}; // Class PropertiesLayout

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  PropertiesLayout& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PropertiesLayout& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PROPERTIES_LAYOUT_H_INCLUDED  defined
