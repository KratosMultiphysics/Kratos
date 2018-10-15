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
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "custom_utilities/table_key_variables.hpp"
#include "custom_utilities/properties_extensions.hpp"

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

  typedef std::pair<std::size_t, double>  VariableKeyArgumentsType;

  typedef std::pair<std::size_t, VariableKeyArgumentsType> ScalarTableArgumentsType;

  typedef std::vector<ScalarTableArgumentsType> ScalarTableArgumentsContainerType;

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
      ,mTableVariables(rOther.mTableVariables)
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
  PropertiesLayout& operator=(PropertiesLayout const& rOther)
  {
    mpData = rOther.mpData;
    mpTables = rOther.mpTables;
    mTableVariables = rOther.mTableVariables;
    mTableArguments = rOther.mTableArguments;

    return *this;
  }

  ///@}
  ///@name Operations
  ///@{

  void RegisterTable(const Variable<double>& rXVariable, const Variable<double>& rYVariable)
  {
    mTableVariables.RegisterTable(rXVariable,rYVariable);
  }

  void Configure(const Properties& rProperties, const GeometryType& rGeometry, const Vector& rShapeFunctions);


  template<class TVariableType>
  void GetValue(const TVariableType& rVariable, typename TVariableType::Type& rValue) const
  {
    typename ScalarTableArgumentsContainerType::const_iterator i;

    if((i = std::find_if(mTableArguments.begin(), mTableArguments.end(), VariableKeyCheck(rVariable.Key()))) != mTableArguments.end())
    {
      rValue = static_cast<const typename TVariableType::Type>(mpTables->at((i->first))[(i->second).second]);
    }
    else
      rValue = mpData->GetValue(rVariable);
  }

  template<class TVariableType>
  void GetValue(const TVariableType& rVariable, typename TVariableType::Type& rValue)
  {
    typename ScalarTableArgumentsContainerType::const_iterator i;

    if((i = std::find_if(mTableArguments.begin(), mTableArguments.end(), VariableKeyCheck(rVariable.Key()))) != mTableArguments.end())
    {
      rValue = static_cast<const typename TVariableType::Type>(mpTables->at((i->first))[(i->second).second]);
    }
    else
      rValue = mpData->GetValue(rVariable);
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

  TableKeyVariables<double,double>  mTableVariables;

  ScalarTableArgumentsContainerType mTableArguments; // all table variable Keys and table arguments considered as "double" type


  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  class VariableKeyCheck
  {
    std::size_t mI;
   public:
    VariableKeyCheck(std::size_t I) : mI(I) {}
    bool operator()(const ScalarTableArgumentsType& I)
    {
      return I.first == mI;
    }
  };

  ///@}
  ///@name Serialization
  ///@{

  friend class Serializer;

  void save(Serializer& rSerializer) const
  {
    rSerializer.save("mTableVariables", mTableVariables);
    rSerializer.save("mTableArguments", mTableArguments);
  }

  void load(Serializer& rSerializer)
  {
    rSerializer.load("mTableVariables", mTableVariables);
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
