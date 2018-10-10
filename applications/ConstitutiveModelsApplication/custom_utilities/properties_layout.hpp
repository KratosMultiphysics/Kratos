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
#include <string>
#include <iostream>
#include <cstddef>
#include <unordered_map>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/data_value_container.h"
#include "includes/process_info.h"
#include "includes/table.h"


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

/// PropertiesLayout encapsulates data shared by different Elements or Conditions
/**
 * PropertiesLayout encapsulates data shared by different Elements or Conditions. It can store any type of data and provides a variable base access to them.
 *  These are all parameters that can be shared between Element. Usually material parameters are common for a set of element, so this category of data is referred as properties.
 * But in general it can be any common parameter for a group of Elements. Sharing these data as properties reduces the memory used by the application and also helps updating them if necessary.
 * As mentioned before PropertiesLayout is a shared data container between Elements or Conditions. In finite element problems there are several parameters which are the same for a set of elements and conditions.
 * Thermal conductivity, elasticity of the material and viscosity of the fluid are examples of these parameters. PropertiesLayout holds these data and is shared by elements or Conditions. This eliminates memory overhead due to redundant copies of these data for each element and Condition. PropertiesLayout also can be used to access nodal data if it is necessary.
 * It is important to mention that accessing the nodal data via PropertiesLayout is not the same as accessing it via Node. When user asks PropertiesLayout for a variable data in a Node, the process starts with finding the variable in the PropertiesLayout data container and if it does not exist then get it from Node.
 * This means that the priority of data is with the one stored in PropertiesLayout and then in Node.
 */
class PropertiesLayout : public IndexedObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PropertiesLayout
    KRATOS_CLASS_POINTER_DEFINITION(PropertiesLayout);

#ifdef  _WIN32 // work around for windows int64_t error
    typedef __int64 int64_t;
#endif
    typedef IndexedObject BaseType;

    /// Type of container used for variables
    typedef DataValueContainer ContainerType;

    typedef Node<3> NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef NodeType::IndexType IndexType;

    typedef Table<double> TableType;

    typedef std::unordered_map<std::size_t, TableType> TablesContainerType; // This is a provisional implementation and should be changed to hash. Pooyan.

    typedef std::pair<const VariableData*, std::size_t>  VariableTableType;

    typedef std::vector<VariableTableType> VariableTableContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PropertiesLayout(IndexType NewId = 0) : BaseType(NewId), mData(), mTables() {}

    /// Copy constructor.
    PropertiesLayout(const PropertiesLayout& rOther) : BaseType(rOther), mData(rOther.mData), mTables(rOther.mTables) {}

    /// Destructor.
    ~PropertiesLayout() override {}


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

    void Configure(const Properties& rProperties, const GeometryType& rGeometry,const Vector& rShapeFunctions)
    {
      mpData = &rProperties.Data();
      mpTables = &rProperties.Tables();

      std::size_t max_size = 0;
      for(auto it = mTables->begin(); it != mTables->.end(); ++it)
        if( it->first > max_size )
          max_size = it->first;

      mTableArguments.resize(max_size);

      for(auto it = mTables->begin(); it != mTables->.end(); ++it)
      {
        mVariableTables.push_back(VariableTableType((it->second)->GetYVariable(),it->first));
        double Variable = 0.0;
        for(std::size_t j=1; j<rShapeFunctions.size(); ++j)
          Variable = rShapeFunctions[j] * rGeometry[j].FastGetSolutionStepValue(pTable->GetYVariable());
        mTableArguments[it->first] = Variable;
      }
    }

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rV)
    {
      if((i = std::find_if(mVariableTables.begin(), mVariableTables.end(), IndexCheck(rV.Key()))) != mVariableTables.end())
      {
        return *static_cast<const TVariableType::Type*>(mpTables[(*i->first)][mTableArguments[]]);
      }

      return mData.GetValue(rV);
    }

    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rV) const
    {
      if((i = std::find_if(mVariableTables.begin(), mVariableTables.end(), IndexCheck(rThisVariable.Key()))) != mVariableTables.end())
      {
        double Variable = 0.0;
        for(std::size_t j=1; j<mpShapeFunctions->size(); ++j)
          Variable = mpShapeFunctions[j] * GetGeometry[j].FastGetSolutionStepValue((i->second)->GetYVariable());

        return *static_cast<const TVariableType::Type*>((*i->second)[Key((i->second)->GetXVariable().Key(), (i->second)->GetYVariable().Key())][Variable]);
      }

      return mData.GetValue(rV);
    }


    bool HasVariables()
    {
        return !mpData->IsEmpty();
    }

    bool HasTables()
    {
        return !mpTables->empty();
    }

    bool IsEmpty()
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
    std::string Info() const override
    {
        return "PropertiesLayout";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream <<  "PropertiesLayout";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        mData.PrintData(rOStream);
        rOStream << "This properties contains " << mTables.size() << " tables";
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

    ContainerType* mpData;

    TablesContainerType* mpTables;

    VariableTablesContainerType mVariableTables;

    Vector mTableArguments; // all table arguments considered as "double" type

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject);
        rSerializer.save("Data", mData);
        rSerializer.save("Tables", mTables);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject);
        rSerializer.load("Data", mData);
        rSerializer.load("Tables", mTables);
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
