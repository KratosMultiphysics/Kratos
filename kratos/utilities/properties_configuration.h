//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_PROPERTIES_CONFIGURATION_H)
#define  KRATOS_PROPERTIES_CONFIGURATION_H

// System includes
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/data_value_container.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

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

/// Short class definition.
/** Detail class definition.
*/
class PropertiesConfiguration
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PropertiesConfiguration
    KRATOS_CLASS_POINTER_DEFINITION(PropertiesConfiguration);

    using NodeType = Node<3>;
    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;
    using SizeType = std::size_t;

    using DoubleVariableType = Variable<double>;
    using ComponentVariableType = VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > >;
    using Array3VariableType = Variable<array_1d<double, 3>>;

    struct ConfigurationParameters
    {
    public:

        ConfigurationParameters(const GeometryType& rGeom,
                                const DataValueContainer& rData) :
                                mpGeometry(&rGeom),
                                mpData(&rData)
                                {}

        void SetShapeFunctionsValues (const Vector& rShapeFunctionsValues) {mpShapeFunctionsValues=&rShapeFunctionsValues;};
        void SetGaussPointIndex      (const int GaussPointIndex)           {mGaussPointIndex=GaussPointIndex;};
        void SetProcessInfo          (ProcessInfo& rProcessInfo)           {mpProcessInfo = &rProcessInfo;};

    private:
        const Vector*               mpShapeFunctionsValues;
        int                         mGaussPointIndex;
        const ProcessInfo*          mpProcessInfo;
        const GeometryType*         mpGeometry;
        const DataValueContainer*   mpData;
    };

    using DoubleFunction = std::function<double(const Properties&,
                                                const ConfigurationParameters&)>;

    using DoubleConfiguration = std::map<DoubleVariableType, DoubleFunction>;

    // std::function<void(int)>

    // using std::unordered_map<DoubleVariableType>

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PropertiesConfiguration()
        : mpParentProperties(nullptr) {}

    /// Destructor.
    virtual ~PropertiesConfiguration() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rV,
                                           const ConfigurationParameters& rParams)
    {
        mpParentProperties->GetValue(rV);
    }


    double GetValue(const DoubleVariableType& rV,
                    const ConfigurationParameters& rParams)
    {
        if(mConfigurationsDouble.find(rV) != mConfigurationsDouble.end())
            return mConfigurationsDouble[rV](*mpParentProperties,
                                             rParams);
        else
            return mpParentProperties->GetValue(rV);
    }

    void Configure(const DoubleVariableType& rV, std::function<double(
                                                const Properties&,
                                                const ConfigurationParameters&)> DoubleFunction)
    {
        mConfigurationsDouble[rV] = DoubleFunction;
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
        std::stringstream buffer;
        buffer << "PropertiesConfiguration" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const {rOStream << "PropertiesConfiguration";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}


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

    Properties* mpParentProperties;
    DoubleConfiguration mConfigurationsDouble;

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

    void save(Serializer& rSerializer) const
    {
        // rSerializer.save("ParentProps", mpParentProperties);
        // TODO check if this works... (=> circular dependencies)
        // If it is not working I could make a "SetAccessorProps"
        // that is called in the load fct of the Properties and passes a "this"
    }

    void load(Serializer& rSerializer)
    {
        // rSerializer.load("ParentProps", mpParentProperties);
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

    // /// Assignment operator.
    // PropertiesConfiguration& operator=(PropertiesConfiguration const& rOther){}

    // /// Copy constructor.
    // PropertiesConfiguration(PropertiesConfiguration const& rOther){}


    ///@}

}; // Class PropertiesConfiguration

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 PropertiesConfiguration& rThis){}

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const PropertiesConfiguration& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PROPERTIES_CONFIGURATION_H  defined
