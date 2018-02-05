//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED )
#define  KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Interpolative Mapper
/** This class implements the Nearest Element Mapping technique.
* Each node on the destination side gets assigned is's closest condition or element (distance to center) 
* on the other side of the interface.
* In the mapping phase every node gets assigned the interpolated value of the condition/element.
* The interpolation is done with the shape funcitons
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
class NearestElementMapper : public Mapper
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestElementMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestElementMapper);

    ///@}
    ///@name Life Cycle
    ///@{

    NearestElementMapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination) : Mapper(
                         rModelPartOrigin, rModelPartDestination) {}

    NearestElementMapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination,
                         Parameters JsonParameters) : Mapper(
                                 rModelPartOrigin, rModelPartDestination, JsonParameters)
    {
        mpMapperCommunicator->InitializeOrigin(MapperUtilities::Condition_Center);
        mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
        mpMapperCommunicator->Initialize();

        mpInverseMapper.reset(); // explicitly specified to be safe
    }

    /// Destructor.
    virtual ~NearestElementMapper() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius) override
    {
        mpMapperCommunicator->UpdateInterface(MappingOptions, SearchRadius);
        if (mpInverseMapper)
        {
            mpInverseMapper->UpdateInterface(MappingOptions, SearchRadius);
        }

        if (MappingOptions.Is(MapperFlags::REMESHED))
        {
            ComputeNumberOfNodesAndConditions();
        }
    }

    /* This function maps from Origin to Destination */
    void Map(const Variable<double>& rOriginVariable,
             const Variable<double>& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        double factor = 1.0f;

        if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            factor = MapperUtilities::ComputeConservativeFactor(
                         mNumNodesOrigin,
                         mNumNodesDestination);
        }

        ProcessMappingOptions(MappingOptions, factor);

        // Creating the function pointers for the InterfaceObjects
        auto function_pointer_origin = std::bind(&GetInterpolatedValueFromGeometryScalar,
                                       std::placeholders::_1,
                                       rOriginVariable,
                                       MappingOptions,
                                       std::placeholders::_2);

        auto function_pointer_destination = std::bind(&SetValueOfNode<double>,
                                            std::placeholders::_1,
                                            std::placeholders::_2,
                                            rDestinationVariable,
                                            MappingOptions,
                                            factor);

        mpMapperCommunicator->TransferVariableData(function_pointer_origin,
                function_pointer_destination,
                rOriginVariable);
    }

    /* This function maps from Origin to Destination */
    void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
             const Variable< array_1d<double, 3> >& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        double factor = 1.0f;

        if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            factor = MapperUtilities::ComputeConservativeFactor(
                         mNumNodesOrigin,
                         mNumNodesDestination);
        }

        ProcessMappingOptions(MappingOptions, factor);

        // Creating the function pointers for the InterfaceObjects
        auto function_pointer_origin = std::bind(&GetInterpolatedValueFromGeometryVector,
                                       std::placeholders::_1,
                                       rOriginVariable,
                                       MappingOptions,
                                       std::placeholders::_2);

        auto function_pointer_destination = std::bind(&SetValueOfNode< array_1d<double, 3> >,
                                            std::placeholders::_1,
                                            std::placeholders::_2,
                                            rDestinationVariable,
                                            MappingOptions,
                                            factor);

        mpMapperCommunicator->TransferVariableData(function_pointer_origin,
                function_pointer_destination,
                rOriginVariable);
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable<double>& rOriginVariable,
                    const Variable<double>& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {
        // Construct the inverse mapper if it hasn't been done before
        // It is constructed with the order of the model_parts changed!
        if (!mpInverseMapper)
        {
            mpInverseMapper = this->Clone(mModelPartDestination,
                                          mModelPartOrigin,
                                          mJsonParameters);
        }
        mpInverseMapper->Map(rDestinationVariable, rOriginVariable, MappingOptions);
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                    const Variable< array_1d<double, 3> >& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {
        // Construct the inverse mapper if it hasn't been done before
        // It is constructed with the order of the model_parts changed!
        if (!mpInverseMapper)
        {
            mpInverseMapper = this->Clone(mModelPartDestination,
                                          mModelPartOrigin,
                                          mJsonParameters);
        }
        mpInverseMapper->Map(rDestinationVariable, rOriginVariable, MappingOptions);
    }

    Mapper::Pointer Clone(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters) override
    {
        return Kratos::make_shared<NearestElementMapper>(rModelPartOrigin,
                                                         rModelPartDestination,
                                                         JsonParameters);
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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "NearestElementMapper" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestElementMapper";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    Mapper::Pointer mpInverseMapper;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static double GetInterpolatedValueFromGeometryScalar(InterfaceObject::Pointer pInterfaceObject, //TODO const
            const Variable<double>& rVariable,
            const Kratos::Flags& rOptions,
            const std::vector<double>& rShapeFunctionValues)
    {
        Geometry<Node<3>>* p_base_geometry = pInterfaceObject->pGetBaseGeometry();

        double interpolated_value = 0.0f;

        for (std::size_t i = 0; i < p_base_geometry->PointsNumber(); ++i)
        {
            interpolated_value += p_base_geometry->GetPoint(i).FastGetSolutionStepValue(rVariable) * rShapeFunctionValues[i];
        }
        return interpolated_value;
    }

    static array_1d<double, 3> GetInterpolatedValueFromGeometryVector(InterfaceObject::Pointer pInterfaceObject, //TODO const
            const Variable< array_1d<double, 3> >& rVariable,
            const Kratos::Flags& rOptions,
            const std::vector<double>& rShapeFunctionValues)
    {
        Geometry<Node<3>>* p_base_geometry = pInterfaceObject->pGetBaseGeometry();

        array_1d<double, 3> interpolated_value;
        interpolated_value[0] = 0.0f;
        interpolated_value[1] = 0.0f;
        interpolated_value[2] = 0.0f;
        for (std::size_t i = 0; i < p_base_geometry->PointsNumber(); ++i)
        {
            interpolated_value[0] += p_base_geometry->GetPoint(i).FastGetSolutionStepValue(rVariable)[0] * rShapeFunctionValues[i];
            interpolated_value[1] += p_base_geometry->GetPoint(i).FastGetSolutionStepValue(rVariable)[1] * rShapeFunctionValues[i];
            interpolated_value[2] += p_base_geometry->GetPoint(i).FastGetSolutionStepValue(rVariable)[2] * rShapeFunctionValues[i];
        }
        return interpolated_value;
    }


    template <typename T>
    static void SetValueOfNode(InterfaceObject::Pointer pInterfaceObject,
                               const T& rValue,
                               const Variable< T >& rVariable,
                               const Kratos::Flags& rOptions,
                               const double Factor)
    {
        Node<3>* p_base_node = pInterfaceObject->pGetBaseNode();

        if (rOptions.Is(MapperFlags::ADD_VALUES))
        {
            p_base_node->FastGetSolutionStepValue(rVariable) += rValue * Factor;
        }
        else
        {
            p_base_node->FastGetSolutionStepValue(rVariable) = rValue * Factor;
        }
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

    /// Assignment operator.
    NearestElementMapper& operator=(NearestElementMapper const& rOther);

    //   /// Copy constructor.
    //   NearestElementMapper(NearestElementMapper const& rOther){}


    ///@}

}; // Class NearestElementMapper

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  NearestElementMapper& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const NearestElementMapper& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED  defined
