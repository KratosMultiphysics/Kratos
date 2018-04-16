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
template<class TSparseSpace, class TDenseSpace>
class NearestElementMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestElementMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestElementMapper);

    ///@}
    ///@name Life Cycle
    ///@{

    NearestElementMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                         : Mapper<TSparseSpace, TDenseSpace>(rModelPartOrigin,
                                  rModelPartDestination) {}

    NearestElementMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters,
                         const bool IsMPIExecution)
                         : Mapper<TSparseSpace, TDenseSpace>(rModelPartOrigin,
                                  rModelPartDestination,
                                  JsonParameters,
                                  IsMPIExecution)
    {
        // mpMapperCommunicator->InitializeOrigin(MapperUtilities::Condition_Center);
        // mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
        // mpMapperCommunicator->Initialize();

        // mpInverseMapper.reset(); // explicitly specified to be safe
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
        // mpMapperCommunicator->UpdateInterface(MappingOptions, SearchRadius);
        // if (mpInverseMapper)
        // {
        //     mpInverseMapper->UpdateInterface(MappingOptions, SearchRadius);
        // }

        // if (MappingOptions.Is(MapperFlags::REMESHED))
        // {
        //     ComputeNumberOfNodesAndConditions();
        // }
    }

    typename Mapper<TSparseSpace, TDenseSpace>::Pointer Clone(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters,
                          const bool IsMPIExecution) override
    {
        return Kratos::make_shared<NearestElementMapper<TSparseSpace, TDenseSpace>>(rModelPartOrigin,
                                                         rModelPartDestination,
                                                         JsonParameters,
                                                         IsMPIExecution);
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
        std::stringstream buffer;
        buffer << "NearestElementMapper" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestElementMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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

    Parameters GetInterfaceParameters() override
    {
        Parameters mapper_interface_parameters = Parameters( R"(
        {
            "mapper_condition_name" : "",
            "use_nodes"      : true
        }  )" );

        return mapper_interface_parameters;
    }

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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // static double GetInterpolatedValueFromGeometryScalar(InterfaceObject::Pointer pInterfaceObject, //TODO const
    //         const Variable<double>& rVariable,
    //         const Kratos::Flags& rOptions,
    //         const std::vector<double>& rShapeFunctionValues)
    // {
    //     Geometry<Node<3>>* p_base_geometry = pInterfaceObject->pGetBaseGeometry();

    //     double interpolated_value = 0.0f;

    //     for (std::size_t i = 0; i < p_base_geometry->PointsNumber(); ++i)
    //     {
    //         interpolated_value += p_base_geometry->GetPoint(i).FastGetSolutionStepValue(rVariable) * rShapeFunctionValues[i];
    //     }
    //     return interpolated_value;
    // }

    // static array_1d<double, 3> GetInterpolatedValueFromGeometryVector(InterfaceObject::Pointer pInterfaceObject, //TODO const
    //         const Variable< array_1d<double, 3> >& rVariable,
    //         const Kratos::Flags& rOptions,
    //         const std::vector<double>& rShapeFunctionValues)
    // {
    //     Geometry<Node<3>>* p_base_geometry = pInterfaceObject->pGetBaseGeometry();

    //     array_1d<double, 3> interpolated_value;
    //     interpolated_value[0] = 0.0f;
    //     interpolated_value[1] = 0.0f;
    //     interpolated_value[2] = 0.0f;
    //     for (std::size_t i = 0; i < p_base_geometry->PointsNumber(); ++i)
    //     {
    //         interpolated_value[0] += p_base_geometry->GetPoint(i).FastGetSolutionStepValue(rVariable)[0] * rShapeFunctionValues[i];
    //         interpolated_value[1] += p_base_geometry->GetPoint(i).FastGetSolutionStepValue(rVariable)[1] * rShapeFunctionValues[i];
    //         interpolated_value[2] += p_base_geometry->GetPoint(i).FastGetSolutionStepValue(rVariable)[2] * rShapeFunctionValues[i];
    //     }
    //     return interpolated_value;
    // }


    // template <typename T>
    // static void SetValueOfNode(InterfaceObject::Pointer pInterfaceObject,
    //                            const T& rValue,
    //                            const Variable< T >& rVariable,
    //                            const Kratos::Flags& rOptions,
    //                            const double Factor)
    // {
    //     Node<3>* p_base_node = pInterfaceObject->pGetBaseNode();

    //     if (rOptions.Is(MapperFlags::ADD_VALUES))
    //     {
    //         p_base_node->FastGetSolutionStepValue(rVariable) += rValue * Factor;
    //     }
    //     else
    //     {
    //         p_base_node->FastGetSolutionStepValue(rVariable) = rValue * Factor;
    //     }
    // }

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
    // NearestElementMapper& operator=(NearestElementMapper const& rOther) {}

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


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED  defined
