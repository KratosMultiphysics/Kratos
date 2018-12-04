//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborators:   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CALCULATE_NODAL_AREA_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_NODAL_AREA_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"

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

/**
 * @brief This struct is used in order to identify when using the hitorical and non historical variables
 */
struct CalculateNodalAreaSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};
    
/** 
 * @class CalculateNodalAreaProcess
 * @ingroup KratosCore 
 * @brief Computes NODAL_AREA
 * @details Calculate the NODAL_AREA for computing the weighted area in each node
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
// template<bool THistorical = true>
// class KRATOS_API(KRATOS_CORE) CalculateNodalAreaProcess
class CalculateNodalAreaProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Index type definition
    typedef std::size_t IndexType;
    
    /// Size type definition
    typedef std::size_t SizeType;
    
    /// Pointer definition of CalculateNodalAreaProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateNodalAreaProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param DomainSize The size of the space, if the value is not provided will compute from the model part
     */
    CalculateNodalAreaProcess(
        ModelPart& rModelPart, 
        const SizeType DomainSize = 0
        ): mrModelPart(rModelPart),
           mDomainSize(DomainSize)
    {
        // In case is not provided we will take from the model part
        if (mDomainSize == 0) {
            const auto& it_element_begin = mrModelPart.ElementsBegin();
            const auto& r_first_element_geometry = it_element_begin->GetGeometry();
            mDomainSize = r_first_element_geometry.WorkingSpaceDimension();
        }
    }

    /// Destructor.
    ~CalculateNodalAreaProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY

        // Set to zero the nodal area        
        const auto it_node_begin = mrModelPart.NodesBegin();
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); ++i) {
            auto it_node = it_node_begin + i;
            it_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        }

        const auto& it_element_begin = mrModelPart.ElementsBegin();
        const auto& r_first_element_geometry = it_element_begin->GetGeometry();
        const std::size_t local_space_dimension = r_first_element_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_first_element_geometry.PointsNumber();
        
        // The integration points
        const auto& integration_method = r_first_element_geometry.GetDefaultIntegrationMethod();
        const auto& integration_points = r_first_element_geometry.IntegrationPoints(integration_method);
        const std::size_t number_of_integration_points = integration_points.size();
        
        Vector N = ZeroVector(number_of_nodes);
        Matrix J0 = ZeroMatrix(mDomainSize, local_space_dimension);
        
        #pragma omp parallel for firstprivate(N, J0)
        for(int i=0; i<static_cast<int>(mrModelPart.Elements().size()); ++i) {
            auto it_elem = it_element_begin + i;
            auto& r_geometry = it_elem->GetGeometry();
            
            // The containers of the shape functions
            const auto& rNcontainer = r_geometry.ShapeFunctionsValues(integration_method);
            
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Getting the shape functions
                noalias(N) = row(rNcontainer, point_number);
                
                // Getting the jacobians and local gradients
                GeometryUtils::JacobianOnInitialConfiguration(r_geometry, integration_points[point_number], J0);
                const double detJ0 = MathUtils<double>::GeneralizedDet(J0);
                const double gauss_point_volume = integration_points[point_number].Weight() * detJ0;
                
                for(std::size_t i_node =0; i_node < number_of_nodes; ++i_node) {
                    #pragma omp atomic 
                    r_geometry[i_node].FastGetSolutionStepValue(NODAL_AREA) += N[i_node] * gauss_point_volume;
                }
            }
        }
    
        mrModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);

        KRATOS_CATCH("");

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
        return "CalculateNodalAreaProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateNodalAreaProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
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
    
    ModelPart& mrModelPart;  /// The model part where the nodal area is computed
    SizeType mDomainSize;    /// The dimension of the space

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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
    CalculateNodalAreaProcess& operator=(CalculateNodalAreaProcess const& rOther);

    /// Copy constructor.
    //CalculateNodalAreaProcess(CalculateNodalAreaProcess const& rOther);


    ///@}

}; // Class CalculateNodalAreaProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateNodalAreaProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateNodalAreaProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_CALCULATE_NODAL_AREA_PROCESS_H_INCLUDED  defined 


