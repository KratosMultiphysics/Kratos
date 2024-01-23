// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef LUMPED_INTEGRATION_UTILITY
#define  LUMPED_INTEGRATION_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/find_conditions_neighbours_process.h"


namespace Kratos
{

class LumpedIntegrationUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef array_1d<double,3> array_3d;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /// Pointer definition of LumpedIntegrationUtility
    KRATOS_CLASS_POINTER_DEFINITION(LumpedIntegrationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    LumpedIntegrationUtility(ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~LumpedIntegrationUtility() = default;

    ///@}
    ///@name Operations
    ///@{

    void CalculateLumpedAreas()
    {
        FindConditionsNeighboursProcess find_conditions_neighbours_process(mrModelPart, mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);
        find_conditions_neighbours_process.Execute();

        auto nodes_begin = mrModelPart.NodesBegin();
        #pragma omp parallel for
        for(int i=0; i<mrModelPart.Nodes().size(); ++i)
        {
            auto& r_node_i = *(nodes_begin + i);

            // Get all neighbour conditions
            const GlobalPointersVector<Condition>& r_conditions = r_node_i.GetValue(NEIGHBOUR_CONDITIONS);

            // loop conditions
            for(unsigned int j=0; j<r_conditions.size(); ++j)
            {
                // Get geometry of current condition
                const Condition& r_condition = r_conditions[j];
                const Condition::GeometryType& r_geom_i = r_condition.GetGeometry();
                r_node_i.GetValue(LUMPED_AREA) += r_geom_i.DomainSize() / r_geom_i.size();
            }
        }
    }

    template <typename T>
    void Integrate(const T& rVariable)
    {
        const auto nodes_begin = mrModelPart.NodesBegin();
        #pragma omp parallel for
        for(int i=0; i<mrModelPart.Nodes().size(); ++i)
        {
            auto& r_node_i = *(nodes_begin + i);
            r_node_i.FastGetSolutionStepValue(rVariable) *= r_node_i.GetValue(LUMPED_AREA);
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const {return "CalculateLumpedAreas";}

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

private:

    ///@name Private Member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}

}; // Class CalculateLumpedAreas

}  // namespace Kratos.

#endif // LUMPED_INTEGRATION_UTILITY_INCLUDED defined
