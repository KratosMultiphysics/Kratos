// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "adjoint_lift_response_function_coordinates_global.h"
// #include "node.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_processes/compute_lift_process.h"

namespace Kratos
{
    AdjointLiftGlobalCoordinatesResponseFunction::AdjointLiftGlobalCoordinatesResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
     : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        // This response function currently only works in 2D!
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const int domain_size = r_current_process_info[DOMAIN_SIZE];
        KRATOS_ERROR_IF(domain_size != 2) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;
    }

    AdjointLiftGlobalCoordinatesResponseFunction::~AdjointLiftGlobalCoordinatesResponseFunction(){}

    void AdjointLiftGlobalCoordinatesResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);
        rResponseGradient.clear();
        ComputeInitialLift();
        auto geom = rAdjointElement.GetGeometry();
        unsigned int NumNodes = geom.PointsNumber();
        double epsilon=1e-6;
        double perturbated_lift=0.0;

        if (rAdjointElement.IsNot(MARKER)){
            for (unsigned int i=0;i<NumNodes;i++)
            {   
                if (geom[i].IsNot(VISITED) && geom[i].IsNot(STRUCTURE)){
                    geom[i].Set(VISITED);
                    Vector resultForce(3);
                    int node_id = geom[i].Id();
                    double unperturbed_potential=geom[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
                    mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(VELOCITY_POTENTIAL) = unperturbed_potential + epsilon ;
                    ComputeLiftProcess(mrModelPart,resultForce).Execute();
                    perturbated_lift=resultForce(1);
                    rResponseGradient(i) = -(mInitialLift-perturbated_lift)/epsilon;
                    mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(VELOCITY_POTENTIAL) = unperturbed_potential;
                }else{
                    rResponseGradient(i) = 0.0;
                }
            }
        }else{
            array_1d<double,3> distances;
            distances=rAdjointElement.GetValue(ELEMENTAL_DISTANCES);
            for (unsigned int i=0;i<NumNodes;i++)
            {   
                if (geom[i].IsNot(VISITED)&& geom[i].IsNot(STRUCTURE)){
                    geom[i].Set(VISITED);
                    Vector resultForce(3);
                    int node_id = geom[i].Id();
                    if(distances[i] > 0){
                        double unperturbed_potential=geom[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(VELOCITY_POTENTIAL) = unperturbed_potential + epsilon ;
                        ComputeLiftProcess(mrModelPart,resultForce).Execute();
                        perturbated_lift=resultForce(1);
                        rResponseGradient(i) = -(mInitialLift-perturbated_lift)/epsilon;
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(VELOCITY_POTENTIAL) = unperturbed_potential;

                    }
                    else{
                        double unperturbed_potential=geom[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = unperturbed_potential + epsilon;
                        ComputeLiftProcess(mrModelPart,resultForce).Execute();
                        perturbated_lift=resultForce(1);
        
                        rResponseGradient(i) = -(mInitialLift-perturbated_lift)/epsilon;
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = unperturbed_potential;
                    }
                    if(distances[i] < 0){
                        double unperturbed_potential=geom[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(VELOCITY_POTENTIAL) = unperturbed_potential + epsilon ;
                        ComputeLiftProcess(mrModelPart,resultForce).Execute();
                        perturbated_lift=resultForce(1);
                        rResponseGradient(i+NumNodes) = -(mInitialLift-perturbated_lift)/epsilon;
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(VELOCITY_POTENTIAL) = unperturbed_potential;

                    }
                    else{
                        double unperturbed_potential=geom[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = unperturbed_potential + epsilon ;
                        ComputeLiftProcess(mrModelPart,resultForce).Execute();
                        perturbated_lift=resultForce(1);
                        rResponseGradient(i+NumNodes) = -(mInitialLift-perturbated_lift)/epsilon;
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = unperturbed_potential;
                    }
                }else{
                    rResponseGradient(i) = 0.0;
                    rResponseGradient(i+NumNodes) = 0.0;
                }
            }
        }
        KRATOS_CATCH("");
    }

    void AdjointLiftGlobalCoordinatesResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();
        KRATOS_CATCH("");
    }
  

    void AdjointLiftGlobalCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("")
    }

    void AdjointLiftGlobalCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("");
    }

    void AdjointLiftGlobalCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("")
    }

    void AdjointLiftGlobalCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        auto geom = rAdjointCondition.GetGeometry();
        unsigned int NumNodes = geom.PointsNumber();
        unsigned int Dim = geom.WorkingSpaceDimension();
        ComputeInitialLift();

        if (rSensitivityGradient.size() != NumNodes*Dim)
            rSensitivityGradient.resize(NumNodes*Dim, false);

        if (rAdjointCondition.IsNot(BOUNDARY)){
            double epsilon=1e-6;
            double perturbated_lift=0.0;
            for (unsigned int i_node=0;i_node<NumNodes;i_node++)
            {
                if (geom[i_node].IsNot(VISITED) && geom[i_node].IsNot(STRUCTURE)){
                    geom[i_node].Set(VISITED); 
                    int node_id = geom[i_node].Id();
                    for (unsigned int i_dim=0;i_dim<Dim;i_dim++)
                    {   
                        Vector resultForce(3);

                        mrModelPart.pGetNode(node_id,0)->GetInitialPosition()[i_dim] += epsilon;
                        mrModelPart.pGetNode(node_id,0)->Coordinates()[i_dim] += epsilon;

                        ComputeLiftProcess(mrModelPart,resultForce).Execute();
                        perturbated_lift=resultForce(1);
                        rSensitivityGradient(i_dim + i_node*Dim) = (mInitialLift-perturbated_lift)/epsilon;
                        
                        mrModelPart.pGetNode(node_id,0)->GetInitialPosition()[i_dim] -= epsilon;
                        mrModelPart.pGetNode(node_id,0)->Coordinates()[i_dim] -= epsilon;
                    }
                }else{
                    for (unsigned int i_dim=0;i_dim<Dim;i_dim++)
                        rSensitivityGradient(i_dim + i_node*Dim) = 0.0;
                }
            }
        }else{
            
            for (unsigned int i_node=0;i_node<NumNodes;i_node++)
            {
                for (unsigned int i_dim=0;i_dim<Dim;i_dim++)
                {   
                    rSensitivityGradient(i_dim + i_node*Dim) = 0.0;
                }
            }
        }
        KRATOS_CATCH("");
    }

    void AdjointLiftGlobalCoordinatesResponseFunction::FinalizeSolutionStep() 
    {
        KRATOS_TRY;
        // #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (mrModelPart.Nodes().size()); ++k)
        {
            auto it_node = mrModelPart.NodesBegin() + k;
            it_node->Set(VISITED,false);
        }

        KRATOS_CATCH("");
    }


    void AdjointLiftGlobalCoordinatesResponseFunction::ComputeInitialLift()
    {
        KRATOS_TRY;
        if (mComputeLift)
        {
            mInitialLift=CalculateValue(mrModelPart);
            mComputeLift=false;
        }
        KRATOS_CATCH("");
    }

    double AdjointLiftGlobalCoordinatesResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;
        Vector resultForce(3);         
        ComputeLiftProcess(rModelPart,resultForce).Execute();
        return resultForce(1);
        KRATOS_CATCH("");
    }
} // namespace Kratos.


