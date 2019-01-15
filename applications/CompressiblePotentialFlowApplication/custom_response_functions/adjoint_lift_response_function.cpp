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
#include "adjoint_lift_response_function.h"
// #include "node.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_processes/compute_lift_level_set_process.h"

namespace Kratos
{
    AdjointLiftResponseFunction::AdjointLiftResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
     : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        // This response function currently only works in 2D!
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const int domain_size = r_current_process_info[DOMAIN_SIZE];
        KRATOS_ERROR_IF(domain_size != 2) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;
    }

    AdjointLiftResponseFunction::~AdjointLiftResponseFunction(){}

    void AdjointLiftResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();
        ComputeInitialLift();
        Vector resultForce(3);
        auto geom = rAdjointElement.GetGeometry();
        unsigned int NumNodes = geom.PointsNumber();
        double epsilon=1e-9;
        double perturbated_lift=0.0;
        if (rAdjointElement.Is(BOUNDARY)){
            if (rAdjointElement.IsNot(MARKER)){
                for (unsigned int i=0;i<NumNodes;i++)
                {   
                    Vector resultForce(3);
                    int node_id = geom[i].Id();
                    double unperturbed_potential=geom[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                    mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential + epsilon ;
                    ComputeLiftLevelSetProcess(mrModelPart,resultForce).Execute();
                    perturbated_lift=resultForce(1);
                    rResponseGradient(i) = (mInitialLift-perturbated_lift)/epsilon;
                    mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential;

                }
            }else{
                array_1d<double,3> distances;
                distances=rAdjointElement.GetValue(ELEMENTAL_DISTANCES);
                for (unsigned int i=0;i<NumNodes;i++)
                {
                    Vector resultForce(3);
                    int node_id = geom[i].Id();
                    if(distances[i] > 0){
                        double unperturbed_potential=geom[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential + epsilon ;
                        ComputeLiftLevelSetProcess(mrModelPart,resultForce).Execute();
                        perturbated_lift=resultForce(1);
                        rResponseGradient(i) = (mInitialLift-perturbated_lift)/epsilon;
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential;

                    }
                    else{
                        double unperturbed_potential=geom[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = unperturbed_potential + epsilon;
                        ComputeLiftLevelSetProcess(mrModelPart,resultForce).Execute();
                        perturbated_lift=resultForce(1);
        
                        rResponseGradient(i) = (mInitialLift-perturbated_lift)/epsilon;
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = unperturbed_potential;
                    }
                }
                for (unsigned int i=0;i<NumNodes;i++)
                {
                    Vector resultForce(3);
                    int node_id = geom[i].Id();
                    if(distances[i] < 0){
                        double unperturbed_potential=geom[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential + epsilon ;
                        ComputeLiftLevelSetProcess(mrModelPart,resultForce).Execute();
                        perturbated_lift=resultForce(1);
                        rResponseGradient(i+NumNodes) = (mInitialLift-perturbated_lift)/epsilon;
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential;

                    }
                    else{
                        double unperturbed_potential=geom[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = unperturbed_potential + epsilon ;
                        ComputeLiftLevelSetProcess(mrModelPart,resultForce).Execute();
                        perturbated_lift=resultForce(1);
                        rResponseGradient(i+NumNodes) = (mInitialLift-perturbated_lift)/epsilon;
                        mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = unperturbed_potential;
                    }
                }
            }

        }else{
                for (unsigned int i=0;i<NumNodes;i++)
                {
                    rResponseGradient(i)=0.0;
                }
                if (rAdjointElement.Is(MARKER)){
                    for (unsigned int i=0;i<NumNodes;i++)
                    {
                        rResponseGradient(i+NumNodes)=0.0;
                    }  
                }
        }

        KRATOS_CATCH("");
    }
    // void AdjointLiftResponseFunction::CalculateGradient(const Element& rAdjointElement,
    //                                const Matrix& rResidualGradient,
    //                                Vector& rResponseGradient,
    //                                const ProcessInfo& rProcessInfo)
    // {
    //     KRATOS_TRY;

    //     if (rResponseGradient.size() != rResidualGradient.size1())
    //         rResponseGradient.resize(rResidualGradient.size1(), false);

    //     rResponseGradient.clear();
    //     ComputeInitialLift();
    //     Vector resultForce(3);
    //     auto geom = rAdjointElement.GetGeometry();
    //     unsigned int NumNodes = geom.PointsNumber();
    //     double epsilon=1e-9;
    //     double perturbated_lift=0.0;

    //     auto geomModelPart = mrModelPart.pGetElement(rAdjointElement.Id(),0)->GetGeometry();
    //     if (rAdjointElement.Is(BOUNDARY)){
    //         if (rAdjointElement.IsNot(MARKER)){
    //             for (unsigned int i=0;i<NumNodes;i++)
    //             {   
    //                 Vector resultForce(3);
    //                 int node_id = geom[i].Id();
    //                 std::vector <double> cp_ini;
    //                 std::vector< array_1d<double,3>> normal_ini;
    //                 mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(NORMAL,normal_ini,mrModelPart.GetProcessInfo());
    //                 mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(PRESSURE,cp_ini,mrModelPart.GetProcessInfo());
    //                 // Vector normal_initial = 
    //                 double unperturbed_potential=geom[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
    //                 mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential + epsilon ;
    //                 std::vector <double> cp;
    //                 std::vector< array_1d<double,3>> normal;
    //                 mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(NORMAL,normal,mrModelPart.GetProcessInfo());
    //                 mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());;
    //                 rResponseGradient(i) = (cp_ini[0]*normal_ini[0][1]-cp[0]*normal[0][1])/epsilon;
    //                 mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential;
    //                 KRATOS_WATCH(normal_ini)
    //                 KRATOS_WATCH(normal)

    //             }
    //         }else{
    //             array_1d<double,3> distances;
    //             distances=rAdjointElement.GetValue(ELEMENTAL_DISTANCES);
    //             for (unsigned int i=0;i<NumNodes;i++)
    //             {
    //                 Vector resultForce(3);
    //                 int node_id = geom[i].Id();
    //                 std::vector <double> cp_ini;
    //                 std::vector< array_1d<double,3>> normal_ini;
    //                 mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(PRESSURE,cp_ini,mrModelPart.GetProcessInfo());
    //                 mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(NORMAL,normal_ini,mrModelPart.GetProcessInfo());

    //                 if(distances[i] > 0){
    //                     double unperturbed_potential=geom[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
    //                     mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential + epsilon ;
    //                     std::vector <double> cp;
    //                     std::vector< array_1d<double,3>> normal;
    //                     mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(NORMAL,normal,mrModelPart.GetProcessInfo());
    //                     mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());;
    //                     rResponseGradient(i) = (cp_ini[0]*normal_ini[0][1]-cp[0]*normal[0][1])/epsilon;
    //                     mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential;

    //                 }
    //                 else{
    //                     double unperturbed_potential=geom[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
    //                     mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = unperturbed_potential + epsilon;
    //                     std::vector <double> cp;
    //                     std::vector< array_1d<double,3>> normal;
    //                     mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(NORMAL,normal,mrModelPart.GetProcessInfo());
    //                     mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());;
    //                     rResponseGradient(i) = (cp_ini[0]*normal_ini[0][1]-cp[0]*normal[0][1])/epsilon;
    //                     mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = unperturbed_potential;
    //                 }
    //             }
    //             for (unsigned int i=0;i<NumNodes;i++)
    //             {
    //                 Vector resultForce(3);
    //                 int node_id = geom[i].Id();
    //                 std::vector <double> cp_ini;
    //                 std::vector< array_1d<double,3>> normal_ini;
    //                 mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(PRESSURE,cp_ini,mrModelPart.GetProcessInfo());
    //                 mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(NORMAL,normal_ini,mrModelPart.GetProcessInfo());

    //                 if(distances[i] < 0){
    //                     double unperturbed_potential=geom[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL);
    //                     mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential + epsilon ;
    //                     std::vector <double> cp;
    //                     std::vector< array_1d<double,3>> normal;
    //                     mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(NORMAL,normal,mrModelPart.GetProcessInfo());
    //                     mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());;
    //                     rResponseGradient(i) = (cp_ini[0]*normal_ini[0][1]-cp[0]*normal[0][1])/epsilon;
    //                     mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(POSITIVE_POTENTIAL) = unperturbed_potential;

    //                 }
    //                 else{
    //                     double unperturbed_potential=geom[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL);
    //                     mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = unperturbed_potential + epsilon;
    //                     std::vector <double> cp;
    //                     std::vector< array_1d<double,3>> normal;
    //                     mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(NORMAL,normal,mrModelPart.GetProcessInfo());
    //                     mrModelPart.pGetElement(rAdjointElement.Id(),0) -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());;
    //                     rResponseGradient(i) = (cp_ini[0]*normal_ini[0][1]-cp[0]*normal[0][1])/epsilon;
    //                     mrModelPart.pGetNode(node_id,0)->FastGetSolutionStepValue(NEGATIVE_POTENTIAL) = unperturbed_potential;
    //                 }
    //             }
    //         }
    //     } else{
    //         for (unsigned int i=0;i<NumNodes;i++)
    //         {
    //             rResponseGradient(i)=0.0;
    //         }
    //         if (rAdjointElement.Is(MARKER)){
    //             for (unsigned int i=0;i<NumNodes;i++)
    //             {
    //                 rResponseGradient(i+NumNodes)=0.0;
    //             }  
    //         }
    //     }

    //     KRATOS_CATCH("");
    // }

    void AdjointLiftResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
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

    void AdjointLiftResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
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

    void AdjointLiftResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
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

    void AdjointLiftResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("");
    }


    void AdjointLiftResponseFunction::ComputeInitialLift()
    {
        if (mComputeLift)
        {
            mInitialLift=CalculateValue(mrModelPart);
            mComputeLift=false;
        }
    }

    double AdjointLiftResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;
        Vector resultForce(3);         
        ComputeLiftLevelSetProcess(rModelPart,resultForce).Execute();
        return resultForce(1);
        KRATOS_CATCH("");
    }
} // namespace Kratos.


