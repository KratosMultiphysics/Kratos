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
#include "adjoint_lift_response_function_coordinates.h"
// #include "node.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_processes/compute_lift_process.h"

namespace Kratos
{
    AdjointLiftCoordinatesResponseFunction::AdjointLiftCoordinatesResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
     : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        // This response function currently only works in 2D!
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const int domain_size = r_current_process_info[DOMAIN_SIZE];
        KRATOS_ERROR_IF(domain_size != 2) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;
    }

    AdjointLiftCoordinatesResponseFunction::~AdjointLiftCoordinatesResponseFunction(){}

    void AdjointLiftCoordinatesResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);
        rResponseGradient.clear();
        auto geom = rAdjointElement.GetGeometry();
        unsigned int NumNodes = geom.PointsNumber();
        double epsilon=1e-6;

        bool elem_is_wing = false;
        unsigned int counter = 0;
        for (unsigned int i=0;i<NumNodes;i++){
            if (geom[i].Is(SOLID))
                counter++;
        }
        if (counter==NumNodes-1)
            elem_is_wing=true;
        
        if (elem_is_wing){
            Vector normal = ComputeNormal(geom);
            std::vector <double> cp_ini;
            Element::Pointer pElem = mrModelPart.pGetElement(rAdjointElement.Id());
            pElem -> GetValueOnIntegrationPoints(PRESSURE,cp_ini,mrModelPart.GetProcessInfo());
            if (rAdjointElement.IsNot(MARKER)){
                for (unsigned int i=0;i<NumNodes;i++)
                {        
                    Vector resultForce(3);
                    pElem -> GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL) += epsilon;
                    std::vector <double> cp;
                    pElem -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());
                    rResponseGradient(i) = normal(1)*(cp_ini[0]-cp[0])/epsilon;
                    pElem -> GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL) -= epsilon;
                }
            }else{
                array_1d<double,3> distances;
                distances=rAdjointElement.GetValue(ELEMENTAL_DISTANCES);
                for (unsigned int i=0;i<NumNodes;i++)
                {   
                    if(distances[i] > 0){
                        pElem -> GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL) += epsilon;
                        std::vector <double> cp;
                        pElem -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());
                        rResponseGradient(i) = normal(1)*(cp_ini[0]-cp[0])/epsilon;
                        pElem -> GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL) -= epsilon;
                    }
                    else{
                        pElem -> GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL) += epsilon;
                        std::vector <double> cp;
                        pElem -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());
                        rResponseGradient(i) = normal(1)*(cp_ini[0]-cp[0])/epsilon;
                        pElem -> GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL) -= epsilon;
                    }
                    if(distances[i] < 0){
                        pElem -> GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL) += epsilon;
                        std::vector <double> cp;
                        pElem -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());
                        rResponseGradient(i) = normal(1)*(cp_ini[0]-cp[0])/epsilon;
                        pElem -> GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_POTENTIAL) -= epsilon;
                    }
                    else{
                        pElem -> GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL) += epsilon;
                        std::vector <double> cp;
                        pElem -> GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());
                        rResponseGradient(i) = normal(1)*(cp_ini[0]-cp[0])/epsilon;
                        pElem -> GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_POTENTIAL) -= epsilon; 
                    }
                }
            }
        }
        KRATOS_CATCH("");
    }

    void AdjointLiftCoordinatesResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
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
  

    void AdjointLiftCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
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

    void AdjointLiftCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
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

    void AdjointLiftCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
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

    void AdjointLiftCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        auto geom = rAdjointCondition.GetGeometry();
        unsigned int NumNodes = geom.PointsNumber();
        unsigned int Dim = geom.WorkingSpaceDimension();

        if (rSensitivityGradient.size() != NumNodes*Dim)
            rSensitivityGradient.resize(NumNodes*Dim, false);

        double epsilon=1e-6;
        
        if (rAdjointCondition.IsNot(BOUNDARY)){
            Vector normal_ini = ComputeNormalCondition(geom);
            std::vector <double> cp_ini;
      
            GeometryType &rGeom = mrModelPart.pGetCondition(rAdjointCondition.Id())->GetGeometry();
            WeakPointerVector<Element> ElementCandidates;
            GetElementCandidates(ElementCandidates, rGeom);

            std::vector<IndexType> NodeIds, ElementNodeIds;
            GetSortedIds(NodeIds, rGeom);
            Element::WeakPointer pElem=FindParentElement(NodeIds, ElementNodeIds, ElementCandidates);

            pElem.lock() -> GetValueOnIntegrationPoints(PRESSURE,cp_ini,mrModelPart.GetProcessInfo());

            for (unsigned int i_node=0;i_node<NumNodes;i_node++)
            {   if (geom[i_node].IsNot(STRUCTURE)){
                    for (unsigned int i_dim=0;i_dim<Dim;i_dim++)
                    {   
                        pElem.lock() -> GetGeometry()[i_node].GetInitialPosition()[i_dim] += epsilon;
                        pElem.lock() -> GetGeometry()[i_node].Coordinates()[i_dim] += epsilon;

                        Vector normal_perturbed = ComputeNormalCondition(pElem.lock()->GetGeometry());
                        std::vector <double> cp_perturbed;
                        pElem.lock() -> GetValueOnIntegrationPoints(PRESSURE,cp_perturbed,mrModelPart.GetProcessInfo());

                        rSensitivityGradient(i_dim + i_node*Dim) = (cp_ini[0]*normal_ini(1)-cp_perturbed[0]*normal_perturbed(1))/epsilon;
                        
                        pElem.lock() -> GetGeometry()[i_node].GetInitialPosition()[i_dim] += epsilon;
                        pElem.lock() -> GetGeometry()[i_node].Coordinates()[i_dim] += epsilon;
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
                    rSensitivityGradient(i_dim + i_node*Dim) = 0.0;                
            }
        }
        KRATOS_WATCH(rSensitivityGradient)

        KRATOS_CATCH("");
    }

    void AdjointLiftCoordinatesResponseFunction::FinalizeSolutionStep() 
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


    void AdjointLiftCoordinatesResponseFunction::ComputeInitialLift()
    {
        KRATOS_TRY;
        if (mComputeLift)
        {
            mInitialLift=CalculateValue(mrModelPart);
            mComputeLift=false;
        }
        KRATOS_CATCH("");
    }
    Vector AdjointLiftCoordinatesResponseFunction::ComputeNormalCondition(Geometry<Node<3>> &rGeom)
    {
        Vector An(3);        
        An[0] = rGeom[1].Y() - rGeom[0].Y();
        An[1] = -(rGeom[1].X() - rGeom[0].X());
        An[2] = 0.00;        
        return An;
    }

    Vector AdjointLiftCoordinatesResponseFunction::ComputeNormal(Geometry<Node<3>> &rGeom)
    {
        Vector An(3);
        Vector x(2),y(2);
        double y_check;
        unsigned int dim=0;
        for (unsigned int i=0;i<rGeom.PointsNumber();i++){
            if (rGeom[i].Is(SOLID)){
                x(dim)=rGeom[i].X();
                y(dim)=rGeom[i].Y();
                dim++;
            }else{
                y_check=rGeom[i].Y();
            }
        }
        An[0] = y(1) - y(0);
        An[1] = -(x(1) - x(0));
        An[2] = 0.00;
        if((An[1]*(y_check-y(0)))>0){
            An = -An;
            std::cout<<"Switching normal facing in adjoint lift coordinates response function!"<<std::endl;
        }
        return An;
    }

    double AdjointLiftCoordinatesResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;
        Vector resultForce(3);         
        ComputeLiftProcess(rModelPart,resultForce).Execute();
        return resultForce(1);
        KRATOS_CATCH("");
    }

    void AdjointLiftCoordinatesResponseFunction::GetElementCandidates(WeakPointerVector<Element> &ElementCandidates, GeometryType &rGeom)
    {
        for (SizeType i = 0; i < rGeom.WorkingSpaceDimension(); i++)
        {
            if (rGeom[i].Has(NEIGHBOUR_ELEMENTS)) {
                WeakPointerVector<Element> &rNodeElementCandidates = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
                for (SizeType j = 0; j < rNodeElementCandidates.size(); j++)
                    ElementCandidates.push_back(rNodeElementCandidates(j));
            } else
                std::cout << "fuck" << std::endl; 
        }
    }

    void AdjointLiftCoordinatesResponseFunction::GetSortedIds(std::vector<IndexType> &Ids,
                        const GeometryType &rGeom)
    {
        Ids.resize(rGeom.PointsNumber());
        for (SizeType i = 0; i < Ids.size(); i++)
            Ids[i] = rGeom[i].Id();
        std::sort(Ids.begin(), Ids.end());
    }

    Element::WeakPointer AdjointLiftCoordinatesResponseFunction::FindParentElement(std::vector<IndexType> &NodeIds,
                            std::vector<IndexType> &ElementNodeIds,
                            WeakPointerVector<Element> ElementCandidates)
    {
        for (SizeType i = 0; i < ElementCandidates.size(); i++)
        {
            
            GeometryType &rElemGeom = ElementCandidates[i].GetGeometry();
            
            GetSortedIds(ElementNodeIds, rElemGeom);

            if (std::includes(ElementNodeIds.begin(), ElementNodeIds.end(), NodeIds.begin(), NodeIds.end()))
            {
                return ElementCandidates(i);
            }
        }
        Element::WeakPointer void_ptr;
        return void_ptr;
    }
} // namespace Kratos.


