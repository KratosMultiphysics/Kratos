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
#include "adjoint_lift_response_function_coordinates_jump.h"
// #include "node.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_processes/compute_lift_process.h"

namespace Kratos
{
    AdjointLiftJumpCoordinatesResponseFunction::AdjointLiftJumpCoordinatesResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
     : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        // This response function currently only works in 2D!
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const int domain_size = r_current_process_info[DOMAIN_SIZE];
        KRATOS_ERROR_IF(domain_size != 2) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;
        // Get id of node where a displacement should be traced
        // const int id_traced_node = this->GetKuttaNodeId();

        // Get pointer to traced node
        // mpTracedNode = rModelPart.pGetNode(id_traced_node);

        this->GetNeighboringElementPointer();
    }

    AdjointLiftJumpCoordinatesResponseFunction::~AdjointLiftJumpCoordinatesResponseFunction(){}

    void AdjointLiftJumpCoordinatesResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);
        rResponseGradient.clear();
        if( rAdjointElement.Id() == mpNeighboringElement->Id() )
        {
            double value= 2.0/10.0;
            unsigned int NumNodes = rAdjointElement.GetGeometry().size();
            for(IndexType i = 0; i < NumNodes; ++i)
            {
                if(rAdjointElement.GetGeometry()[i].Is(STRUCTURE))
                {
                    rResponseGradient[i] = value;
                    rResponseGradient[i+NumNodes] = -value;
                }
            }
        }
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
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
  

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
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

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
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

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
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

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
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

    void AdjointLiftJumpCoordinatesResponseFunction::ComputeInitialLift()
    {
        KRATOS_TRY;
        if (mComputeLift)
        {
            mInitialLift=CalculateValue(mrModelPart);
            mComputeLift=false;
        }
        KRATOS_CATCH("");
    }
    Vector AdjointLiftJumpCoordinatesResponseFunction::ComputeNormalCondition(Geometry<Node<3>> &rGeom)
    {
        Vector An(3);        
        An[0] = rGeom[1].Y() - rGeom[0].Y();
        An[1] = -(rGeom[1].X() - rGeom[0].X());
        An[2] = 0.00;        
        return An;
    }


    Vector AdjointLiftJumpCoordinatesResponseFunction::ComputeNormal(Geometry<Node<3>> &rGeom)
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
        }
        return An;
    }

    Vector AdjointLiftJumpCoordinatesResponseFunction::ComputeNormalPerturbed(Geometry<Node<3>> &rGeom,Vector normal_ini)
    {
        Vector An(3);
        Vector x(2),y(2);

        unsigned int dim=0;
        for (unsigned int i=0;i<rGeom.PointsNumber();i++){
            if (rGeom[i].Is(SOLID)){
                x(dim)=rGeom[i].X();
                y(dim)=rGeom[i].Y();
                dim++;
            }
        }
        An(0) = y(1) - y(0);
        An(1) = -(x(1) - x(0));
        An(2) = 0.00;

        for (std::size_t i_dim=0;i_dim<rGeom.WorkingSpaceDimension();i_dim++){
            if((An[i_dim]*(normal_ini(i_dim)))<0)
                An(i_dim) = -An(i_dim);
        }
        return An;
    }

    double AdjointLiftJumpCoordinatesResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;
        Vector resultForce(3);         
        ComputeLiftProcess(rModelPart,resultForce).Execute();
        return resultForce(1);
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::GetElementCandidates(WeakPointerVector<Element> &ElementCandidates, GeometryType &rGeom)
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

    void AdjointLiftJumpCoordinatesResponseFunction::GetSortedIds(std::vector<IndexType> &Ids,
                        const GeometryType &rGeom)
    {
        Ids.resize(rGeom.PointsNumber());
        for (SizeType i = 0; i < Ids.size(); i++)
            Ids[i] = rGeom[i].Id();
        std::sort(Ids.begin(), Ids.end());
    }

    Element::WeakPointer AdjointLiftJumpCoordinatesResponseFunction::FindParentElement(std::vector<IndexType> &NodeIds,
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
    void AdjointLiftJumpCoordinatesResponseFunction::GetNeighboringElementPointer()
    {
        KRATOS_TRY;

        for (auto elem_it = mrModelPart.Elements().ptr_begin(); elem_it != mrModelPart.Elements().ptr_end(); ++elem_it)
        {   
            if ((*elem_it)->Is(MARKER) && (*elem_it)->Is(STRUCTURE)){
                mpNeighboringElement = (*elem_it);
                return;
            }
        }
        KRATOS_ERROR << "No neighboring element is available for the traced node." << std::endl;

        KRATOS_CATCH("");
    }
} // namespace Kratos.


