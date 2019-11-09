//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nuñez, based on Martin Fusseder work, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "adjoint_lift_response_function_coordinates_jump.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{
    AdjointLiftJumpCoordinatesResponseFunction::AdjointLiftJumpCoordinatesResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
     : AdjointPotentialResponseFunction(rModelPart, ResponseSettings)
    {
        // This response function currently only works in 2D!
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const int domain_size = r_current_process_info[DOMAIN_SIZE];
        KRATOS_ERROR_IF(domain_size != 2) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;

        // Reading the reference chord from the parameters
        mReferenceChord = ResponseSettings["reference_chord"].GetDouble();
        const double eps = std::numeric_limits<double>::epsilon();
        KRATOS_ERROR_IF(mReferenceChord < eps)
            << "The reference chord should be larger than 0." << mReferenceChord << std::endl;

    }

    AdjointLiftJumpCoordinatesResponseFunction::~AdjointLiftJumpCoordinatesResponseFunction(){}

    void AdjointLiftJumpCoordinatesResponseFunction::InitializeSolutionStep() {
        // Get pointer to element that contains the traced node
        this->GetNeighboringElementPointer();
    }

    double AdjointLiftJumpCoordinatesResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        Element tracedElem = rModelPart.GetElement(mpNeighboringElement->Id());
        const array_1d<double, 3> v_inf = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_VELOCITY);
        double free_stream_velocity_norm = norm_2(v_inf);
        unsigned int NumNodes = tracedElem.GetGeometry().size();
        double lift_coefficient=0.0;
        for(IndexType i = 0; i < NumNodes; ++i)
        {
            if(tracedElem.GetGeometry()[i].GetValue(TRAILING_EDGE))
            {
                double potential = tracedElem.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
                double aux_potential = tracedElem.GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
                double potential_jump = std::abs(potential - aux_potential);

                lift_coefficient = 2.0 * potential_jump / (free_stream_velocity_norm * mReferenceChord);
            }
        }

        return lift_coefficient;

        KRATOS_CATCH("");
    }

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
            const array_1d<double, 3> v_inf = rProcessInfo.GetValue(FREE_STREAM_VELOCITY);
            double v_norm = norm_2(v_inf);
            double derivative = 2.0 / (v_norm * mReferenceChord);
            unsigned int NumNodes = rAdjointElement.GetGeometry().size();
            for(IndexType i = 0; i < NumNodes; ++i)
            {
                if(rAdjointElement.GetGeometry()[i].GetValue(TRAILING_EDGE))
                {
                    rResponseGradient[i] = derivative;
                    rResponseGradient[i+NumNodes] = -derivative;
                    return;
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
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }


    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::GetNeighboringElementPointer()
    {
        KRATOS_TRY;

        for (auto elem_it = mrModelPart.Elements().ptr_begin(); elem_it != mrModelPart.Elements().ptr_end(); ++elem_it)
        {
            if ((*elem_it)->Is(STRUCTURE)){
                mpNeighboringElement = (*elem_it);
                return;
            }
        }
        KRATOS_ERROR << "No neighboring element is available for the traced node." << std::endl;

        KRATOS_CATCH("");
    }
} // namespace Kratos.


