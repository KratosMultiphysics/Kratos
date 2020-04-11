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
#include "adjoint_local_stress_response_function.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_base_element.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
    AdjointLocalStressResponseFunction::AdjointLocalStressResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        // Get traced element
        const int id_of_traced_element = ResponseSettings["traced_element_id"].GetInt();
        mpTracedElement = rModelPart.pGetElement(id_of_traced_element);

        // Tell traced element the stress type
        mTracedStressType = StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"].GetString());
        mpTracedElement->SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType) );

        // Get info how and where to treat the stress
        mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment( ResponseSettings["stress_treatment"].GetString() );

        if(mStressTreatment == StressTreatment::GaussPoint || mStressTreatment == StressTreatment::Node) {
            mIdOfLocation = ResponseSettings["stress_location"].GetInt();
            KRATOS_ERROR_IF(mIdOfLocation < 1) << "Chose a 'stress_location' > 0. Specified 'stress_location': " << mIdOfLocation << std::endl;
        }

        if(ResponseSettings.Has("add_particular_solution")) {
            mAddParticularSolution = ResponseSettings["add_particular_solution"].GetBool();
        }
    }

    AdjointLocalStressResponseFunction::~AdjointLocalStressResponseFunction(){}

    void AdjointLocalStressResponseFunction::FinalizeSolutionStep()
    {
        KRATOS_TRY;

        if(mAddParticularSolution) {
            this->CalculateParticularSolution();
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if(rAdjointElement.Id() == mpTracedElement->Id())
        {
            Matrix stress_displacement_derivative;

            if(mStressTreatment == StressTreatment::Mean)
            {
                mpTracedElement->Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
                this->ExtractMeanStressDerivative(stress_displacement_derivative, rResponseGradient);
            }
            else if(mStressTreatment == StressTreatment::GaussPoint)
            {
                mpTracedElement->Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
                this->ExtractGaussPointStressDerivative(stress_displacement_derivative, rResponseGradient);
            }
            else if(mStressTreatment == StressTreatment::Node)
            {
                mpTracedElement->Calculate(STRESS_DISP_DERIV_ON_NODE, stress_displacement_derivative, rProcessInfo);
                this->ExtractNodeStressDerivative(stress_displacement_derivative, rResponseGradient);
            }

            KRATOS_ERROR_IF(rResponseGradient.size() != rResidualGradient.size1())
                 << "Size of stress displacement derivative does not fit!" << std::endl;

            rResponseGradient *= (-1);
        }
        else
        {
            if(rResponseGradient.size() != rResidualGradient.size1())
                rResponseGradient.resize(rResidualGradient.size1(), false);

            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if(rAdjointElement.Id() == mpTracedElement->Id())
        {
            ProcessInfo process_info = rProcessInfo;
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix,
                                                                    rSensitivityGradient, process_info);
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("")
    }

    void AdjointLocalStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if(rAdjointElement.Id() == mpTracedElement->Id())
        {
            ProcessInfo process_info = rProcessInfo;
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix,
                                                                    rSensitivityGradient, process_info);
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    double AdjointLocalStressResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        double stress_value = 0.0;

        if(mStressTreatment == StressTreatment::Mean)
            stress_value = CalculateMeanElementStress(rModelPart);
        else if (mStressTreatment == StressTreatment::GaussPoint)
            stress_value = CalculateGaussPointStress(rModelPart);
        else if (mStressTreatment == StressTreatment::Node)
            stress_value = CalculateNodeStress(rModelPart);
        return stress_value;

        KRATOS_CATCH("");
    }

    double AdjointLocalStressResponseFunction::CalculateMeanElementStress(ModelPart& rModelPart)
    {
        double stress_value = 0.0;

        Vector element_stress;

        Element& r_element = rModelPart.GetElement(mpTracedElement->Id());
        StressCalculation::CalculateStressOnGP(r_element, mTracedStressType, element_stress, rModelPart.GetProcessInfo());

        const SizeType stress_vec_size = element_stress.size();

        for(IndexType i = 0; i < stress_vec_size; ++i)
            stress_value += element_stress[i];

        stress_value /= stress_vec_size;

        return stress_value;
    }

    double AdjointLocalStressResponseFunction::CalculateGaussPointStress(ModelPart& rModelPart)
    {
        Vector element_stress;

        Element& r_element = rModelPart.GetElement(mpTracedElement->Id());
        StressCalculation::CalculateStressOnGP(r_element, mTracedStressType, element_stress, rModelPart.GetProcessInfo());

        const SizeType stress_vec_size = element_stress.size();

        if(stress_vec_size >= mIdOfLocation)
            return element_stress[mIdOfLocation - 1];

        KRATOS_ERROR << "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                        stress_vec_size  << "!"<< std::endl;
    }

    double AdjointLocalStressResponseFunction::CalculateNodeStress(ModelPart& rModelPart)
    {
        Vector element_stress;

        Element& r_element = rModelPart.GetElement(mpTracedElement->Id());
        StressCalculation::CalculateStressOnNode(r_element, mTracedStressType, element_stress, rModelPart.GetProcessInfo());

        const SizeType num_ele_nodes = mpTracedElement->GetGeometry().PointsNumber();

        if(num_ele_nodes >= mIdOfLocation)
            return element_stress[mIdOfLocation - 1];

        KRATOS_ERROR << "Chosen Node is not available. The element has only " <<
                        num_ele_nodes  << " nodes."<< std::endl;

    }

    void AdjointLocalStressResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rAdjointElement,
                                      const std::string& rVariableName,
                                      const Matrix& rSensitivityMatrix,
                                      Vector& rSensitivityGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, rVariableName);

        Matrix stress_design_variable_derivative;

        if(mStressTreatment == StressTreatment::Mean)
        {
            rAdjointElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            if (stress_design_variable_derivative.size1() == 0)
                rSensitivityGradient = ZeroVector(0);
            else
                this->ExtractMeanStressDerivative(stress_design_variable_derivative, rSensitivityGradient);
        }
        else if(mStressTreatment == StressTreatment::GaussPoint)
        {
            rAdjointElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            if (stress_design_variable_derivative.size1() == 0)
                rSensitivityGradient = ZeroVector(0);
            else
                this->ExtractGaussPointStressDerivative(stress_design_variable_derivative, rSensitivityGradient);
        }
        else if(mStressTreatment == StressTreatment::Node)
        {
            rAdjointElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_NODE, stress_design_variable_derivative, rProcessInfo);
            if (stress_design_variable_derivative.size1() == 0)
                rSensitivityGradient = ZeroVector(0);
            else
                this->ExtractNodeStressDerivative(stress_design_variable_derivative, rSensitivityGradient);
        }

        KRATOS_ERROR_IF(rSensitivityGradient.size() != rSensitivityMatrix.size1())
             << "Size of partial stress design variable derivative does not fit!" << std::endl;

        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, "");

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
    {
        KRATOS_TRY;

        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();
        double stress_derivative_value = 0.0;

        if(rResponseGradient.size() != num_of_derivatives_per_stress)
            rResponseGradient.resize(num_of_derivatives_per_stress, false);

        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
        {
            for(IndexType stress_it = 0; stress_it < num_of_stress_positions; ++stress_it)
                stress_derivative_value += rStressDerivativesMatrix(deriv_it, stress_it);

            stress_derivative_value /= num_of_stress_positions;

            rResponseGradient[deriv_it] = stress_derivative_value;
            stress_derivative_value = 0.0;
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::ExtractNodeStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
    {
        KRATOS_TRY;

        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();

        if(rResponseGradient.size() != num_of_derivatives_per_stress)
            rResponseGradient.resize(num_of_derivatives_per_stress, false);

        KRATOS_ERROR_IF_NOT(num_of_stress_positions >= mIdOfLocation ) << "Chosen node is not available. The element has only " <<
                            num_of_stress_positions << " nodes."<< std::endl;

        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
            rResponseGradient[deriv_it] = rStressDerivativesMatrix(deriv_it, (mIdOfLocation-1));

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::ExtractGaussPointStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
    {
        KRATOS_TRY;

        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();

        if(rResponseGradient.size() != num_of_derivatives_per_stress)
            rResponseGradient.resize(num_of_derivatives_per_stress, false);

        KRATOS_ERROR_IF_NOT(num_of_stress_positions >= mIdOfLocation ) <<
                "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                                num_of_stress_positions  << "!"<< std::endl;

        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
            rResponseGradient[deriv_it] = rStressDerivativesMatrix(deriv_it, (mIdOfLocation-1));

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateParticularSolution() const
    {
        KRATOS_TRY;

        std::string element_name;
        CompareElementsAndConditionsUtility::GetRegisteredName(*mpTracedElement, element_name);

        if(element_name == "AdjointFiniteDifferenceCrBeamElementLinear3D2N" || element_name == "AdjointFiniteDifferenceTrussLinearElement3D2N") {
            Vector particular_solution;
            // delivers particular solution of influence function in local coordinates
            this->CalculateParticularSolutionLinearElement2N(particular_solution);
            // transform particular solution into global coordinates
            Matrix transformation_matrix = ZeroMatrix(particular_solution.size());
            Matrix local_element_orientation = ZeroMatrix(3);
            mpTracedElement->Calculate(LOCAL_ELEMENT_ORIENTATION, local_element_orientation, mrModelPart.GetProcessInfo());

            const SizeType dimension(3);
            KRATOS_ERROR_IF_NOT((particular_solution.size()%dimension)==0) << "Size of particular solution does not fit!" << std::endl;
            const SizeType check_a(particular_solution.size()/dimension);

            SizeType iterator_count(0);
            for (SizeType i=0;  i<check_a;++i){
                iterator_count = i*dimension;
                project(transformation_matrix, range(iterator_count,iterator_count+dimension),range(iterator_count,iterator_count+dimension)) += local_element_orientation;
            }

            particular_solution = prod(transformation_matrix, particular_solution);
            // set particular solution as non-historical result
            mpTracedElement->SetValue(ADJOINT_PARTICULAR_DISPLACEMENT, particular_solution);
        } else {
            KRATOS_ERROR << "CalculateParticularSolution not available for " << element_name << "!" << std::endl;
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateParticularSolutionLinearElement2N(Vector& rResult) const
    {
        KRATOS_TRY;

        DofsVectorType dofs_of_element;
        mpTracedElement->GetDofList(dofs_of_element, mrModelPart.GetProcessInfo());
        rResult = ZeroVector(dofs_of_element.size());

        Array1DComponentsPointerType p_traced_dof;
        this->FindVariableComponent(p_traced_dof);

        if(mStressTreatment == StressTreatment::Mean) {
            this->CalculateMeanParticularSolutionLinearElement2N(rResult, dofs_of_element, p_traced_dof);
        } else if(mStressTreatment == StressTreatment::GaussPoint) {
            this->CalculateGPParticularSolutionLinearElement2N(rResult, dofs_of_element, p_traced_dof);
        } else if(mStressTreatment == StressTreatment::Node) {
            this->CalculateNodeParticularSolutionLinearElement2N(rResult, dofs_of_element, p_traced_dof);
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateMeanParticularSolutionLinearElement2N(Vector& rResult,
                        DofsVectorType &rElementalDofList, const Array1DComponentsPointerType TracedDof) const
    {
        KRATOS_TRY;

        if(rResult.size() != rElementalDofList.size()) {
            rResult.resize(rElementalDofList.size());
        }

        const unsigned int num_GP = mpTracedElement->GetGeometry().IntegrationPointsNumber(mpTracedElement->GetIntegrationMethod());
        const double prefactor = 1.0 / (1.0 + num_GP);

        const IndexType id_node_1 = mpTracedElement->GetGeometry()[0].Id();
        const IndexType id_node_2 = mpTracedElement->GetGeometry()[1].Id();

        for(IndexType gp_it = 0; gp_it < num_GP; ++gp_it) {
            for(IndexType i = 0; i < rElementalDofList.size(); ++i) {
                if (rElementalDofList[i]->GetVariable() == *TracedDof) {
                    if (rElementalDofList[i]->Id() == id_node_1) {
                        rResult[i] += prefactor * (num_GP - gp_it);
                    } else if (rElementalDofList[i]->Id() == id_node_2) {
                        rResult[i] += -prefactor * (gp_it + 1);
                    }
                }
            }
        }
        rResult /= num_GP;

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateGPParticularSolutionLinearElement2N(Vector& rResult,
                            DofsVectorType &rElementalDofList, const Array1DComponentsPointerType TracedDof) const
    {
        KRATOS_TRY;

        if(rResult.size() != rElementalDofList.size()) {
            rResult.resize(rElementalDofList.size());
        }

        const unsigned int num_GP = mpTracedElement->GetGeometry().IntegrationPointsNumber(mpTracedElement->GetIntegrationMethod());
        const double prefactor = 1.0 / (1.0 + num_GP);

        const IndexType id_node_1 = mpTracedElement->GetGeometry()[0].Id();
        const IndexType id_node_2 = mpTracedElement->GetGeometry()[1].Id();

        for(IndexType i = 0; i < rElementalDofList.size(); ++i) {
            if (rElementalDofList[i]->GetVariable() == *TracedDof) {
                if (rElementalDofList[i]->Id() == id_node_1) {
                    rResult[i] = prefactor * (num_GP + 1 - mIdOfLocation);
                } else if (rElementalDofList[i]->Id() == id_node_2) {
                    rResult[i] = -prefactor * mIdOfLocation;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateNodeParticularSolutionLinearElement2N(Vector& rResult,
                     DofsVectorType &rElementalDofList, const Array1DComponentsPointerType TracedDof) const
    {
        KRATOS_TRY;

        if(rResult.size() != rElementalDofList.size()) {
            rResult.resize(rElementalDofList.size());
        }

        for(IndexType i = 0; i < rElementalDofList.size(); ++i) {
            if (rElementalDofList[i]->Id() == mpTracedElement->GetGeometry()[mIdOfLocation-1].Id() &&
                rElementalDofList[i]->GetVariable() == *TracedDof) {
                if (mIdOfLocation == 1) {
                    rResult[i] = 1.0;
                } else if (mIdOfLocation == 2) {
                    rResult[i] = -1.0;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::FindVariableComponent(Array1DComponentsPointerType& rTracedDof) const
    {
        KRATOS_TRY;

        switch (mTracedStressType)
        {
            case TracedStressType::MX:
            {
                rTracedDof = &ADJOINT_ROTATION_X;
                break;
            }
            case TracedStressType::MY:
            {
                rTracedDof = &ADJOINT_ROTATION_Y;
                break;
            }
            case TracedStressType::MZ:
            {
                rTracedDof = &ADJOINT_ROTATION_Z;
                break;
            }
            case TracedStressType::FX:
            {
                rTracedDof = &ADJOINT_DISPLACEMENT_X;
                break;
            }
            case TracedStressType::FY:
            {
                rTracedDof = &ADJOINT_DISPLACEMENT_Y;
                break;
            }
            case TracedStressType::FZ:
            {
                rTracedDof = &ADJOINT_DISPLACEMENT_Z;
                break;
            }
            default:
                KRATOS_ERROR << "Invalid stress type! Stress type not supported for particular solution!" << std::endl;
        }

        KRATOS_CATCH("");
    }

} // namespace Kratos.

