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

namespace Kratos
{
    AdjointLocalStressResponseFunction::AdjointLocalStressResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        ResponseData stress_response_data;

        // Get traced element
        const int id_of_traced_element = ResponseSettings["traced_element_id"].GetInt();
        mpTracedElement = rModelPart.pGetElement(id_of_traced_element);

        // Tell traced element the stress type
        TracedStressType traced_stress_type = stress_response_data.ConvertStressType(ResponseSettings["stress_type"].GetString());
        mpTracedElement->SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress_type) );

        // Get info how and where to treat the stress
        mStressTreatment = stress_response_data.ConvertStressTreatment( ResponseSettings["stress_treatment"].GetString() );
        KRATOS_ERROR_IF(mStressTreatment == StressTreatment::StressTreatmentNotAvailable) << "Chosen option for stress treatmeant is not available! Chose 'GP','node' or 'mean'!" << std::endl;

        if(mStressTreatment == StressTreatment::GaussPoint || mStressTreatment == StressTreatment::Node)
        {
            mIdOfLocation = ResponseSettings["stress_location"].GetInt();
            KRATOS_ERROR_IF(mIdOfLocation < 1) << "Chose a 'stress_location' > 0. Specified 'stress_location': " << mIdOfLocation << std::endl;
        }
    }

    AdjointLocalStressResponseFunction::~AdjointLocalStressResponseFunction(){}


    double AdjointLocalStressResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        double stress_value = 0.0;

        // Working variables
        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
        Vector element_stress;

        if(mStressTreatment == StressTreatment::Mean || mStressTreatment == StressTreatment::GaussPoint)
            mpTracedElement->Calculate(STRESS_ON_GP, element_stress, r_current_process_info);
        else
            mpTracedElement->Calculate(STRESS_ON_NODE, element_stress, r_current_process_info);

        const SizeType stress_vec_size = element_stress.size();

        if(mStressTreatment == StressTreatment::Mean)
        {
            for(IndexType i = 0; i < stress_vec_size; ++i)
                stress_value += element_stress[i];

            stress_value /= stress_vec_size;
        }
        else if(mStressTreatment == StressTreatment::GaussPoint)
        {
            if(stress_vec_size >= mIdOfLocation)
                stress_value = element_stress[mIdOfLocation - 1];
            else
                KRATOS_ERROR << "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                                stress_vec_size  << "!"<< std::endl;
        }
        else if(mStressTreatment == StressTreatment::Node)
        {
            const SizeType num_ele_nodes = mpTracedElement->GetGeometry().PointsNumber();
            if(num_ele_nodes >= mIdOfLocation)
                stress_value = element_stress[mIdOfLocation - 1];
            else
                KRATOS_ERROR << "Chosen Node is not available. The element has only " <<
                                num_ele_nodes  << " nodes."<< std::endl;

        }

        return stress_value;

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateGradient(const Element& rAdjointElem, const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo)
    {
        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        if(rAdjointElem.Id() == mpTracedElement->Id())
        {
            Matrix stress_displacement_derivative;
            if(mStressTreatment == StressTreatment::Mean || mStressTreatment == StressTreatment::GaussPoint)
                mpTracedElement->Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            else
                mpTracedElement->Calculate(STRESS_DISP_DERIV_ON_NODE, stress_displacement_derivative, rProcessInfo);

            const SizeType num_of_dofs = stress_displacement_derivative.size1();
            const SizeType num_of_deriv = stress_displacement_derivative.size2();
            double stress_displ_deriv_value = 0.0;

            KRATOS_ERROR_IF(rResponseGradient.size() != stress_displacement_derivative.size1())
                 << "Size of stress displacement derivative does not fit!" << std::endl;

            for (IndexType dof_it = 0 ; dof_it < num_of_dofs; ++dof_it)
            {
                if(mStressTreatment == StressTreatment::Mean)
                {
                    for(IndexType GP_it = 0; GP_it < num_of_deriv; ++GP_it)
                        stress_displ_deriv_value += stress_displacement_derivative(dof_it, GP_it);

                    stress_displ_deriv_value /= num_of_deriv;
                }
                else if(mStressTreatment == StressTreatment::GaussPoint)
                {
                    if(num_of_deriv >= mIdOfLocation)
                        stress_displ_deriv_value = stress_displacement_derivative(dof_it, (mIdOfLocation-1));
                    else
                        KRATOS_ERROR << "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                                    num_of_deriv  << "!"<< std::endl;
                }
                else if(mStressTreatment == StressTreatment::Node)
                {
                    if(num_of_deriv >= mIdOfLocation)
                        stress_displ_deriv_value = stress_displacement_derivative(dof_it, (mIdOfLocation-1));
                    else
                        KRATOS_ERROR << "Chosen node is not available. The element has only " <<
                                    num_of_deriv  << " nodes."<< std::endl;

                }
                rResponseGradient[dof_it] = (-1) * stress_displ_deriv_value;
                stress_displ_deriv_value = 0.0;
            }
        }
    }

    void AdjointLocalStressResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<double>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
          KRATOS_TRY


        if(rAdjointElem.Id() == mpTracedElement->Id())
        {
            rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, rVariable.Name());

            Matrix stress_design_variable_derivative;
            if(mStressTreatment == StressTreatment::Mean || mStressTreatment == StressTreatment::GaussPoint)
                rAdjointElem.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            else
                rAdjointElem.Calculate(STRESS_DESIGN_DERIVATIVE_ON_NODE, stress_design_variable_derivative, rProcessInfo);

            const SizeType num_of_DV = stress_design_variable_derivative.size1();
            const SizeType num_of_deriv = stress_design_variable_derivative.size2();
            double stress_DV_deriv_value = 0.0;

            if(rResponseGradient.size() != stress_design_variable_derivative.size1())
                rResponseGradient.resize(stress_design_variable_derivative.size1(), false);
            KRATOS_ERROR_IF(rResponseGradient.size() != rDerivativesMatrix.size1())
                 << "Size of partial stress design variable derivative does not fit!" << std::endl;

            for (IndexType dv_it = 0 ; dv_it < num_of_DV; ++dv_it)
            {
                if(mStressTreatment == StressTreatment::Mean)
                {
                    for(IndexType GP_it = 0; GP_it < num_of_deriv; ++GP_it)
                        stress_DV_deriv_value += stress_design_variable_derivative(dv_it, GP_it);

                    stress_DV_deriv_value /= num_of_deriv;
                }
                else if(mStressTreatment == StressTreatment::GaussPoint)
                {
                    if(num_of_deriv >= mIdOfLocation)
                        stress_DV_deriv_value = stress_design_variable_derivative(dv_it, (mIdOfLocation-1));
                    else
                        KRATOS_ERROR << "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                                    num_of_deriv  << "!"<< std::endl;
                }
                else if(mStressTreatment == StressTreatment::Node)
                {
                    if(num_of_deriv >= mIdOfLocation)
                        stress_DV_deriv_value = stress_design_variable_derivative(dv_it, (mIdOfLocation-1));
                    else
                        KRATOS_ERROR << "Chosen node is not available. The element has only " <<
                                    num_of_deriv  << " nodes."<< std::endl;
                }
                rResponseGradient[dv_it] =  stress_DV_deriv_value;
                stress_DV_deriv_value = 0.0;
            }

            rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, "");
        }
        else
        {
            if (rResponseGradient.size() != rDerivativesMatrix.size1())
                      rResponseGradient.resize(rDerivativesMatrix.size1(), false);
            rResponseGradient.clear();
        }

        KRATOS_CATCH("")
    }

    void AdjointLocalStressResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                     const Variable<double>& rVariable,
                                     const Matrix& rDerivativesMatrix,
                                     Vector& rResponseGradient,
                                     ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
                  rResponseGradient.resize(rDerivativesMatrix.size1(), false);
        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if(rAdjointElem.Id() == mpTracedElement->Id())
        {
            rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, rVariable.Name());

            Matrix stress_design_variable_derivative;
            if(mStressTreatment == StressTreatment::Mean || mStressTreatment == StressTreatment::GaussPoint)
                rAdjointElem.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            else
                rAdjointElem.Calculate(STRESS_DESIGN_DERIVATIVE_ON_NODE, stress_design_variable_derivative, rProcessInfo);

            const SizeType num_of_DV = stress_design_variable_derivative.size1();
            const SizeType  num_of_deriv = stress_design_variable_derivative.size2();
            double stress_DV_deriv_value = 0.0;

            if(rResponseGradient.size() != stress_design_variable_derivative.size1())
                rResponseGradient.resize(stress_design_variable_derivative.size1(), false);
            KRATOS_ERROR_IF(rResponseGradient.size() != rDerivativesMatrix.size1())
                << "Size of partial stress design variable derivative does not fit!" << std::endl;

            for (IndexType dv_it = 0 ; dv_it < num_of_DV; ++dv_it)
            {
                if(mStressTreatment == StressTreatment::Mean)
                {
                    for(IndexType GP_it = 0; GP_it < num_of_deriv; ++GP_it)
                        stress_DV_deriv_value += stress_design_variable_derivative(dv_it, GP_it);

                    stress_DV_deriv_value /= num_of_deriv;
                }
                else if(mStressTreatment == StressTreatment::GaussPoint)
                {
                    if(num_of_deriv >= mIdOfLocation)
                        stress_DV_deriv_value = stress_design_variable_derivative(dv_it, (mIdOfLocation-1));
                    else
                        KRATOS_ERROR << "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                                    num_of_deriv  << "!"<< std::endl;
                }
                else if(mStressTreatment == StressTreatment::Node)
                {
                    if(num_of_deriv >= mIdOfLocation)
                        stress_DV_deriv_value = stress_design_variable_derivative(dv_it, (mIdOfLocation-1));
                    else
                        KRATOS_ERROR << "Chosen node is not available. The element has only " <<
                                    num_of_deriv  << " nodes."<< std::endl;
                }
                rResponseGradient[dv_it] = stress_DV_deriv_value;
                stress_DV_deriv_value = 0.0;
            }

            rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, "");
        }
        else
        {
            if (rResponseGradient.size() != rDerivativesMatrix.size1())
                      rResponseGradient.resize(rDerivativesMatrix.size1(), false);
            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if(rResponseGradient.size() != rDerivativesMatrix.size1())
              rResponseGradient.resize(rDerivativesMatrix.size1(), false);
        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

} // namespace Kratos.

