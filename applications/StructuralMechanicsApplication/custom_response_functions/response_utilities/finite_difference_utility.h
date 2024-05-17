// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/** \brief FiniteDifferenceUtility
 *
 * This class calculates the derivatives of different element quantities (e.g. RHS, LHS, mass-matrix, ...)
 * with respect to a design variable (e.g. nodal-coordinate, property).
 */


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) FiniteDifferenceUtility
{
public:

    typedef Variable<double> array_1d_component_type;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    static void CalculateRightHandSideDerivative(Element& rElement,
                                                const Vector& rRHS,
                                                const Variable<double>& rDesignVariable,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo);

    template <typename TElementType>
    static void CalculateRightHandSideDerivative(TElementType& rElement,
                                                const Vector& rRHS,
                                                const array_1d_component_type& rDesignVariable,
                                                Node& rNode,
                                                const double& rPertubationSize,
                                                Vector& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if( rDesignVariable == SHAPE_SENSITIVITY_X || rDesignVariable == SHAPE_SENSITIVITY_Y || rDesignVariable == SHAPE_SENSITIVITY_Z )
        {
            const IndexType coord_dir =
                FiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

            // define working variables
            Vector RHS_perturbed;

            if (rOutput.size() != rRHS.size())
                rOutput.resize(rRHS.size(), false);

            // perturb the design variable
            rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
            rNode.Coordinates()[coord_dir] += rPertubationSize;

            // compute LHS after perturbation
            rElement.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);

            // compute derivative of RHS w.r.t. design variable with finite differences
            noalias(rOutput) = (RHS_perturbed - rRHS) / rPertubationSize;

            // unperturb the design variable
            rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
            rNode.Coordinates()[coord_dir] -= rPertubationSize;
        }
        if( rDesignVariable == TEMPERATURE )
        {
            //const IndexType coord_dir =
            //    FiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

            // define working variables
            Vector RHS_perturbed;

            if (rOutput.size() != rRHS.size())
                rOutput.resize(rRHS.size(), false);
            //std::cout << " In finite difference utility header" << std::endl;
            //std::cout << "design variable is " << rDesignVariable << std::endl;
            //std::cout << "unperturbed rhs is " << rRHS << std::endl;
            //std::cout << "unperturbed temp is " << rNode.FastGetSolutionStepValue(rDesignVariable) << std::endl;
            // perturb the design variable
            rNode.FastGetSolutionStepValue(rDesignVariable) += rPertubationSize;
            //double& temp = rNode.FastGetSolutionStepValue(rDesignVariable);
            //rNode.SetValue(rDesignVariable, (temp+rPertubationSize));
            //rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
            //rNode.Coordinates()[coord_dir] += rPertubationSize;
            //std::cout << "perturbed temp is " << rNode.FastGetSolutionStepValue(rDesignVariable) << std::endl;
            // compute LHS after perturbation
            rElement.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);

            //std::cout << "perturbed rhs is " << RHS_perturbed << std::endl;
            
            // compute derivative of RHS w.r.t. design variable with finite differences
            noalias(rOutput) = (RHS_perturbed - rRHS) / rPertubationSize;
            //std::cout << "difference of rhs is " << (RHS_perturbed - rRHS) << std::endl;
            // unperturb the design variable
            rNode.FastGetSolutionStepValue(rDesignVariable) -= rPertubationSize;
            //std::cout << "reset temp is " << rNode.FastGetSolutionStepValue(rDesignVariable) << std::endl;
            //rNode.SetValue(rDesignVariable, (temp-rPertubationSize));
            //rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
            //rNode.Coordinates()[coord_dir] -= rPertubationSize;
        }
        else
        {
            KRATOS_WARNING("FiniteDifferenceUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size() != 0) )
                rOutput.resize(0,false);
        }

        KRATOS_CATCH("");
    }

    static void CalculateLeftHandSideDerivative(Element& rElement,
                                                const Matrix& rLHS,
                                                const array_1d_component_type& rDesignVariable,
                                                Node& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo);

    static void CalculateMassMatrixDerivative(Element& rElement,
                                                const Matrix& rMassMatrix,
                                                const array_1d_component_type& rDesignVariable,
                                                Node& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo);

private:

    static std::size_t GetCoordinateDirection(const array_1d_component_type& rDesignVariable);

}; // class FiniteDifferenceUtility.



}  // namespace Kratos.


