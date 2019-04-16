//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klauss B Sautter
//                   Philipp Bucher
//

#if !defined(KRATOS_STATIC_CONDENSATION_UTILITY_H_INCLUDED )
#define  KRATOS_STATIC_CONDENSATION_UTILITY_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"


namespace Kratos
{

    /**
     * @namespace StaticCondensationUtility
     *
     * @brief This utilitiy condenses given degrees of freedom from any element stiffness matrix to model e.g. hinges
     *
     * @author Klaus B Sautter, Philipp Bucher
     */



	namespace StaticCondensationUtility
	{
		typedef Element ElementType;
		typedef std::size_t SizeType;
		typedef Matrix MatrixType;

        /**
         * @brief This function is the main operation of this utility. It sorts the reference matrix w.r.t. the given dofs and condenses the reference matrix by using the following inputs:
         * @param rTheElement The current element
         * @param rLeftHandSideMatrix The matrix which will be condensed
         * @param rDofList The list containing the dofs to be condensed
         */

		void CondenseLeftHandSide(
			ElementType& rTheElement,
			MatrixType& rLeftHandSideMatrix,
			const std::vector<int> & rDofList);

        /**
         * @brief This function calculates the 4 schur-complements linking the dofs to be condensed to the dofs to remain by using the following inputs:
         * @param rTheElement The current element
         * @param rLeftHandSideMatrix The matrix which will be condensed
         * @param rDofList The list containing the dofs to be condensed
         */

		std::vector<MatrixType> CalculateSchurComplements(
			ElementType& rTheElement,
			const MatrixType& rLeftHandSideMatrix,
			const std::vector<int> & rDofList);

        /**
         * @brief This function creates a list containing all dofs to remain by using the following inputs:
         * @param rTheElement The current element
         * @param rDofList The list containing the dofs to be condensed
         */

		std::vector<int> CreateRemainingDofList(
			ElementType& rTheElement,
			const std::vector<int> & rDofList);

        /**
         * @brief This function creates the single schur-complements, called by CalculateSchurComplements, by using the following inputs:
         * @param Submatrix The current submatrix to be filled (schur-complement i)
         * @param rLeftHandSideMatrix The matrix which will be condensed
         * @param rVecA RemainingDofs or CondensedDof (according to schur-complement i)
         * @param rVecB RemainingDofs or CondensedDof (according to schur-complement i)
         * @param rSizeA Size of RemainingDofs or CondensedDof (according to schur-complement i)
         * @param rSizeB Size of RemainingDofs or CondensedDof (according to schur-complement i)
         */

		void FillSchurComplements(
			MatrixType& Submatrix,
			const MatrixType& rLeftHandSideMatrix,
			const std::vector<int>& rVecA,
			const std::vector<int>& rVecB,
			const SizeType& rSizeA,
			const SizeType& rSizeB); //maybe inline

        /**
         * @brief This function re-calculates the condensed degree of freedom in relation to the remaining dofs by using the following inputs:
         * @param rTheElement The current element
         * @param rLocalizedDofVector The localized remaining dof values
         * @rValues rValues The complete localized dof values after re-calculation
         * @rDofList rValues The list containing the dofs to be condensed
         * @rLeftHandSideMatrix The matrix which will be condensed
         */

		void ConvertingCondensation(
			ElementType& rTheElement,
			Vector& rLocalizedDofVector,
			Vector& rValues,
			const std::vector<int>& rDofList,
			const MatrixType& rLeftHandSideMatrix);

        /**
         * @brief This function returns the number of dofs of the respective element by using the following input:
         * @param rTheElement The current element
         */

		SizeType GetNumDofsElement(ElementType& rTheElement);

	}  // namespace StaticCondensationUtility

}  // namespace Kratos.

#endif // KRATOS_STATIC_CONDENSATION_UTILITY_H_INCLUDED  defined


