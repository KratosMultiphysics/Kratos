//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Klauss B Sautter
//                   Philipp Bucher
//                   Vahid Galavi

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{

/**
 * @class GeoStaticCondensationUtility
 *
 * @brief This utility condenses given degrees of freedom from any element stiffness matrix to model e.g. hinges
 *
 */

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoStaticCondensationUtility
{
public:
    /**
     * @brief This function is the main operation of this utility. It sorts the reference matrix w.r.t. the given dofs and condenses the reference matrix by using the following inputs:
     * @param rTheElement The current element
     * @param rLeftHandSideMatrix The matrix which will be condensed
     * @param rDofList The list containing the dofs to be condensed
     */

    static void CondenseLeftHandSide(Element&                rTheElement,
                                     Matrix&                 rLeftHandSideMatrix,
                                     const std::vector<int>& rDofList);

    /**
     * @brief This function calculates the 4 schur-complements linking the dofs to be condensed to the dofs to remain by using the following inputs:
     * @param rTheElement The current element
     * @param rLeftHandSideMatrix The matrix which will be condensed
     * @param rDofList The list containing the dofs to be condensed
     */

    static std::vector<Matrix> CalculateSchurComplements(Element&      rTheElement,
                                                         const Matrix& rLeftHandSideMatrix,
                                                         const std::vector<int>& rDofList);

    /**
     * @brief This function creates a list containing all dofs to remain by using the following inputs:
     * @param rTheElement The current element
     * @param rDofList The list containing the dofs to be condensed
     */

    static std::vector<int> CreateRemainingDofList(Element& rTheElement, const std::vector<int>& rDofList);

    /**
     * @brief This function creates the single schur-complements, called by CalculateSchurComplements, by using the following inputs:
     * @param Submatrix The current submatrix to be filled (schur-complement i)
     * @param rLeftHandSideMatrix The matrix which will be condensed
     * @param rVecA RemainingDofs or CondensedDof (according to schur-complement i)
     * @param rVecB RemainingDofs or CondensedDof (according to schur-complement i)
     * @param rSizeA Size of RemainingDofs or CondensedDof (according to schur-complement i)
     * @param rSizeB Size of RemainingDofs or CondensedDof (according to schur-complement i)
     */

    static void FillSchurComplements(Matrix&                 Submatrix,
                                     const Matrix&           rLeftHandSideMatrix,
                                     const std::vector<int>& rVecA,
                                     const std::vector<int>& rVecB,
                                     const std::size_t&      rSizeA,
                                     const std::size_t&      rSizeB);

    /**
     * @brief This function re-calculates the condensed degree of freedom in relation to the remaining dofs by using the following inputs:
     * @param rTheElement The current element
     * @param rLocalizedDofVector The localized remaining dof values
     * @rValues rValues The complete localized dof values after re-calculation
     * @rDofList rValues The list containing the dofs to be condensed
     * @rLeftHandSideMatrix The matrix which will be condensed
     */

    static void ConvertingCondensation(Element&                rTheElement,
                                       Vector&                 rLocalizedDofVector,
                                       Vector&                 rValues,
                                       const std::vector<int>& rDofList,
                                       const Matrix&           rLeftHandSideMatrix);

    /**
     * @brief This function returns the number of dofs of the respective element by using the following input:
     * @param rTheElement The current element
     */

    static std::size_t GetNumDofsElement(Element& rTheElement);

}; // class StaticCondensationUtility

} // namespace Kratos.
