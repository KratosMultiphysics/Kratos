// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//
#pragma once

#include "element_variables.h"
#include "geo_mechanics_application_constants.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
class DrainagePolicy
{
public:
    virtual ~DrainagePolicy()                                                            = default;
    [[nodiscard]] virtual std::unique_ptr<DrainagePolicy<TDim, TNumNodes>> Clone() const = 0;

protected:
    using VectorType       = Vector;
    using MatrixType       = UPwSmallStrain::MatrixType;
    using ElementVariables = UPwSmallStrain::ElementVariables<TDim, TNumNodes>;

    /*virtual void CalculateAndAddLHS(MatrixType&       rLeftHandSideMatrix,
                                    ElementVariables& rVariables,
                                    void (*rCalculateAndAddStiffnessMatrix)(Matrix&, ElementVariables&),
                                    void (*rCalculateAndAddCompressibilityMatrix)(Matrix&, ElementVariables&),
                                    void (*rCalculateAndAddCouplingMatrix)(Matrix&, ElementVariables&),
                                    void (*rCalculateAndAddPermeabilityMatrix)(Matrix&, ElementVariables&)) const = 0;*/

    virtual void CalculateAndAddLHS(
        MatrixType&                                      rLeftHandSideMatrix,
        ElementVariables&                                rVariables,
        std::function<void(Matrix&, ElementVariables&)>& rCalculateAndAddStiffnessMatrix,
        std::function<void(Matrix&, ElementVariables&)>& rCalculateAndAddCompressibilityMatrix,
        std::function<void(Matrix&, ElementVariables&)>& rCalculateAndAddCouplingMatrix,
        std::function<void(Matrix&, ElementVariables&)>& rCalculateAndAddPermeabilityMatrix) const = 0;

    virtual void CalculateAndAddRHS(
        VectorType&       rRightHandSideVector,
        ElementVariables& rVariables,
        unsigned int      GPoint,
        void (*CalculateAndAddStiffnessForce)(VectorType&, ElementVariables&, unsigned int),
        void (*CalculateAndAddMixBodyForce)(VectorType&, ElementVariables&),
        void (*CalculateAndAddCouplingTerms)(VectorType&, ElementVariables&),
        void (*CalculateAndAddCompressibilityFlow)(VectorType&, ElementVariables&),
        void (*CalculateAndAddPermeabilityFlow)(VectorType&, ElementVariables&),
        void (*CalculateAndAddFluidBodyFlow)(VectorType&, ElementVariables&)) = 0;
};

} // namespace Kratos