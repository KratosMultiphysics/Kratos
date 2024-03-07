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

#include "drainage_policy.h"
#include "element_variables.h"
#include "geo_mechanics_application_constants.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
class UndrainedCase : public DrainagePolicy<TDim, TNumNodes>
{
public:
    std::unique_ptr<DrainagePolicy<TDim, TNumNodes>> Clone() const
    {
        return std::make_unique<UndrainedCase<TDim, TNumNodes>>();
    }

protected:
    using VectorType       = Vector;
    using MatrixType       = Matrix;
    using ElementVariables = UPwSmallStrain::ElementVariables<TDim, TNumNodes>;

    void CalculateAndAddLHS(
        MatrixType&       rLeftHandSideMatrix,
        ElementVariables& rVariables,
        /*void              (*rCalculateAndAddStiffnessMatrix)(Matrix&, ElementVariables&),
        void              (*rCalculateAndAddCompressibilityMatrix)(Matrix&, ElementVariables&),
        void              (*rCalculateAndAddCouplingMatrix)(Matrix&, ElementVariables&),
        void              (*rCalculateAndAddPermeabilityMatrix)(Matrix&, ElementVariables&)*/
        std::function<void(MatrixType&, ElementVariables&)>& rCalculateAndAddStiffnessMatrix,
        std::function<void(MatrixType&, ElementVariables&)>& rCalculateAndAddCompressibilityMatrix,
        std::function<void(MatrixType&, ElementVariables&)>& rCalculateAndAddCouplingMatrix,
        std::function<void(MatrixType&, ElementVariables&)>& rCalculateAndAddPermeabilityMatrix) const override
    {
        KRATOS_TRY;
        rCalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);
        rCalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);
        rCalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);
        KRATOS_CATCH("");
    }

    void CalculateAndAddRHS(VectorType&       rRightHandSideVector,
                            ElementVariables& rVariables,
                            unsigned int      GPoint,
                            void (*CalculateAndAddStiffnessForce)(VectorType&, ElementVariables&, unsigned int),
                            void (*CalculateAndAddMixBodyForce)(VectorType&, ElementVariables&),
                            void (*CalculateAndAddCouplingTerms)(VectorType&, ElementVariables&),
                            void (*CalculateAndAddCompressibilityFlow)(VectorType&, ElementVariables&),
                            void (*CalculateAndAddPermeabilityFlow)(VectorType&, ElementVariables&),
                            void (*CalculateAndAddFluidBodyFlow)(VectorType&, ElementVariables&)) override
    {
        KRATOS_TRY;
        CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables, GPoint);
        CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);
        CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);
        CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);
        KRATOS_CATCH("");
    }
};

} // namespace Kratos