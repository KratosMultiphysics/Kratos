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
class SteadyStatePwDrainageCase : public DrainagePolicy<TDim, TNumNodes>
{
public:
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

        rCalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);

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

        CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);
        CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

        KRATOS_CATCH("");
    }

    std::unique_ptr<DrainagePolicy> Clone() const
    {
        return std::make_unique<SteadyStatePwDrainageCase>();
    }
};

} // namespace Kratos