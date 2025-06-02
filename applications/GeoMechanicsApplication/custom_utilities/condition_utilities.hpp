// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "geo_mechanics_application_variables.h"
#include "includes/element.h"

namespace Kratos
{

class ConditionUtilities
{
public:
    template <unsigned int TDim, unsigned int TNumNodes>
    static inline void CalculateNuMatrix(BoundedMatrix<double, TDim, TDim * TNumNodes>& rNu,
                                         const Matrix&                                  NContainer,
                                         const unsigned int&                            GPoint)
    {
        for (unsigned int i = 0; i < TDim; ++i) {
            unsigned int index = i - TDim;
            for (unsigned int j = 0; j < TNumNodes; ++j) {
                index += TDim;
                rNu(i, index) = NContainer(GPoint, j);
            }
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static inline void InterpolateVariableWithComponents(array_1d<double, TDim>& rVector,
                                                         const Matrix&           Ncontainer,
                                                         const array_1d<double, TDim * TNumNodes>& VariableWithComponents,
                                                         const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(TDim);

        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            for (unsigned int idim = 0; idim < TDim; ++idim) {
                rVector[idim] += Ncontainer(GPoint, i) * VariableWithComponents[index++];
            }
        }
    }

    static inline void GetDisplacementsVector(array_1d<double, 4>&         rDisplacementVector,
                                              const Element::GeometryType& rGeom)
    {
        // Line_2d_2
        for (unsigned int node = 0; node < 2; ++node) {
            std::copy_n(rGeom[node].FastGetSolutionStepValue(DISPLACEMENT).begin(), 2,
                        rDisplacementVector.begin() + node * 2);
        }
    }

    static inline void GetDisplacementsVector(array_1d<double, 12>&        rDisplacementVector,
                                              const Element::GeometryType& rGeom)
    {
        // Quadrilateral_3d_4
        for (unsigned int node = 0; node < 4; ++node) {
            std::copy_n(rGeom[node].FastGetSolutionStepValue(DISPLACEMENT).begin(), 3,
                        rDisplacementVector.begin() + node * 3);
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static inline void GetFaceLoadVector(array_1d<double, TDim * TNumNodes>& rFaceLoadVector,
                                         const Element::GeometryType&        rGeom)
    {
        const auto& variable = TDim == 2u ? LINE_LOAD : SURFACE_LOAD;

        for (unsigned int node = 0; node < TNumNodes; ++node) {
            std::copy_n(rGeom[node].FastGetSolutionStepValue(variable).begin(), TDim,
                        rFaceLoadVector.begin() + node * TDim);
        }
    }

    static double CalculateIntegrationCoefficient(const Matrix& rJacobian, double Weight)
    {
        auto vector = Vector{rJacobian.size1(), 0.0};
        if (rJacobian.size1() == 2) {
            vector = column(rJacobian, 0);
        } else if (rJacobian.size1() == 3) {
            MathUtils<>::CrossProduct(vector, column(rJacobian, 0), column(rJacobian, 1));
        }
        return Weight * MathUtils<>::Norm(vector);
    }

}; /* Class ConditionUtilities*/
} /* namespace Kratos.*/
