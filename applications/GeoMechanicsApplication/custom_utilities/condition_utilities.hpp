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

#if !defined(KRATOS_CONDITION_UTILITIES )
#define  KRATOS_CONDITION_UTILITIES

// Project includes
#include "includes/element.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos
{

class ConditionUtilities
{

public:

    //----------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void CalculateNuMatrix(BoundedMatrix<double, TDim, TDim*TNumNodes>& rNu,
                                         const Matrix& NContainer,
                                         const unsigned int& GPoint)
    {
        for (unsigned int i=0; i < TDim; ++i) {
            unsigned int index = i - TDim;
            for (unsigned int j=0; j < TNumNodes; ++j) {
                index += TDim;
                rNu(i, index) = NContainer(GPoint, j);
            }
        }
    }

    //----------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    static inline void InterpolateVariableWithComponents(array_1d<double,TDim>& rVector,
                                                         const Matrix& Ncontainer,
                                                         const array_1d<double,TDim*TNumNodes>& VariableWithComponents,
                                                         const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(TDim);

        unsigned int index = 0;
        for (unsigned int i=0; i<TNumNodes; ++i) {
            for (unsigned int idim=0; idim<TDim; ++idim) {
                rVector[idim] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            }
        }
    }

    //----------------------------------------------------------------------------------------
    static inline void GetDisplacementsVector(array_1d<double,4>& rDisplacementVector,
                                              const Element::GeometryType& Geom)
    {
        //Line_2d_2
        array_1d<double,3> DisplacementAux;
        unsigned int index = 0;
        for (unsigned int i=0; i<2; ++i) {
            noalias(DisplacementAux) = Geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            rDisplacementVector[index++] = DisplacementAux[0];
            rDisplacementVector[index++] = DisplacementAux[1];
        }
    }

    //----------------------------------------------------------------------------------------
    static inline void GetDisplacementsVector(array_1d<double,12>& rDisplacementVector,
                                              const Element::GeometryType& Geom)
    {
        //Quadrilateral_3d_4
        array_1d<double,3> DisplacementAux;
        unsigned int index = 0;
        for (unsigned int i=0; i<4; ++i) {
            noalias(DisplacementAux) = Geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            rDisplacementVector[index++] = DisplacementAux[0];
            rDisplacementVector[index++] = DisplacementAux[1];
            rDisplacementVector[index++] = DisplacementAux[2];
        }
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template<unsigned int TNumNodes >
    static inline void GetFaceLoadVector(array_1d<double,3*TNumNodes>& rFaceLoadVector,
                                         const Element::GeometryType& Geom)
    {

        // for 3D geometry
        const unsigned int TDim = 3;
        array_1d<double,3> FaceLoadAux;
        unsigned int index = 0;
        for (unsigned int i=0; i<TNumNodes; ++i) {
            noalias(FaceLoadAux) = Geom[i].FastGetSolutionStepValue(SURFACE_LOAD);
            for (unsigned int idim=0; idim<TDim; ++idim) {
                rFaceLoadVector[index++] = FaceLoadAux[idim];
            }
        }
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    template<unsigned int TNumNodes >
    static inline void GetFaceLoadVector(array_1d<double,2*TNumNodes>& rFaceLoadVector,
                                         const Element::GeometryType& Geom)
    {

        // for 2D geometry
        const unsigned int TDim = 2;
        array_1d<double,3> FaceLoadAux;
        unsigned int index = 0;
        for (unsigned int i=0; i<TNumNodes; ++i) {
            noalias(FaceLoadAux) = Geom[i].FastGetSolutionStepValue(LINE_LOAD);
            for (unsigned int idim=0; idim < TDim; ++idim) {
                rFaceLoadVector[index++] = FaceLoadAux[idim];
            }
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static double CalculateIntegrationCoefficient(const Matrix& rJacobian, double Weight)
    {
        auto normal_vector = Vector{TDim, 0.0};

        if constexpr (TDim == 2)
        {
            normal_vector = column(rJacobian, 0);
        }
        else if constexpr (TDim == 3)
        {
            MathUtils<>::CrossProduct(normal_vector, column(rJacobian, 0),
                                      column(rJacobian, 1));
        }
        return Weight * MathUtils<>::Norm(normal_vector);
    }

}; /* Class ConditionUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_CONDITION_UTILITIES defined */
