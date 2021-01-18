//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

#if !defined(KRATOS_ELEMENT_UTILITIES )
#define  KRATOS_ELEMENT_UTILITIES

// Project includes
#include "includes/define.h"

// System includes
//#include <cmath>

// Project includes
//#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes

namespace Kratos
{

class ElementUtilities
{


public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void InterpolateVariableWithComponents(array_1d<double,2>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,3>, 3>& VariableWithComponents,const unsigned int& GPoint)
    {
    KRATOS_TRY

        noalias(rVector) = ZeroVector(2);

        for(unsigned int i=0; i<3; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];
        }

    KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void InterpolateVariableWithComponents(array_1d<double,2>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,3>, 4>& VariableWithComponents,const unsigned int& GPoint)
    {
    KRATOS_TRY

        noalias(rVector) = ZeroVector(2);

        for(unsigned int i=0; i<4; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];
        }

    KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,3>, 4>& VariableWithComponents,const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(3);

        for(unsigned int i=0; i<4; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[i][2];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,3>, 8>& VariableWithComponents,const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(3);

        for(unsigned int i=0; i<8; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[i][2];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void FillArray1dOutput(array_1d<double,3>& rOutputValue, const array_1d<double,2>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue[0];
        rOutputValue[1] = ComputedValue[1];
        rOutputValue[2] = 0.0;
    }

    //----------------------------------------------------------------------------------------

    static inline void FillArray1dOutput(array_1d<double,3>& rOutputValue, const array_1d<double,3>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue[0];
        rOutputValue[1] = ComputedValue[1];
        rOutputValue[2] = ComputedValue[2];
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// CalculateExtrapolationMatrix
    /// The matrix contains the shape functions at each GP evaluated at each node.
    /// Rows: nodes
    /// Columns: GP

    static inline void Calculate2DExtrapolationMatrix(BoundedMatrix<double,3,3>& rExtrapolationMatrix)
    {
        //Triangle_2d_3
        //GI_GAUSS_2

        rExtrapolationMatrix(0,0) = 1.6666666666666666666; rExtrapolationMatrix(0,1) = -0.33333333333333333333; rExtrapolationMatrix(0,2) = -0.33333333333333333333;
        rExtrapolationMatrix(1,0) = -0.33333333333333333333; rExtrapolationMatrix(1,1) = 1.6666666666666666666; rExtrapolationMatrix(1,2) = -0.33333333333333333333;
        rExtrapolationMatrix(2,0) = -0.33333333333333333333; rExtrapolationMatrix(2,1) = -0.33333333333333333333; rExtrapolationMatrix(2,2) = 1.6666666666666666666;
    }

    //----------------------------------------------------------------------------------------

    static inline void Calculate2DExtrapolationMatrix(BoundedMatrix<double,4,4>& rExtrapolationMatrix)
    {
        //Quadrilateral_2d_4
        //GI_GAUSS_2

        rExtrapolationMatrix(0,0) = 1.8660254037844386; rExtrapolationMatrix(0,1) = -0.5; rExtrapolationMatrix(0,2) = 0.13397459621556132; rExtrapolationMatrix(0,3) = -0.5;
        rExtrapolationMatrix(1,0) = -0.5; rExtrapolationMatrix(1,1) = 1.8660254037844386; rExtrapolationMatrix(1,2) = -0.5; rExtrapolationMatrix(1,3) = 0.13397459621556132;
        rExtrapolationMatrix(2,0) = 0.13397459621556132; rExtrapolationMatrix(2,1) = -0.5; rExtrapolationMatrix(2,2) = 1.8660254037844386; rExtrapolationMatrix(2,3) = -0.5;
        rExtrapolationMatrix(3,0) = -0.5; rExtrapolationMatrix(3,1) = 0.13397459621556132; rExtrapolationMatrix(3,2) = -0.5; rExtrapolationMatrix(3,3) = 1.8660254037844386;
    }

    //----------------------------------------------------------------------------------------

    static inline void Calculate3DExtrapolationMatrix(BoundedMatrix<double,4,4>& rExtrapolationMatrix)
    {
        //Tetrahedra_3d_4
        //GI_GAUSS_2

        rExtrapolationMatrix(0,0) = -0.309016988749894905; rExtrapolationMatrix(0,1) = -0.3090169887498949046; rExtrapolationMatrix(0,2) = -0.309016988749894905; rExtrapolationMatrix(0,3) = 1.9270509662496847144;
        rExtrapolationMatrix(1,0) = 1.9270509662496847144; rExtrapolationMatrix(1,1) = -0.30901698874989490481; rExtrapolationMatrix(1,2) = -0.3090169887498949049; rExtrapolationMatrix(1,3) = -0.30901698874989490481;
        rExtrapolationMatrix(2,0) = -0.30901698874989490473; rExtrapolationMatrix(2,1) = 1.9270509662496847143; rExtrapolationMatrix(2,2) = -0.3090169887498949049; rExtrapolationMatrix(2,3) = -0.30901698874989490481;
        rExtrapolationMatrix(3,0) = -0.3090169887498949048; rExtrapolationMatrix(3,1) = -0.30901698874989490471; rExtrapolationMatrix(3,2) = 1.9270509662496847143; rExtrapolationMatrix(3,3) = -0.30901698874989490481;
    }

    //----------------------------------------------------------------------------------------

    static inline void Calculate3DExtrapolationMatrix(BoundedMatrix<double,8,8>& rExtrapolationMatrix)
    {
        //Hexahedra_3d_8
        //GI_GAUSS_2

        rExtrapolationMatrix(0,0) = 2.549038105676658; rExtrapolationMatrix(0,1) = -0.6830127018922192; rExtrapolationMatrix(0,2) = 0.18301270189221927; rExtrapolationMatrix(0,3) = -0.6830127018922192;
        rExtrapolationMatrix(0,4) = -0.6830127018922192; rExtrapolationMatrix(0,5) = 0.18301270189221927; rExtrapolationMatrix(0,6) = -0.04903810567665795; rExtrapolationMatrix(0,7) = 0.18301270189221927;

        rExtrapolationMatrix(1,0) = -0.6830127018922192; rExtrapolationMatrix(1,1) = 2.549038105676658; rExtrapolationMatrix(1,2) = -0.6830127018922192; rExtrapolationMatrix(1,3) = 0.18301270189221927;
        rExtrapolationMatrix(1,4) = 0.18301270189221927; rExtrapolationMatrix(1,5) = -0.6830127018922192; rExtrapolationMatrix(1,6) = 0.18301270189221927; rExtrapolationMatrix(1,7) = -0.04903810567665795;

        rExtrapolationMatrix(2,0) = 0.18301270189221927; rExtrapolationMatrix(2,1) = -0.6830127018922192; rExtrapolationMatrix(2,2) = 2.549038105676658; rExtrapolationMatrix(2,3) = -0.6830127018922192;
        rExtrapolationMatrix(2,4) = -0.04903810567665795; rExtrapolationMatrix(2,5) = 0.18301270189221927; rExtrapolationMatrix(2,6) = -0.6830127018922192; rExtrapolationMatrix(2,7) = 0.18301270189221927;

        rExtrapolationMatrix(3,0) = -0.6830127018922192; rExtrapolationMatrix(3,1) = 0.18301270189221927; rExtrapolationMatrix(3,2) = -0.6830127018922192; rExtrapolationMatrix(3,3) = 2.549038105676658;
        rExtrapolationMatrix(3,4) = 0.18301270189221927; rExtrapolationMatrix(3,5) = -0.04903810567665795; rExtrapolationMatrix(3,6) = 0.18301270189221927; rExtrapolationMatrix(3,7) = -0.6830127018922192;

        rExtrapolationMatrix(4,0) = -0.6830127018922192; rExtrapolationMatrix(4,1) = 0.18301270189221927; rExtrapolationMatrix(4,2) = -0.04903810567665795; rExtrapolationMatrix(4,3) = 0.18301270189221927;
        rExtrapolationMatrix(4,4) = 2.549038105676658; rExtrapolationMatrix(4,5) = -0.6830127018922192; rExtrapolationMatrix(4,6) = 0.18301270189221927; rExtrapolationMatrix(4,7) = -0.6830127018922192;

        rExtrapolationMatrix(5,0) = 0.18301270189221927; rExtrapolationMatrix(5,1) = -0.6830127018922192; rExtrapolationMatrix(5,2) = 0.18301270189221927; rExtrapolationMatrix(5,3) = -0.04903810567665795;
        rExtrapolationMatrix(5,4) = -0.6830127018922192; rExtrapolationMatrix(5,5) = 2.549038105676658; rExtrapolationMatrix(5,6) = -0.6830127018922192; rExtrapolationMatrix(5,7) = 0.18301270189221927;

        rExtrapolationMatrix(6,0) = -0.04903810567665795; rExtrapolationMatrix(6,1) = 0.18301270189221927; rExtrapolationMatrix(6,2) = -0.6830127018922192; rExtrapolationMatrix(6,3) = 0.18301270189221927;
        rExtrapolationMatrix(6,4) = 0.18301270189221927; rExtrapolationMatrix(6,5) = -0.6830127018922192; rExtrapolationMatrix(6,6) = 2.549038105676658; rExtrapolationMatrix(6,7) = -0.6830127018922192;

        rExtrapolationMatrix(7,0) = 0.18301270189221927; rExtrapolationMatrix(7,1) = -0.04903810567665795; rExtrapolationMatrix(7,2) = 0.18301270189221927; rExtrapolationMatrix(7,3) = -0.6830127018922192;
        rExtrapolationMatrix(7,4) = -0.6830127018922192; rExtrapolationMatrix(7,5) = 0.18301270189221927; rExtrapolationMatrix(7,6) = -0.6830127018922192; rExtrapolationMatrix(7,7) = 2.549038105676658;
}

}; /* Class ElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_ELEMENT_UTILITIES defined */
