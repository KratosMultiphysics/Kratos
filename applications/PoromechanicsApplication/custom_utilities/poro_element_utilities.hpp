//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_PORO_ELEMENT_UTILITIES )
#define  KRATOS_PORO_ELEMENT_UTILITIES

// System includes
//#include <cmath>

// Project includes
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

class PoroElementUtilities
{

typedef std::size_t IndexType;

public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,2,6>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Triangle_2d_3
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,2) = Ncontainer(GPoint,1); rNu(0,4) = Ncontainer(GPoint,2);
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,3) = Ncontainer(GPoint,1); rNu(1,5) = Ncontainer(GPoint,2);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,2,8>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Quadrilateral_2d_4
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,2) = Ncontainer(GPoint,1); rNu(0,4) = Ncontainer(GPoint,2); rNu(0,6) = Ncontainer(GPoint,3);
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,3) = Ncontainer(GPoint,1); rNu(1,5) = Ncontainer(GPoint,2); rNu(1,7) = Ncontainer(GPoint,3);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,3,12>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Tetrahedra_3d_4
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,3) = Ncontainer(GPoint,1); rNu(0,6) = Ncontainer(GPoint,2); rNu(0,9)  = Ncontainer(GPoint,3);
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,4) = Ncontainer(GPoint,1); rNu(1,7) = Ncontainer(GPoint,2); rNu(1,10) = Ncontainer(GPoint,3);
        rNu(2,2) = Ncontainer(GPoint,0); rNu(2,5) = Ncontainer(GPoint,1); rNu(2,8) = Ncontainer(GPoint,2); rNu(2,11) = Ncontainer(GPoint,3);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,3,18>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Prism_3d_6
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,3) = Ncontainer(GPoint,1); rNu(0,6) = Ncontainer(GPoint,2);
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,4) = Ncontainer(GPoint,1); rNu(1,7) = Ncontainer(GPoint,2);
        rNu(2,2) = Ncontainer(GPoint,0); rNu(2,5) = Ncontainer(GPoint,1); rNu(2,8) = Ncontainer(GPoint,2);

        rNu(0,9) = Ncontainer(GPoint,3); rNu(0,12) = Ncontainer(GPoint,4); rNu(0,15) = Ncontainer(GPoint,5);
        rNu(1,10) = Ncontainer(GPoint,3); rNu(1,13) = Ncontainer(GPoint,4); rNu(1,16) = Ncontainer(GPoint,5);
        rNu(2,11) = Ncontainer(GPoint,3); rNu(2,14) = Ncontainer(GPoint,4); rNu(2,17) = Ncontainer(GPoint,5);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,3,24>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Hexahedron_3d_8
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,3) = Ncontainer(GPoint,1); rNu(0,6) = Ncontainer(GPoint,2); rNu(0,9)  = Ncontainer(GPoint,3);
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,4) = Ncontainer(GPoint,1); rNu(1,7) = Ncontainer(GPoint,2); rNu(1,10) = Ncontainer(GPoint,3);
        rNu(2,2) = Ncontainer(GPoint,0); rNu(2,5) = Ncontainer(GPoint,1); rNu(2,8) = Ncontainer(GPoint,2); rNu(2,11) = Ncontainer(GPoint,3);

        rNu(0,12) = Ncontainer(GPoint,4); rNu(0,15) = Ncontainer(GPoint,5); rNu(0,18) = Ncontainer(GPoint,6); rNu(0,21) = Ncontainer(GPoint,7);
        rNu(1,13) = Ncontainer(GPoint,4); rNu(1,16) = Ncontainer(GPoint,5); rNu(1,19) = Ncontainer(GPoint,6); rNu(1,22) = Ncontainer(GPoint,7);
        rNu(2,14) = Ncontainer(GPoint,4); rNu(2,17) = Ncontainer(GPoint,5); rNu(2,20) = Ncontainer(GPoint,6); rNu(2,23) = Ncontainer(GPoint,7);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculateNuElementMatrix(BoundedMatrix<double,3,9>& rNut, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Triangle_2d_3
        rNut(0,0) = Ncontainer(GPoint,0); rNut(0,3) = Ncontainer(GPoint,1); rNut(0,6) = Ncontainer(GPoint,2);
        rNut(1,1) = Ncontainer(GPoint,0); rNut(1,4) = Ncontainer(GPoint,1); rNut(1,7) = Ncontainer(GPoint,2);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuElementMatrix(BoundedMatrix<double,3,12>& rNut, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Quadrilateral_2d_4
        rNut(0,0) = Ncontainer(GPoint,0); rNut(0,3) = Ncontainer(GPoint,1); rNut(0,6) = Ncontainer(GPoint,2); rNut(0,9) = Ncontainer(GPoint,3);
        rNut(1,1) = Ncontainer(GPoint,0); rNut(1,4) = Ncontainer(GPoint,1); rNut(1,7) = Ncontainer(GPoint,2); rNut(1,10) = Ncontainer(GPoint,3);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuElementMatrix(BoundedMatrix<double,4,16>& rNut, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Tetrahedra_3d_4
        rNut(0,0) = Ncontainer(GPoint,0); rNut(0,4) = Ncontainer(GPoint,1); rNut(0,8)  = Ncontainer(GPoint,2); rNut(0,12) = Ncontainer(GPoint,3);
        rNut(1,1) = Ncontainer(GPoint,0); rNut(1,5) = Ncontainer(GPoint,1); rNut(1,9)  = Ncontainer(GPoint,2); rNut(1,13) = Ncontainer(GPoint,3);
        rNut(2,2) = Ncontainer(GPoint,0); rNut(2,6) = Ncontainer(GPoint,1); rNut(2,10) = Ncontainer(GPoint,2); rNut(2,14) = Ncontainer(GPoint,3);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuElementMatrix(BoundedMatrix<double,4,24>& rNut, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Prism_3d_6
        rNut(0,0) = Ncontainer(GPoint,0); rNut(0,4) = Ncontainer(GPoint,1); rNut(0,8) = Ncontainer(GPoint,2);
        rNut(1,1) = Ncontainer(GPoint,0); rNut(1,5) = Ncontainer(GPoint,1); rNut(1,9) = Ncontainer(GPoint,2);
        rNut(2,2) = Ncontainer(GPoint,0); rNut(2,6) = Ncontainer(GPoint,1); rNut(2,10) = Ncontainer(GPoint,2);

        rNut(0,12) = Ncontainer(GPoint,3); rNut(0,16) = Ncontainer(GPoint,4); rNut(0,20) = Ncontainer(GPoint,5);
        rNut(1,13) = Ncontainer(GPoint,3); rNut(1,17) = Ncontainer(GPoint,4); rNut(1,21) = Ncontainer(GPoint,5);
        rNut(2,14) = Ncontainer(GPoint,3); rNut(2,18) = Ncontainer(GPoint,4); rNut(2,22) = Ncontainer(GPoint,5);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuElementMatrix(BoundedMatrix<double,4,32>& rNut, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Hexahedron_3d_8
        rNut(0,0) = Ncontainer(GPoint,0); rNut(0,4) = Ncontainer(GPoint,1); rNut(0,8) = Ncontainer(GPoint,2); rNut(0,12) = Ncontainer(GPoint,3);
        rNut(1,1) = Ncontainer(GPoint,0); rNut(1,5) = Ncontainer(GPoint,1); rNut(1,9) = Ncontainer(GPoint,2); rNut(1,13) = Ncontainer(GPoint,3);
        rNut(2,2) = Ncontainer(GPoint,0); rNut(2,6) = Ncontainer(GPoint,1); rNut(2,10) = Ncontainer(GPoint,2); rNut(2,14) = Ncontainer(GPoint,3);

        rNut(0,16) = Ncontainer(GPoint,4); rNut(0,20) = Ncontainer(GPoint,5); rNut(0,24) = Ncontainer(GPoint,6); rNut(0,28) = Ncontainer(GPoint,7);
        rNut(1,17) = Ncontainer(GPoint,4); rNut(1,21) = Ncontainer(GPoint,5); rNut(1,25) = Ncontainer(GPoint,6); rNut(1,29) = Ncontainer(GPoint,7);
        rNut(2,18) = Ncontainer(GPoint,4); rNut(2,22) = Ncontainer(GPoint,5); rNut(2,26) = Ncontainer(GPoint,6); rNut(2,30) = Ncontainer(GPoint,7);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,2>& rVector,const Matrix& Ncontainer,
                                                        const array_1d<double,6>& VariableWithComponents,const unsigned int& GPoint)
    {
        //Triangle_2d_3
        noalias(rVector) = ZeroVector(2);

        unsigned int index = 0;
        for(unsigned int i=0; i<3; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,2>& rVector,const Matrix& Ncontainer,
                                                        const array_1d<double,8>& VariableWithComponents,const unsigned int& GPoint)
    {
        //Quadrilateral_2d_4
        noalias(rVector) = ZeroVector(2);

        unsigned int index = 0;
        for(unsigned int i=0; i<4; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer,
                                                        const array_1d<double,12>& VariableWithComponents,const unsigned int& GPoint)
    {
        //Tetrahedra_3d_4
        noalias(rVector) = ZeroVector(3);

        unsigned int index = 0;
        for(unsigned int i=0; i<4; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer,
                                                        const array_1d<double,18>& VariableWithComponents,const unsigned int& GPoint)
    {
        //Prism_3d_6
        noalias(rVector) = ZeroVector(3);

        unsigned int index = 0;
        for(unsigned int i=0; i<6; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer,
                                                        const array_1d<double,24>& VariableWithComponents,const unsigned int& GPoint)
    {
        //Hexahedra_3d_8
        noalias(rVector) = ZeroVector(3);

        unsigned int index = 0;
        for(unsigned int i=0; i<8; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
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

    static inline void GetNodalVariableVector(array_1d<double,6>& rNodalVariableVector, const Element::GeometryType& Geom,
                                            const Variable<array_1d<double,3>>& Variable, IndexType SolutionStepIndex = 0)
    {
        //Triangle_2d_3
        array_1d<double,3> NodalVariableAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<3; i++)
        {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable,SolutionStepIndex);
            rNodalVariableVector[index++] = NodalVariableAux[0];
            rNodalVariableVector[index++] = NodalVariableAux[1];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetNodalVariableVector(array_1d<double,8>& rNodalVariableVector, const Element::GeometryType& Geom,
                                            const Variable<array_1d<double,3>>& Variable, IndexType SolutionStepIndex = 0)
    {
        //Quadrilateral_2d_4
        array_1d<double,3> NodalVariableAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<4; i++)
        {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable,SolutionStepIndex);
            rNodalVariableVector[index++] = NodalVariableAux[0];
            rNodalVariableVector[index++] = NodalVariableAux[1];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetNodalVariableVector(array_1d<double,12>& rNodalVariableVector, const Element::GeometryType& Geom,
                                            const Variable<array_1d<double,3>>& Variable, IndexType SolutionStepIndex = 0)
    {
        //Tetrahedra_3d_4
        array_1d<double,3> NodalVariableAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<4; i++)
        {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable,SolutionStepIndex);
            rNodalVariableVector[index++] = NodalVariableAux[0];
            rNodalVariableVector[index++] = NodalVariableAux[1];
            rNodalVariableVector[index++] = NodalVariableAux[2];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetNodalVariableVector(array_1d<double,18>& rNodalVariableVector, const Element::GeometryType& Geom,
                                            const Variable<array_1d<double,3>>& Variable, IndexType SolutionStepIndex = 0)
    {
        //Prism_3d_6
        array_1d<double,3> NodalVariableAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<6; i++)
        {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable,SolutionStepIndex);
            rNodalVariableVector[index++] = NodalVariableAux[0];
            rNodalVariableVector[index++] = NodalVariableAux[1];
            rNodalVariableVector[index++] = NodalVariableAux[2];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetNodalVariableVector(array_1d<double,24>& rNodalVariableVector, const Element::GeometryType& Geom,
                                            const Variable<array_1d<double,3>>& Variable, IndexType SolutionStepIndex = 0)
    {
        //Hexahedron_3d_8
        array_1d<double,3> NodalVariableAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<8; i++)
        {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable,SolutionStepIndex);
            rNodalVariableVector[index++] = NodalVariableAux[0];
            rNodalVariableVector[index++] = NodalVariableAux[1];
            rNodalVariableVector[index++] = NodalVariableAux[2];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculatePermeabilityMatrix(BoundedMatrix<double,2,2>& rPermeabilityMatrix,
                                                    const Element::PropertiesType& Prop)
    {
        //2D
        rPermeabilityMatrix(0,0) = Prop[PERMEABILITY_XX];
        rPermeabilityMatrix(1,1) = Prop[PERMEABILITY_YY];

        rPermeabilityMatrix(0,1) = Prop[PERMEABILITY_XY];
        rPermeabilityMatrix(1,0) = rPermeabilityMatrix(0,1);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculatePermeabilityMatrix(BoundedMatrix<double,3,3>& rPermeabilityMatrix,
                                                    const Element::PropertiesType& Prop)
    {
        //3D
        rPermeabilityMatrix(0,0) = Prop[PERMEABILITY_XX];
        rPermeabilityMatrix(1,1) = Prop[PERMEABILITY_YY];
        rPermeabilityMatrix(2,2) = Prop[PERMEABILITY_ZZ];

        rPermeabilityMatrix(0,1) = Prop[PERMEABILITY_XY];
        rPermeabilityMatrix(1,0) = rPermeabilityMatrix(0,1);

        rPermeabilityMatrix(1,2) = Prop[PERMEABILITY_YZ];
        rPermeabilityMatrix(2,1) = rPermeabilityMatrix(1,2);

        rPermeabilityMatrix(2,0) = Prop[PERMEABILITY_ZX];
        rPermeabilityMatrix(0,2) = rPermeabilityMatrix(2,0);
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void InvertMatrix2(BoundedMatrix<double,2,2>& rInvertedMatrix,
                                    const BoundedMatrix<double,2,2>& InputMatrix)
    {
        double InputMatrixDet = InputMatrix(0,0)*InputMatrix(1,1)-InputMatrix(0,1)*InputMatrix(1,0);

        rInvertedMatrix(0,0) =  InputMatrix(1,1)/InputMatrixDet;
        rInvertedMatrix(0,1) = -InputMatrix(0,1)/InputMatrixDet;
        rInvertedMatrix(1,0) = -InputMatrix(1,0)/InputMatrixDet;
        rInvertedMatrix(1,1) =  InputMatrix(0,0)/InputMatrixDet;
    }

    //----------------------------------------------------------------------------------------
/*
	template<class T>
	bool InvertMatrix(const T& input, T& inverse)
	{
		typedef permutation_matrix<std::size_t> pmatrix;

		// create a working copy of the input
		T A(input);

		// create a permutation matrix for the LU-factorization
		pmatrix pm(A.size1());

		// perform LU-factorization
		int res = lu_factorize(A, pm);
		if (res != 0)
			return false;

		// create identity matrix of "inverse"
		inverse.assign(identity_matrix<double> (A.size1()));

		// backsubstitute to get the inverse
		lu_substitute(A, pm, inverse);

		return true;
	}
*/

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void AssembleUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,6,6>& UBlockMatrix)
    {
        //Triangle_2d_3
        unsigned int Global_i, Global_j, Local_i, Local_j;

        for(unsigned int i = 0; i < 3; i++)
        {
            Global_i = i * (2 + 1);
            Local_i = i * 2;

            for(unsigned int j = 0; j < 3; j++)
            {
                Global_j = j * (2 + 1);
                Local_j = j * 2;

                rLeftHandSideMatrix(Global_i,Global_j)     += UBlockMatrix(Local_i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1)   += UBlockMatrix(Local_i,Local_j+1);
                rLeftHandSideMatrix(Global_i+1,Global_j)   += UBlockMatrix(Local_i+1,Local_j);
                rLeftHandSideMatrix(Global_i+1,Global_j+1) += UBlockMatrix(Local_i+1,Local_j+1);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,8,8>& UBlockMatrix)
    {
        //Quadrilateral_2d_4
        unsigned int Global_i, Global_j, Local_i, Local_j;

        for(unsigned int i = 0; i < 4; i++)
        {
            Global_i = i * (2 + 1);
            Local_i = i * 2;

            for(unsigned int j = 0; j < 4; j++)
            {
                Global_j = j * (2 + 1);
                Local_j = j * 2;

                rLeftHandSideMatrix(Global_i,Global_j)     += UBlockMatrix(Local_i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1)   += UBlockMatrix(Local_i,Local_j+1);
                rLeftHandSideMatrix(Global_i+1,Global_j)   += UBlockMatrix(Local_i+1,Local_j);
                rLeftHandSideMatrix(Global_i+1,Global_j+1) += UBlockMatrix(Local_i+1,Local_j+1);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,12,12>& UBlockMatrix)
    {
        //Tetrahedra_3d_4
        unsigned int Global_i, Global_j, Local_i, Local_j;

        for(unsigned int i = 0; i < 4; i++)
        {
            Global_i = i * (3 + 1);
            Local_i = i * 3;

            for(unsigned int j = 0; j < 4; j++)
            {
                Global_j = j * (3 + 1);
                Local_j = j * 3;

                rLeftHandSideMatrix(Global_i,Global_j)     += UBlockMatrix(Local_i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1)   += UBlockMatrix(Local_i,Local_j+1);
                rLeftHandSideMatrix(Global_i+1,Global_j)   += UBlockMatrix(Local_i+1,Local_j);
                rLeftHandSideMatrix(Global_i+1,Global_j+1) += UBlockMatrix(Local_i+1,Local_j+1);

                rLeftHandSideMatrix(Global_i,Global_j+2)   += UBlockMatrix(Local_i,Local_j+2);
                rLeftHandSideMatrix(Global_i+1,Global_j+2) += UBlockMatrix(Local_i+1,Local_j+2);
                rLeftHandSideMatrix(Global_i+2,Global_j+1) += UBlockMatrix(Local_i+2,Local_j+1);
                rLeftHandSideMatrix(Global_i+2,Global_j)   += UBlockMatrix(Local_i+2,Local_j);
                rLeftHandSideMatrix(Global_i+2,Global_j+2) += UBlockMatrix(Local_i+2,Local_j+2);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,18,18>& UBlockMatrix)
    {
        //Prism_3d_6
        unsigned int Global_i, Global_j, Local_i, Local_j;

        for(unsigned int i = 0; i < 6; i++)
        {
            Global_i = i * (3 + 1);
            Local_i = i * 3;

            for(unsigned int j = 0; j < 6; j++)
            {
                Global_j = j * (3 + 1);
                Local_j = j * 3;

                rLeftHandSideMatrix(Global_i,Global_j)     += UBlockMatrix(Local_i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1)   += UBlockMatrix(Local_i,Local_j+1);
                rLeftHandSideMatrix(Global_i+1,Global_j)   += UBlockMatrix(Local_i+1,Local_j);
                rLeftHandSideMatrix(Global_i+1,Global_j+1) += UBlockMatrix(Local_i+1,Local_j+1);

                rLeftHandSideMatrix(Global_i,Global_j+2)   += UBlockMatrix(Local_i,Local_j+2);
                rLeftHandSideMatrix(Global_i+1,Global_j+2) += UBlockMatrix(Local_i+1,Local_j+2);
                rLeftHandSideMatrix(Global_i+2,Global_j+1) += UBlockMatrix(Local_i+2,Local_j+1);
                rLeftHandSideMatrix(Global_i+2,Global_j)   += UBlockMatrix(Local_i+2,Local_j);
                rLeftHandSideMatrix(Global_i+2,Global_j+2) += UBlockMatrix(Local_i+2,Local_j+2);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,24,24>& UBlockMatrix)
    {
        //Hexahedra_3d_8
        unsigned int Global_i, Global_j, Local_i, Local_j;

        for(unsigned int i = 0; i < 8; i++)
        {
            Global_i = i * (3 + 1);
            Local_i = i * 3;

            for(unsigned int j = 0; j < 8; j++)
            {
                Global_j = j * (3 + 1);
                Local_j = j * 3;

                rLeftHandSideMatrix(Global_i,Global_j)     += UBlockMatrix(Local_i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1)   += UBlockMatrix(Local_i,Local_j+1);
                rLeftHandSideMatrix(Global_i+1,Global_j)   += UBlockMatrix(Local_i+1,Local_j);
                rLeftHandSideMatrix(Global_i+1,Global_j+1) += UBlockMatrix(Local_i+1,Local_j+1);

                rLeftHandSideMatrix(Global_i,Global_j+2)   += UBlockMatrix(Local_i,Local_j+2);
                rLeftHandSideMatrix(Global_i+1,Global_j+2) += UBlockMatrix(Local_i+1,Local_j+2);
                rLeftHandSideMatrix(Global_i+2,Global_j+1) += UBlockMatrix(Local_i+2,Local_j+1);
                rLeftHandSideMatrix(Global_i+2,Global_j)   += UBlockMatrix(Local_i+2,Local_j);
                rLeftHandSideMatrix(Global_i+2,Global_j+2) += UBlockMatrix(Local_i+2,Local_j+2);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUPBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,6,3>& UPBlockMatrix)
    {
        //Triangle_2d_3
        unsigned int Global_i, Global_j, Local_i;

        for(unsigned int i = 0; i < 3; i++)
        {
            Global_i = i * (2 + 1);
            Local_i = i * 2;

            for(unsigned int j = 0; j < 3; j++)
            {
                Global_j = j * (2 + 1) + 2;

                rLeftHandSideMatrix(Global_i,Global_j)   += UPBlockMatrix(Local_i,j);
                rLeftHandSideMatrix(Global_i+1,Global_j) += UPBlockMatrix(Local_i+1,j);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUPBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,8,4>& UPBlockMatrix)
    {
        //Quadrilateral_2d_4
        unsigned int Global_i, Global_j, Local_i;

        for(unsigned int i = 0; i < 4; i++)
        {
            Global_i = i * (2 + 1);
            Local_i = i * 2;

            for(unsigned int j = 0; j < 4; j++)
            {
                Global_j = j * (2 + 1) + 2;

                rLeftHandSideMatrix(Global_i,Global_j)   += UPBlockMatrix(Local_i,j);
                rLeftHandSideMatrix(Global_i+1,Global_j) += UPBlockMatrix(Local_i+1,j);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUPBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,12,4>& UPBlockMatrix)
    {
        //Tetrahedra_3d_4
        unsigned int Global_i, Global_j, Local_i;

        for(unsigned int i = 0; i < 4; i++)
        {
            Global_i = i * (3 + 1);
            Local_i = i * 3;

            for(unsigned int j = 0; j < 4; j++)
            {
                Global_j = j * (3 + 1) + 3;

                rLeftHandSideMatrix(Global_i,Global_j)   += UPBlockMatrix(Local_i,j);
                rLeftHandSideMatrix(Global_i+1,Global_j) += UPBlockMatrix(Local_i+1,j);
                rLeftHandSideMatrix(Global_i+2,Global_j) += UPBlockMatrix(Local_i+2,j);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUPBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,18,6>& UPBlockMatrix)
    {
        //Prism_3d_6
        unsigned int Global_i, Global_j, Local_i;

        for(unsigned int i = 0; i < 6; i++)
        {
            Global_i = i * (3 + 1);
            Local_i = i * 3;

            for(unsigned int j = 0; j < 6; j++)
            {
                Global_j = j * (3 + 1) + 3;

                rLeftHandSideMatrix(Global_i,Global_j)   += UPBlockMatrix(Local_i,j);
                rLeftHandSideMatrix(Global_i+1,Global_j) += UPBlockMatrix(Local_i+1,j);
                rLeftHandSideMatrix(Global_i+2,Global_j) += UPBlockMatrix(Local_i+2,j);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUPBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,24,8>& UPBlockMatrix)
    {
        //Hexahedra_3d_8
        unsigned int Global_i, Global_j, Local_i;

        for(unsigned int i = 0; i < 8; i++)
        {
            Global_i = i * (3 + 1);
            Local_i = i * 3;

            for(unsigned int j = 0; j < 8; j++)
            {
                Global_j = j * (3 + 1) + 3;

                rLeftHandSideMatrix(Global_i,Global_j)   += UPBlockMatrix(Local_i,j);
                rLeftHandSideMatrix(Global_i+1,Global_j) += UPBlockMatrix(Local_i+1,j);
                rLeftHandSideMatrix(Global_i+2,Global_j) += UPBlockMatrix(Local_i+2,j);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssemblePUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,3,6>& PUBlockMatrix)
    {
        //Triangle_2d_3
        unsigned int Global_i, Global_j, Local_j;

        for(unsigned int i = 0; i < 3; i++)
        {
            Global_i = i * (2 + 1) + 2;

            for(unsigned int j = 0; j < 3; j++)
            {
                Global_j = j * (2 + 1);
                Local_j = j * 2;

                rLeftHandSideMatrix(Global_i,Global_j)   += PUBlockMatrix(i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1) += PUBlockMatrix(i,Local_j+1);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssemblePUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,4,8>& PUBlockMatrix)
    {
        //Quadrilateral_2d_4
        unsigned int Global_i, Global_j, Local_j;

        for(unsigned int i = 0; i < 4; i++)
        {
            Global_i = i * (2 + 1) + 2;

            for(unsigned int j = 0; j < 4; j++)
            {
                Global_j = j * (2 + 1);
                Local_j = j * 2;

                rLeftHandSideMatrix(Global_i,Global_j)   += PUBlockMatrix(i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1) += PUBlockMatrix(i,Local_j+1);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssemblePUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,4,12>& PUBlockMatrix)
    {
        //Tetrahedra_3d_4
        unsigned int Global_i, Global_j, Local_j;

        for(unsigned int i = 0; i < 4; i++)
        {
            Global_i = i * (3 + 1) + 3;

            for(unsigned int j = 0; j < 4; j++)
            {
                Global_j = j * (3 + 1);
                Local_j = j * 3;

                rLeftHandSideMatrix(Global_i,Global_j)   += PUBlockMatrix(i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1) += PUBlockMatrix(i,Local_j+1);
                rLeftHandSideMatrix(Global_i,Global_j+2) += PUBlockMatrix(i,Local_j+2);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssemblePUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,6,18>& PUBlockMatrix)
    {
        //Prism_3d_6
        unsigned int Global_i, Global_j, Local_j;

        for(unsigned int i = 0; i < 6; i++)
        {
            Global_i = i * (3 + 1) + 3;

            for(unsigned int j = 0; j < 6; j++)
            {
                Global_j = j * (3 + 1);
                Local_j = j * 3;

                rLeftHandSideMatrix(Global_i,Global_j)   += PUBlockMatrix(i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1) += PUBlockMatrix(i,Local_j+1);
                rLeftHandSideMatrix(Global_i,Global_j+2) += PUBlockMatrix(i,Local_j+2);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssemblePUBlockMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,8,24>& PUBlockMatrix)
    {
        //Hexahedra_3d_8
        unsigned int Global_i, Global_j, Local_j;

        for(unsigned int i = 0; i < 8; i++)
        {
            Global_i = i * (3 + 1) + 3;

            for(unsigned int j = 0; j < 8; j++)
            {
                Global_j = j * (3 + 1);
                Local_j = j * 3;

                rLeftHandSideMatrix(Global_i,Global_j)   += PUBlockMatrix(i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1) += PUBlockMatrix(i,Local_j+1);
                rLeftHandSideMatrix(Global_i,Global_j+2) += PUBlockMatrix(i,Local_j+2);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    template< class TMatrixType >
    static inline void AssemblePBlockMatrix(Matrix& rLeftHandSideMatrix,const TMatrixType& PBlockMatrix, const unsigned int& Dim, const unsigned int& NumNodes)
    {
        unsigned int Global_i, Global_j;

        for(unsigned int i = 0; i < NumNodes; i++)
        {
            Global_i = i * (Dim + 1) + Dim;

            for(unsigned int j = 0; j < NumNodes; j++)
            {
                Global_j = j * (Dim + 1) + Dim;

                rLeftHandSideMatrix(Global_i,Global_j) += PBlockMatrix(i,j);
            }
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockVector(Vector& rRightHandSideVector, const array_1d<double,6>& UBlockVector)
    {
        //Triangle_2d_3
        unsigned int Global_i, Local_i;

        for(unsigned int i = 0; i < 3; i++)
        {
            Global_i = i * (2 + 1);
            Local_i  = i * 2;

            rRightHandSideVector[Global_i]   += UBlockVector[Local_i];
            rRightHandSideVector[Global_i+1] += UBlockVector[Local_i+1];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockVector(Vector& rRightHandSideVector, const array_1d<double,8>& UBlockVector)
    {
        //Quadrilateral_2d_4
        unsigned int Global_i, Local_i;

        for(unsigned int i = 0; i < 4; i++)
        {
            Global_i = i * (2 + 1);
            Local_i  = i * 2;

            rRightHandSideVector[Global_i]   += UBlockVector[Local_i];
            rRightHandSideVector[Global_i+1] += UBlockVector[Local_i+1];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockVector(Vector& rRightHandSideVector, const array_1d<double,12>& UBlockVector)
    {
        //Tetrahedra_3d_4
        unsigned int Global_i, Local_i;

        for(unsigned int i = 0; i < 4; i++)
        {
            Global_i = i * (3 + 1);
            Local_i  = i * 3;

            rRightHandSideVector[Global_i]   += UBlockVector[Local_i];
            rRightHandSideVector[Global_i+1] += UBlockVector[Local_i+1];
            rRightHandSideVector[Global_i+2] += UBlockVector[Local_i+2];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockVector(Vector& rRightHandSideVector, const array_1d<double,18>& UBlockVector)
    {
        //Prism_3d_6
        unsigned int Global_i, Local_i;

        for(unsigned int i = 0; i < 6; i++)
        {
            Global_i = i * (3 + 1);
            Local_i  = i * 3;

            rRightHandSideVector[Global_i]   += UBlockVector[Local_i];
            rRightHandSideVector[Global_i+1] += UBlockVector[Local_i+1];
            rRightHandSideVector[Global_i+2] += UBlockVector[Local_i+2];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockVector(Vector& rRightHandSideVector, const array_1d<double,24>& UBlockVector)
    {
        //Hexahedra_3d_8
        unsigned int Global_i, Local_i;

        for(unsigned int i = 0; i < 8; i++)
        {
            Global_i = i * (3 + 1);
            Local_i  = i * 3;

            rRightHandSideVector[Global_i]   += UBlockVector[Local_i];
            rRightHandSideVector[Global_i+1] += UBlockVector[Local_i+1];
            rRightHandSideVector[Global_i+2] += UBlockVector[Local_i+2];
        }
    }

    //----------------------------------------------------------------------------------------

    template< class TVectorType >
    static inline void AssemblePBlockVector(Vector& rRightHandSideVector,const TVectorType& PBlockVector, const unsigned int& Dim, const unsigned int& NumNodes)
    {
        unsigned int Global_i;

        for(unsigned int i = 0; i < NumNodes; i++)
        {
            Global_i = i * (Dim + 1) + Dim;

            rRightHandSideVector[Global_i] += PBlockVector[i];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// CalculateExtrapolationMatrix
    /// The matrix contains the shape functions at each GP evaluated at each node.
    /// Rows: nodes
    /// Columns: GP

    static inline void CalculateExtrapolationMatrix(BoundedMatrix<double,4,4>& rExtrapolationMatrix)
    {
        //Quadrilateral_2d_4
        //GI_GAUSS_2

        rExtrapolationMatrix(0,0) = 1.8660254037844386; rExtrapolationMatrix(0,1) = -0.5; rExtrapolationMatrix(0,2) = 0.13397459621556132; rExtrapolationMatrix(0,3) = -0.5;
        rExtrapolationMatrix(1,0) = -0.5; rExtrapolationMatrix(1,1) = 1.8660254037844386; rExtrapolationMatrix(1,2) = -0.5; rExtrapolationMatrix(1,3) = 0.13397459621556132;
        rExtrapolationMatrix(2,0) = 0.13397459621556132; rExtrapolationMatrix(2,1) = -0.5; rExtrapolationMatrix(2,2) = 1.8660254037844386; rExtrapolationMatrix(2,3) = -0.5;
        rExtrapolationMatrix(3,0) = -0.5; rExtrapolationMatrix(3,1) = 0.13397459621556132; rExtrapolationMatrix(3,2) = -0.5; rExtrapolationMatrix(3,3) = 1.8660254037844386;
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateExtrapolationMatrix(BoundedMatrix<double,8,8>& rExtrapolationMatrix)
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

}; /* Class PoroElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_PORO_ELEMENT_UTILITIES defined */
