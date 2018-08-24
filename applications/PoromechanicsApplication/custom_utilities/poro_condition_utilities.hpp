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


#if !defined(KRATOS_PORO_CONDITION_UTILITIES )
#define  KRATOS_PORO_CONDITION_UTILITIES

// Project includes
#include "includes/element.h"

namespace Kratos
{

class PoroConditionUtilities
{
    
public:
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,2,4>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Line_2d_2
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,2) = Ncontainer(GPoint,1); 
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,3) = Ncontainer(GPoint,1);
    }
    
    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,3,9>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Triangle_3d_3
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,3) = Ncontainer(GPoint,1); rNu(0,6) = Ncontainer(GPoint,2);
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,4) = Ncontainer(GPoint,1); rNu(1,7) = Ncontainer(GPoint,2);
        rNu(2,2) = Ncontainer(GPoint,0); rNu(2,5) = Ncontainer(GPoint,1); rNu(2,8) = Ncontainer(GPoint,2);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,3,12>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Quadrilateral_3d_4
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,3) = Ncontainer(GPoint,1); rNu(0,6) = Ncontainer(GPoint,2); rNu(0,9) = Ncontainer(GPoint,3);
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,4) = Ncontainer(GPoint,1); rNu(1,7) = Ncontainer(GPoint,2); rNu(1,10) = Ncontainer(GPoint,3);
        rNu(2,2) = Ncontainer(GPoint,0); rNu(2,5) = Ncontainer(GPoint,1); rNu(2,8) = Ncontainer(GPoint,2); rNu(2,11) = Ncontainer(GPoint,3);
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,2>& rVector,const Matrix& Ncontainer, 
                                                        const array_1d<double,4>& VariableWithComponents,const unsigned int& GPoint)
    {        
        //Line_2d_2
        noalias(rVector) = ZeroVector(2);
        
        unsigned int index = 0;
        for(unsigned int i=0; i<2; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }
        
    //----------------------------------------------------------------------------------------
    
    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer, 
                                                        const array_1d<double,9>& VariableWithComponents,const unsigned int& GPoint)
    {
        //Triangle_3d_3
        noalias(rVector) = ZeroVector(3);
        
        unsigned int index = 0;
        for(unsigned int i=0; i<3; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }

    //----------------------------------------------------------------------------------------
    
    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer, 
                                                        const array_1d<double,12>& VariableWithComponents,const unsigned int& GPoint)
    {
        //Quadrilateral_3d_4
        noalias(rVector) = ZeroVector(3);
        
        unsigned int index = 0;
        for(unsigned int i=0; i<4; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[index++];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void GetNodalVariableVector(array_1d<double,4>& rNodalVariableVector, const Element::GeometryType& Geom, 
                                            const Variable<array_1d<double,3>>& Variable)
    {
        //Line_2d_2
        array_1d<double,3> NodalVariableAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<2; i++)
        {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable);
            rNodalVariableVector[index++] = NodalVariableAux[0];
            rNodalVariableVector[index++] = NodalVariableAux[1];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetNodalVariableVector(array_1d<double,9>& rNodalVariableVector, const Element::GeometryType& Geom,
                                            const Variable<array_1d<double,3>>& Variable)
    {
        //Triangle_3d_3
        array_1d<double,3> NodalVariableAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<3; i++)
        {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable);
            rNodalVariableVector[index++] = NodalVariableAux[0];
            rNodalVariableVector[index++] = NodalVariableAux[1];
            rNodalVariableVector[index++] = NodalVariableAux[2];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetNodalVariableVector(array_1d<double,12>& rNodalVariableVector, const Element::GeometryType& Geom,
                                            const Variable<array_1d<double,3>>& Variable)
    {
        //Quadrilateral_3d_4
        array_1d<double,3> NodalVariableAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<4; i++)
        {
            noalias(NodalVariableAux) = Geom[i].FastGetSolutionStepValue(Variable);
            rNodalVariableVector[index++] = NodalVariableAux[0];
            rNodalVariableVector[index++] = NodalVariableAux[1];
            rNodalVariableVector[index++] = NodalVariableAux[2];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void AssembleUBlockVector(Vector& rRightHandSideVector, const array_1d<double,4>& UBlockVector)
    {
        //Line_2d_2  
        unsigned int Global_i, Local_i;

        for(unsigned int i = 0; i < 2; i++)
        {
            Global_i = i * (2 + 1);
            Local_i  = i * 2;

            rRightHandSideVector[Global_i]   += UBlockVector[Local_i];
            rRightHandSideVector[Global_i+1] += UBlockVector[Local_i+1];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockVector(Vector& rRightHandSideVector, const array_1d<double,9>& UBlockVector)
    {      
        //Triangle_3d_3  
        unsigned int Global_i, Local_i;

        for(unsigned int i = 0; i < 3; i++)
        {
            Global_i = i * (3 + 1);
            Local_i  = i * 3;

            rRightHandSideVector[Global_i]   += UBlockVector[Local_i];
            rRightHandSideVector[Global_i+1] += UBlockVector[Local_i+1];
            rRightHandSideVector[Global_i+2] += UBlockVector[Local_i+2];
        }
    }
    
    //----------------------------------------------------------------------------------------

    static inline void AssembleUBlockVector(Vector& rRightHandSideVector, const array_1d<double,12>& UBlockVector)
    {        
        //Quadrilateral_3d_4
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

	static inline void AssembleUPMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,4,2>& UPBlockMatrix)
    {        
        //Line_2d_2
        unsigned int Global_i, Global_j, Local_i;

        for(unsigned int i = 0; i < 2; i++)
        {
            Global_i = i * (2 + 1);
            Local_i = i * 2;

            for(unsigned int j = 0; j < 2; j++)
            {
                Global_j = j * (2 + 1) + 2;

                rLeftHandSideMatrix(Global_i,Global_j)   += UPBlockMatrix(Local_i,j);
                rLeftHandSideMatrix(Global_i+1,Global_j) += UPBlockMatrix(Local_i+1,j);
            }
        }
    }
    
//----------------------------------------------------------------------------------------
 
	static inline void AssembleUPMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,9,3>& UPBlockMatrix)
    {
        //Triangle_3d_3  
        unsigned int Global_i, Global_j, Local_i;

        for(unsigned int i = 0; i < 3; i++)
        {
            Global_i = i * (3 + 1);
            Local_i = i * 3;

            for(unsigned int j = 0; j < 3; j++)
            {
                Global_j = j * (3 + 1) + 3;

                rLeftHandSideMatrix(Global_i,Global_j)   += UPBlockMatrix(Local_i,j);
                rLeftHandSideMatrix(Global_i+1,Global_j) += UPBlockMatrix(Local_i+1,j);
                rLeftHandSideMatrix(Global_i+2,Global_j) += UPBlockMatrix(Local_i+2,j);
            }
        }
    }
    
//----------------------------------------------------------------------------------------
 
	static inline void AssembleUPMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,12,4>& UPBlockMatrix)
    {
        //Quadrilateral_3d_4
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

	static inline void AssemblePUMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,2,4>& PUBlockMatrix)
    {        
        //Line_2d_2
        unsigned int Global_i, Global_j, Local_j;

        for(unsigned int i = 0; i < 2; i++)
        {
            Global_i = i * (2 + 1) + 2;

            for(unsigned int j = 0; j < 2; j++)
            {
                Global_j = j * (2 + 1);
                Local_j = j * 2;

                rLeftHandSideMatrix(Global_i,Global_j)   += PUBlockMatrix(i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1) += PUBlockMatrix(i,Local_j+1);
            }
        }
    }

//----------------------------------------------------------------------------------------

	static inline void AssemblePUMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,3,9>& PUBlockMatrix)
    {
        //Triangle_3d_3 
        unsigned int Global_i, Global_j, Local_j;

        for(unsigned int i = 0; i < 3; i++)
        {
            Global_i = i * (3 + 1) + 3;

            for(unsigned int j = 0; j < 3; j++)
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

	static inline void AssemblePUMatrix(Matrix& rLeftHandSideMatrix, const BoundedMatrix<double,4,12>& PUBlockMatrix)
    {
        //Quadrilateral_3d_4
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
    
}; /* Class PoroConditionUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_PORO_CONDITION_UTILITIES defined */
