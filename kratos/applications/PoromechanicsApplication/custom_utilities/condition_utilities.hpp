//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_CONDITION_UTILITIES )
#define  KRATOS_CONDITION_UTILITIES

// Project includes
#include "includes/element.h"

namespace Kratos
{

class ConditionUtilities
{
    
public:
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(boost::numeric::ublas::bounded_matrix<double,2,4>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Line_2d_2
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,2) = Ncontainer(GPoint,1); 
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,3) = Ncontainer(GPoint,1);
    }
    
    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(boost::numeric::ublas::bounded_matrix<double,3,9>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Triangle_3d_3
        rNu(0,0) = Ncontainer(GPoint,0); rNu(0,3) = Ncontainer(GPoint,1); rNu(0,6) = Ncontainer(GPoint,2);
        rNu(1,1) = Ncontainer(GPoint,0); rNu(1,4) = Ncontainer(GPoint,1); rNu(1,7) = Ncontainer(GPoint,2);
        rNu(2,2) = Ncontainer(GPoint,0); rNu(2,5) = Ncontainer(GPoint,1); rNu(2,8) = Ncontainer(GPoint,2);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(boost::numeric::ublas::bounded_matrix<double,3,12>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
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

    static inline void GetDisplacementsVector(array_1d<double,4>& rDisplacementVector, const Element::GeometryType& Geom)
    {
        //Line_2d_2
        array_1d<double,3> DisplacementAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<2; i++)
        {
            noalias(DisplacementAux) = Geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            rDisplacementVector[index++] = DisplacementAux[0];
            rDisplacementVector[index++] = DisplacementAux[1];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetDisplacementsVector(array_1d<double,12>& rDisplacementVector, const Element::GeometryType& Geom)
    {
        //Quadrilateral_3d_4
        array_1d<double,3> DisplacementAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<4; i++)
        {
            noalias(DisplacementAux) = Geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            rDisplacementVector[index++] = DisplacementAux[0];
            rDisplacementVector[index++] = DisplacementAux[1];
            rDisplacementVector[index++] = DisplacementAux[2];
        }
    }
    
    //----------------------------------------------------------------------------------------

    static inline void GetFaceLoadVector(array_1d<double,4>& rFaceLoadVector, const Element::GeometryType& Geom)
    {
        //Line_2d_2
        array_1d<double,3> FaceLoadAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<2; i++)
        {
            noalias(FaceLoadAux) = Geom[i].FastGetSolutionStepValue(FACE_LOAD);
            rFaceLoadVector[index++] = FaceLoadAux[0];
            rFaceLoadVector[index++] = FaceLoadAux[1];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetFaceLoadVector(array_1d<double,9>& rFaceLoadVector, const Element::GeometryType& Geom)
    {
        //Triangle_3d_3
        array_1d<double,3> FaceLoadAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<3; i++)
        {
            noalias(FaceLoadAux) = Geom[i].FastGetSolutionStepValue(FACE_LOAD);
            rFaceLoadVector[index++] = FaceLoadAux[0];
            rFaceLoadVector[index++] = FaceLoadAux[1];
            rFaceLoadVector[index++] = FaceLoadAux[2];
        }
    }

    //----------------------------------------------------------------------------------------

    static inline void GetFaceLoadVector(array_1d<double,12>& rFaceLoadVector, const Element::GeometryType& Geom)
    {
        //Quadrilateral_3d_4
        array_1d<double,3> FaceLoadAux;
        unsigned int index = 0;
        for(unsigned int i=0; i<4; i++)
        {
            noalias(FaceLoadAux) = Geom[i].FastGetSolutionStepValue(FACE_LOAD);
            rFaceLoadVector[index++] = FaceLoadAux[0];
            rFaceLoadVector[index++] = FaceLoadAux[1];
            rFaceLoadVector[index++] = FaceLoadAux[2];
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
    
}; /* Class ConditionUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_CONDITION_UTILITIES defined */
