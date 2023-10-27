
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


#if !defined(KRATOS_INTERFACE_ELEMENT_UTILITIES )
#define  KRATOS_INTERFACE_ELEMENT_UTILITIES

// Project includes
#include "includes/element.h"

namespace Kratos
{

class InterfaceElementUtilities
{
    
public:
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // NOTE: 
    // Taking into account that we evaluate the shape functions associated to displacements (Nu) at the lobatto integration points of order 1 (at the mid plane of the interface),
    // and that we need to reduce the dimension of those shape functions in order to interpolate the actual relative displacement at the interface, we must multiply Nu by 2.0. In essence:
    // 2D-quadrilateral: the displacement is interpolated as in a line, but the pore pressure is interpolated as in a quadrilateral
    // 3D-wedge: the displacement is interpolated as in a triangle, but the pore pressure is interpolated as in a wedge
    // 3D-hexahedra: the displacement is interpolated as in a quadrilateral, but the pore pressure is interpolated as in a hexahedra

    static inline void CalculateNuMatrix(BoundedMatrix<double,2,4>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Line_interface_2d_2
        rNu(0,0) = -2.0*Ncontainer(GPoint,0); rNu(0,2) = 2.0*Ncontainer(GPoint,1); 
        rNu(1,1) = -2.0*Ncontainer(GPoint,0); rNu(1,3) = 2.0*Ncontainer(GPoint,1);
    }
    
    //----------------------------------------------------------------------------------------
    
    static inline void CalculateNuMatrix(BoundedMatrix<double,2,8>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Quadrilateral_interface_2d_4
        rNu(0,0) = -2.0*Ncontainer(GPoint,0); rNu(0,2) = -2.0*Ncontainer(GPoint,1); 
        rNu(1,1) = -2.0*Ncontainer(GPoint,0); rNu(1,3) = -2.0*Ncontainer(GPoint,1);
        
        rNu(0,4) =  2.0*Ncontainer(GPoint,2); rNu(0,6) =  2.0*Ncontainer(GPoint,3);
        rNu(1,5) =  2.0*Ncontainer(GPoint,2); rNu(1,7) =  2.0*Ncontainer(GPoint,3);        
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,3,12>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Quadrilateral_interface_3d_4
        rNu(0,0) = -2.0*Ncontainer(GPoint,0); rNu(0,3) = -2.0*Ncontainer(GPoint,1); 
        rNu(1,1) = -2.0*Ncontainer(GPoint,0); rNu(1,4) = -2.0*Ncontainer(GPoint,1); 
        rNu(2,2) = -2.0*Ncontainer(GPoint,0); rNu(2,5) = -2.0*Ncontainer(GPoint,1); 
        
        rNu(0,6) = 2.0*Ncontainer(GPoint,2); rNu(0,9)  = 2.0*Ncontainer(GPoint,3);
        rNu(1,7) = 2.0*Ncontainer(GPoint,2); rNu(1,10) = 2.0*Ncontainer(GPoint,3);
        rNu(2,8) = 2.0*Ncontainer(GPoint,2); rNu(2,11) = 2.0*Ncontainer(GPoint,3);
    }
    
    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,3,18>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Prism_interface_3d_6
        rNu(0,0) = -2.0*Ncontainer(GPoint,0); rNu(0,3) = -2.0*Ncontainer(GPoint,1); rNu(0,6) = -2.0*Ncontainer(GPoint,2); 
        rNu(1,1) = -2.0*Ncontainer(GPoint,0); rNu(1,4) = -2.0*Ncontainer(GPoint,1); rNu(1,7) = -2.0*Ncontainer(GPoint,2); 
        rNu(2,2) = -2.0*Ncontainer(GPoint,0); rNu(2,5) = -2.0*Ncontainer(GPoint,1); rNu(2,8) = -2.0*Ncontainer(GPoint,2); 
        
        rNu(0,9)  = 2.0*Ncontainer(GPoint,3); rNu(0,12) = 2.0*Ncontainer(GPoint,4); rNu(0,15) = 2.0*Ncontainer(GPoint,5);
        rNu(1,10) = 2.0*Ncontainer(GPoint,3); rNu(1,13) = 2.0*Ncontainer(GPoint,4); rNu(1,16) = 2.0*Ncontainer(GPoint,5);
        rNu(2,11) = 2.0*Ncontainer(GPoint,3); rNu(2,14) = 2.0*Ncontainer(GPoint,4); rNu(2,17) = 2.0*Ncontainer(GPoint,5);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuMatrix(BoundedMatrix<double,3,24>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Hexahedral_interface_3d_8
        rNu(0,0) = -2.0*Ncontainer(GPoint,0); rNu(0,3) = -2.0*Ncontainer(GPoint,1); rNu(0,6) = -2.0*Ncontainer(GPoint,2); rNu(0,9)  = -2.0*Ncontainer(GPoint,3);
        rNu(1,1) = -2.0*Ncontainer(GPoint,0); rNu(1,4) = -2.0*Ncontainer(GPoint,1); rNu(1,7) = -2.0*Ncontainer(GPoint,2); rNu(1,10) = -2.0*Ncontainer(GPoint,3);
        rNu(2,2) = -2.0*Ncontainer(GPoint,0); rNu(2,5) = -2.0*Ncontainer(GPoint,1); rNu(2,8) = -2.0*Ncontainer(GPoint,2); rNu(2,11) = -2.0*Ncontainer(GPoint,3);

        rNu(0,12) = 2.0*Ncontainer(GPoint,4); rNu(0,15) = 2.0*Ncontainer(GPoint,5); rNu(0,18) = 2.0*Ncontainer(GPoint,6); rNu(0,21) = 2.0*Ncontainer(GPoint,7);
        rNu(1,13) = 2.0*Ncontainer(GPoint,4); rNu(1,16) = 2.0*Ncontainer(GPoint,5); rNu(1,19) = 2.0*Ncontainer(GPoint,6); rNu(1,22) = 2.0*Ncontainer(GPoint,7);
        rNu(2,14) = 2.0*Ncontainer(GPoint,4); rNu(2,17) = 2.0*Ncontainer(GPoint,5); rNu(2,20) = 2.0*Ncontainer(GPoint,6); rNu(2,23) = 2.0*Ncontainer(GPoint,7);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculateNuElementMatrix(BoundedMatrix<double,3,12>& rNut, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Quadrilateral_interface_2d_4
        rNut(0,0) = -2.0*Ncontainer(GPoint,0); rNut(0,3) = -2.0*Ncontainer(GPoint,1); 
        rNut(1,1) = -2.0*Ncontainer(GPoint,0); rNut(1,4) = -2.0*Ncontainer(GPoint,1); 
        
        rNut(0,6) = 2.0*Ncontainer(GPoint,2); rNut(0,9)  = 2.0*Ncontainer(GPoint,3);
        rNut(1,7) = 2.0*Ncontainer(GPoint,2); rNut(1,10) = 2.0*Ncontainer(GPoint,3);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuElementMatrix(BoundedMatrix<double,4,24>& rNut, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Prism_interface_3d_6
        rNut(0,0) = -2.0*Ncontainer(GPoint,0); rNut(0,4) = -2.0*Ncontainer(GPoint,1); rNut(0,8)  = -2.0*Ncontainer(GPoint,2); 
        rNut(1,1) = -2.0*Ncontainer(GPoint,0); rNut(1,5) = -2.0*Ncontainer(GPoint,1); rNut(1,9)  = -2.0*Ncontainer(GPoint,2); 
        rNut(2,2) = -2.0*Ncontainer(GPoint,0); rNut(2,6) = -2.0*Ncontainer(GPoint,1); rNut(2,10) = -2.0*Ncontainer(GPoint,2); 
        
        rNut(0,12) = 2.0*Ncontainer(GPoint,3); rNut(0,16) = 2.0*Ncontainer(GPoint,4); rNut(0,20)  = 2.0*Ncontainer(GPoint,5);
        rNut(1,13) = 2.0*Ncontainer(GPoint,3); rNut(1,17) = 2.0*Ncontainer(GPoint,4); rNut(1,21)  = 2.0*Ncontainer(GPoint,5);
        rNut(2,14) = 2.0*Ncontainer(GPoint,3); rNut(2,18) = 2.0*Ncontainer(GPoint,4); rNut(2,22)  = 2.0*Ncontainer(GPoint,5);
    }

    //----------------------------------------------------------------------------------------

    static inline void CalculateNuElementMatrix(BoundedMatrix<double,4,32>& rNut, const Matrix& Ncontainer, const unsigned int& GPoint)
    {
        //Hexahedral_interface_3d_8
        rNut(0,0) = -2.0*Ncontainer(GPoint,0); rNut(0,4) = -2.0*Ncontainer(GPoint,1); rNut(0,8)  = -2.0*Ncontainer(GPoint,2); rNut(0,12) = -2.0*Ncontainer(GPoint,3);
        rNut(1,1) = -2.0*Ncontainer(GPoint,0); rNut(1,5) = -2.0*Ncontainer(GPoint,1); rNut(1,9)  = -2.0*Ncontainer(GPoint,2); rNut(1,13) = -2.0*Ncontainer(GPoint,3);
        rNut(2,2) = -2.0*Ncontainer(GPoint,0); rNut(2,6) = -2.0*Ncontainer(GPoint,1); rNut(2,10) = -2.0*Ncontainer(GPoint,2); rNut(2,14) = -2.0*Ncontainer(GPoint,3);

        rNut(0,16) = 2.0*Ncontainer(GPoint,4); rNut(0,20) = 2.0*Ncontainer(GPoint,5); rNut(0,24)  = 2.0*Ncontainer(GPoint,6); rNut(0,28) = 2.0*Ncontainer(GPoint,7);
        rNut(1,17) = 2.0*Ncontainer(GPoint,4); rNut(1,21) = 2.0*Ncontainer(GPoint,5); rNut(1,25)  = 2.0*Ncontainer(GPoint,6); rNut(1,29) = 2.0*Ncontainer(GPoint,7);
        rNut(2,18) = 2.0*Ncontainer(GPoint,4); rNut(2,22) = 2.0*Ncontainer(GPoint,5); rNut(2,26)  = 2.0*Ncontainer(GPoint,6); rNut(2,30) = 2.0*Ncontainer(GPoint,7);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculatePermeabilityMatrix(BoundedMatrix<double,2,2>& rPermeabilityMatrix,
                                                    const double& JointWidth, const double& Transversal_Permeability_Coeff)
    {
        //Quadrilateral_interface_2d_4
        rPermeabilityMatrix(0,0) = JointWidth*JointWidth/12.0;
        rPermeabilityMatrix(1,1) = Transversal_Permeability_Coeff;
    }
    
    //----------------------------------------------------------------------------------------
    
    static inline void CalculatePermeabilityMatrix(BoundedMatrix<double,3,3>& rPermeabilityMatrix,
                                                    const double& JointWidth, const double& Transversal_Permeability_Coeff)
    {
        //Prism_interface_3d_6 & Hexahedral_interface_3d_8
        rPermeabilityMatrix(0,0) = JointWidth*JointWidth/12.0;
        rPermeabilityMatrix(1,1) = JointWidth*JointWidth/12.0;
        rPermeabilityMatrix(2,2) = Transversal_Permeability_Coeff;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculateVoigtVector(array_1d<double,2>& rVoigtVector)
    {
        //Quadrilateral_interface_2d_4
        rVoigtVector[0] = 0.0;
        rVoigtVector[1] = 1.0;
    }
    
    //----------------------------------------------------------------------------------------
    
    static inline void CalculateVoigtVector(array_1d<double,3>& rVoigtVector)
    {
        //Prism_interface_3d_6 & Hexahedral_interface_3d_8
        rVoigtVector[0] = 0.0;
        rVoigtVector[1] = 0.0;
        rVoigtVector[2] = 1.0;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void CalculateLinkPermeabilityMatrix(BoundedMatrix<double,2,2>& rPermeabilityMatrix,
                                                    const double& JointWidth)
    {
        //Quadrilateral_interface_2d_4
        rPermeabilityMatrix(0,0) = JointWidth*JointWidth/12.0;
        rPermeabilityMatrix(1,1) = JointWidth*JointWidth/12.0;
    }
    
    //----------------------------------------------------------------------------------------
    
    static inline void CalculateLinkPermeabilityMatrix(BoundedMatrix<double,3,3>& rPermeabilityMatrix,
                                                    const double& JointWidth)
    {
        //Prism_interface_3d_6 & Hexahedral_interface_3d_8
        rPermeabilityMatrix(0,0) = JointWidth*JointWidth/12.0;
        rPermeabilityMatrix(1,1) = JointWidth*JointWidth/12.0;
        rPermeabilityMatrix(2,2) = JointWidth*JointWidth/12.0;
    }

}; /* Class InterfaceElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_INTERFACE_ELEMENT_UTILITIES defined */
