
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
#include "includes/serializer.h"
#include "includes/constitutive_law.h"

// Application includes
#include "poromechanics_application_variables.h"

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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void AddInitialInterfaceStresses2D(Vector& rStressVector,
                                                        ConstitutiveLaw::Parameters& rValues,
                                                        const Element::GeometryType& geometry)
    {
        // Obtain necessary variables
        const Vector& N = rValues.GetShapeFunctionsValues();
        const unsigned int number_of_nodes = geometry.size();
        const unsigned int dimension = 2;
        const unsigned int strainSize = 3;   

        // Create components to initialise the stresses
        Vector nodal_initial_stress_vector(strainSize);
        Matrix nodal_initial_stress_tensor(dimension, dimension); 
        Vector gp_initial_stress_vector(strainSize);
        noalias(gp_initial_stress_vector) = ZeroVector(strainSize);
        Vector LocalInitialStress(dimension);
        noalias(LocalInitialStress) = ZeroVector(dimension);

        // Contribution of each of the nodal values to the corresponding GP
        for (unsigned int i=0; i < number_of_nodes; i++) {
            const Matrix& r_initial_stress_tensor = geometry[i].GetSolutionStepValue(INITIAL_STRESS_TENSOR);
            for(unsigned int j=0; j < dimension; j++) {
                for(unsigned int k=0; k < dimension; k++) {
                    nodal_initial_stress_tensor(j,k) = r_initial_stress_tensor(j,k);
                }
            }
            noalias(nodal_initial_stress_vector) = MathUtils<double>::StressTensorToVector(nodal_initial_stress_tensor);

            for(unsigned int j=0; j < strainSize; j++){
                gp_initial_stress_vector[j] += N[i] * nodal_initial_stress_vector[j];
            }
        }

        // Define mid-plane points for quadrilateral_interface_2d_4
        array_1d<double, 3> pmid0;
        array_1d<double, 3> pmid1;
        noalias(pmid0) = 0.5 * (geometry.GetPoint( 0 ) + geometry.GetPoint( 3 ));
        noalias(pmid1) = 0.5 * (geometry.GetPoint( 1 ) + geometry.GetPoint( 2 ));

        // Unitary vector in local x direction
        array_1d<double, 3> Vx;
        noalias(Vx) = pmid1 - pmid0;
        double inv_norm_x = 1.0/norm_2(Vx);
        Vx[0] *= inv_norm_x; // cos
        Vx[1] *= inv_norm_x; // sin

        // Define and assign the rotation matrix
        BoundedMatrix<double,dimension,strainSize> RotationInterface;
        noalias(RotationInterface) = ZeroMatrix(dimension,strainSize);
        RotationInterface(0,0) = -Vx[1]*Vx[0];
        RotationInterface(0,1) = Vx[1]*Vx[0];
        RotationInterface(0,2) = Vx[0]*Vx[0] - Vx[1]*Vx[1];
        RotationInterface(1,0) = Vx[1]*Vx[1];
        RotationInterface(1,1) = Vx[0]*Vx[0];
        RotationInterface(1,2) = -2*Vx[1]*Vx[0];
    
        // Local stress vector
        noalias(LocalInitialStress) = prod(RotationInterface, trans(gp_initial_stress_vector));
    
        // Update the stress vector
        rStressVector[0] += LocalInitialStress[0];
        rStressVector[1] += LocalInitialStress[1];

    }

    //----------------------------------------------------------------------------------------

    static inline void AddInitialInterfaceStresses3D(Vector& rStressVector,
                                                        ConstitutiveLaw::Parameters& rValues,
                                                        const Element::GeometryType& geometry)
    {
        // Obtain necessary variables
        const Vector& N = rValues.GetShapeFunctionsValues();
        const unsigned int number_of_nodes = geometry.size();
        const unsigned int dimension = 3;
        const unsigned int strainSize = 6;

        // Create components to initialise the stresses
        Vector nodal_initial_stress_vector(strainSize);
        Matrix nodal_initial_stress_tensor(dimension, dimension);
        Vector gp_initial_stress_vector(strainSize);
        noalias(gp_initial_stress_vector) = ZeroVector(strainSize);
        Vector LocalInitialStress(strainSize);
        noalias(LocalInitialStress) = ZeroVector(strainSize);

        // Contribution of each of the nodal values to the corresponding GP
        for (unsigned int i=0; i < number_of_nodes; i++) {
            const Matrix& r_initial_stress_tensor = geometry[i].GetSolutionStepValue(INITIAL_STRESS_TENSOR);
            for (unsigned int j=0; j < dimension; j++) {
                for (unsigned int k=0; k < dimension; k++) {
                    nodal_initial_stress_tensor(j,k) = r_initial_stress_tensor(j,k);
                }
            }
            noalias(nodal_initial_stress_vector) = MathUtils<double>::StressTensorToVector(nodal_initial_stress_tensor);

            for (unsigned int j=0; j < strainSize; j++) {
                gp_initial_stress_vector[j] += N[i] * nodal_initial_stress_vector[j];
            }
        }

        // Create and obtain the rotation matrix
        BoundedMatrix<double,strainSize,strainSize> RotationMatrixInterface;
        noalias(RotationMatrixInterface) = ZeroMatrix(strainSize, strainSize);
        ObtainRotationMatrix(RotationMatrixInterface, geometry);

        // Rotate the stresses by multiplying the nodal stress by the rotation matrix
        noalias(LocalInitialStress) = prod(RotationMatrixInterface, trans(gp_initial_stress_vector));
        
        // Assign the corresponding values to the stress vector
        rStressVector[0] += LocalInitialStress[0];
        rStressVector[1] += LocalInitialStress[1];
        rStressVector[2] += LocalInitialStress[2];

    }

    //----------------------------------------------------------------------------------------

    // Utility to compute the 6x6 rotation matrix to rotate the stresses from global to local
    // The derivation is taken from Cook - Concepts and applications of FEA - Page 274
    static inline void ObtainRotationMatrix(BoundedMatrix<double,6,6>& RotationMatrixInterface, const Element::GeometryType& Geom)
    {
        // Define mid-plane points
        array_1d<double, 3> pmid0;
        array_1d<double, 3> pmid1;
        array_1d<double, 3> pmid2;

        // Differentiate between prisms and heaxhedra and define the coordinates of the mid-plane
        const unsigned int number_of_nodes = Geom.size();
        if (number_of_nodes == 6) {
            noalias(pmid0) = 0.5 * (Geom.GetPoint( 0 ) + Geom.GetPoint( 3 ));
            noalias(pmid1) = 0.5 * (Geom.GetPoint( 1 ) + Geom.GetPoint( 4 ));
            noalias(pmid2) = 0.5 * (Geom.GetPoint( 2 ) + Geom.GetPoint( 5 ));
        } else if (number_of_nodes == 8) {
            noalias(pmid0) = 0.5 * (Geom.GetPoint( 0 ) + Geom.GetPoint( 4 ));
            noalias(pmid1) = 0.5 * (Geom.GetPoint( 1 ) + Geom.GetPoint( 5 ));
            noalias(pmid2) = 0.5 * (Geom.GetPoint( 2 ) + Geom.GetPoint( 6 ));
        }

        // Unitary vector in local x direction
        array_1d<double, 3> Vx;
        noalias(Vx) = pmid1 - pmid0;
        double inv_norm = 1.0/norm_2(Vx);
        Vx[0] *= inv_norm; // l1
        Vx[1] *= inv_norm; // l2
        Vx[2] *= inv_norm; // l3

        // Unitary vector in local z direction
        array_1d<double, 3> Vy;
        noalias(Vy) = pmid2 - pmid0;
        array_1d<double, 3> Vz;
        MathUtils<double>::CrossProduct(Vz, Vx, Vy);
        inv_norm = 1.0/norm_2(Vz);
        Vz[0] *= inv_norm; // n1
        Vz[1] *= inv_norm; // n2
        Vz[2] *= inv_norm; // n3

        // Unitary vector in local y direction
        MathUtils<double>::CrossProduct(Vy, Vz, Vx); // Vy = [m1, m2, m3]

        // Create the inverse of the tranpose matrix as described in Eq. (8.2-8)
        BoundedMatrix<double,3,3> T11;
        T11(0,0) = Vx[0]*Vx[0];
        T11(1,0) = Vx[1]*Vx[1];
        T11(2,0) = Vx[2]*Vx[2];
        T11(0,1) = Vy[0]*Vy[0];
        T11(1,1) = Vy[1]*Vy[1];
        T11(2,1) = Vy[2]*Vy[2];
        T11(0,2) = Vz[0]*Vz[0];
        T11(1,2) = Vz[1]*Vz[1];
        T11(2,2) = Vz[2]*Vz[2];

        BoundedMatrix<double,3,3> T12;
        T12(0,0) = Vx[0]*Vy[0];
        T12(1,0) = Vx[1]*Vy[1];
        T12(2,0) = Vx[2]*Vy[2];
        T12(0,1) = Vy[0]*Vz[0];
        T12(1,1) = Vy[1]*Vz[1];
        T12(2,1) = Vy[2]*Vz[2];
        T12(0,2) = Vz[0]*Vx[0];
        T12(1,2) = Vz[1]*Vx[1];
        T12(2,2) = Vz[2]*Vx[2];

        BoundedMatrix<double,3,3> T21;
        T21(0,0) = 2*Vx[0]*Vx[1];
        T21(1,0) = 2*Vx[1]*Vx[2];
        T21(2,0) = 2*Vx[2]*Vx[0];
        T21(0,1) = 2*Vy[0]*Vy[1];
        T21(1,1) = 2*Vy[1]*Vy[2];
        T21(2,1) = 2*Vy[2]*Vy[0];
        T21(0,2) = 2*Vz[0]*Vz[1];
        T21(1,2) = 2*Vz[1]*Vz[2];
        T21(2,2) = 2*Vz[2]*Vz[0];

        BoundedMatrix<double,3,3> T22;
        T22(0,0) = Vx[0]*Vy[1] + Vx[1]*Vy[0];
        T22(1,0) = Vx[1]*Vy[2] + Vx[2]*Vy[1];
        T22(2,0) = Vx[2]*Vy[0] + Vx[0]*Vy[2];
        T22(0,1) = Vy[0]*Vz[1] + Vy[1]*Vz[0];
        T22(1,1) = Vy[1]*Vz[2] + Vy[2]*Vz[1];
        T22(2,1) = Vy[2]*Vz[0] + Vy[0]*Vz[2];
        T22(0,2) = Vz[0]*Vx[1] + Vz[1]*Vx[0];
        T22(1,2) = Vz[1]*Vx[2] + Vz[2]*Vx[1];
        T22(2,2) = Vz[2]*Vx[0] + Vz[0]*Vx[2];

        // Assign the values to the final Rotation Matrix
        RotationMatrixInterface(0,0) = T11(0,0);
        RotationMatrixInterface(1,0) = T11(1,0);
        RotationMatrixInterface(2,0) = T11(2,0);
        RotationMatrixInterface(0,2) = T11(0,1);
        RotationMatrixInterface(1,1) = T11(1,1);
        RotationMatrixInterface(2,2) = T11(2,1);
        RotationMatrixInterface(0,2) = T11(0,2);
        RotationMatrixInterface(2,2) = T11(1,2);
        RotationMatrixInterface(2,2) = T11(2,2);

        RotationMatrixInterface(0,3) = 2*T12(0,0);
        RotationMatrixInterface(1,3) = 2*T12(1,0);
        RotationMatrixInterface(2,3) = 2*T12(2,0);
        RotationMatrixInterface(0,4) = 2*T12(0,1);
        RotationMatrixInterface(1,4) = 2*T12(1,1);
        RotationMatrixInterface(2,4) = 2*T12(2,1);
        RotationMatrixInterface(0,5) = 2*T12(0,2);
        RotationMatrixInterface(1,5) = 2*T12(1,2);
        RotationMatrixInterface(2,5) = 2*T12(2,2);

        RotationMatrixInterface(3,0) = 0.5*T21(0,0);
        RotationMatrixInterface(4,0) = 0.5*T21(1,0);
        RotationMatrixInterface(5,0) = 0.5*T21(2,0);
        RotationMatrixInterface(3,1) = 0.5*T21(0,1);
        RotationMatrixInterface(4,1) = 0.5*T21(1,1);
        RotationMatrixInterface(5,1) = 0.5*T21(2,1);
        RotationMatrixInterface(3,2) = 0.5*T21(0,2);
        RotationMatrixInterface(4,2) = 0.5*T21(1,2);
        RotationMatrixInterface(5,2) = 0.5*T21(2,2);

        RotationMatrixInterface(3,3) = T22(0,0);
        RotationMatrixInterface(4,3) = T22(1,0);
        RotationMatrixInterface(5,3) = T22(2,0);
        RotationMatrixInterface(3,4) = T22(0,1);
        RotationMatrixInterface(4,4) = T22(1,1);
        RotationMatrixInterface(5,4) = T22(2,1);
        RotationMatrixInterface(3,5) = T22(0,2);
        RotationMatrixInterface(4,5) = T22(1,2);
        RotationMatrixInterface(5,5) = T22(2,2);

    }

}; /* Class InterfaceElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_INTERFACE_ELEMENT_UTILITIES defined */
