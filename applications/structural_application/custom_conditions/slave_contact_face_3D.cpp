/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
/* *********************************************************   
*          
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2007-07-17 16:15:08 $
*   Revision:            $Revision: 1.5 $
*
* ***********************************************************/


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/slave_contact_face_3D.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    SlaveContactFace3D::SlaveContactFace3D( IndexType NewId, 
                                  GeometryType::Pointer pGeometry) : 
            Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    SlaveContactFace3D::SlaveContactFace3D( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties) : 
            Condition( NewId, pGeometry, pProperties )
    {
//         for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
//         {
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_X));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Y));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Z));
//         }
        mpMasterElements = ContactMasterContainerType::Pointer( new ContactMasterContainerType() );
        GetValue( LAMBDAS ).resize(GetGeometry().IntegrationPoints().size());
        noalias(GetValue( LAMBDAS )) = ZeroVector( GetGeometry().IntegrationPoints().size());
        GetValue( LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2 );
        noalias(GetValue( LAMBDAS_T )) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( GAPS ).resize( GetGeometry().IntegrationPoints().size());
        noalias(GetValue( GAPS )) = ZeroVector( GetGeometry().IntegrationPoints().size());
        GetValue( DELTA_LAMBDAS ).resize( GetGeometry().IntegrationPoints().size() );
        noalias(GetValue( DELTA_LAMBDAS )) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( DELTA_LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2 );
        noalias(GetValue( DELTA_LAMBDAS_T )) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( PENALTY ).resize( GetGeometry().IntegrationPoints().size() );
        noalias(GetValue( PENALTY )) =  ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( PENALTY_T ).resize( GetGeometry().IntegrationPoints().size() );
		noalias(GetValue( PENALTY_T )) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( IS_CONTACT_SLAVE ) = 1;
		GetValue( IS_CONTACT_MASTER ) = 0;    }
    
    Condition::Pointer SlaveContactFace3D::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new SlaveContactFace3D(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    /**
     * Destructor. Never to be called manually
     */
    SlaveContactFace3D::~SlaveContactFace3D()
    {
    }
    
    //************************************************************************************
    //************************************************************************************
    
    /**
     * returns the tangential vectors of the current surface in an arbitrary point
     */
//     Matrix SlaveContactFace3D::TangentialVectors( GeometryType::CoordinatesArrayType& rPoint )
//     {
// //         setting up result matrix
//         Matrix T = ZeroMatrix( 2, 3 );
// //         shape function gradients
//         Matrix DN = ZeroMatrix(GetGeometry().PointsNumber(),2);
//         GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
//         calculating tangential vectors
//         for( int n=0; n<GetGeometry().PointsNumber(); n++ )
//         {
//             T(0,0) += GetGeometry().GetPoint(n).X()*DN(n,0);
//             T(0,1) += GetGeometry().GetPoint(n).Y()*DN(n,0);
//             T(0,2) += GetGeometry().GetPoint(n).Z()*DN(n,0);
//             T(1,0) += GetGeometry().GetPoint(n).X()*DN(n,1);
//             T(1,1) += GetGeometry().GetPoint(n).Y()*DN(n,1);
//             T(1,2) += GetGeometry().GetPoint(n).Z()*DN(n,1);
//         }
//         return( T );
//     }
    //************************************************************************************
    //************************************************************************************
    
    /**
     * returns normal vector in arbitrary point 
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector 
     */
//     Vector SlaveContactFace3D::NormalVector( GeometryType::CoordinatesArrayType& rPoint )
//     {
//         Vector Result = ZeroVector(3);
// //         getting tangential vectors
//         Matrix T = TangentialVectors( rPoint );
//         calculating normal vector
//         Result[0] = T(0,1)*T(1,2)-T(0,2)*T(1,1);
//         Result[1] = T(0,2)*T(1,0)-T(0,0)*T(1,2);
//         Result[2] = T(0,0)*T(1,1)-T(0,1)*T(1,0);
//         SD_MathUtils<double>::Normalize( Result );
//         return( Result );
//     }
  
    //************************************************************************************
    //************************************************************************************
   
   
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void SlaveContactFace3D::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
    {
//         calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType matrix = Matrix();
        CalculateAll( matrix, rRightHandSideVector, 
                      rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, 
                      CalculateResidualVectorFlag);
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * calculates this contact element's local contributions
     */
    void SlaveContactFace3D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                              VectorType& rRightHandSideVector, 
                                              ProcessInfo& rCurrentProcessInfo)
    {/*
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);*/
    }
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED!
     */
//     void SlaveContactFace3D::CalculateCrossElementarySystemContributions( SecondaryConditionContainerType& SecondaryConditions, EquationIdVectorType& PrimaryEquationId, ProcessInfo& rCurrentProcessInfo )
//     {
    /*
        KRATOS_TRY
        
        SecondaryConditions.clear();
        EquationIdVector( PrimaryEquationId, rCurrentProcessInfo );
        
        //checking for existent master surface
        if( mpMasterElements->size() != GetGeometry().IntegrationPoints().size() ) 
        {
            KRATOS_ERROR( std::invalid_argument, "Not enough master elements for each integration point found", "" );
        }
        //manually setting dimension to 3:
        unsigned int dim = 3;
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points 
                = GetGeometry().IntegrationPoints();
        const GeometryType::ShapeFunctionsGradientsType& DN_De 
                = GetGeometry().ShapeFunctionsLocalGradients();
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();
                
        //calculating current jacobian
        GeometryType::JacobiansType J;
        J = GetGeometry().Jacobian(J);
        //calculating determinant of jacobian
        Vector<double> DetJ( integration_points.size() );
        GetGeometry().DeterminantOfJacobian( DetJ );
        
        //loop over all integration points
        unsigned int PointNumber = 0;
        for( ContactMasterContainerType::ptr_iterator it
             = mpMasterElements->ptr_begin();
             it != mpMasterElements->ptr_end(); ++it )
        {
            MasterContactFace3D current = **it;
            //storing master element's equation ID vector
            SecondaryCondition CrossCondition;
            EquationIdVectorType CurrentSecondaryIdVector;
            current.EquationIdVector( CurrentSecondaryIdVector, rCurrentProcessInfo );
            CrossCondition.SecondaryEquationIdVector = CurrentSecondaryIdVector;
            //LHS sizes
            unsigned int PrimarySize = dim*GetGeometry().size();
            unsigned int SecondarySize = dim*current.GetGeometry().size();
            //setting up LHS contribution
            LHS_ContributionType CurrentLHS = ZeroMatrix( PrimarySize, SecondarySize );
            
            //calculating normal vector on slave elements
            Vector vSlave = ZeroVector(3);
            Vector t1 = ZeroVector(3);//first tangential vector
            Vector t2 = ZeroVector(3);//second tangential vector
            for( int n=0; n<GetGeometry().PointsNumber(); n++ )
            {
                t1[0] += GetGeometry().GetPoint(n).X()*DN_De[PointNumber][n][0];
                t1[1] += GetGeometry().GetPoint(n).Y()*DN_De[PointNumber][n][0];
                t1[2] += GetGeometry().GetPoint(n).Z()*DN_De[PointNumber][n][0];
                t2[0] += GetGeometry().GetPoint(n).X()*DN_De[PointNumber][n][1];
                t2[1] += GetGeometry().GetPoint(n).Y()*DN_De[PointNumber][n][1];
                t2[2] += GetGeometry().GetPoint(n).Z()*DN_De[PointNumber][n][1];
            }
            vSlave[0] = t1[1]*t2[2]-t1[2]*t2[1];
            vSlave[1] = t1[2]*t2[0]-t1[0]*t2[2];
            vSlave[2] = t1[0]*t2[1]-t1[1]*t2[0];
            //calculating dA
            double dA = MathUtils<double>::Norm3( vSlave );
            //getting normal vector on master element
            Vector coords = ZeroVector(3);
            //calculating global coordinates of current integration point 
            //loop over all nodes
            for( int n=0; n<GetGeometry().PointsNumber(); n++ )
            {
                coords[0] += GetGeometry().GetPoint(n).X()*Ncontainer[PointNumber][n];
                coords[1] += GetGeometry().GetPoint(n).Y()*Ncontainer[PointNumber][n];
                coords[2] += GetGeometry().GetPoint(n).Z()*Ncontainer[PointNumber][n];
            }
            GeometryType::PointType currentPoint(3, coords[0], coords[1], coords[2] );
            Vector v = current.NormalVector( mpContactPartnersLocal[PointNumber] );
            Vector SecondaryN = ZeroVector( current.GetGeometry().size() );
            for( int NodeNumber = 0; NodeNumber < SecondaryN.size(); NodeNumber++ )
            {
                SecondaryN[NodeNumber] = current.GetGeometry().ShapeFunctionValue(
                        NodeNumber, mpContactPartnersLocal[PointNumber] );
            }
            //updating gap function
            double Gap = -(v[0]*(currentPoint.X()-mpContactPartnersGlobal[PointNumber].X())
                       +v[1]*(currentPoint.Y()-mpContactPartnersGlobal[PointNumber].Y())
                        +v[2]*(currentPoint.Z()-mpContactPartnersGlobal[PointNumber].Z()));
            
            //calculating normal contact stress
            double normalStress = mLambdas[PointNumber]+penalty*mGaps[PointNumber];
            if( normalStress < 0.0 ) normalStress = 0.0;
            //Getting Integration Weight for current quadrature Pointer
            double IntegrationWeight = integration_points[PointNumber].Weight();
            
            if( mGaps[PointNumber] > 0.0 )
            {
                //adding cross elemental contribution
                for( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
                {
                    for( unsigned int i=0; i<dim; i++ )
                    {
                        for( unsigned int sec=0; sec < current.GetGeometry().size(); sec++ )
                        {
                            for( unsigned int j=0; j<dim; j++ )
                            {
                                CurrentLHS[prim*dim+i][sec*dim+j] -=
                                        Ncontainer[PointNumber][prim]*SecondaryN[sec]*IntegrationWeight
                                        *dA*penalty*v[i]*v[j];
                            }
                        }
                    }
                }
            }
            CrossCondition.LHS_Contribution = CurrentLHS;
            SecondaryConditions.push_back( CrossCondition );
            PointNumber++;
        }
        KRATOS_CATCH("")
    */
//     }//CalculateCrossElementarySystemContributions
           
    //************************************************************************************
    //************************************************************************************
    /**
     * Calculation of stiffness matrix contribution due to linearization of normal vector.
     * In this mehtod the components referring to the slave element as primary index are
     * calculated
     * Hence there is a loop over all registered Master contact elements serving as secondary indices.
     * The result is stored in a contribution container returned to the builder
     * @param SecondaryConditions the contribution container the resulting contributions are stored in
     * @param PrimaryEquationId the equationID vector of the primary element (here: the current slave
     * element)
     * @param rCurrentProcessInfo the current system's process info
     * 
     * (REMOVED)
     */
//     void SlaveContactFace3D::CalculateNormalLinearizationElementarySystemContributions(
//             SecondaryConditionContainerType& SecondaryConditions, 
//     EquationIdVectorType& PrimaryEquationId, ProcessInfo& rCurrentProcessInfo )
//     {
    /*
        SecondaryConditions.clear();
        EquationIdVector( PrimaryEquationId, rCurrentProcessInfo );
        
        //checking for existent master surface
        if( mpMasterElements->size() != GetGeometry().IntegrationPoints().size() ) 
        {
            KRATOS_ERROR( std::invalid_argument, "Not enough master elements for each integration point found", "" );
        }
        //integration points
        const GeometryType::IntegrationPointsArrayType& integration_points 
                = GetGeometry().IntegrationPoints();
        //shape function gradients on slave element 
        const GeometryType::ShapeFunctionsGradientsType& DN_De 
                = GetGeometry().ShapeFunctionsLocalGradients();
        
        unsigned int PointNumber = 0;
        //loop over all master conditions (or integration points)
        for( ContactMasterContainerType::ptr_iterator it= mpMasterElements->ptr_begin();
             it != mpMasterElements->ptr_end(); it++ )
        {
            
            MasterContactFace3D current = **it;
            
            //setting up variables
            EquationIdVectorType PrimaryEquationId;
            EquationIdVector( PrimaryEquationId, rCurrentProcessInfo );
            EquationIdVectorType SecondaryEquationId;
            current.EquationIdVector( SecondaryEquationId, rCurrentProcessInfo );
            
            //auxiliary information
            const unsigned int dim = 3;
            
            GeometryType::PointType currentPoint(3);
            currentPoint.X() = GetGeometry().IntegrationPoints()[PointNumber].X();
            currentPoint.Y() = GetGeometry().IntegrationPoints()[PointNumber].Y();
            currentPoint.Z() = GetGeometry().IntegrationPoints()[PointNumber].Z();
            
            //calculating normal vector on slave elements
            Vector vSlave = ZeroVector(3);
            Vector t1 = ZeroVector(3);//first tangential vector
            Vector t2 = ZeroVector(3);//second tangential vector
            for( int n=0; n<GetGeometry().PointsNumber(); n++ )
            {
                t1[0] += GetGeometry().GetPoint(n).X()*DN_De[PointNumber][n][0];
                t1[1] += GetGeometry().GetPoint(n).Y()*DN_De[PointNumber][n][0];
                t1[2] += GetGeometry().GetPoint(n).Z()*DN_De[PointNumber][n][0];
                t2[0] += GetGeometry().GetPoint(n).X()*DN_De[PointNumber][n][1];
                t2[1] += GetGeometry().GetPoint(n).Y()*DN_De[PointNumber][n][1];
                t2[2] += GetGeometry().GetPoint(n).Z()*DN_De[PointNumber][n][1];
            }
            vSlave[0] = t1[1]*t2[2]-t1[2]*t2[1];
            vSlave[1] = t1[2]*t2[0]-t1[0]*t2[2];
            vSlave[2] = t1[0]*t2[1]-t1[1]*t2[0];
            //calculating dA
            double dA = MathUtils<double>::Norm3( vSlave );
            double Gap = mGaps[PointNumber];
            double aStress = mLambdas[PointNumber]+penalty*Gap;
            if( aStress < 0.0 ) aStress = 0.0;
//             else std::cout << "contact stress: " << normalStress << std::endl;
            //Getting Integration Weight for current quadrature Pointer
            double Weight = integration_points[PointNumber].Weight();
            aStress = aStress*Weight*dA;
            
            //first order derivatives of shape function
            Matrix DN = ZeroMatrix(GetGeometry().PointsNumber(),2);
            current.GetGeometry().ShapeFunctionsGradients( DN, mpContactPartnersLocal[PointNumber] );
            
            //second order derivatives of shape functions
            GeometryType::ShapeFunctionsGradientsType D2N; 
            GetGeometry().SecondOrderShapeFunctionDerivatives( D2N, mpContactPartnersLocal[PointNumber] );
            
                    
            if( Gap != 0.0 )
            {
                
                Vector normalVector = current.NormalVector( mpContactPartnersLocal[PointNumber] );
                Matrix T = current.TangentialVectors( mpContactPartnersLocal[PointNumber] );
            
            //LHS sizes
                unsigned int PrimarySize = dim*GetGeometry().size();
                unsigned int SecondarySize = dim*current.GetGeometry().size();
            //setting up result
                LHS_ContributionType CurrentLHS;
                
            //linearization of normal vector
                Matrix m = ZeroMatrix(2,2);
            
                for( int i=0; i<2; i++ )
                {
                    for( int j=0; j<2; j++ )
                    {
                        m[i][j] = T[i][0]*T[j][0]+T[i][1]*T[j][1]+T[i][2]*T[j][2];
                    }
                }
            
            //setting up Phi hypermatrix (2,2,3)
                vector<Matrix> Phi(2);
                for( int i=0; i<2; i++ )
                {
                    Phi[i] = ZeroMatrix(2,3);
                    for( int j=0; j<2; j++ )
                    {
                        for( int n=0; n<current.GetGeometry().size(); n++ )
                        {
                            Phi[i][j][0] += current.GetGeometry().GetPoint(n).X()*D2N[n][i][j];
                            Phi[i][j][1] += current.GetGeometry().GetPoint(n).Y()*D2N[n][i][j];
                            Phi[i][j][2] += current.GetGeometry().GetPoint(n).Z()*D2N[n][i][j];
                        }
                    }
                }
            // Kappa matrix
                Matrix Kappa = ZeroMatrix(2,2);
            
                for( int i=0; i<2; i++ )
                {
                    for( int j=0; j<2; j++ )
                    {
                        Kappa[i][j] = Phi[i][j][0]*normalVector[0]
                                + Phi[i][j][1]*normalVector[1]
                                + Phi[i][j][2]*normalVector[2];
                    }
                }
            
            //A matrix
                Matrix A = ZeroMatrix(2,2);
                for( int i=0; i<2; i++ )
                {
                    for( int j=0; j<2; j++ )
                    {
                        A[i][j] = m[i][j]+Gap*Kappa[i][j];
                    }
                }
                
                //calculating contributions Slave-Master 
                CurrentLHS = ZeroMatrix( PrimarySize, SecondarySize );
                //loop over all nodes on primary element 
                for( int prim = 0; prim < GetGeometry().size(); prim++ )
                {
                    double N = GetGeometry().ShapeFunctionValue( 
                                prim, currentPoint );
                    //loop over all dimensions 
                    for( int i=0; i<dim; i++ )
                    {
                        
                        double PsiPrim1 = N*T[0][i];
                        double PsiPrim2 = N*T[1][i];
                        //loop over all nodes of secondary element
                        for( int sec = 0; sec < current.GetGeometry().size(); sec++ )
                        {
                            double Nsec = current.GetGeometry().ShapeFunctionValue(
                                        sec, mpContactPartnersLocal[PointNumber] );
                            //loop over all dimensions
                            for( int j=0; j<dim; j++ )
                            {
                                double PsiSec1 = -Nsec*T[0][j]-Gap*normalVector[j]*DN[sec][0];
                                double PsiSec2 = -Nsec*T[1][j]-Gap*normalVector[j]*DN[sec][1];
                                //loop over gamma
                                for( int gamma = 0; gamma < 2; gamma++ )
                                {
                                    //loop over beta
                                    for( int beta=0; beta < 2; beta++ )
                                    {               
                                        CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                                *((Gap*normalVector[j]*DN[sec][gamma]
                                                +(Kappa[0][gamma]/A[0][0]+Kappa[1][gamma]/A[0][1])
                                                *PsiSec1+(Kappa[0][gamma]/A[1][0]+Kappa[1][gamma]/A[1][1])
                                                *PsiSec2)*m[gamma][beta]
                                                *((Kappa[0][beta]/A[0][0]+Kappa[1][beta]/A[0][1])
                                                *PsiPrim1+(Kappa[0][beta]/A[1][0]+Kappa[1][beta]/A[1][1])
                                                *PsiPrim2));
                                        
                                        CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                                *(Kappa[gamma][beta]*(1.0/A[0][beta]*PsiSec1+1.0/A[1][beta]*PsiSec2)
                                                *(1.0/A[0][gamma]*PsiPrim1+1.0/A[1][gamma]*PsiPrim2));
                                        
                                    }//loop over beta
                                    
                                    CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                            *((1.0/A[0][gamma]*PsiPrim1+1.0/A[1][gamma]*PsiPrim2)
                                            *normalVector[j]*DN[sec][gamma]);
                                }//loop over gamma
                            }//loop over j
                        }//loop over sec
                    }//loop over i
                }//loop over prim 
                
                //adding contribution
                SecondaryCondition FirstContribution;
                FirstContribution.SecondaryEquationIdVector = SecondaryEquationId;
                FirstContribution.PrimaryEquationIdVector = PrimaryEquationId;
                FirstContribution.LHS_Contribution = CurrentLHS;
                SecondaryConditions.push_back( FirstContribution );
                
                //calculating contributions Slave-Slave  
                CurrentLHS = ZeroMatrix( PrimarySize, PrimarySize );
                //loop over all nodes on primary element 
                for( int prim = 0; prim < GetGeometry().size(); prim++ )
                {
                    double N = GetGeometry().ShapeFunctionValue( 
                                prim, currentPoint );
                    //loop over all dimensions 
                    for( int i=0; i<dim; i++ )
                    {
                        
                        double PsiPrim1 = N*T[0][i];
                        double PsiPrim2 = N*T[1][i];
                        //loop over all nodes of secondary element
                        for( int sec = 0; sec < GetGeometry().size(); sec++ )
                        {
                            double Nsec = GetGeometry().ShapeFunctionValue(
                                        sec, currentPoint );
                            //loop over all dimensions
                            for( int j=0; j<dim; j++ )
                            {
                                double PsiSec1 = Nsec*T[0][j];
                                double PsiSec2 = Nsec*T[1][j];
                                //loop over gamma
                                for( int gamma = 0; gamma < 2; gamma++ )
                                {
                                    //loop over beta
                                    for( int beta=0; beta < 2; beta++ )
                                    {
                                        CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                                *(Gap*((Kappa[0][gamma]/A[0][0]+Kappa[1][gamma]/A[0][1])
                                                *PsiSec1+(Kappa[0][gamma]/A[1][0]+Kappa[1][gamma]/A[1][1])
                                                *PsiSec2)*m[gamma][beta]
                                                *((Kappa[0][beta]/A[0][0]+Kappa[1][beta]/A[0][1])
                                                *PsiPrim1+(Kappa[0][beta]/A[1][0]+Kappa[1][beta]/A[1][1])
                                                *PsiPrim2));
                                        CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                                *(Kappa[gamma][beta]*(1.0/A[0][beta]*PsiSec1+1.0/A[1][beta]
                                                *PsiSec2)*(1.0/A[0][gamma]*PsiPrim1+1.0/A[1][gamma]*PsiPrim2));
                                    }//loop over beta
                                }//loop over gamma
                            }//loop over j
                        }//loop over sec
                    }//loop over i
                }//loop over prim 
                
                
                
                //adding contribution
                SecondaryCondition SecondContribution;
                SecondContribution.SecondaryEquationIdVector = PrimaryEquationId;
                SecondContribution.PrimaryEquationIdVector = PrimaryEquationId;
                SecondContribution.LHS_Contribution = CurrentLHS;
                SecondaryConditions.push_back( SecondContribution );
                
                
            }//if Gap != 0
            //updating integration point index
            PointNumber++;       
            
        }//loop over all master conditions
    */
//     }//CalculateNormalLinearizationElementarySystemContributions
    //************************************************************************************
    //************************************************************************************
    /**
     * calculates the contact related contributions to the system
     * does nothing as assembling is switched to link condition
     */
    void SlaveContactFace3D::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                      VectorType& rRightHandSideVector, 
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
        /*
        KRATOS_TRY
        //checking for existent master surface
        if( mpMasterElements->size() != GetGeometry().IntegrationPoints().size() ) 
        {
            KRATOS_ERROR( std::invalid_argument, "Not enough master elements for each integration point found", "" );
        }
        
        unsigned int number_of_nodes = GetGeometry().size();
        
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        //manually setting dimension to 3:
        dim = 3;
        //resizing as needed the LHS
        int MatSize=number_of_nodes*dim;
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            rLeftHandSideMatrix.resize(MatSize,MatSize);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize); //resetting LHS
        }
        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            rRightHandSideVector.resize(MatSize);
            rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
        }
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points 
                = GetGeometry().IntegrationPoints();
        const GeometryType::ShapeFunctionsGradientsType& DN_De 
                = GetGeometry().ShapeFunctionsLocalGradients();
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();
                
        //calculating current jacobian
        GeometryType::JacobiansType J;
        J = GetGeometry().Jacobian(J);
        //calculating determinant of jacobian
        Vector<double> DetJ( integration_points.size() );
        GetGeometry().DeterminantOfJacobian( DetJ );
                
        //auxiliary terms
        //Vector BodyForce;	
        //Matrix B;
        //Matrix F(dim,dim);
        //Matrix D;
        //Matrix C(dim,dim);
        //Vector StrainVector;
        //Vector StressVector;
        
        ////sizing work matrices
        //MatrixType DN_DX(DN_De[0].size1(),DN_De[0].size2());
        
        //loop over all integration points
        unsigned int PointNumber = 0;
        for( ContactMasterContainerType::ptr_iterator it
             = mpMasterElements->ptr_begin();
             it != mpMasterElements->ptr_end(); ++it )
        {
            //calculating normal vector on slave elements
            Vector vSlave = ZeroVector(3);
            Vector t1 = ZeroVector(3);//first tangential vector
            Vector t2 = ZeroVector(3);//second tangential vector
            for( int n=0; n<GetGeometry().PointsNumber(); n++ )
            {
                t1[0] += GetGeometry().GetPoint(n).X()*DN_De[PointNumber][n][0];
                t1[1] += GetGeometry().GetPoint(n).Y()*DN_De[PointNumber][n][0];
                t1[2] += GetGeometry().GetPoint(n).Z()*DN_De[PointNumber][n][0];
                t2[0] += GetGeometry().GetPoint(n).X()*DN_De[PointNumber][n][1];
                t2[1] += GetGeometry().GetPoint(n).Y()*DN_De[PointNumber][n][1];
                t2[2] += GetGeometry().GetPoint(n).Z()*DN_De[PointNumber][n][1];
            }
            vSlave[0] = t1[1]*t2[2]-t1[2]*t2[1];
            vSlave[1] = t1[2]*t2[0]-t1[0]*t2[2];
            vSlave[2] = t1[0]*t2[1]-t1[1]*t2[0];
            //calculating dA
            double dA = MathUtils<double>::Norm3( vSlave );
            //getting normal vector on master element
            MasterContactFace3D current = **it;
            Vector coords = ZeroVector(3);
            //calculating global coordinates of current integration point 
            //loop over all nodes
            for( int n=0; n<GetGeometry().PointsNumber(); n++ )
            {
                coords[0] += GetGeometry().GetPoint(n).X()*Ncontainer[PointNumber][n];
                coords[1] += GetGeometry().GetPoint(n).Y()*Ncontainer[PointNumber][n];
                coords[2] += GetGeometry().GetPoint(n).Z()*Ncontainer[PointNumber][n];
            }
            GeometryType::PointType currentPoint(3, coords[0], coords[1], coords[2] );
            Vector v = current.NormalVector( mpContactPartnersLocal[PointNumber] );
            //updating contact partner coordinates
            mpContactPartnersGlobal[PointNumber].X() = 0.0;
            mpContactPartnersGlobal[PointNumber].Y() = 0.0;
            mpContactPartnersGlobal[PointNumber].Z() = 0.0;
            for( int n=0; n<current.GetGeometry().PointsNumber(); n++ )
            {
                double N = current.GetGeometry().ShapeFunctionValue( n, mpContactPartnersLocal[PointNumber] );
                mpContactPartnersGlobal[PointNumber].X() += current.GetGeometry().GetPoint(n).X()*N;
                mpContactPartnersGlobal[PointNumber].Y() += current.GetGeometry().GetPoint(n).Y()*N;
                mpContactPartnersGlobal[PointNumber].Z() += current.GetGeometry().GetPoint(n).Z()*N;
            }
            //updating gap function
            double Gap = -(v[0]*(currentPoint.X()-mpContactPartnersGlobal[PointNumber].X())
                        +v[1]*(currentPoint.Y()-mpContactPartnersGlobal[PointNumber].Y())
                        +v[2]*(currentPoint.Z()-mpContactPartnersGlobal[PointNumber].Z()));
            mGaps[PointNumber] = Gap;
//             std::cout << "penetration: " << Gap << std::endl;
            if( Gap > 0 )
            {
//                 std::cout << "!!!!!! CONTACT !!!!!!!" << std::endl;
                
            }
            //calculating normal contact stress
            double normalStress = mLambdas[PointNumber]+penalty*mGaps[PointNumber];
            if( normalStress < 0.0 ) normalStress = 0.0;
//             else std::cout << "contact stress: " << normalStress << std::endl;
            //Getting Integration Weight for current quadrature Pointer
            double IntegrationWeight = integration_points[PointNumber].Weight();
            //adding contributions to the residual vector
            if (CalculateResidualVectorFlag == true)
            {
                if (normalStress != 0.0)
                {
                    CalculateAndAdd_PressureForce( rRightHandSideVector,
                        Ncontainer[PointNumber],
                        v, normalStress, IntegrationWeight, dA );
                    Vector LocalCoords = ZeroVector(3);
                    LocalCoords[0] = mpContactPartnersLocal[PointNumber].X();
                    LocalCoords[1] = mpContactPartnersLocal[PointNumber].Y();
                    LocalCoords[2] = mpContactPartnersLocal[PointNumber].Z();
                    EquationIdVectorType EqId;
                    EquationIdVector( EqId, rCurrentProcessInfo );
                    (*it)->AddContactStress( LocalCoords, normalStress, IntegrationWeight, dA, v, penalty, Gap, EqId, Ncontainer[PointNumber] );
                }
            }
            if( CalculateStiffnessMatrixFlag == true )
            {
                if( mGaps[PointNumber] > 0.0 )
                {
                    CalculateAndAddKc( rLeftHandSideMatrix, Ncontainer[PointNumber],
                                       IntegrationWeight, dA, v );
                }
            }
            PointNumber++;
        }
        KRATOS_CATCH("")
        */
    } // CalculateAll
    
    //***********************************************************************
    //***********************************************************************
    /**
     * System matrix contribution due to contact energy
     * TODO: implement mixed elementary contribution
     */
    void SlaveContactFace3D::CalculateAndAddKc( Matrix& K,
            const Vector& N,
            double weight,
            double dA,
            Vector v )
    {/*
        KRATOS_TRY
        //adding pure elemental contribution
        for( unsigned int prim=0; prim < GetGeometry().size(); prim++ )
        {
            for( unsigned int i=0; i < 3; i++ )
            {
                for( unsigned int sec=0; sec < GetGeometry().size(); sec++ )
                {
                    for( unsigned int j=0; j<3; j++ )
                    {
                        K(prim*3+i,sec*3+j) += N[prim]*N[sec]*v[i]*v[j]*weight*dA*penalty;
                    }
                }
            }
            for( unsigned int i=0; i < GetGeometry().size(); i++ )
            {
                
            }
        }
        KRATOS_CATCH("")*/
    }
    
    //***********************************************************************
    //***********************************************************************
    /**
     * TO BE TESTED!!!
     */
    void SlaveContactFace3D::CalculateAndAdd_PressureForce( Vector& residualvector,
            const Vector& N,
            Vector& v3,
            double pressure,
            double weight,
            double DetJ )
    {/*
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = 3;
        for (unsigned int i=0;i<number_of_nodes;i++)
        {
            int index = dim*i;
            double coeff = pressure * N[i] * weight * DetJ;
            residualvector[index]   += coeff*v3[0];
            residualvector[index+1] += coeff*v3[1];
            residualvector[index+2] += coeff*v3[2];
        }*/
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
    void SlaveContactFace3D::EquationIdVector( EquationIdVectorType& rResult, 
                                          ProcessInfo& CurrentProcessInfo )
    {
        /*
        int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int index;
        unsigned int dim = 3;
        rResult.resize(number_of_nodes*dim);
        for (int i=0;i<number_of_nodes;i++)
        {
            index = i*dim;
            rResult[index]   = (GetGeometry()[i].GetDof(DISPLACEMENT_X)).EquationId();
            rResult[index+1] = (GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());
            rResult[index+2] = (GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId());
        }
        */
    }
    //************************************************************************************
    //************************************************************************************
    void SlaveContactFace3D::MasterElementsEquationIdVectors(EquationIdVectorContainerType& rResult,
            ProcessInfo& rCurrentProcessInfo )
    {/*
        rResult.clear();
        EquationIdVectorType EquationIdVector; 
        for( ContactMasterContainerType::ptr_iterator it = mpMasterElements->ptr_begin();
             it != mpMasterElements->ptr_end(); ++it )
        {
            (*it)->EquationIdVector( EquationIdVector, rCurrentProcessInfo );
            rResult.push_back( EquationIdVector );
        }*/
    } 
    
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
    void SlaveContactFace3D::GetDofList( DofsVectorType& ConditionalDofList,
                                    ProcessInfo& CurrentProcessInfo)
    {
        /*
        unsigned int dim = 3;
        ConditionalDofList.resize(GetGeometry().size()*dim);
        unsigned int index;
        for (unsigned int i=0;i<GetGeometry().size();i++)
        {
            index = i*dim;
            ConditionalDofList[index]   = (GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ConditionalDofList[index+1] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            ConditionalDofList[index+2] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
        */
    }
} // Namespace Kratos
