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
#include "custom_conditions/master_contact_face_3D.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    MasterContactFace3D::MasterContactFace3D( IndexType NewId, 
                                  GeometryType::Pointer pGeometry) : 
            Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    MasterContactFace3D::MasterContactFace3D( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties) : 
            Condition( NewId, pGeometry, pProperties )
    {
//         for( unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
//         {
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_X));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Y));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Z));
//         }
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
		GetValue( IS_CONTACT_MASTER ) = 1;
		GetValue( IS_CONTACT_SLAVE ) = 0;
    }
    
    Condition::Pointer MasterContactFace3D::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new MasterContactFace3D(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    /**
     * Destructor. Never to be called manually
     */
    MasterContactFace3D::~MasterContactFace3D()
    {
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * returns condition type info
     */
//     int MasterContactFace3D::IsContactType()
//     {
//         return( 2 );
//     }
    //************************************************************************************
    //************************************************************************************
    /**
     * returns the tangential vectors of the current surface in an arbitrary point
     * TODO: TO BE REIMPLEMENTED!!!
     */
//     Matrix MasterContactFace3D::TangentialVectors( GeometryType::CoordinatesArrayType& rPoint )
//     {/*
//         //setting up result matrix
//         Matrix T = ZeroMatrix( 2, 3 );
//         //shape function gradients
//         Matrix DN = ZeroMatrix(GetGeometry().PointsNumber(),2);
//         DN = GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
//         //calculating tangential vectors
//         for( int n=0; n<GetGeometry().PointsNumber(); n++ )
//         {
//             T(0,0) += GetGeometry().GetPoint(n).X()*DN(n,0);
//             T(0,1) += GetGeometry().GetPoint(n).Y()*DN(n,0);
//             T(0,2) += GetGeometry().GetPoint(n).Z()*DN(n,0);
//             T(1,0) += GetGeometry().GetPoint(n).X()*DN(n,1);
//             T(1,1) += GetGeometry().GetPoint(n).Y()*DN(n,1);
//             T(1,2) += GetGeometry().GetPoint(n).Z()*DN(n,1);
//         }
//         return( T );*/
//     }
    //************************************************************************************
    //************************************************************************************
    
    /**
     * returns normal vector in arbitrary point 
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector 
     */
//     Vector MasterContactFace3D::NormalVector( GeometryType::CoordinatesArrayType& rPoint )
//     {/*
//         Vector Result = ZeroVector(3);
//         //getting tangential vectors
//         Matrix T = TangentialVectors( rPoint );
//         //calculating normal vector
//         Result[0] = T(0,1)*T(1,2)-T(0,2)*T(1,1);
//         Result[1] = T(0,2)*T(1,0)-T(0,0)*T(1,2);
//         Result[2] = T(0,0)*T(1,1)-T(0,1)*T(1,0);
//         SD_MathUtils<double>::Normalize( Result );
//         return( Result );*/
//     }
    
    /**
     * returns closest point on current condition element with regard to given point in global coordinates
     * @param rResult a Point in global coordinates being overwritten by the desired information
     * @param rLocalResult a Point in local coordinates being overwritten by the desired information
     * @param rCandidate a Node of the current surface in global coordinates which lies closest to rPoint
     * @param rPoint the point in global coordinates the closest point on the current condition element is to
     * be calculated for
     * @return true if an orthogonal projection of the given point lies within the boundaries of the current 
     * condition element
     */
     /**
     * returns closest point on current condition element with regard to given point in global coordinates
     * @param rResultGlobal a Point in global coordinates being overwritten by the desired information
     * @param rResultLocal a Point in global coordinates being overwritten by the desired information
     * @param rSlaveContactGlobalPoint the point in global coordinates the closest point on the current condition element is to
     * @param rCandidateGlobal the closest node to rSlaveContactGlobalPoint on current
     * surface
     * be calculated for
     * @return true if an orthogonal projection of the given point lies within the boundaries of the current 
     * condition element
      */
    bool MasterContactFace3D::ClosestPoint( GeometryType::CoordinatesArrayType& rResultGlobal, 
                                            GeometryType::CoordinatesArrayType& rResultLocal,
                                            const GeometryType::CoordinatesArrayType& rSlaveContactGlobalPoint,
                                            const GeometryType::CoordinatesArrayType& rCandidateGlobal
                                          )
    {/*
        double Xi1 = 0.0;
        double Xi2 = 0.0;
        double deltaXi1 = 0.0;
        double deltaXi2 = 0.0;
        Matrix localCoords;
        localCoords = GetGeometry().PointsLocalCoordinates( localCoords );
        //determining local coordinates for rResult
        for( int n=0; n<GetGeometry().PointsNumber(); n++ )
        {
            if(    fabs(rCandidateGlobal[0]-GetGeometry().GetPoint(n).X()) < 1e-7
                   && fabs(rCandidateGlobal[1]-GetGeometry().GetPoint(n).Y()) < 1e-7
                   && fabs(rCandidateGlobal[2]-GetGeometry().GetPoint(n).Z()) < 1e-7 )
            {
                Xi1 = localCoords(n,0);
                Xi2 = localCoords(n,1);
                break;
            }
        }
        //setting up LocalPoint
        rResultLocal[0] = Xi1;
        rResultLocal[1] = Xi2;
        rResultLocal[2] = 0.0;
        //setting up rResult
        rResultGlobal = rCandidateGlobal;
        //searching for orthogonal projection
        for( int k=0; k<1000; k++ )
        {
            //setting up tangential vectors
            Vector t1 = ZeroVector(3);//first tangential vector
            Vector t2 = ZeroVector(3);//second tangential vector
            //derivatives of tangential vectors
            Vector dt11 = ZeroVector(3);
            Vector dt12 = ZeroVector(3);
            Vector dt21 = ZeroVector(3);
            Vector dt22 = ZeroVector(3);
            
            //retrieving first order derivatives in current solution point 
            Matrix DN = ZeroMatrix(GetGeometry().PointsNumber(),2);
            GetGeometry().ShapeFunctionsLocalGradients( DN, rResultLocal );
            //retrieving second order derivatives in current solution point 
            GeometryType::ShapeFunctionsSecondDerivativesType D2N;
            GetGeometry().ShapeFunctionsSecondDerivatives( D2N, rResultLocal );
            for( int n=0; n<GetGeometry().PointsNumber(); n++ )
            {
                //contribution to tangential vectors
                t1[0] += GetGeometry().GetPoint(n).X()*row(DN,n)[0];
                t1[1] += GetGeometry().GetPoint(n).Y()*row(DN,n)[0];
                t1[2] += GetGeometry().GetPoint(n).Z()*row(DN,n)[0];
                t2[0] += GetGeometry().GetPoint(n).X()*row(DN,n)[1];
                t2[1] += GetGeometry().GetPoint(n).Y()*row(DN,n)[1];
                t2[2] += GetGeometry().GetPoint(n).Z()*row(DN,n)[1];
                //contribution to derivatives of tangential vectors
                dt11[0] += GetGeometry().GetPoint(n).X()*D2N[n](0,0);
                dt11[1] += GetGeometry().GetPoint(n).Y()*D2N[n](0,0);
                dt11[2] += GetGeometry().GetPoint(n).Z()*D2N[n](0,0);
                dt12[0] += GetGeometry().GetPoint(n).X()*D2N[n](0,1);
                dt12[1] += GetGeometry().GetPoint(n).Y()*D2N[n](0,1);
                dt12[2] += GetGeometry().GetPoint(n).Z()*D2N[n](0,1);
                dt21[0] += GetGeometry().GetPoint(n).X()*D2N[n](1,0);
                dt21[1] += GetGeometry().GetPoint(n).Y()*D2N[n](1,0);
                dt21[2] += GetGeometry().GetPoint(n).Z()*D2N[n](1,0);
                dt22[0] += GetGeometry().GetPoint(n).X()*D2N[n](1,1);
                dt22[1] += GetGeometry().GetPoint(n).Y()*D2N[n](1,1);
                dt22[2] += GetGeometry().GetPoint(n).Z()*D2N[n](1,1);
            }
            //defining auxiliary terms
            double A1 = ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*t1[0])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*t1[1])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*t1[2]);
            double A2 = ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*t2[0])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*t2[1])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*t2[2]);
            double B11 = (-t1[0]*t1[0]-t1[1]*t1[1]-t1[2]*t1[2])
                        + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt11[0])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt11[1])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt11[2]);
            double B12 = (-t2[0]*t1[0]-t2[1]*t1[1]-t2[2]*t1[2])
                        + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt12[0])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt12[1])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt12[2]);
            double B21 = (-t1[0]*t2[0]-t1[1]*t2[1]-t1[2]*t2[2])
                        + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt21[0])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt21[1])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt21[2]);
            double B22 = (-t2[0]*t2[0]-t2[1]*t2[1]-t2[2]*t2[2])
                        + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt22[0])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt22[1])
                        +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt22[2]);
            //calculating update for Xi
            deltaXi1 = -A1*B22/(B11*B22-B12*B21)+A2*B12/(B11*B22-B12*B21);
            deltaXi2 =  A2*B21/(B11*B22-B12*B21)-A2*B11/(B11*B22-B12*B21);
            //updating Xi
            Xi1 += deltaXi1;
            Xi2 += deltaXi2;
            //updating LocalPoint
            rResultLocal[0] = Xi1;
            rResultLocal[1] = Xi2;
            //updating rResult
            rResultGlobal = ZeroVector( 3 );
            rResultGlobal = GetGeometry().GlobalCoordinates( rResultGlobal, rResultLocal );
            
            if( fabs(deltaXi1) < 1e-7 && fabs(deltaXi2) < 1e-7 )
            {
                //check whether contact point lies within elementary boundaries
                if( fabs(Xi1) <= 1.0 && fabs(Xi2) <= 1.0 )
                {
//                     std::cout << "found matching point after " << k << " iteratons" << std::endl;
                    return true;
                }
                else
                {
//                     std::cout << "no matching point inside this element" << std::endl;
                    return false;
                }
            }
        }*/
        return false;
    }
    //************************************************************************************
    //************************************************************************************
    /**
     * applies the contact stress from a dedicated slave condition.
     * @param coords the coordinates of the slave condition's partner point on the
     * current master condition in local coordinates
     * @param Stress the value of the current contact stress in current point, augmented by
     * the current integration weight and differential area in slave element
     * @param NormalDirection the normal vector on current master surface
     * 
     * REMOVED: assembling is switched to linking objects
     */
    /*
    void MasterContactFace3D::AddContactStress( Vector coords,
            const double Stress, const double Weight, const double dA, const Vector NormalDirection, double Penalty, double Gap, EquationIdVectorType SlaveEquationId, Vector SlaveNcontainer )
    {
        StressConditionType condition;
        condition.NormalDirection = NormalDirection;
        condition.Stress = Stress;
        condition.Weight = Weight;
        condition.dA = dA;
        condition.coords = coords;
        condition.Penalty = Penalty;
        condition.Gap = Gap;
        condition.SlaveEquationId = SlaveEquationId;
        condition.SlaveNcontainer = SlaveNcontainer;
        mpStressConditions->push_back( condition );
    }
    */
            
           
    
    //************************************************************************************
    //************************************************************************************
    
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void MasterContactFace3D::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
    {/*
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType matrix = Matrix();
        CalculateAll( matrix, rRightHandSideVector, 
                      rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, 
                      CalculateResidualVectorFlag);*/
    }
    
    //************************************************************************************
    //************************************************************************************
    
    /**
     * calculates the cross-diagonal terms of the system stiffness matrix
     * resulting from interaction between contact elements
     * 
     * Does nothing as assembling is to be switched to linking objects
     * (REMOVED)
     */
//     void MasterContactFace3D::CalculateCrossElementarySystemContributions( SecondaryConditionContainerType& SecondaryConditions, EquationIdVectorType& PrimaryEquationId, ProcessInfo& rCurrentProcessInfo )
//     {
        /*
        KRATOS_TRY

        
        SecondaryConditions.clear();
        EquationIdVector( PrimaryEquationId, rCurrentProcessInfo );
        
        //manually setting dimension to 3:
        unsigned int dim = 3;
        //loop over all stress conditions
        unsigned int PointNumber = 0;
        for( StressConditionContainerType::ptr_iterator it = mpStressConditions->ptr_begin();
             it != mpStressConditions->ptr_end(); ++it )
        {
            StressConditionType current = **it;
            //storing master element's equation ID vector
            SecondaryCondition CrossCondition;
            CrossCondition.SecondaryEquationIdVector = current.SlaveEquationId;
            //LHS sizes
            unsigned int PrimarySize = dim*GetGeometry().size();
            unsigned int SecondarySize = dim*current.SlaveNcontainer.size();
            //setting up LHS contribution
            LHS_ContributionType CurrentLHS = ZeroMatrix( PrimarySize, SecondarySize );
            
            Vector PrimaryN = ZeroVector( GetGeometry().size() );
            for( int NodeNumber = 0; NodeNumber < PrimaryN.size(); NodeNumber++ )
            {
                GeometryType::PointType currentPoint(3);
                currentPoint.X() = current.coords[0];
                currentPoint.Y() = current.coords[1];
                currentPoint.Z() = current.coords[2];
                PrimaryN[NodeNumber] = GetGeometry().ShapeFunctionValue(
                        NodeNumber, currentPoint );
            }
                    
            //adding cross elemental contribution
            for( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
            {
                for( int i=0; i<dim; i++ )
                {
                    for( unsigned int sec = 0; sec < current.SlaveNcontainer.size(); sec++ )
                    {
                        for( unsigned int j = 0; j<dim; j++ )
                        {
                            CurrentLHS[prim*dim+i][sec*dim+j] -= PrimaryN[prim]*current.SlaveNcontainer[sec]
                                    *current.Weight*current.dA*current.Penalty*current.NormalDirection[i]
                                    *current.NormalDirection[j];
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
//     }
    //************************************************************************************
    //************************************************************************************
    /**
     * Calculation of stiffness matrix contribution due to linearization of normal vector.
     * In this mehtod the components referring to the master element as primary index are
     * calculated
     * Hence there is a loop over all registered slave contact points serving as secondary indices.
     * The result is stored in a contribution container returned to the builder
     * @param SecondaryConditions the contribution container the resulting contributions are stored in
     * @param PrimaryEquationId the equationID vector of the primary element (here: the current master
     * element)
     * @param rCurrentProcessInfo the current system's process info
     * 
     * Does nothing as assembling is to be switched to linking objects
     * (REMOVED)
     */
//     void MasterContactFace3D::CalculateNormalLinearizationElementarySystemContributions(
//             SecondaryConditionContainerType& SecondaryConditions, 
//             EquationIdVectorType& PrimaryEquationId, ProcessInfo& rCurrentProcessInfo )
//     {
        /*
        SecondaryConditions.clear();
        EquationIdVector( PrimaryEquationId, rCurrentProcessInfo );
        
        //loop over all applied contact conditions
        for( StressConditionContainerType::ptr_iterator it= mpStressConditions->ptr_begin();
             it != mpStressConditions->ptr_end(); it++ )
        {
            StressConditionType current = **it;
            
            //setting up variables
            EquationIdVectorType PrimaryEquationId;
            EquationIdVector( PrimaryEquationId, rCurrentProcessInfo );
            EquationIdVectorType SecondaryEquationId = current.SlaveEquationId;
            
            //auxiliary information
            const unsigned int dim = 3;
            
            GeometryType::PointType currentPoint(3);
            currentPoint.X() = current.coords[0];
            currentPoint.Y() = current.coords[1];
            currentPoint.Z() = current.coords[2];
            
            double aStress = current.Stress*current.dA*current.Weight;
            
            //first order derivatives of shape function
            Matrix DN = ZeroMatrix(GetGeometry().PointsNumber(),2);
            GetGeometry().ShapeFunctionsGradients( DN, currentPoint );
            
            //second order derivatives of shape functions
            GeometryType::ShapeFunctionsGradientsType D2N; 
            GetGeometry().SecondOrderShapeFunctionDerivatives( D2N, currentPoint );
            
            double Gap = current.Gap;
            if( Gap != 0.0 )
            {
                Vector normalVector = current.NormalDirection;
                Matrix T = TangentialVectors( currentPoint );
            
            //LHS sizes
                unsigned int PrimarySize = dim*GetGeometry().size();
                unsigned int SecondarySize = dim*current.SlaveNcontainer.size();
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
                        for( int n=0; n<GetGeometry().size(); n++ )
                        {
                            Phi[i][j][0] += GetGeometry().GetPoint(n).X()*D2N[n][i][j];
                            Phi[i][j][1] += GetGeometry().GetPoint(n).Y()*D2N[n][i][j];
                            Phi[i][j][2] += GetGeometry().GetPoint(n).Z()*D2N[n][i][j];
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
                
                //calculating contributions Master-Slave 
                CurrentLHS = ZeroMatrix( PrimarySize, SecondarySize );
                //loop over all nodes on primary element 
                for( int prim = 0; prim < GetGeometry().size(); prim++ )
                {
                    double N = GetGeometry().ShapeFunctionValue( 
                                prim, currentPoint );
                    //loop over all dimensions 
                    for( int i=0; i<dim; i++ )
                    {
                        
                        double PsiPrim1 = -N*T[0][i]-Gap*normalVector[i]*DN[prim][0];
                        double PsiPrim2 = -N*T[1][i]-Gap*normalVector[i]*DN[prim][1];
                        //loop over all nodes of secondary element
                        for( int sec = 0; sec < current.SlaveNcontainer.size(); sec++ )
                        {
                            //loop over all dimensions
                            for( int j=0; j<dim; j++ )
                            {
                                double PsiSec1 = current.SlaveNcontainer[sec]*T[0][j];
                                double PsiSec2 = current.SlaveNcontainer[sec]*T[1][j];
                                //loop over gamma
                                for( int gamma = 0; gamma < 2; gamma++ )
                                {
                                    //loop over beta
                                    for( int beta=0; beta < 2; beta++ )
                                    {
                                        CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                                *(Gap*((Kappa[0][gamma]/A[0][0]+Kappa[1][gamma]/A[0][1])
                                                *PsiSec1+(Kappa[0][gamma]/A[1][0]+Kappa[1][gamma]/A[1][1])
                                                *PsiPrim2)*m[gamma][beta]*(normalVector[j]*DN[prim][beta]
                                                +(Kappa[0][beta]/A[0][0]+Kappa[1][beta]/A[0][1])*PsiPrim1
                                                +(Kappa[0][beta]/A[1][0]+Kappa[1][beta]/A[1][1])*PsiPrim2));
                                        CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                                *(Kappa[gamma][beta]*(1.0/A[0][beta]*PsiSec1+1.0/A[1][beta]*PsiSec2)
                                                *(1.0/A[0][gamma]*PsiPrim1+1.0/A[1][gamma]*PsiPrim2));
                                    }//loop over beta
                                    CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                            *((1.0/A[0][gamma]*PsiSec1+1.0/A[1][gamma]*PsiSec2)*normalVector[i]
                                            *DN[prim][gamma]);
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
                
                //calculating contributions Master-Master 
                CurrentLHS = ZeroMatrix( PrimarySize, PrimarySize );
                //loop over all nodes on primary element 
                for( int prim = 0; prim < GetGeometry().size(); prim++ )
                {
                    double N = GetGeometry().ShapeFunctionValue( 
                                prim, currentPoint );
                    //loop over all dimensions 
                    for( int i=0; i<dim; i++ )
                    {
                        
                        double PsiPrim1 = -N*T[0][i]-Gap*normalVector[i]*DN[prim][0];
                        double PsiPrim2 = -N*T[1][i]-Gap*normalVector[i]*DN[prim][1];
                        //loop over all nodes of secondary element
                        for( int sec = 0; sec < GetGeometry().size(); sec++ )
                        {
                            double Nsec = GetGeometry().ShapeFunctionValue(
                                        sec, currentPoint );
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
                                                *(Gap*(normalVector[j]*DN[sec][gamma]
                                                +(Kappa[0][gamma]/A[0][0]+Kappa[1][gamma]/A[0][1])
                                                *PsiSec1+(Kappa[0][gamma]/A[1][0]+Kappa[1][gamma]/A[1][1])
                                                *PsiSec2)*m[gamma][beta]*(normalVector[j]*DN[prim][beta]
                                                +(Kappa[0][beta]/A[0][0]+Kappa[1][beta]/A[0][1])*PsiPrim1
                                                +(Kappa[0][beta]/A[1][0]+Kappa[1][beta]/A[1][1])*PsiPrim2));
                                        CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                                *(Kappa[gamma][beta]*(1.0/A[0][beta]*PsiSec1+1.0/A[1][beta]*PsiSec2)
                                                *(1.0/A[0][gamma]*PsiPrim1+1.0/A[1][gamma]*PsiPrim2));
                                                
                                    }//loop over beta
                                    CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                            *((1.0/A[0][gamma]*PsiSec1+1.0/A[1][gamma]*PsiSec2)
                                            *normalVector[i]*DN[prim][gamma]);
                                    CurrentLHS[prim*dim+i][sec*dim+j] += aStress
                                            *((1.0/A[0][gamma]*PsiPrim1+1.0/A[1][gamma]*PsiPrim2)
                                            *normalVector[j]*DN[sec][gamma]);
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
                    
        }//loop over all applied conditions
        
        */
                    
//     }//CalculateNormalLinearizationElementarySystemContributions
    
    
    //************************************************************************************
    //************************************************************************************
    
    
    /**
     * calculates this contact element's local contributions
     */
    void MasterContactFace3D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
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
     * calculates the contact related contributions to the system
     * Does nothing as assembling is to be switched to linking objects
     */
    void MasterContactFace3D::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                      VectorType& rRightHandSideVector, 
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
        /*
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        //unsigned int dim = GetGeometry().WorkingSpaceDimension();
        //manually setting dimension to 3:
        unsigned int dim = 3;
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
                        
        //calculating current jacobian
        GeometryType::JacobiansType J;
        J = GetGeometry().Jacobian(J);
        //calculating determinant of jacobian
        Vector<double> DetJ( integration_points.size() );
        GetGeometry().DeterminantOfJacobian( DetJ );
        
        //loop over all applied contact conditions
        if( mpStressConditions->size() != 0 )
        {
            for( StressConditionContainerType::ptr_iterator it= mpStressConditions->ptr_begin();
                 it != mpStressConditions->ptr_end(); it++ )
            {
                double Stress = (*it)->Stress;
//                 std::cout << "Stress on corresponding master point: " << Stress << std::endl;
                GeometryType::PointType currentPoint(3);
                currentPoint.X() = (*it)->coords[0];
                currentPoint.Y() = (*it)->coords[1];
                currentPoint.Z() = (*it)->coords[2];
                Vector Ncontainer = ZeroVector( GetGeometry().PointsNumber() );
                for( int n=0; n<GetGeometry().PointsNumber(); n++ )
                {
                    Ncontainer[n] = GetGeometry().ShapeFunctionValue( n, currentPoint );
                }
                if( CalculateResidualVectorFlag == true )
                {
                    CalculateAndAdd_PressureForce( rRightHandSideVector,
                                               Ncontainer,
                                               (*it)->NormalDirection, (*it)->Stress, (*it)->Weight, (*it)->dA );
                }
                
                if( CalculateStiffnessMatrixFlag == true )
                {
                    //linearization of gap
                    CalculateAndAddKc( rLeftHandSideMatrix,
                                       Ncontainer,
                                       (*it)->Weight,
                                       (*it)->dA,
                                       (*it)->Penalty,
                                       (*it)->NormalDirection
                                     );
                }
            }
        }
        mpStressConditions->clear();
        
        KRATOS_CATCH("")
        */
    } // CalculateAll
    
    //***********************************************************************
    //***********************************************************************
    /**
     * System matrix contribution due to contact energy
     * TO BE IMPLEMENTED
     */
    void MasterContactFace3D::CalculateAndAddKc( Matrix& K,
            const Vector& N,
            double weight,
            double dA,
            double penalty,
            Vector v
                                               )
    {
		/*
       KRATOS_TRY
        //adding pure elemental contribution
        for( unsigned int prim=0; prim < GetGeometry().size(); prim++ )
        {
            {
                for( unsigned int i=0; i<3; i++ )
                {
                    for( unsigned int sec=0; sec < GetGeometry().size(); sec++ )
                    {
                        for( unsigned int j=0; j<3; j++ )
                        {
                            K(prim*3+i,sec*3+j) += N[prim]*N[sec]*v[i]*v[j]*weight*dA*penalty;
                        }
                    }
                }
            }
        }
        KRATOS_CATCH("")
		*/
    }
    
    //***********************************************************************
    //***********************************************************************
    /**
     *
     */
    void MasterContactFace3D::CalculateAndAdd_PressureForce( Vector& residualvector,
            const Vector& N,
            Vector& v3,
            double pressure,
            double weight, double dA 
                                                     )
    {/*
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = 3;
        for (unsigned int i=0;i<number_of_nodes;i++)
        {
            int index = dim*i;
            double coeff = -pressure * N[i] * weight * dA;
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
   void MasterContactFace3D::EquationIdVector( EquationIdVectorType& rResult, 
                                         ProcessInfo& CurrentProcessInfo 
                                       )
   {
	   rResult.resize(0);
//        int number_of_nodes = GetGeometry().PointsNumber();
//        unsigned int index;
//        unsigned int dim = 3;
//        rResult.resize(number_of_nodes*dim);
//        for (int i=0;i<number_of_nodes;i++)
//        {
//            index = i*dim;
//            rResult[index]   = (GetGeometry()[i].GetDof(DISPLACEMENT_X)).EquationId();
//            rResult[index+1] = (GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());
//            rResult[index+2] = (GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId());
//        }
   }
       
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
   void MasterContactFace3D::GetDofList( DofsVectorType& ConditionalDofList,
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
