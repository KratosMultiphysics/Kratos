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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-20 08:57:46 $
//   Revision:            $Revision: 1.4 $
//
//
// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/contact_link_3D_lagrange_tying.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    Contact_Link_3D_Lagrange_Tying::Contact_Link_3D_Lagrange_Tying( IndexType NewId,
                                                                   GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    Contact_Link_3D_Lagrange_Tying::Contact_Link_3D_Lagrange_Tying( IndexType NewId, GeometryType::Pointer pGeometry,
                                                                   PropertiesType::Pointer pProperties
                                                                   )
    : Condition( NewId, pGeometry, pProperties )
    {
    }
    
    /**
     * Constructor of the contact link 3D, creates a link between a quadrature point on the
     * slave surface and its closest point projection on the master surface
     * @param pGeometry the link does not have a special geometry
     * @param pProperties
     * @param Master MasterContactFace3D the closest point projection is located on
     * @param Slave SlaveContactFace3D the quadrature point is located on
     * @param MasterContactLocalPoint local coordinates in Master of the closest point projection
     * @param SlaveContactLocalPoint local coordinates in Slave of the quadrature point
     * @param SlaveIntegrationPointIndex integration point index of the quadrature point
     */
    Contact_Link_3D_Lagrange_Tying::Contact_Link_3D_Lagrange_Tying( IndexType NewId, GeometryType::Pointer pGeometry,
                                                                   PropertiesType::Pointer pProperties,
                                                                   Condition::Pointer Master,
                                                                   Condition::Pointer Slave,
                                                                  Point& MasterContactLocalPoint,
                                                                  Point& SlaveContactLocalPoint,
                                                                   int SlaveIntegrationPointIndex
                                                                   )
    : Condition( NewId, pGeometry, pProperties )
    {
        GetValue( CONTACT_LINK_MASTER ) = Master;
        GetValue( CONTACT_LINK_SLAVE ) = Slave;
        //
        GetValue( MASTER_CONTACT_LOCAL_POINT ) = MasterContactLocalPoint;
        GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT ) = MasterContactLocalPoint;
        GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT ) = MasterContactLocalPoint;
        //
        GetValue( SLAVE_CONTACT_LOCAL_POINT ) = SlaveContactLocalPoint;
        //Test for calculating coordinates at time step midpoint
        GetValue( MASTER_CONTACT_GLOBAL_POINT ) = GlobalCoordinates( GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_GLOBAL_POINT ), GetValue( MASTER_CONTACT_LOCAL_POINT ) );
        //Test for calculating coordinates at time step midpoint
        GetValue( SLAVE_CONTACT_GLOBAL_POINT ) = GlobalCoordinates( GetValue( CONTACT_LINK_SLAVE ), GetValue( SLAVE_CONTACT_GLOBAL_POINT ), GetValue( SLAVE_CONTACT_LOCAL_POINT ) );
        //
        GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ) = SlaveIntegrationPointIndex;
        
        GetValue( CONTACT_LINK_M ).resize( 2, 2, false );
        noalias( GetValue( CONTACT_LINK_M ) ) = ZeroMatrix( 2, 2 );
        
        mvMaster.resize( 3 );
        mTMaster.resize( 2, 3 );
    }
    
    //********************************************************
    //**** Operations ****************************************
    //********************************************************
    
    
    Condition::Pointer Contact_Link_3D_Lagrange_Tying::Create( IndexType NewId,
                                                              NodesArrayType const& ThisNodes,
                                                              PropertiesType::Pointer pProperties ) const
    {
        return Condition::Pointer( new Contact_Link_3D_Lagrange_Tying( NewId, GetGeometry().Create( ThisNodes ),
                                                                      pProperties ) );
    }
    
    /**
     * Destructor. Never to be called manually
     */
    Contact_Link_3D_Lagrange_Tying::~Contact_Link_3D_Lagrange_Tying()
    {
    }
    
    
    void Contact_Link_3D_Lagrange_Tying::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        
        noalias( mvMaster ) =  NormalVector( GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_LOCAL_POINT ) );
        
        noalias( mTMaster ) = TangentialVectors( GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_LOCAL_POINT ) );
    }
    
    /**
     * returns the tangential vectors of the current surface in an arbitrary point due to current configuration
     * @param Surface Surface Condition for that the Tangential Vector should be calculated
     * @param rPoint local coordinates of the point in Surface for that the Tangential Vector should be calculated
     * @return Matrix(2,3) of the two tangential vectors
     */
    Matrix Contact_Link_3D_Lagrange_Tying::TangentialVectors( Condition::Pointer Surface,
                                                             const GeometryType::CoordinatesArrayType& rPoint )
    {
        //setting up result matrix
        Matrix T( 2, 3 );
        noalias( T ) = ZeroMatrix( 2, 3 );
        //shape function gradients
        Matrix DN = ZeroMatrix( Surface->GetGeometry().PointsNumber(), 2 );
        Surface->GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
        //calculating tangential vectors
        
        for ( unsigned int n = 0; n < Surface->GetGeometry().PointsNumber(); n++ )
        {
            T( 0, 0 ) += ( Surface->GetGeometry().GetPoint( n ).X0()
                          + Surface->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
            * DN( n, 0 );
            T( 0, 1 ) += ( Surface->GetGeometry().GetPoint( n ).Y0()
                          + Surface->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
            * DN( n, 0 );
            T( 0, 2 ) += ( Surface->GetGeometry().GetPoint( n ).Z0()
                          + Surface->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
            * DN( n, 0 );
            T( 1, 0 ) += ( Surface->GetGeometry().GetPoint( n ).X0()
                          + Surface->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
            * DN( n, 1 );
            T( 1, 1 ) += ( Surface->GetGeometry().GetPoint( n ).Y0()
                          + Surface->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
            * DN( n, 1 );
            T( 1, 2 ) += ( Surface->GetGeometry().GetPoint( n ).Z()
                          + Surface->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
            * DN( n, 1 );
        }
        
        return( T );
    }
    
    /**
     * returns the tangential vectors of the current surface in an arbitrary point due to reference configuration
     * @param Surface Surface Condition for that the Tangential Vector should be calculated
     * @param rPoint local coordinates of the point in Surface for that the Tangential Vector should be calculated
     * @return Matrix(2,3) of the two tangential vectors
     */
    Matrix Contact_Link_3D_Lagrange_Tying::TangentialVectors_inOrigin( Condition::Pointer Surface,
                                                                      const GeometryType::CoordinatesArrayType& rPoint )
    {
        //setting up result matrix
        Matrix T( 2, 3 );
        noalias( T ) = ZeroMatrix( 2, 3 );
        //shape function gradients
        Matrix DN = ZeroMatrix( Surface->GetGeometry().PointsNumber(), 2 );
        Surface->GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
        //calculating tangential vectors
        
        for ( unsigned int n = 0; n < Surface->GetGeometry().PointsNumber(); n++ )
        {
            T( 0, 0 ) += Surface->GetGeometry().GetPoint( n ).X0() * DN( n, 0 );
            T( 0, 1 ) += Surface->GetGeometry().GetPoint( n ).Y0() * DN( n, 0 );
            T( 0, 2 ) += Surface->GetGeometry().GetPoint( n ).Z0() * DN( n, 0 );
            T( 1, 0 ) += Surface->GetGeometry().GetPoint( n ).X0() * DN( n, 1 );
            T( 1, 1 ) += Surface->GetGeometry().GetPoint( n ).Y0() * DN( n, 1 );
            T( 1, 2 ) += Surface->GetGeometry().GetPoint( n ).Z0() * DN( n, 1 );
        }
        
        return( T );
    }
    
    /**
     * returns the normalized tangential vectors of the current surface in an arbitrary point due to current configuration
     * @param Surface Surface Condition for that the Tangential Vector should be calculated
     * @param rPoint local coordinates of the point in Surface for that the Tangential Vector should be calculated
     * @return Matrix(2,3) of the two tangential vectors
     */
    Matrix Contact_Link_3D_Lagrange_Tying::TangentialVectorsGlobal( Condition::Pointer Surface,
                                                                   const GeometryType::CoordinatesArrayType& rPoint )
    {
        //setting up result matrix
        Matrix T( 2, 3 );
        noalias( T ) = TangentialVectors( Surface, rPoint );
        
        //shape function gradients
        double normT1 = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );
        double normT2 = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );
        T( 0, 0 ) = T( 0, 0 ) / normT1;
        T( 0, 1 ) = T( 0, 0 ) / normT1;
        T( 0, 2 ) = T( 0, 0 ) / normT1;
        T( 1, 0 ) = T( 0, 0 ) / normT2;
        T( 1, 1 ) = T( 0, 0 ) / normT2;
        T( 1, 2 ) = T( 0, 0 ) / normT2;
        return( T );
    }
    
    /**
     * returns the normal vector in arbitrary point
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector
     */
    Vector Contact_Link_3D_Lagrange_Tying::NormalVector( Condition::Pointer Surface,
                                                        const GeometryType::CoordinatesArrayType& rPoint )
    {
        Vector Result( 3 );
        noalias( Result ) = ZeroVector( 3 );
        //getting tangential vectors
        Matrix T( 2, 3 );
        noalias( T ) = TangentialVectors( Surface, rPoint );
        //calculating normal vector
        Result[0] = T( 0, 1 ) * T( 1, 2 ) - T( 0, 2 ) * T( 1, 1 );
        Result[1] = T( 0, 2 ) * T( 1, 0 ) - T( 0, 0 ) * T( 1, 2 );
        Result[2] = T( 0, 0 ) * T( 1, 1 ) - T( 0, 1 ) * T( 1, 0 );
        
        SD_MathUtils<double>::Normalize( Result );
        
        return( Result );
    }
    
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void Contact_Link_3D_Lagrange_Tying::CalculateRightHandSide( VectorType& rRightHandSideVector,
                                                                ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType matrix = Matrix();
        CalculateAll( matrix, rRightHandSideVector,
                     rCurrentProcessInfo,
                     CalculateStiffnessMatrixFlag,
                     CalculateResidualVectorFlag );
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * calculates this contact element's local contributions
     */
    void Contact_Link_3D_Lagrange_Tying::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                                              VectorType& rRightHandSideVector,
                                                              ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                     CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }
    
    /**
     * This function calculates all system contributions due to the contact problem
     * with regard to the current master and slave partners.
     * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node
     */
    void Contact_Link_3D_Lagrange_Tying::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo,
                                                      bool CalculateStiffnessMatrixFlag,
                                                      bool CalculateResidualVectorFlag )
    {
        KRATOS_TRY
        
        //**********************************************
        //setting up the dimensions of the contributions
        //**********************************************
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        unsigned int dim = 3;
        
        //resizing as needed the LHS
        int MatSize = ( MasterNN + SlaveNN + 1 ) * dim;
        
        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
            noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
        }
        
        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            rRightHandSideVector.resize( MatSize, false );
            noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
        }
        
        //fixed normal vector on master element
        Vector vMaster( 3 );
        
        noalias( vMaster ) = mvMaster;
        
        //calculating and updating gap
        double Gap = inner_prod(
                                vMaster, ( GetValue( MASTER_CONTACT_GLOBAL_POINT ) - GetValue( SLAVE_CONTACT_GLOBAL_POINT ) ) );
        
        GetValue( CONTACT_LINK_SLAVE )->GetValue( GAPS )[GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )] = Gap;
        
        if( Gap > 1.0e-10 )
        {
            //*******************************************************************************
            //****************** adding all contributions ***********************************
            //*******************************************************************************
        
            //contribution to RHS
            //subtracting relDisp from reaction vector
        
//             for ( unsigned int node = 0; node < MasterNN; node++ )
//             {
//                rRightHandSideVector[dim*node] -= vMaster[0]*Gap;
//                rRightHandSideVector[dim*node+1] -= vMaster[1]*Gap;
//                rRightHandSideVector[dim*node+2] -= vMaster[2]*Gap;
//             }
        
            //for ( unsigned int node = 0; node < SlaveNN; node++ )
            //{
            //    rRightHandSideVector[MasterNN*dim+dim*node] += vMaster[0]*Gap;
            //    rRightHandSideVector[MasterNN*dim+dim*node+1] += vMaster[1]*Gap;
            //    rRightHandSideVector[MasterNN*dim+dim*node+2] += vMaster[2]*Gap;
            //}
            
            for ( unsigned int i = (MasterNN+SlaveNN)*dim; i < (MasterNN+SlaveNN+1)*dim; i++ )
            {
               rRightHandSideVector[i] -= vMaster[0]*Gap;
            }
        

            //contribution to LHS
            for ( unsigned int i = 0; i < dim; i++ )
            {
                for ( unsigned int slave_node = 0; slave_node < SlaveNN; slave_node++ )
                {
                    for ( unsigned int master_node = 0; master_node < MasterNN; master_node++ )
                    {
                        rLeftHandSideMatrix( master_node * dim + i, MasterNN * dim + SlaveNN * dim + i ) -= vMaster[i] * GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionValue( master_node, GetValue( MASTER_CONTACT_LOCAL_POINT ) );
                        rLeftHandSideMatrix( MasterNN * dim + SlaveNN * dim + i, master_node * dim + i ) -= vMaster[i] * GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionValue( master_node, GetValue( MASTER_CONTACT_LOCAL_POINT ) );
                        rLeftHandSideMatrix( MasterNN * dim + slave_node * dim + i, MasterNN * dim + SlaveNN * dim + i ) += vMaster[i] * GetValue( CONTACT_LINK_SLAVE )->GetGeometry().ShapeFunctionValue( slave_node, GetValue( SLAVE_CONTACT_LOCAL_POINT ) );
                        rLeftHandSideMatrix( MasterNN * dim + SlaveNN * dim + i, MasterNN * dim + slave_node * dim + i ) += vMaster[i] * GetValue( CONTACT_LINK_SLAVE )->GetGeometry().ShapeFunctionValue( slave_node, GetValue( SLAVE_CONTACT_LOCAL_POINT ) );
                    }
                }
            }
        }
        KRATOS_CATCH( "" )
    } // CalculateAll
    
    /**
     * This function calculates the system contributions to the global damp matrix due to the contact problem
     * with regard to the current master and slave partners.
     * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node
     */
    void Contact_Link_3D_Lagrange_Tying::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        
        KRATOS_TRY
        
        KRATOS_CATCH( "" )
    }
    
    //***********************************************************************
    //***********************************************************************
    /**
     * Calculates the residual contribution due to contact stresses on
     * both the master and slave conditions. The contributions are stored
     * into the residual vector.
     * @param residualvector the right hand side vector of the current
     * linking condition
     * @param NMaster the shape function values of the master element in
     * the current contact point
     * @param NSlave the shape function values of the slave element in
     * the current contact point
     * @param vMaster the normal vector to the master surface in
     * current contact point
     * @param normalStress the value of contact stress in current contact
     * point
     * @param SlaveIntegrationWeight the integration weight in current slave
     * integration point
     * @param dASlave the differential area element of the current slave
     * integration point
     */
    void Contact_Link_3D_Lagrange_Tying::CalculateAndAdd_RHS( Vector& residualvector,
                                                             const Vector& NMaster,
                                                             const Vector& NSlave,
                                                             const Vector& vMaster,
                                                             const Matrix& T,
                                                             const Vector& tangentialStresses,
                                                             double Gap,
                                                             double normalStress,
                                                             double SlaveIntegrationWeight,
                                                             double dASlave
                                                             )
    {
        //**********************************************
        //BEGIN OF MASTER: CalculateAndAdd_PressureForce
        //**********************************************
        
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int dim = 3;
        
        for ( unsigned int i = 0; i < MasterNN; i++ )
        {
            int index = dim * i;
            double coeff = -normalStress * NMaster[i] * SlaveIntegrationWeight * dASlave;
            residualvector[index]   = coeff * vMaster[0];
            residualvector[index+1] = coeff * vMaster[1];
            residualvector[index+2] = coeff * vMaster[2];
        }
        
        //********************************************
        //END OF MASTER: CalculateAndAdd_PressureForce
        //********************************************
        
        //*********************************************
        //BEGIN OF SLAVE: CalculateAndAdd_PressureForce
        //*********************************************
        
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        
        for ( unsigned int i = 0; i < SlaveNN; i++ )
        {
            int index = MasterNN * dim + dim * i;
            double coeff = normalStress * NSlave[i] * SlaveIntegrationWeight * dASlave;
            residualvector[index]   = coeff * vMaster[0];
            residualvector[index+1] = coeff * vMaster[1];
            residualvector[index+2] = coeff * vMaster[2];
        }
        
        //*******************************************
        //END OF SLAVE: CalculateAndAdd_PressureForce
        //*******************************************
        //         KRATOS_WATCH( residualvector );
        
        
        //*******************************************
        //CalculateAndAdd RHS-Vector due to frictional contact
        //*******************************************
        
        
        //**********************************************
        //BEGIN OF MASTER: CalculateAndAdd_PressureForce
        //**********************************************
        
        if ( GetValue( CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT] > 0.0 )
        {
            
            Vector norm_T( 2 );
            
            norm_T( 0 ) = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );
            norm_T( 1 ) = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );
            
            for ( unsigned int i = 0; i < MasterNN; i++ )
            {
                for ( unsigned int j = 0 ; j < dim ; j++ )
                {
                    Vector Xi = ZeroVector( 2 );
                    
                    Xi[0] = -NMaster[i] * T( 0, j ) / norm_T( 0 );
                    
                    Xi[1] = -NMaster[i] * T( 1, j ) / norm_T( 1 );
                    
                    residualvector[ i*dim+j] -= ( tangentialStresses[0] * Xi[0] +
                                                 tangentialStresses[1] * Xi[1] ) * SlaveIntegrationWeight * dASlave;
                }
            }
            
            //         **********************************************
            //         END OF MASTER: CalculateAndAdd_PressureForce
            //         **********************************************
            
            //         **********************************************
            //         BEGIN OF SLAVE: CalculateAndAdd_PressureForce
            //         **********************************************
            
            
            for ( unsigned int i = 0; i < SlaveNN; i++ )
            {
                for ( unsigned int j = 0; j < dim; j++ )
                {
                    Vector Xi = ZeroVector( 2 );
                    
                    Xi[0] = NSlave[i] * T( 0, j ) / norm_T( 0 );
                    
                    Xi[1] = NSlave[i] * T( 1, j ) / norm_T( 1 );
                    
                    residualvector[MasterNN*dim +i*dim+j] -= ( tangentialStresses[0] * Xi[0] +
                                                              tangentialStresses[1] * Xi[1] ) * SlaveIntegrationWeight * dASlave;
                }
            }
            
            //         **********************************************
            //         END OF SLAVE: CalculateAndAdd_PressureForce
            //         **********************************************
            
        }
        
    }
    
    /**
     * This function calculates updates the local and global coordinates
     * of the master contact partner in order to follow the movement of
     * the slave surface along the master surface
     */
    void Contact_Link_3D_Lagrange_Tying::UpdateMasterLocalPoint( )
    {
        double Xi1 = GetValue( MASTER_CONTACT_LOCAL_POINT )[0];
        double Xi2 = GetValue( MASTER_CONTACT_LOCAL_POINT )[1];
        double deltaXi1 = 0.0;
        double deltaXi2 = 0.0;
        
        for ( int k = 0; k < 1000; k++ )
        {
            // KRATOS_WATCH( GetValue( MASTER_CONTACT_GLOBAL_POINT ) );
            //setting up tangential vectors
            Vector t1 = ZeroVector( 3 );//first tangential vector
            Vector t2 = ZeroVector( 3 );//second tangential vector
            //derivatives of tangential vectors
            Vector dt11 = ZeroVector( 3 );
            Vector dt12 = ZeroVector( 3 );
            Vector dt21 = ZeroVector( 3 );
            Vector dt22 = ZeroVector( 3 );
            
            //retrieving first order derivatives in current solution point
            Matrix DN = ZeroMatrix( GetValue( CONTACT_LINK_MASTER )->GetGeometry().PointsNumber(), 2 );
            GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionsLocalGradients( DN, GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT ) );
            //retrieving second order derivatives in current solution point
            GeometryType::ShapeFunctionsSecondDerivativesType D2N;
            GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionsSecondDerivatives( D2N, GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT ) );
            
            for ( unsigned  int n = 0; n < GetValue( CONTACT_LINK_MASTER )->GetGeometry().PointsNumber(); n++ )
            {
                //contribution to tangential vectors
                t1[0] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).X0()
                          + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                * DN( n, 0 );
                t1[1] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Y0()
                          + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                * DN( n, 0 );
                t1[2] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Z0()
                          + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                * DN( n, 0 );
                t2[0] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).X0()
                          + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                * DN( n, 1 );
                t2[1] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Y0()
                          + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                * DN( n, 1 );
                t2[2] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Z0()
                          + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                * DN( n, 1 );
                //contribution to derivatives of tangential vectors
                dt11[0] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).X0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                * D2N[n]( 0, 0 );
                dt11[1] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Y0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                * D2N[n]( 0, 0 );
                dt11[2] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Z0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                * D2N[n]( 0, 0 );
                dt12[0] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).X0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                * D2N[n]( 0, 1 );
                dt12[1] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Y0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                * D2N[n]( 0, 1 );
                dt12[2] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Z0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                * D2N[n]( 0, 1 );
                dt21[0] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).X0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                * D2N[n]( 1, 0 );
                dt21[1] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Y0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                * D2N[n]( 1, 0 );
                dt21[2] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Z0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                * D2N[n]( 1, 0 );
                dt22[0] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).X0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) )
                * D2N[n]( 1, 1 );
                dt22[1] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Y0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) )
                * D2N[n]( 1, 1 );
                dt22[2] += ( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).Z0()
                            + GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) )
                * D2N[n]( 1, 1 );
            }
            
            //defining auxiliary terms
            double A1 = (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0] ) * t1[0] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1] ) * t1[1] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2] ) * t1[2] );
            
            double A2 = (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0] ) * t2[0] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1] ) * t2[1] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2] ) * t2[2] );
            
            double B11 = ( -t1[0] * t1[0] - t1[1] * t1[1] - t1[2] * t1[2] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0] ) * dt11[0] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1] ) * dt11[1] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2] ) * dt11[2] );
            
            double B12 = ( -t2[0] * t1[0] - t2[1] * t1[1] - t2[2] * t1[2] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0] ) * dt12[0] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1] ) * dt12[1] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2] ) * dt12[2] );
            
            double B21 = ( -t1[0] * t2[0] - t1[1] * t2[1] - t1[2] * t2[2] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0] ) * dt21[0] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1] ) * dt21[1] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2] ) * dt21[2] );
            
            double B22 = ( -t2[0] * t2[0] - t2[1] * t2[1] - t2[2] * t2[2] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0] ) * dt22[0] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1] ) * dt22[1] )
            + (( GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2] - GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2] ) * dt22[2] );
            
            //calculating update for Xi
            deltaXi1 = -A1 * B22 / ( B11 * B22 - B12 * B21 ) + A2 * B12 / ( B11 * B22 - B12 * B21 );
            
            deltaXi2 =  A2 * B21 / ( B11 * B22 - B12 * B21 ) - A2 * B11 / ( B11 * B22 - B12 * B21 );
            
            //updating Xi
            Xi1 += deltaXi1;
            
            Xi2 += deltaXi2;
            
            //updating LocalPoint
            GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT )[0] = Xi1;
            
            GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT )[1] = Xi2;
            
            //updating rResult
            GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT ) = ZeroVector( 3 );
            
            GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT ) = GlobalCoordinates( GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT ), GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT ) );
            
            if ( fabs( deltaXi1 ) < 1e-7 && fabs( deltaXi2 ) < 1e-7 )
            {
                return;
            }
        }
        
        std::cout << "******** ATTENTION: NO MAPPING TO MASTER SURFACE FOUND ************" << std::endl;
        
        return;
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * Setting up the EquationIdVector for the current partners.
     * All conditions are assumed to be defined in 3D space with 3 DOFs per node.
     * All Equation IDs are given Master first, Slave second
     */
    void Contact_Link_3D_Lagrange_Tying::EquationIdVector( EquationIdVectorType& rResult,
                                                          ProcessInfo& CurrentProcessInfo )
    {
        //determining size of DOF list
        //dimension of space
        unsigned int dim = 3;
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        
        unsigned int index;
        rResult.resize(( MasterNN + SlaveNN + 1 )*dim, false );
        
        for ( unsigned int i = 0; i < MasterNN; i++ )
        {
            index = i * dim;
            rResult[index]   = ( GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].GetDof( DISPLACEMENT_X ) ).EquationId();
            rResult[index+1] = ( GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId() );
            rResult[index+2] = ( GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId() );
        }
        
        for ( unsigned  int i = 0; i < SlaveNN; i++ )
        {
            index = MasterNN * dim + i * dim;
            rResult[index]   = ( GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].GetDof( DISPLACEMENT_X ) ).EquationId();
            rResult[index+1] = ( GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId() );
            rResult[index+2] = ( GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId() );
        }
        
        index = MasterNN * dim + SlaveNN * dim;
        
        rResult[index] = ( GetGeometry()[0].GetDof( LAGRANGE_DISPLACEMENT_X ) ).EquationId();
        rResult[index+1] = ( GetGeometry()[0].GetDof( LAGRANGE_DISPLACEMENT_Y ) ).EquationId();
        rResult[index+2] = ( GetGeometry()[0].GetDof( LAGRANGE_DISPLACEMENT_Z ) ).EquationId();
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * Setting up the DOF list for the current partners.
     * All conditions are assumed to be defined in 3D space with 3 DOFs per Node.
     * All DOF are given Master first, Slave second
     */
    void Contact_Link_3D_Lagrange_Tying::GetDofList( DofsVectorType& ConditionalDofList,
                                                    ProcessInfo& CurrentProcessInfo )
    {
        //determining size of DOF list
        //dimension of space
        unsigned int dim = 3;
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        
        ConditionalDofList.resize(( MasterNN + SlaveNN + 1 )*dim );
        unsigned int index;
        //setting up master DOFs
        
        for ( unsigned int i = 0; i < MasterNN; i++ )
        {
            index = i * dim;
            ConditionalDofList[index]   = ( GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ConditionalDofList[index+1] = ( GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            ConditionalDofList[index+2] = ( GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
        
        //setting up slave DOFs
        for ( unsigned int i = 0; i < SlaveNN; i++ )
        {
            index = MasterNN * dim + i * dim;
            ConditionalDofList[index]   = ( GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ConditionalDofList[index+1] = ( GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            ConditionalDofList[index+2] = ( GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
        
        index = MasterNN * dim + SlaveNN * dim;
        
        ConditionalDofList[index] = ( GetGeometry()[0].pGetDof( LAGRANGE_DISPLACEMENT_X ) );
        ConditionalDofList[index+1] = ( GetGeometry()[0].pGetDof( LAGRANGE_DISPLACEMENT_Y ) );
        ConditionalDofList[index+2] = ( GetGeometry()[0].pGetDof( LAGRANGE_DISPLACEMENT_Z ) );
    }
    
    //new functions includes
    
   Point& Contact_Link_3D_Lagrange_Tying::GlobalCoordinates( Condition::Pointer Surface,Point& rResult,Point const& LocalCoordinates )
    {
        noalias( rResult ) = ZeroVector( 3 );
        
        for ( IndexType i = 0 ; i < Surface->GetGeometry().size() ; i++ )
        {
            double shape_func = Surface->GetGeometry().ShapeFunctionValue( i, LocalCoordinates );
            
            rResult( 0 ) += shape_func *
            (( Surface->GetGeometry()[i] ).X0()
             + ( Surface->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_X ) );
            
            rResult( 1 ) += shape_func *
            (( Surface->GetGeometry()[i] ).Y0()
             + ( Surface->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_Y ) );
            
            rResult( 2 ) += shape_func *
            (( Surface->GetGeometry()[i] ).Z0()
             + ( Surface->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_Z ) );
        }
        
        return rResult;
    }
    
    /**
     * returns the relative tangential velocity between the quadrature point on the slave surface and its
     * closest point projection
     * @param T Matrix of the Tanegntial Vectors on the Master Surface in the current configuration
     * @return tangential velocity
     */
    Vector Contact_Link_3D_Lagrange_Tying::GetRelativTangentialVelocity( Matrix& T )
    {
        Vector result( 2 );
        
        Vector slave_velo( 3 );
        Vector master_velo( 3 );
        
        noalias( slave_velo ) = ZeroVector( 3 );
        noalias( master_velo ) = ZeroVector( 3 );
        
        for ( IndexType i = 0 ; i < GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size() ; i++ )
        {
            double shape_func = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().ShapeFunctionValue( i, GetValue( SLAVE_CONTACT_LOCAL_POINT ) );
            slave_velo += (( GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) ) *
            shape_func;
        }
        
        for ( IndexType i = 0 ; i < GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() ; i++ )
        {
            double shape_func = GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionValue( i, GetValue( MASTER_CONTACT_LOCAL_POINT ) );
            master_velo += (( GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) ) *
            shape_func;
        }
        
        Vector norm_T( 2 );
        
        norm_T( 0 ) = sqrt( T( 0, 0 ) * T( 0, 0 ) + T( 0, 1 ) * T( 0, 1 ) + T( 0, 2 ) * T( 0, 2 ) );
        
        norm_T( 1 ) = sqrt( T( 1, 0 ) * T( 1, 0 ) + T( 1, 1 ) * T( 1, 1 ) + T( 1, 2 ) * T( 1, 2 ) );
        
        result( 0 ) = (( slave_velo( 0 ) - master_velo( 0 ) ) * T( 0, 0 ) + ( slave_velo( 1 ) - master_velo( 1 ) ) * T( 0, 1 )
                       + ( slave_velo( 2 ) - master_velo( 2 ) ) * T( 0, 2 ) ) / norm_T( 0 );
        
        result( 1 ) = (( slave_velo( 0 ) - master_velo( 0 ) ) * T( 1, 0 ) + ( slave_velo( 1 ) - master_velo( 1 ) ) * T( 1, 1 )
                       + ( slave_velo( 2 ) - master_velo( 2 ) ) * T( 1, 2 ) ) / norm_T( 1 );
        
        return result;
    }
    
    /**
     * returns the relative velocity between the quadrature point on the slave surface and its
     * closest point projection
     * @return relative velocity
     */
    Vector Contact_Link_3D_Lagrange_Tying::GetRelativVelocity()
    {
        Vector result( 3 );
        
        Vector slave_velo( 3 );
        Vector master_velo( 3 );
        
        noalias( slave_velo ) = ZeroVector( 3 );
        noalias( master_velo ) = ZeroVector( 3 );
        
        for ( IndexType i = 0 ; i < GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size() ; i++ )
        {
            double shape_func = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().ShapeFunctionValue( i, GetValue( SLAVE_CONTACT_LOCAL_POINT ) );
            slave_velo += (( GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) ) *
            shape_func;
        }
        
        for ( IndexType i = 0 ; i < GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() ; i++ )
        {
            double shape_func = GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionValue( i, GetValue( MASTER_CONTACT_LOCAL_POINT ) );
            master_velo += (( GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) ) *
            shape_func;
        }
        
        result( 0 ) = ( slave_velo( 0 ) - master_velo( 0 ) );
        
        result( 1 ) = ( slave_velo( 1 ) - master_velo( 1 ) );
        
        result( 2 ) = ( slave_velo( 2 ) - master_velo( 2 ) );
        
        return result;
    }
    
    /// Print information about this object.
    void Contact_Link_3D_Lagrange_Tying::PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "Condition #" << Id();
    }
    
    /// Print object's data.
    void Contact_Link_3D_Lagrange_Tying::PrintData( std::ostream& rOStream ) const
    {
        rOStream << "Contact_Link_3D_Lagrange_Tying" << std::endl;
    }
} // Namespace Kratos
