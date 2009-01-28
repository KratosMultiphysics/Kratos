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
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-05-06 13:54:06 $
//   Revision:            $Revision: 1.6 $
//
//
// System includes 

// External includes 

// Project includes 
#include "includes/define.h"
#include "custom_conditions/contact_link_3D.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    ContactLink3D::ContactLink3D( IndexType NewId, 
                                  GeometryType::Pointer pGeometry)
    : Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    ContactLink3D::ContactLink3D( IndexType NewId, GeometryType::Pointer pGeometry,  
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
    ContactLink3D::ContactLink3D( IndexType NewId, GeometryType::Pointer pGeometry,  
                                  PropertiesType::Pointer pProperties,
                                  Condition::Pointer Master, 
                                  Condition::Pointer Slave,
                                  Point<3>& MasterContactLocalPoint,
                                  Point<3>& SlaveContactLocalPoint,
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
        GetValue( MASTER_CONTACT_GLOBAL_POINT ) = GlobalCoordinates(GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_GLOBAL_POINT ), GetValue( MASTER_CONTACT_LOCAL_POINT ) );
//Test for calculating coordinates at time step midpoint
        GetValue( SLAVE_CONTACT_GLOBAL_POINT ) = GlobalCoordinates(GetValue( CONTACT_LINK_SLAVE ), GetValue( SLAVE_CONTACT_GLOBAL_POINT ), GetValue( SLAVE_CONTACT_LOCAL_POINT ) );
//
        GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ) = SlaveIntegrationPointIndex;

        GetValue(CONTACT_LINK_M ).resize(2,2,false);
        noalias(GetValue(CONTACT_LINK_M )) = ZeroMatrix(2,2);
    }
    
    //********************************************************
    //**** Operations ****************************************
    //********************************************************
            
    
    Condition::Pointer ContactLink3D::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new ContactLink3D(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    /**
     * Destructor. Never to be called manually
     */
    ContactLink3D::~ContactLink3D()
    {
    }

    /**
     * returns the tangential vectors of the current surface in an arbitrary point due to current configuration
     * @param Surface Surface Condition for that the Tangential Vector should be calculated
     * @param rPoint local coordinates of the point in Surface for that the Tangential Vector should be calculated
     * @return Matrix(2,3) of the two tangential vectors
     */
    Matrix ContactLink3D::TangentialVectors( Condition::Pointer Surface, 
                                             const GeometryType::CoordinatesArrayType& rPoint )
    {
        //setting up result matrix
		Matrix T( 2, 3 );
        noalias(T) = ZeroMatrix( 2, 3 );
        //shape function gradients
        Matrix DN = ZeroMatrix( Surface->GetGeometry().PointsNumber(),2);
        Surface->GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
        //calculating tangential vectors
        for( unsigned int n=0; n<Surface->GetGeometry().PointsNumber(); n++ )
        {
            T(0,0) += (Surface->GetGeometry().GetPoint(n).X0()
		+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
		*DN(n,0);
            T(0,1) += (Surface->GetGeometry().GetPoint(n).Y0()
		+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
		*DN(n,0);
            T(0,2) += (Surface->GetGeometry().GetPoint(n).Z0()
		+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
		*DN(n,0);
            T(1,0) += (Surface->GetGeometry().GetPoint(n).X0()
		+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
		*DN(n,1);
            T(1,1) += (Surface->GetGeometry().GetPoint(n).Y0()
		+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
		*DN(n,1);
            T(1,2) += (Surface->GetGeometry().GetPoint(n).Z()
		+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
		*DN(n,1);
        }
        return( T );
    }
    /**
     * returns the tangential vectors of the current surface in an arbitrary point due to reference configuration
     * @param Surface Surface Condition for that the Tangential Vector should be calculated
     * @param rPoint local coordinates of the point in Surface for that the Tangential Vector should be calculated
     * @return Matrix(2,3) of the two tangential vectors
     */
    Matrix ContactLink3D::TangentialVectors_inOrigin( Condition::Pointer Surface, 
                                             const GeometryType::CoordinatesArrayType& rPoint )
    {
        //setting up result matrix
		Matrix T( 2, 3 );
        noalias(T) = ZeroMatrix( 2, 3 );
        //shape function gradients
        Matrix DN = ZeroMatrix( Surface->GetGeometry().PointsNumber(),2);
        Surface->GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
        //calculating tangential vectors
        for( unsigned int n=0; n<Surface->GetGeometry().PointsNumber(); n++ )
        {
            T(0,0) += Surface->GetGeometry().GetPoint(n).X0()*DN(n,0);
            T(0,1) += Surface->GetGeometry().GetPoint(n).Y0()*DN(n,0);
            T(0,2) += Surface->GetGeometry().GetPoint(n).Z0()*DN(n,0);
            T(1,0) += Surface->GetGeometry().GetPoint(n).X0()*DN(n,1);
            T(1,1) += Surface->GetGeometry().GetPoint(n).Y0()*DN(n,1);
            T(1,2) += Surface->GetGeometry().GetPoint(n).Z0()*DN(n,1);
        }
        return( T );
    }
    /**
     * returns the normalized tangential vectors of the current surface in an arbitrary point due to current configuration
     * @param Surface Surface Condition for that the Tangential Vector should be calculated
     * @param rPoint local coordinates of the point in Surface for that the Tangential Vector should be calculated
     * @return Matrix(2,3) of the two tangential vectors
     */
    Matrix ContactLink3D::TangentialVectorsGlobal( Condition::Pointer Surface, 
                                             const GeometryType::CoordinatesArrayType& rPoint )
    {
        //setting up result matrix
        Matrix T( 2, 3);
        noalias(T) = TangentialVectors( Surface, rPoint );

        //shape function gradients
        double normT1= sqrt(T(0,0)*T(0,0)+T(0,1)*T(0,1)+T(0,2)*T(0,2));
        double normT2=sqrt(T(1,0)*T(1,0)+T(1,1)*T(1,1)+T(1,2)*T(1,2));
        T(0,0) = T(0,0)/normT1;
        T(0,1) = T(0,0)/normT1;
        T(0,2) = T(0,0)/normT1;
        T(1,0) = T(0,0)/normT2;
        T(1,1) = T(0,0)/normT2;
        T(1,2) = T(0,0)/normT2;
        return( T );
    }
    /**
     * returns the normal vector in arbitrary point 
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector 
     */
    Vector ContactLink3D::NormalVector( Condition::Pointer Surface,
                                        const GeometryType::CoordinatesArrayType& rPoint )
    {
        Vector Result(3);
        noalias(Result) = ZeroVector(3);
        //getting tangential vectors
        Matrix T(2,3);
        noalias(T) = TangentialVectors( Surface, rPoint );
        //calculating normal vector
        Result[0] = T(0,1)*T(1,2)-T(0,2)*T(1,1);
        Result[1] = T(0,2)*T(1,0)-T(0,0)*T(1,2);
        Result[2] = T(0,0)*T(1,1)-T(0,1)*T(1,0);
        
        SD_MathUtils<double>::Normalize( Result );

        return( Result );
    }
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void ContactLink3D::CalculateRightHandSide( VectorType& rRightHandSideVector, 
            ProcessInfo& rCurrentProcessInfo)
    {
        //calculation flags
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
    void ContactLink3D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                              VectorType& rRightHandSideVector, 
                                              ProcessInfo& rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }
    /**
     * This function calculates all system contributions due to the contact problem
     * with regard to the current master and slave partners.
     * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node 
     */
    void ContactLink3D::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                      VectorType& rRightHandSideVector, 
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        //**********************************************
        //setting up the dimensions of the contributions
        //**********************************************
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        unsigned int dim = 3;
        
        //resizing as needed the LHS
        int MatSize=(MasterNN+SlaveNN)*dim;
        
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            rLeftHandSideMatrix.resize(MatSize,MatSize,false);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize); //resetting LHS
        }
        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            rRightHandSideVector.resize(MatSize,false);
            noalias(rRightHandSideVector) = ZeroVector(MatSize); //resetting RHS
        }
        //*******************************
        //Calculating general information
        //*******************************
        
        //calculating shape function values for current slave element
        Vector SlaveShapeFunctionValues( GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size() );
		noalias(SlaveShapeFunctionValues) = ZeroVector( GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size() );
        for( IndexType PointNumber = 0; 
             PointNumber<GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size(); PointNumber++ )
        {
            SlaveShapeFunctionValues[PointNumber] 
                    = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().ShapeFunctionValue( PointNumber, 
            GetValue( SLAVE_CONTACT_LOCAL_POINT ) );
        }
        //calculating shape function values for current master element
        Vector MasterShapeFunctionValues( GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );
        noalias(MasterShapeFunctionValues) = ZeroVector( GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );
        for( IndexType PointNumber = 0;
             PointNumber < GetValue( CONTACT_LINK_MASTER )->GetGeometry().size(); PointNumber++ )
        {
            MasterShapeFunctionValues[PointNumber]
                    = GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionValue( PointNumber,
            GetValue( MASTER_CONTACT_LOCAL_POINT ) );
        }
        
        //getting first order derivatives of master element shape functions
        Matrix MasterDN = ZeroMatrix( MasterNN, dim-1 );
        noalias(MasterDN) = GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionsLocalGradients( MasterDN,
        GetValue( MASTER_CONTACT_LOCAL_POINT ) );

        //calculating normal vector on slave element
		Matrix TSlave(2,3);
		noalias(TSlave) = TangentialVectors_inOrigin( GetValue( CONTACT_LINK_SLAVE ), GetValue( SLAVE_CONTACT_LOCAL_POINT ) );

        Vector vSlaveNonNormalized(3); 
		noalias(vSlaveNonNormalized)= ZeroVector(3);

        vSlaveNonNormalized[0] = TSlave(0,1)*TSlave(1,2)-TSlave(0,2)*TSlave(1,1);
        vSlaveNonNormalized[1] = TSlave(0,2)*TSlave(1,0)-TSlave(0,0)*TSlave(1,2);
        vSlaveNonNormalized[2] = TSlave(0,0)*TSlave(1,1)-TSlave(0,1)*TSlave(1,0);

        double dASlave = MathUtils<double>::Norm3( vSlaveNonNormalized );

        //calculating normal vector on master element
        Vector vMaster(3);
        noalias(vMaster) = NormalVector( GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_LOCAL_POINT ) );
        //updating contact point coordinates
        noalias(GetValue( MASTER_CONTACT_GLOBAL_POINT ) )= ZeroVector(3);
        noalias(GetValue( SLAVE_CONTACT_GLOBAL_POINT )) = ZeroVector(3);

        noalias(GetValue( MASTER_CONTACT_GLOBAL_POINT )) = GlobalCoordinates(GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_GLOBAL_POINT ), GetValue( MASTER_CONTACT_LOCAL_POINT ));

        noalias(GetValue( SLAVE_CONTACT_GLOBAL_POINT )) = GlobalCoordinates(GetValue( CONTACT_LINK_SLAVE ), GetValue( SLAVE_CONTACT_GLOBAL_POINT ), GetValue( SLAVE_CONTACT_LOCAL_POINT ));

        //calculating and updating gap
        double Gap = inner_prod( 
                    vMaster, (GetValue( MASTER_CONTACT_GLOBAL_POINT )-GetValue( SLAVE_CONTACT_GLOBAL_POINT )) );

        GetValue( CONTACT_LINK_SLAVE )->GetValue( GAPS )[GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )] = Gap;

		Vector distVec(3);

		noalias(distVec)= (GetValue( MASTER_CONTACT_GLOBAL_POINT )-GetValue( SLAVE_CONTACT_GLOBAL_POINT ));

        //retrieving penalty value 
        double Penalty = GetValue( CONTACT_LINK_SLAVE )->GetValue(PENALTY)[GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )];
// 		KRATOS_WATCH( Penalty );

        double Penalty_T = GetValue( CONTACT_LINK_SLAVE )->GetValue(PENALTY_T)[GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )];
 
        //calculating normal contact stress
        double normalStress = ( GetValue( CONTACT_LINK_SLAVE )->GetValue(LAMBDAS)[GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )] )
                    +Penalty * Gap;

        if( normalStress <= 0.0 )
        {
            normalStress = 0.0;
  	    GetValue( CONTACT_LINK_SLAVE )->GetValue(NORMAL_STRESS)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ))= 0.0;
            return;
        } 

  		GetValue( CONTACT_LINK_SLAVE )->GetValue(NORMAL_STRESS)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ))= normalStress;

        //tangential vector on master surface
        Matrix T(2,3);
        noalias(T) = TangentialVectors( GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_LOCAL_POINT ) );
        
		Vector norm_T(2);
	
		norm_T(0)= sqrt(T(0,0)*T(0,0)+T(0,1)*T(0,1)+T(0,2)*T(0,2));

		norm_T(1)= sqrt(T(1,0)*T(1,0)+T(1,1)*T(1,1)+T(1,2)*T(1,2));

		Vector vMasterNonNormalized = ZeroVector(3);

		vMasterNonNormalized[0] = T(0,1)*T(1,2)-T(0,2)*T(1,1);
        vMasterNonNormalized[1] = T(0,2)*T(1,0)-T(0,0)*T(1,2);
        vMasterNonNormalized[2] = T(0,0)*T(1,1)-T(0,1)*T(1,0);

		double NormvMaster= MathUtils<double>::Norm3(vMasterNonNormalized);
        //
        Matrix m(2,2);
        noalias(m) = ZeroMatrix(2,2);//m[alpha][beta]=T[alpha]*T[beta]
        for( int i=0; i<2; i++ )
        {
            for( int j=0; j<2; j++ )
            {
                m(i,j) = (T(i,0)*T(j,0)+T(i,1)*T(j,1)+T(i,2)*T(j,2))/(norm_T(i)*norm_T(j));
            }
        }
        GetValue(CONTACT_LINK_M)=m; 

        Vector tangentialStresses(2);
        noalias(tangentialStresses) = ZeroVector(2);
        Vector tangentialStresses_trial(2);
        noalias(tangentialStresses_trial) = ZeroVector(2);
		Vector relativTangentialVelocity(2);
		noalias(relativTangentialVelocity)= ZeroVector(2);
		Vector relativVelocity(3);
		noalias(relativVelocity)= ZeroVector(3);
        double normTangentialStresses_trial = 0.0;
                    
        bool Stick = true;
                
        if( GetValue( CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT] > 0.0 )
        {
		//TEST Calculation of relative velocity between Master and Slave surface
		noalias(relativTangentialVelocity)= GetRelativTangentialVelocity(T);

		noalias(relativVelocity)= GetRelativVelocity();

		tangentialStresses_trial[0] =  
			( GetValue( CONTACT_LINK_SLAVE )->GetValue(LAMBDAS_T)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ),0) )
                    	+relativTangentialVelocity(0)*Penalty_T;

		tangentialStresses_trial[1] =  
			( GetValue( CONTACT_LINK_SLAVE )->GetValue(LAMBDAS_T)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ),1) )
                    	+relativTangentialVelocity(1)*Penalty_T;

            	normTangentialStresses_trial =  sqrt(tangentialStresses_trial[0]*m(0,0)*tangentialStresses_trial[0]
                        + tangentialStresses_trial[0]*m(0,1)*tangentialStresses_trial[1] 
                        + tangentialStresses_trial[1]*m(1,0)*tangentialStresses_trial[0]
                        + tangentialStresses_trial[1]*m(1,1)*tangentialStresses_trial[1]);
        	if( normTangentialStresses_trial > GetValue( 
                CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]*normalStress )//Slip
            	{
                	tangentialStresses[0] = GetValue( 
                        	CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        	*normalStress *tangentialStresses_trial[0]/normTangentialStresses_trial;//
                	tangentialStresses[1] = GetValue( 
                        	CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        	*normalStress *tangentialStresses_trial[1]/normTangentialStresses_trial;//
                	Stick = false;

  			GetValue( CONTACT_LINK_SLAVE )->GetValue(TANGENTIAL_STRESS)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ))= GetValue( 
                	CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]*normalStress;

  			GetValue( CONTACT_LINK_SLAVE )->GetValue(STICK)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ))= 0.0;
            	} 
            	else//Stick
            	{
  			GetValue( CONTACT_LINK_SLAVE )->GetValue(TANGENTIAL_STRESS)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ))= normTangentialStresses_trial;

  			GetValue( CONTACT_LINK_SLAVE )->GetValue(STICK)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ))= 1.0;

                	tangentialStresses[0]=tangentialStresses_trial[0];
                	tangentialStresses[1]=tangentialStresses_trial[1];
                	Stick = true;   
            	}
        }
        //Calculating slave element's current integration weight
        double SlaveIntegrationWeight = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().IntegrationPoints()[
                    GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )].Weight();

         //BEGIN OF ADDING RESIDUAL CONTRIBUTIONS
        if (CalculateResidualVectorFlag == true)
        {
            if (normalStress > 0.0)
            {
                CalculateAndAdd_RHS( rRightHandSideVector,
                                     MasterShapeFunctionValues,
                                     SlaveShapeFunctionValues,
                                     vMaster,
                                     T,
                                     tangentialStresses,
                                     Gap,
                                     normalStress,
                                     SlaveIntegrationWeight,
                                     dASlave );
            }
        }
        //END OF ADDING RESIDUAL CONTRIBUTIONS
        //************************************
        
        if( CalculateStiffnessMatrixFlag == false )
        {
            return;
        }

    //*******************************************************************************
    //****************** adding all contributions ***********************************
    //*******************************************************************************

     //*************************************************************
     //BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF CONTACT STRESS 
     //*************************************************************

     //BEGIN OF CONTRIBUTION DUE TO DISPLACEMENTS ON CONTACT SURFACES

     //BEGIN OF CONTRIBUTION: MASTER-MASTER
     	for(unsigned int prim=0; prim < MasterNN;  prim++ )
     	{
		for( unsigned int i=0; i<dim; i++ )
              	{
                  	for( unsigned int sec=0; sec < MasterNN; sec++ )
                  	{
                      		for( unsigned int j=0; j<dim; j++ )
                      		{
                            		rLeftHandSideMatrix(prim*dim+i,sec*dim+j) 
                                    		+= MasterShapeFunctionValues[prim]*vMaster[i]
                                    		* MasterShapeFunctionValues[sec]*vMaster[j]
                                    		* Penalty* SlaveIntegrationWeight*dASlave;
                      		}
                    	}
              	}
       	}
       //END OF CONTRIBUTION: MASTER-MASTER
       //BEGIN OF CONTRIBUTION: MASTER-SLAVE
        for( unsigned int prim = 0; prim < MasterNN; prim++ )
        {
            	for( unsigned int i=0; i<dim; i++ )
            	{
                	for( unsigned int sec = 0; sec < SlaveNN; sec++ )
                	{
                    		for( unsigned int j = 0; j<dim; j++ )
                    		{
                        		rLeftHandSideMatrix(prim*dim+i,sec*dim+j+MasterNN*dim) 
                                		-= MasterShapeFunctionValues[prim]* vMaster[i]
                                		*SlaveShapeFunctionValues[sec]*vMaster[j]*Penalty
                                		*SlaveIntegrationWeight*dASlave;
                    		}
                	}
            	}
        }
        //END OF CONTRIBUTION: MASTER-SLAVE
        //BEGIN OF CONTRIBUTION: SLAVE-MASTER
        for( unsigned int prim = 0; prim < SlaveNN; prim++ )
        {
            	for( unsigned int i=0; i<dim; i++ )
            	{
                	for( unsigned int sec=0; sec < MasterNN; sec++ )
                	{
                    		for( unsigned int j=0; j<dim; j++ )
                    		{
                        		rLeftHandSideMatrix(prim*dim+i+MasterNN*dim,sec*dim+j) 
                                		-= SlaveShapeFunctionValues[prim]*vMaster[i]
                                		*MasterShapeFunctionValues[sec]*vMaster[j]*Penalty
                                		*SlaveIntegrationWeight*dASlave;
                    		}
                	}
            	}
        } 
       //END OF CONTRIBUTION: SLAVE-MASTER
       //BEGIN OF CONTRIBUTION: SLAVE-SLAVE
        for(unsigned int prim=0; prim < SlaveNN; prim++ )
        {
            	for( unsigned int i=0; i < dim; i++ )
            	{
                	for( unsigned int sec=0;  sec < SlaveNN;  sec++ )
                	{
                    		for( unsigned int j=0; j<dim; j++ )
                    		{
                        	      rLeftHandSideMatrix(prim*dim+i+MasterNN*dim,sec*dim+j+MasterNN*dim) 
                                	+= SlaveShapeFunctionValues[prim]*vMaster[i]
                                	*SlaveShapeFunctionValues[sec]*vMaster[j]*Penalty
                                	* SlaveIntegrationWeight*dASlave;
                    		}
                	}
            	}
        }
        //END OF CONTRIBUTION: SLAVE-SLAVE
        //END OF CONTRIBUTION DUE TO DISPLACEMENTS ON CONTACT SURFACES

	//BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF NORMALVECTOR (vMaster)
	//BEGIN OF CONTRIBUTION: MASTER-MASTER
        for( unsigned int prim=0; prim < MasterNN; prim++ )
        {
                for( unsigned int i=0; i<dim; i++ )
                {
                  	for( unsigned int sec=0; sec < MasterNN;  sec++ )
                    	{
                        	for( unsigned int j=0; j<dim; j++ )
                        	{
					std::vector<int> coeff(3);
					std::vector<int> index(3);
					double NormvMaster_DU= 0.0;
					switch(j)
					{
						case 0:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 1; index[1]= 2;
							coeff[2]= -1; index[2]= 1;
							break;
						case 1:
							coeff[0]= -1; index[0]= 2;
							coeff[1]= 0; index[1]= 1;
							coeff[2]= 1; index[2]= 0;
							break;
						case 2:
							coeff[0]= 1; index[0]= 1;
							coeff[1]= -1; index[1]= 0;
							coeff[2]= 0; index[2]= 2;
							break;
						default:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 0; index[1]= 0;
							coeff[2]= 0; index[2]= 0;;
					}


					for( unsigned int k=0; k<3; k++)
					{
						NormvMaster_DU +=
							coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))
							*vMasterNonNormalized(k)/NormvMaster;
					}

					for( unsigned int k=0; k<3; k++)
					{
					   	rLeftHandSideMatrix(prim*dim+i,sec*dim+j) 
                      		              		+= MasterShapeFunctionValues[prim]*vMaster[i]*
							coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))/NormvMaster
							*distVec(k)*Penalty
							* SlaveIntegrationWeight*dASlave;

					   	rLeftHandSideMatrix(prim*dim+i,sec*dim+j) 
                	                    		+= (-1)*MasterShapeFunctionValues[prim]*vMaster[i]
							*vMasterNonNormalized(k)
							/(NormvMaster*NormvMaster)
							*NormvMaster_DU
							*distVec(k)*Penalty
							* SlaveIntegrationWeight*dASlave;
					}
                   		}
                    	}
                }
        }
	//END OF CONTRIBUTION: MASTER-MASTER
	//BEGIN OF CONTRIBUTION: SLAVE-MASTER
        for( unsigned int prim = 0; prim < SlaveNN; prim++ )
        {
            	for( unsigned int i=0; i<dim; i++ )
            	{
                	for( unsigned int sec=0; sec < MasterNN; sec++ )
                	{
                    		for( unsigned int j=0; j<dim; j++ )
                    		{				
					std::vector<int> coeff(3);
					std::vector<int> index(3);
					double NormvMaster_DU= 0.0;
					switch(j)
					{
						case 0:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 1; index[1]= 2;
							coeff[2]= -1; index[2]= 1;
							break;
						case 1:
							coeff[0]= -1; index[0]= 2;
							coeff[1]= 0; index[1]= 1;
							coeff[2]= 1; index[2]= 0;
							break;
						case 2:
							coeff[0]= 1; index[0]= 1;
							coeff[1]= -1; index[1]= 0;
							coeff[2]= 0; index[2]= 2;
							break;
						default:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 0; index[1]= 0;
							coeff[2]= 0; index[2]= 0;
					}

					for( unsigned int k=0; k<3; k++)
					{
						NormvMaster_DU +=
							coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))
							*vMasterNonNormalized(k)/NormvMaster;

					}
	
					for( unsigned int k=0; k<3; k++)
					{
      		                 		rLeftHandSideMatrix(prim*dim+i+MasterNN*dim,sec*dim+j) 
      			                          	-= SlaveShapeFunctionValues[prim]*vMaster[i]*
							coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))
							*distVec(k)*Penalty/NormvMaster
                      		        		*SlaveIntegrationWeight*dASlave;

                        			rLeftHandSideMatrix(prim*dim+i+MasterNN*dim,sec*dim+j) 
                       			         	-= (-1)*SlaveShapeFunctionValues[prim]*vMaster[i]
							*vMasterNonNormalized(k)
							/(NormvMaster*NormvMaster)*NormvMaster_DU
							*distVec(k)*Penalty
                               				*SlaveIntegrationWeight*dASlave;
					}
                    		}
                	}
            	}
        }
	//END OF CONTRIBUTION: SLAVE-MASTER
	//END OF CONTRIBUTION DUE TO LINEARIZATION OF NORMALVECTOR (vMaster)
	//BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF NORMALVECTOR IN VARIATION(vMaster)
	//BEGIN OF CONTRIBUTION: MASTER-MASTER
        for( unsigned int prim=0; prim < MasterNN; prim++ )
        {
                for( unsigned int i=0; i<dim; i++ )
                {
                    	for( unsigned int sec=0; sec < MasterNN; sec++ )
                    	{
                        	for( unsigned int j=0; j<dim; j++ )
                        	{
					std::vector<int> coeff(3);
					std::vector<int> index(3);
					double NormvMaster_DU= 0.0;
					switch(j)
					{
						case 0:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 1; index[1]= 2;
							coeff[2]= -1; index[2]= 1;
							break;
						case 1:
							coeff[0]= -1; index[0]= 2;
							coeff[1]= 0; index[1]= 1;
							coeff[2]= 1; index[2]= 0;
							break;
						case 2:
							coeff[0]= 1; index[0]= 1;
							coeff[1]= -1; index[1]= 0;
							coeff[2]= 0; index[2]= 2;
							break;
						default:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 0; index[1]= 0;
							coeff[2]= 0; index[2]= 0;;
					}


					for( unsigned int k=0; k<3; k++)
					{
						NormvMaster_DU +=
							coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))
							*vMasterNonNormalized(k)/NormvMaster;
					}

				        rLeftHandSideMatrix(prim*dim+i,sec*dim+j) 
               		                     	+= MasterShapeFunctionValues[prim]*normalStress
						*1/NormvMaster*coeff[i]*(MasterDN(sec,1)*T(0,index[j])
						-MasterDN(sec,0)*T(1,index[j]))
						* SlaveIntegrationWeight*dASlave;

				   	rLeftHandSideMatrix(prim*dim+i,sec*dim+j) 
                                    		+= (-1)*MasterShapeFunctionValues[prim]*normalStress
						*vMasterNonNormalized(i)
						/(NormvMaster*NormvMaster)
						*NormvMaster_DU
						* SlaveIntegrationWeight*dASlave;
                        	}
                    	}
                }
        }
	//END OF CONTRIBUTION: MASTER-MASTER
	//BEGIN OF CONTRIBUTION: SLAVE-MASTER
        for( unsigned int prim = 0; prim < SlaveNN; prim++ )
        {
            	for( unsigned int i=0; i<dim; i++ )
            	{
                	for( unsigned int sec=0; sec < MasterNN; sec++ )
                	{
                    		for( unsigned int j=0; j<dim; j++ )
                    		{				
					std::vector<int> coeff(3);
					std::vector<int> index(3);
					double NormvMaster_DU= 0.0;
					switch(j)
					{
						case 0:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 1; index[1]= 2;
							coeff[2]= -1; index[2]= 1;
							break;
						case 1:
							coeff[0]= -1; index[0]= 2;
							coeff[1]= 0; index[1]= 1;
							coeff[2]= 1; index[2]= 0;
							break;
						case 2:
							coeff[0]= 1; index[0]= 1;
							coeff[1]= -1; index[1]= 0;
							coeff[2]= 0; index[2]= 2;
							break;
						default:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 0; index[1]= 0;
							coeff[2]= 0; index[2]= 0;
					}
					for( unsigned int k=0; k<3; k++)
					{
						NormvMaster_DU +=
							coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))
							*vMasterNonNormalized(k)/NormvMaster;
					}

				        rLeftHandSideMatrix(prim*dim+i+MasterNN*dim,sec*dim+j) 
                                    		-= SlaveShapeFunctionValues[prim]*normalStress
						*1/NormvMaster*coeff[i]*(MasterDN(sec,1)*T(0,index[j])
						-MasterDN(sec,0)*T(1,index[j]))
						* SlaveIntegrationWeight*dASlave;

				   	rLeftHandSideMatrix(prim*dim+i+MasterNN*dim,sec*dim+j) 
                                    		+= SlaveShapeFunctionValues[prim]*normalStress
						*vMasterNonNormalized(i)
						/(NormvMaster*NormvMaster)
						*NormvMaster_DU
						* SlaveIntegrationWeight*dASlave;
				}
                	}
            	}
        }
        //***********************************************************
	//END OF CONTRIBUTION: SLAVE-MASTER
        //***********************************************************
	//END OF CONTRIBUTION DUE TO LINEARIZATION OF NORMALVECTOR IN VARIATION(vMaster)
        //***********************************************************
        //END OF CONTRIBUTION DUE TO LINEARIZATION OF CONTACT STRESS 
        //***********************************************************
        if( GetValue( CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT] < 0.0 )
		return;
        //****************************************************************
        //BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF FRICTIONAL STRESS 
        //****************************************************************

	//BEGIN MASTER-MASTER
 	for( unsigned int prim=0; prim < MasterNN; prim++ )
        {
        	for( unsigned int i=0; i<dim; i++ )
                {
                	Vector XiPrim = ZeroVector(2);

			XiPrim[0]=-MasterShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			XiPrim[1]=-MasterShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < MasterNN; sec++ )
                        {
                            	for( unsigned int j=0; j<dim; j++ )
                            	{

					Vector tangentialStresses_DU(2);

					Vector tangentialStresses_trial_DU(2);		

					tangentialStresses_trial_DU(0)=
						Penalty_T*(MasterDN(sec,0)/norm_T(0)*relativVelocity(j)
						-1/(norm_T(0)*norm_T(0)*norm_T(0))
						*(relativVelocity(0)*T(0,0)+relativVelocity(1)*T(0,1)
						+relativVelocity(2)*T(0,2))*T(0,j)*MasterDN(sec,0));

					tangentialStresses_trial_DU(1)=
						Penalty_T*(MasterDN(sec,1)/norm_T(1)*relativVelocity(j)
						-1/(norm_T(1)*norm_T(1)*norm_T(1))
						*(relativVelocity(0)*T(1,0)+relativVelocity(1)*T(1,1)
						+relativVelocity(2)*T(1,2))*T(1,j)*MasterDN(sec,1));

					if(Stick)
					{
						tangentialStresses_DU(0)=
							tangentialStresses_trial_DU(0);

						tangentialStresses_DU(1)=
							tangentialStresses_trial_DU(1);
					}
					else
					{
				   		Vector normTangentialStresses_trial_tangStresses(2);

				   		normTangentialStresses_trial_tangStresses(0)=
							1/(2*normTangentialStresses_trial)*
							(2*m(0,0)*tangentialStresses_trial[0]
                				 	+(m(0,1)+m(1,0))*tangentialStresses_trial[1]);

				   		normTangentialStresses_trial_tangStresses(1)=
							1/(2*normTangentialStresses_trial)*
							(2*m(1,1)*tangentialStresses_trial[1]
                			 		+(m(0,1)+m(1,0))*tangentialStresses_trial[0]);

				   		double normTangentialStresses_trial_DU=
							normTangentialStresses_trial_tangStresses(0)*
							tangentialStresses_trial_DU(0)+
							normTangentialStresses_trial_tangStresses(1)*
							tangentialStresses_trial_DU(1);

				  		 double normalStress_DU= MasterShapeFunctionValues[sec]
                                			*vMaster[j]*Penalty;

				   		tangentialStresses_DU(0)=
							GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
							*normalStress*(
							tangentialStresses_trial_DU(0)/normTangentialStresses_trial-tangentialStresses_trial(0)
							/(normTangentialStresses_trial*normTangentialStresses_trial)
							*normTangentialStresses_trial_DU)
							+GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
							*normalStress_DU*tangentialStresses_trial[0]/normTangentialStresses_trial;

				   		tangentialStresses_DU(1)=
							GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
							*normalStress*(
							tangentialStresses_trial_DU(1)/normTangentialStresses_trial-tangentialStresses_trial(1)
							/(normTangentialStresses_trial*normTangentialStresses_trial)
							*normTangentialStresses_trial_DU)
							+GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
							*normalStress_DU*tangentialStresses_trial[1]/normTangentialStresses_trial;
					}

                                	rLeftHandSideMatrix(prim*dim+i,sec*dim+j) +=
						(tangentialStresses_DU[0]*XiPrim[0]+
                        			tangentialStresses_DU[1] * XiPrim[1]) 
						* SlaveIntegrationWeight* dASlave;
                            	}
			}
      		}
     	}
       	//END OF CONTRIBUTION: MASTER-MASTER
        //BEGIN OF CONTRIBUTION: SLAVE -MASTER
            
        for( unsigned int prim=0; prim < SlaveNN;  prim++ )
        {
          	for( unsigned int i=0; i<dim; i++ )
                {
                	Vector XiPrim = ZeroVector(2);
	
			XiPrim[0]=SlaveShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			XiPrim[1]=SlaveShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < MasterNN; sec++ )
                        {
                            	for( unsigned int j=0; j<dim; j++ )
                            	{
					Vector tangentialStresses_DU(2);

					Vector tangentialStresses_trial_DU(2);		

					tangentialStresses_trial_DU(0)=
						Penalty_T*(MasterDN(sec,0)/norm_T(0)*relativVelocity(j)
						-1/(norm_T(0)*norm_T(0)*norm_T(0))*(relativVelocity(0)*T(0,0)+relativVelocity(1)*T(0,1)+relativVelocity(2)*T(0,2))*T(0,j)*MasterDN(sec,0));

					tangentialStresses_trial_DU(1)=
						Penalty_T*(MasterDN(sec,1)/norm_T(1)*relativVelocity(j)
						-1/(norm_T(1)*norm_T(1)*norm_T(1))*(relativVelocity(0)*T(1,0)+relativVelocity(1)*T(1,1)+relativVelocity(2)*T(1,2))*T(1,j)*MasterDN(sec,1));

					if(Stick)
					{
						tangentialStresses_DU(0)=
							tangentialStresses_trial_DU(0);

						tangentialStresses_DU(1)=
							tangentialStresses_trial_DU(1);
					}
					else
					{
					   Vector normTangentialStresses_trial_tangStresses(2);

					   normTangentialStresses_trial_tangStresses(0)=
						1/(2*normTangentialStresses_trial)*
						(2*m(0,0)*tangentialStresses_trial[0]
                				 +(m(0,1)+m(1,0))*tangentialStresses_trial[1]);

					   normTangentialStresses_trial_tangStresses(1)=
						1/(2*normTangentialStresses_trial)*
						(2*m(1,1)*tangentialStresses_trial[1]
                				 +(m(0,1)+m(1,0))*tangentialStresses_trial[0]);

					   double normTangentialStresses_trial_DU=
						normTangentialStresses_trial_tangStresses(0)*
						tangentialStresses_trial_DU(0)+
						normTangentialStresses_trial_tangStresses(1)*
						tangentialStresses_trial_DU(1);

					   double normalStress_DU= MasterShapeFunctionValues[sec]
                          	 	     	*vMaster[j]
                          	     	 	*Penalty;

					   tangentialStresses_DU(0)=
						GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
						*normalStress*(
						tangentialStresses_trial_DU(0)/normTangentialStresses_trial-tangentialStresses_trial(0)
						/(normTangentialStresses_trial*normTangentialStresses_trial)
						*normTangentialStresses_trial_DU)
						+GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]*normalStress_DU*tangentialStresses_trial[0]/normTangentialStresses_trial;

				   	tangentialStresses_DU(1)=
						GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
						*normalStress*(
						tangentialStresses_trial_DU(1)/normTangentialStresses_trial-tangentialStresses_trial(1)
						/(normTangentialStresses_trial*normTangentialStresses_trial)
						*normTangentialStresses_trial_DU)
						+GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]*normalStress_DU*tangentialStresses_trial[1]/normTangentialStresses_trial;
					}

                             	rLeftHandSideMatrix(prim*dim+i+MasterNN*dim,sec*dim+j) +=
					(tangentialStresses_DU[0]*XiPrim[0]+
                        		tangentialStresses_DU[1] * XiPrim[1]) 
					* SlaveIntegrationWeight* dASlave;
                            	}
                        }
               	}
   	}
       //END OF CONTRIBUTION: SLAVE-MASTER 
        //****************************************************************
        //BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF VARIATION 
        //****************************************************************
	   Matrix kronecker(3,3);
	   noalias(kronecker)= ZeroMatrix(3,3);
	   for(int i=0; i<3; i++)
		kronecker(i,i)= 1.0;

 	for( unsigned int prim=0; prim < MasterNN; prim++ )
        {
                for( unsigned int i=0; i<dim; i++ )
                {
                        for( unsigned int sec=0; sec < MasterNN; sec++ )
                        {
                            	for( unsigned int j=0; j<dim; j++ )
                            	{
                			Vector XiPrim_DU= ZeroVector(2);

                			XiPrim_DU[0]=-MasterShapeFunctionValues[prim]*
					(kronecker(i,j)*MasterDN(sec,0)/norm_T(0)
					-T(0,i)/(norm_T(0)*norm_T(0)*norm_T(0))*T(0,j)*MasterDN(sec,0));

                			XiPrim_DU[1]=-MasterShapeFunctionValues[prim]*
					(kronecker(i,j)*MasterDN(sec,1)/norm_T(1)
					-T(1,i)/(norm_T(1)*norm_T(1)*norm_T(1))*T(1,j)*MasterDN(sec,1));

					rLeftHandSideMatrix(prim*dim+i,sec*dim+j) +=
						(tangentialStresses[0]*XiPrim_DU[0]+
                        			tangentialStresses[1]*XiPrim_DU[1]) 
						* SlaveIntegrationWeight* dASlave;
				}
                        }
               }
      }
            //END OF CONTRIBUTION: MASTER-MASTER
 
            //BEGIN OF CONTRIBUTION: SLAVE -MASTER
      for( unsigned int prim=0; prim < SlaveNN; prim++ )
       {
         	for( unsigned int i=0; i<dim; i++ )
             	{
                        for( unsigned int sec=0; sec < MasterNN;  sec++ )
                        {
                            	for( unsigned int j=0; j<dim; j++ )
                            	{
                			Vector XiPrim_DU= ZeroVector(2);

                			XiPrim_DU[0]=SlaveShapeFunctionValues[prim]*
					(kronecker(i,j)*MasterDN(sec,0)/norm_T(0)
					-T(0,i)/(norm_T(0)*norm_T(0)*norm_T(0))*T(0,j)*MasterDN(sec,0));

                			XiPrim_DU[1]=SlaveShapeFunctionValues[prim]*
					(kronecker(i,j)*MasterDN(sec,1)/norm_T(1)
					-T(1,i)/(norm_T(1)*norm_T(1)*norm_T(1))*T(1,j)*MasterDN(sec,1));

                         	       rLeftHandSideMatrix(MasterNN*dim+prim*dim+i,sec*dim+j) +=
						(tangentialStresses[0]*XiPrim_DU[0]+
                        			tangentialStresses[1]*XiPrim_DU[1]) 
						*SlaveIntegrationWeight* dASlave;
				}
                        }
                }
        }
        //****************************************************************
        //END OF CONTRIBUTION DUE TO LINEARIZATION OF VARIATION 
        //****************************************************************

	if(Stick)  return;
//             //***************************************************************** 
//             //BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF FRICTIONAL STRESS 
//		// ONLY IF STICK IS FALSE
//             //*****************************************************************
        //BEGIN OF CONTRIBUTION MASTER-SLAVE
        for( unsigned int prim=0; prim < MasterNN; prim++ )
        {
   		for( unsigned int i=0; i<dim; i++ )
     		{
    	            	Vector XiPrim = ZeroVector(2);

			XiPrim[0]=-MasterShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			XiPrim[1]=-MasterShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < SlaveNN; sec++ )
                        {
				for( unsigned int j=0; j<dim; j++ )
                            	{
					double normalStress_DU= (-1)*SlaveShapeFunctionValues[sec]
                                		*vMaster[j]*Penalty;

					Vector tangentialStresses_DU(2);		

					tangentialStresses_DU(0)=GetValue( 
                        			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        			*normalStress_DU *tangentialStresses_trial[0]/normTangentialStresses_trial;

					tangentialStresses_DU(1)=GetValue( 
                        			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        			*normalStress_DU *tangentialStresses_trial[1]/normTangentialStresses_trial;

                                	rLeftHandSideMatrix(prim*dim+i,MasterNN*dim+sec*dim+j) +=
						(tangentialStresses_DU[0]*XiPrim[0]+
                        			tangentialStresses_DU[1]*XiPrim[1]) 
						* SlaveIntegrationWeight* dASlave;
                            	}
                        }
                }
	}
       //END OF CONTRIBUTION: MASTER-SLAVE
       //BEGIN OF CONTRIBUTION: SLAVE -SLAVE
        for( unsigned int prim=0; prim < SlaveNN; prim++ )
        {
        	for( unsigned int i=0; i<dim; i++ )
              	{
  	              	Vector XiPrim = ZeroVector(2);

			XiPrim[0]=SlaveShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			XiPrim[1]=SlaveShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < SlaveNN; sec++ )
                        {
                        	for( unsigned int j=0; j<dim; j++ )
                        	{
					double normalStress_DU= (-1)*SlaveShapeFunctionValues[sec]
                                		*vMaster[j]*Penalty;

					Vector tangentialStresses_DU(2);		

					tangentialStresses_DU(0)=GetValue( 
                        			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        			*normalStress_DU *tangentialStresses_trial[0]/normTangentialStresses_trial;

					tangentialStresses_DU(1)=GetValue( 
                        			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
						*normalStress_DU *tangentialStresses_trial[1]/normTangentialStresses_trial;

                                	rLeftHandSideMatrix(MasterNN*dim+prim*dim+i,MasterNN*dim+sec*dim+j) 	+=
						(tangentialStresses_DU[0]*XiPrim[0]+
                        			tangentialStresses_DU[1]*XiPrim[1]) 
						* SlaveIntegrationWeight* dASlave;
			    	}
                        }
              	}
   	}
       //BEGIN OF CONTRIBUTION: SLAVE -SLAVE
//             //***************************************************************** 
//             //END OF CONTRIBUTION DUE TO LINEARIZATION OF FRICTIONAL STRESS 
//		// ONLY IF STICK IS FALSE
//             //*****************************************************************

//             //***************************************************************** 
//             //BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF NORMALVECTOR 
//		// ONLY IF STICK IS FALSE
//             //*****************************************************************
 	for( unsigned int prim=0; prim < MasterNN; prim++ )
        {
             	for( unsigned int i=0; i<dim; i++ )
          	{
               	 	Vector XiPrim = ZeroVector(2);

			XiPrim[0]=-MasterShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			XiPrim[1]=-MasterShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < MasterNN; sec++ )
                        {
                            	for( unsigned int j=0; j<dim; j++ )
                           	{
					std::vector<int> coeff(3);
					std::vector<int> index(3);
					double NormvMaster_DU= 0.0;
					switch(j)
					{
						case 0:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 1; index[1]= 2;
							coeff[2]= -1; index[2]= 1;
							break;
						case 1:
							coeff[0]= -1; index[0]= 2;
							coeff[1]= 0; index[1]= 1;
							coeff[2]= 1; index[2]= 0;
							break;
						case 2:
							coeff[0]= 1; index[0]= 1;
							coeff[1]= -1; index[1]= 0;
							coeff[2]= 0; index[2]= 2;
							break;
						default:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 0; index[1]= 0;
							coeff[2]= 0; index[2]= 0;;
					}
					for( unsigned int k=0; k<3; k++)
					{
						NormvMaster_DU +=
							coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))
							*vMasterNonNormalized(k)/NormvMaster;
					}
					double normalStress_DU= 0.0;

					for( unsigned int k=0; k<3; k++)
					{
					   	normalStress_DU
                  		                  	+= coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))/NormvMaster
							*distVec(k)*Penalty;

					   	normalStress_DU
              		                      		+= (-1)*vMasterNonNormalized(k)
							/(NormvMaster*NormvMaster)
							*NormvMaster_DU
							*distVec(k)*Penalty;
	
					}
					Vector tangentialStresses_DU(2);		

					tangentialStresses_DU(0)=GetValue( 
                        			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        			*(normalStress_DU *tangentialStresses_trial[0])/normTangentialStresses_trial
						;

					tangentialStresses_DU(1)=GetValue( 
                        			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        			*(normalStress_DU *tangentialStresses_trial[1])/normTangentialStresses_trial
						;

					rLeftHandSideMatrix(prim*dim+i,sec*dim+j) +=
						(tangentialStresses_DU[0]*XiPrim[0]+
                        			tangentialStresses_DU[1]*XiPrim[1]) 
						* SlaveIntegrationWeight* dASlave;
				}
                       	}
               	}
 	}
        //END OF CONTRIBUTION: MASTER-MASTER
        //BEGIN OF CONTRIBUTION: SLAVE -MASTER
        for( unsigned int prim=0; prim < SlaveNN;  prim++ )
        {
          	for( unsigned int i=0; i<dim; i++ )
       		{
      	  	        Vector XiPrim = ZeroVector(2);

			XiPrim[0]=SlaveShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			XiPrim[1]=SlaveShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < MasterNN; sec++ )
                        {
                            	for( unsigned int j=0; j<dim; j++ )
                            	{
					std::vector<int> coeff(3);
					std::vector<int> index(3);
					double NormvMaster_DU= 0.0;
					switch(j)
					{
						case 0:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 1; index[1]= 2;
							coeff[2]= -1; index[2]= 1;
							break;
						case 1:
							coeff[0]= -1; index[0]= 2;
							coeff[1]= 0; index[1]= 1;
							coeff[2]= 1; index[2]= 0;
							break;
						case 2:
							coeff[0]= 1; index[0]= 1;
							coeff[1]= -1; index[1]= 0;
							coeff[2]= 0; index[2]= 2;
							break;
						default:
							coeff[0]= 0; index[0]= 0;
							coeff[1]= 0; index[1]= 0;
							coeff[2]= 0; index[2]= 0;;
					}

	
					for( unsigned int k=0; k<3; k++)
					{
						NormvMaster_DU +=
							coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))
							*vMasterNonNormalized(k)/NormvMaster;
					}

					double normalStress_DU= 0.0;

					for( unsigned int k=0; k<3; k++)
					{
					   	normalStress_DU
                  		                  	+= coeff[k]*(T(0,index[k])*MasterDN(sec,1)
							-MasterDN(sec,0)*T(1,index[k]))/NormvMaster
							*distVec(k)*Penalty;

					   	normalStress_DU
              		                      		+= (-1)*vMasterNonNormalized(k)
							/(NormvMaster*NormvMaster)
							*NormvMaster_DU
							*distVec(k)*Penalty;
					}

					Vector tangentialStresses_DU(2);		

					tangentialStresses_DU(0)=GetValue( 
                        			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        			*(normalStress_DU *tangentialStresses_trial[0])/normTangentialStresses_trial
						;

					tangentialStresses_DU(1)=GetValue( 
                        			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        			*(normalStress_DU *tangentialStresses_trial[1])/normTangentialStresses_trial
						;

                                	rLeftHandSideMatrix(MasterNN*dim+prim*dim+i,sec*dim+j) +=
						(tangentialStresses_DU[0]*XiPrim[0]+
                        			tangentialStresses_DU[1]*XiPrim[1]) 
						*SlaveIntegrationWeight* dASlave;
				}
                        }
                }
        }
        //END OF CONTRIBUTION: SLAVE-MASTER 
        //***************************************************************** 
        //END OF CONTRIBUTION DUE TO LINEARIZATION OF NORMALVECTOR 
        //*****************************************************************  
//             //END OF CONTRIBUTION DUE TO LINEARIZATION OF NORMALVECTOR 
//		// ONLY IF STICK IS FALSE
//             //*****************************************************************
//             //*****************************************************************
//             //***************************************************************** 
//             //BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF m 
//		// ONLY IF STICK IS FALSE
//             //*****************************************************************
 	for( unsigned int prim=0; prim < MasterNN; prim++ )
  	{
        	for( unsigned int i=0; i<dim; i++ )
          	{
                	Vector XiPrim = ZeroVector(2);

			XiPrim[0]=-MasterShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			XiPrim[1]=-MasterShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < MasterNN; sec++ )
                        {
                            for( unsigned int j=0; j<dim; j++ )
                            {
				Vector norm_T_DU(2);
				noalias(norm_T_DU)= ZeroVector(2);
				
				norm_T_DU(0)=1/norm_T(0)*T(0,j)*MasterDN(sec,0);
				norm_T_DU(1)=1/norm_T(1)*T(1,j)*MasterDN(sec,1);
					//LinearisationOfm
				Matrix m_DU(2,2);
				noalias(m_DU)=ZeroMatrix(2,2);
				for(int alpha=0; alpha<2; alpha++)
				{
					for(int beta=0; beta<2; beta++)
					{
					m_DU(alpha,beta)=	
					(MasterDN(sec,alpha)*T(beta,j)+
					T(alpha,j)*MasterDN(sec,beta))/(norm_T(alpha)*norm_T(beta))
					-(T(alpha,j)*T(beta,j))/(norm_T(alpha)*norm_T(alpha)*norm_T(beta))
					*norm_T_DU(alpha)
					-(T(alpha,j)*T(beta,j))/(norm_T(alpha)*norm_T(beta)*norm_T(beta))
					*norm_T_DU(beta);	
					}
				}					

				double normTangentialStresses_trial_DU=
					1/(2*normTangentialStresses_trial)*(
					tangentialStresses_trial[0]*m_DU(0,0)
					*tangentialStresses_trial[0]+
					tangentialStresses_trial[1]*m_DU(1,0)
					*tangentialStresses_trial[0]+
					tangentialStresses_trial[0]*m_DU(0,1)
					*tangentialStresses_trial[1]+
					tangentialStresses_trial[1]*m_DU(1,1)
					*tangentialStresses_trial[1]);

				Vector tangentialStresses_DU(2);		

				tangentialStresses_DU(0)=
					GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
					*normalStress*((-1)*tangentialStresses_trial(0)
					/(normTangentialStresses_trial*normTangentialStresses_trial)
					*normTangentialStresses_trial_DU);

			   	tangentialStresses_DU(1)=
					GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
					*normalStress*((-1)*tangentialStresses_trial(1)
					/(normTangentialStresses_trial*normTangentialStresses_trial)
					*normTangentialStresses_trial_DU);

				rLeftHandSideMatrix(prim*dim+i,sec*dim+j) +=
					(tangentialStresses_DU[0]*XiPrim[0]+
                        		tangentialStresses_DU[1]*XiPrim[1]) 
					* SlaveIntegrationWeight* dASlave;
				}
                       	}
                }
    	}
        //END OF CONTRIBUTION: MASTER-MASTER
        //BEGIN OF CONTRIBUTION: SLAVE -MASTER
       	for( unsigned int prim=0; prim < SlaveNN; prim++ )
     	{
        	for( unsigned int i=0; i<dim; i++ )
   		{
                	Vector XiPrim = ZeroVector(2);

			XiPrim[0]=SlaveShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			XiPrim[1]=SlaveShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < MasterNN; sec++ )
                        {
                            	for( unsigned int j=0; j<dim; j++ )
                            	{
					Vector norm_T_DU(2);
					noalias(norm_T_DU)= ZeroVector(2);
				
					norm_T_DU(0)=1/norm_T(0)*T(0,j)*MasterDN(sec,0);
					norm_T_DU(1)=1/norm_T(1)*T(1,j)*MasterDN(sec,1);
					//LinearisationOfm
					Matrix m_DU(2,2);
					noalias(m_DU)=ZeroMatrix(2,2);
					for(int alpha=0; alpha<2; alpha++)
					{
						for(int beta=0; beta<2; beta++)
						{
						m_DU(alpha,beta)=	
						(MasterDN(sec,alpha)*T(beta,j)+
						T(alpha,j)*MasterDN(sec,beta))/(norm_T(alpha)*norm_T(beta))
						-(T(alpha,j)*T(beta,j))/(norm_T(alpha)*norm_T(alpha)*norm_T(beta))
						*norm_T_DU(alpha)
						-(T(alpha,j)*T(beta,j))/(norm_T(alpha)*norm_T(beta)*norm_T(beta))
						*norm_T_DU(beta);	
						}
					}					

				double normTangentialStresses_trial_DU=
					1/(2*normTangentialStresses_trial)*(
					tangentialStresses_trial[0]*m_DU(0,0)
					*tangentialStresses_trial[0]+
					tangentialStresses_trial[1]*m_DU(1,0)
					*tangentialStresses_trial[0]+
					tangentialStresses_trial[0]*m_DU(0,1)
					*tangentialStresses_trial[1]+
					tangentialStresses_trial[1]*m_DU(1,1)
					*tangentialStresses_trial[1]);

				Vector tangentialStresses_DU(2);		

				tangentialStresses_DU(0)=
					GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
					*normalStress*((-1)*tangentialStresses_trial(0)
					/(normTangentialStresses_trial*normTangentialStresses_trial)
					*normTangentialStresses_trial_DU);

			   	tangentialStresses_DU(1)=
					GetValue(CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
					*normalStress*((-1)*tangentialStresses_trial(1)
					/(normTangentialStresses_trial*normTangentialStresses_trial)
					*normTangentialStresses_trial_DU);

                                rLeftHandSideMatrix(MasterNN*dim+prim*dim+i,sec*dim+j) +=
					(tangentialStresses_DU[0]*XiPrim[0]+
                        		tangentialStresses_DU[1]*XiPrim[1]) 
					*SlaveIntegrationWeight* dASlave;
				}
                       	}
             	}
	}
//             //***************************************************************** 
//             //END OF CONTRIBUTION DUE TO LINEARIZATION OF m
//             //***************************************************************** 
//             //BEGIN OF CONTRIBUTION DUE TO LINEARIZATION OF NORMALVECTOR 
//		// ONLY IF STICK IS FALSE
//             //*****************************************************************
//             //*****************************************************************
        //***********************************************************
        //***********************************************************
        //END OF CONTRIBUTION DUE TO FRICTIONAL SRESS
        //***********************************************************
        //***********************************************************
        KRATOS_CATCH("")
    } // CalculateAll
    
    /**
     * This function calculates the system contributions to the global damp matrix due to the contact problem
     * with regard to the current master and slave partners.
     * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node 
     */
    void ContactLink3D::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
    {

	KRATOS_TRY

	if( !(GetValue( CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT] > 0.0 ))
        {
		return;
	}
        //**********************************************
        //setting up the dimensions of the contributions
        //**********************************************
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        unsigned int dim = 3;
        
        //resizing as needed the LHS
        int MatSize=(MasterNN+SlaveNN)*dim;
        
        rDampMatrix.resize(MatSize,MatSize,false);
        noalias(rDampMatrix) = ZeroMatrix(MatSize,MatSize); //resetting LHS

        //calculating shape function values for current slave element
        Vector SlaveShapeFunctionValues(GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size() );
	noalias(SlaveShapeFunctionValues) = ZeroVector( GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size() );
        for( IndexType PointNumber = 0; 
             PointNumber<GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size(); PointNumber++ )
        {
            SlaveShapeFunctionValues[PointNumber] 
                    = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().ShapeFunctionValue( PointNumber, 
            GetValue( SLAVE_CONTACT_LOCAL_POINT ) );
        }
        
        //calculating shape function values for current master element
        Vector MasterShapeFunctionValues(GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );
	noalias(MasterShapeFunctionValues) = ZeroVector( GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );

        for( IndexType PointNumber = 0;
             PointNumber < GetValue( CONTACT_LINK_MASTER )->GetGeometry().size(); PointNumber++ )
        {
            MasterShapeFunctionValues[PointNumber]
                    = GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionValue( PointNumber,
            GetValue( MASTER_CONTACT_LOCAL_POINT ) );
        }
        //calculating areas increment
	Matrix TSlave(2,3);
	noalias(TSlave) = TangentialVectors_inOrigin( GetValue( CONTACT_LINK_SLAVE ), GetValue( SLAVE_CONTACT_LOCAL_POINT ) );

        Vector vSlaveNonNormalized(3);
	noalias(vSlaveNonNormalized) = ZeroVector(3);

        vSlaveNonNormalized(0) = TSlave(0,1)*TSlave(1,2)-TSlave(0,2)*TSlave(1,1);
        vSlaveNonNormalized(1) = TSlave(0,2)*TSlave(1,0)-TSlave(0,0)*TSlave(1,2);
        vSlaveNonNormalized(2) = TSlave(0,0)*TSlave(1,1)-TSlave(0,1)*TSlave(1,0);

        double dASlave = MathUtils<double>::Norm3( vSlaveNonNormalized );

        //calculating normal vector on master element
        Vector vMaster(3);
	noalias(vMaster) = NormalVector( GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_LOCAL_POINT ) );
        //updating contact point coordinates
        noalias(GetValue( MASTER_CONTACT_GLOBAL_POINT )) = ZeroVector(3);
        noalias(GetValue( SLAVE_CONTACT_GLOBAL_POINT )) = ZeroVector(3);

        noalias(GetValue( MASTER_CONTACT_GLOBAL_POINT )) = GlobalCoordinates(GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_GLOBAL_POINT ), GetValue( MASTER_CONTACT_LOCAL_POINT ));

        noalias(GetValue( SLAVE_CONTACT_GLOBAL_POINT )) = GlobalCoordinates(GetValue( CONTACT_LINK_SLAVE ), GetValue( SLAVE_CONTACT_GLOBAL_POINT ), GetValue( SLAVE_CONTACT_LOCAL_POINT ));

        //calculating and updating gap
        double Gap = inner_prod( 
                    vMaster, (GetValue( MASTER_CONTACT_GLOBAL_POINT )-GetValue( SLAVE_CONTACT_GLOBAL_POINT )) );

        GetValue( CONTACT_LINK_SLAVE )->GetValue( GAPS )[GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )] = Gap;
        
        //retrieving penalty value 
        double Penalty = GetValue( CONTACT_LINK_SLAVE )->GetValue(PENALTY)[GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )];

        double Penalty_T = GetValue( CONTACT_LINK_SLAVE )->GetValue(PENALTY_T)[GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )];
 
        //calculating normal contact stress
        double normalStress = ( GetValue( CONTACT_LINK_SLAVE )->GetValue(LAMBDAS)[GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )] )
                    +Penalty * Gap;
        
        if( normalStress < 0.0)
        {
            normalStress = 0.0;
            return;
        }  
        //tangential vector on master surface
        Matrix T(2,3);
	noalias(T) = TangentialVectors( GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_LOCAL_POINT ) );
       
	Vector norm_T(2);
                
	norm_T(0)= sqrt(T(0,0)*T(0,0)+T(0,1)*T(0,1)+T(0,2)*T(0,2));

	norm_T(1)= sqrt(T(1,0)*T(1,0)+T(1,1)*T(1,1)+T(1,2)*T(1,2));
        //
        Matrix m(2,2);
	noalias(m) = ZeroMatrix(2,2);//m[alpha][beta]=T[alpha]*T[beta]
        for( int i=0; i<2; i++ )
        {
            for( int j=0; j<2; j++ )
            {
                m(i,j) = (T(i,0)*T(j,0)+T(i,1)*T(j,1)+T(i,2)*T(j,2))/(norm_T(i)*norm_T(j));
            }
        }

        Vector tangentialStresses(2);
	noalias(tangentialStresses)= ZeroVector(2);
        Vector tangentialStresses_trial(2);
	noalias(tangentialStresses_trial) = ZeroVector(2);
	Vector tangentialVelocity(2);
	noalias(tangentialVelocity)= ZeroVector(2);
        double normTangentialStresses_trial;
                    
        bool Stick = true;

	//TEST Calculation of relative velocity between Master and Slave surface
	noalias(tangentialVelocity)= GetRelativTangentialVelocity(T);
	
	tangentialStresses_trial[0] =  
		( GetValue( CONTACT_LINK_SLAVE )->GetValue(LAMBDAS_T)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ),0) )
                    +tangentialVelocity(0)*Penalty_T;

	tangentialStresses_trial[1] =  
		( GetValue( CONTACT_LINK_SLAVE )->GetValue(LAMBDAS_T)(GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ),1) )
                +tangentialVelocity(1)*Penalty_T;

         normTangentialStresses_trial =  sqrt(tangentialStresses_trial[0]*m(0,0)*tangentialStresses_trial[0]
                 + tangentialStresses_trial[0]*m(0,1)*tangentialStresses_trial[1] 
                 + tangentialStresses_trial[1]*m(1,0)*tangentialStresses_trial[0]
                 + tangentialStresses_trial[1]*m(1,1)*tangentialStresses_trial[1]);
            
         if( normTangentialStresses_trial > GetValue( 
               CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]*normalStress )//Slip
         {
                tangentialStresses[0] = GetValue( 
                        CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        *normalStress *tangentialStresses_trial[0]/normTangentialStresses_trial;//
                tangentialStresses[1] = GetValue( 
                        CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        *normalStress *tangentialStresses_trial[1]/normTangentialStresses_trial;//	

                Stick = false;
         }
         else//Stick
         {
                tangentialStresses[0]=tangentialStresses_trial[0];
                tangentialStresses[1]=tangentialStresses_trial[1];
                Stick = true;   
         }

        //Calculating slave element's current integration weight
        double SlaveIntegrationWeight = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().IntegrationPoints()[
                    GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )].Weight();

		//BEGIN OF CONTRIBUTION: MASTER-MASTER
		// std::cout << "MASTER-MASTER" << std::endl;

            for( unsigned int prim=0; prim < MasterNN; prim++ )
            {
                    for( unsigned int i=0; i<dim; i++ )
                    {
                	Vector Xi = ZeroVector(2);

			Xi[0]=-MasterShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			Xi[1]=-MasterShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < MasterNN; sec++ )
                        {
                            for( unsigned int j=0; j<dim; j++ )
                            {
				//Derivative RelativeTangentialVelocity/MasterVelocity
				Vector tangentialStresses_DV(2);

				Vector tangentialStresses_trial_DV(2);
				
				tangentialStresses_trial_DV(0)=
					Penalty_T*T(0,j)/norm_T(0)*(-1)*MasterShapeFunctionValues[sec];

				tangentialStresses_trial_DV(1)=
					Penalty_T*T(1,j)/norm_T(1)*(-1)*MasterShapeFunctionValues[sec];

				if(Stick)
				{
				  tangentialStresses_DV(0)=
					tangentialStresses_trial_DV(0);

				   tangentialStresses_DV(1)=
					tangentialStresses_trial_DV(1);
				}
				else
				{
				   Vector normTangentialStresses_trial_tangStresses(2);

				   double normTangentialStresses_trial_DV;

				   normTangentialStresses_trial_tangStresses(0)=
					1/(2*normTangentialStresses_trial)*
					(2*m(0,0)*tangentialStresses_trial[0]
                			 +(m(0,1)+m(1,0))*tangentialStresses_trial[1]);

				   normTangentialStresses_trial_tangStresses(1)=
					1/(2*normTangentialStresses_trial)*
					(2*m(1,1)*tangentialStresses_trial[1]
                			 +(m(0,1)+m(1,0))*tangentialStresses_trial[0]);

				   normTangentialStresses_trial_DV=  
					normTangentialStresses_trial_tangStresses(0)
					*tangentialStresses_trial_DV(0)
					+
					normTangentialStresses_trial_tangStresses(1)
					*tangentialStresses_trial_DV(1);
					
				   tangentialStresses_DV(0)=
					GetValue( 
                        		CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        		*normalStress*
					(tangentialStresses_trial_DV(0)/normTangentialStresses_trial
					-tangentialStresses_trial(0)/pow(normTangentialStresses_trial,2)*normTangentialStresses_trial_DV);

				   tangentialStresses_DV(1)=
					GetValue( 
                        		CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        		*normalStress*
					(tangentialStresses_trial_DV(1)/normTangentialStresses_trial
					-tangentialStresses_trial(1)/pow(normTangentialStresses_trial,2)*normTangentialStresses_trial_DV);
				  }

				rDampMatrix(prim*dim+i,sec*dim+j) 
                                        +=(tangentialStresses_DV[0]*Xi[0]+tangentialStresses_DV[1]*Xi[1])
					* SlaveIntegrationWeight* dASlave;
                            }
                      }
                 }
            }


            //END OF CONTRIBUTION: MASTER-MASTER
 
 
            //BEGIN OF CONTRIBUTION MASTER-SLAVE
            for( unsigned int prim=0; prim < MasterNN; prim++ )
            {
                    for( unsigned int i=0; i<dim; i++ )
                    {

                	Vector Xi = ZeroVector(2);

			Xi[0]=-MasterShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			Xi[1]=-MasterShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < SlaveNN; sec++ )
                        {
			
                            	for( unsigned int j=0; j<dim; j++ )
                            	{
					Vector tangentialStresses_DV(2);

					Vector tangentialStresses_trial_DV(2);
				
					tangentialStresses_trial_DV(0)=
						Penalty_T*T(0,j)/norm_T(0)*SlaveShapeFunctionValues[sec];

					tangentialStresses_trial_DV(1)=
						Penalty_T*T(1,j)/norm_T(1)*SlaveShapeFunctionValues[sec];

					if(Stick)
					{
					  tangentialStresses_DV(0)=
						tangentialStresses_trial_DV(0);

				 	  tangentialStresses_DV(1)=
						tangentialStresses_trial_DV(1);
					}
					else
					{
					   Vector normTangentialStresses_trial_tangStresses(2);

					   double normTangentialStresses_trial_DV;
	
					   normTangentialStresses_trial_tangStresses(0)=
						1/(2*normTangentialStresses_trial)*
						(2*m(0,0)*tangentialStresses_trial[0]
                				 +(m(0,1)+m(1,0))*tangentialStresses_trial[1]);

					   normTangentialStresses_trial_tangStresses(1)=
						1/(2*normTangentialStresses_trial)*
						(2*m(1,1)*tangentialStresses_trial[1]
                				 +(m(0,1)+m(1,0))*tangentialStresses_trial[0]);

					   normTangentialStresses_trial_DV=  
						normTangentialStresses_trial_tangStresses(0)
						*tangentialStresses_trial_DV(0)
						+
						normTangentialStresses_trial_tangStresses(1)
						*tangentialStresses_trial_DV(1);
					
					   tangentialStresses_DV(0)=
						GetValue( 
                       	 			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        			*normalStress*
						(tangentialStresses_trial_DV(0)/normTangentialStresses_trial
						-tangentialStresses_trial(0)/pow(normTangentialStresses_trial,2)*normTangentialStresses_trial_DV);

					   tangentialStresses_DV(1)=
						GetValue( 
                        			CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        			*normalStress*
						(tangentialStresses_trial_DV(1)/normTangentialStresses_trial
						-tangentialStresses_trial(1)/pow(normTangentialStresses_trial,2)*normTangentialStresses_trial_DV);
				  	}

				rDampMatrix(prim*dim+i,MasterNN*dim+sec*dim+j) 
                                        +=(tangentialStresses_DV[0]*Xi[0]+tangentialStresses_DV[1]*Xi[1])
					* SlaveIntegrationWeight* dASlave;
                        	}
                        }
                }
       }
        //END OF CONTRIBUTION: MASTER-SLAVE
        //BEGIN OF CONTRIBUTION: SLAVE -MASTER
            for( unsigned int prim=0; prim < SlaveNN; prim++ )
            {
                    for( unsigned int i=0; i<dim; i++ )
                    {

                	Vector Xi = ZeroVector(2);

			Xi[0]=SlaveShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			Xi[1]=SlaveShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < MasterNN;  sec++ )
                        {

                            for( unsigned int j=0; j<dim; j++ )
                            {
				Vector tangentialStresses_DV(2);

				Vector tangentialStresses_trial_DV(2);
				
				tangentialStresses_trial_DV(0)=
					Penalty_T*T(0,j)/norm_T(0)*(-1)*MasterShapeFunctionValues[sec];

				tangentialStresses_trial_DV(1)=
					Penalty_T*T(1,j)/norm_T(1)*(-1)*MasterShapeFunctionValues[sec];

				if(Stick)
				{
				  tangentialStresses_DV(0)=
					tangentialStresses_trial_DV(0);

				   tangentialStresses_DV(1)=
					tangentialStresses_trial_DV(1);
				}
				else
				{
				   Vector normTangentialStresses_trial_tangStresses(2);

				   double normTangentialStresses_trial_DV;


				   normTangentialStresses_trial_tangStresses(0)=
					1/(2*normTangentialStresses_trial)*
					(2*m(0,0)*tangentialStresses_trial[0]
                			 +(m(0,1)+m(1,0))*tangentialStresses_trial[1]);

				   normTangentialStresses_trial_tangStresses(1)=
					1/(2*normTangentialStresses_trial)*
					(2*m(1,1)*tangentialStresses_trial[1]
                			 +(m(0,1)+m(1,0))*tangentialStresses_trial[0]);

				   normTangentialStresses_trial_DV=  
					normTangentialStresses_trial_tangStresses(0)
					*tangentialStresses_trial_DV(0)
					+
					normTangentialStresses_trial_tangStresses(1)
					*tangentialStresses_trial_DV(1);
					
				   tangentialStresses_DV(0)=
					GetValue( 
                        		CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        		*normalStress*
					(tangentialStresses_trial_DV(0)/normTangentialStresses_trial
					-tangentialStresses_trial(0)/pow(normTangentialStresses_trial,2)*normTangentialStresses_trial_DV);

				   tangentialStresses_DV(1)=
					GetValue( 
                        		CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        		*normalStress*
					(tangentialStresses_trial_DV(1)/normTangentialStresses_trial
					-tangentialStresses_trial(1)/pow(normTangentialStresses_trial,2)*normTangentialStresses_trial_DV);
				  }

				rDampMatrix(MasterNN*dim+prim*dim+i,sec*dim+j) 
                                        +=(tangentialStresses_DV[0]*Xi[0]+tangentialStresses_DV[1]*Xi[1])
					* SlaveIntegrationWeight* dASlave;
                            }
                        }
                    }
            }
            //END OF CONTRIBUTION: SLAVE-MASTER 
            //BEGIN OF CONTRIBUTION: SLAVE -SLAVE
            for( unsigned int prim=0;  prim < SlaveNN; prim++ )
            {
                    for( unsigned int i=0; i<dim; i++ )
                    {
                	Vector Xi = ZeroVector(2);

			Xi[0]=SlaveShapeFunctionValues[prim] * T(0,i)/norm_T(0);

			Xi[1]=SlaveShapeFunctionValues[prim] * T(1,i)/norm_T(1);

                        for( unsigned int sec=0; sec < SlaveNN; sec++ )
                        {
                            for( unsigned int j=0; j<dim; j++ )
                            {
				Vector tangentialStresses_DV(2);

				Vector tangentialStresses_trial_DV(2);
				
				tangentialStresses_trial_DV(0)=
					Penalty_T*T(0,j)/norm_T(0)*SlaveShapeFunctionValues[sec];

				tangentialStresses_trial_DV(1)=
					Penalty_T*T(1,j)/norm_T(1)*SlaveShapeFunctionValues[sec];

				if(Stick)
				{
				  tangentialStresses_DV(0)=
					tangentialStresses_trial_DV(0);

				   tangentialStresses_DV(1)=
					tangentialStresses_trial_DV(1);
				}
				else
				{
				   Vector normTangentialStresses_trial_tangStresses(2);

				   double normTangentialStresses_trial_DV;


				   normTangentialStresses_trial_tangStresses(0)=
					1/(2*normTangentialStresses_trial)*
					(2*m(0,0)*tangentialStresses_trial[0]
                			 +(m(0,1)+m(1,0))*tangentialStresses_trial[1]);


				   normTangentialStresses_trial_tangStresses(1)=
					1/(2*normTangentialStresses_trial)*
					(2*m(1,1)*tangentialStresses_trial[1]
                			 +(m(0,1)+m(1,0))*tangentialStresses_trial[0]);

				   normTangentialStresses_trial_DV=  
					normTangentialStresses_trial_tangStresses(0)
					*tangentialStresses_trial_DV(0)
					+
					normTangentialStresses_trial_tangStresses(1)
					*tangentialStresses_trial_DV(1);
					
				   tangentialStresses_DV(0)=
					GetValue( 
                        		CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        		*normalStress*
					(tangentialStresses_trial_DV(0)/normTangentialStresses_trial
					-tangentialStresses_trial(0)/pow(normTangentialStresses_trial,2)*normTangentialStresses_trial_DV);

				   tangentialStresses_DV(1)=
					GetValue( 
                        		CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT]
                        		*normalStress*
					(tangentialStresses_trial_DV(1)/normTangentialStresses_trial
					-tangentialStresses_trial(1)/pow(normTangentialStresses_trial,2)*normTangentialStresses_trial_DV);
				  }

				rDampMatrix(MasterNN*dim+prim*dim+i,MasterNN*dim+sec*dim+j) 
                                        +=(tangentialStresses_DV[0]*Xi[0]+tangentialStresses_DV[1]*Xi[1])
					* SlaveIntegrationWeight* dASlave;
			    }
                        }
                   }
               }
        KRATOS_CATCH("")
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
    void ContactLink3D::CalculateAndAdd_RHS( Vector& residualvector,
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

        for( unsigned int i=0; i<MasterNN; i++ )
        {
            int index = dim*i;
            double coeff = -normalStress * NMaster[i] * SlaveIntegrationWeight * dASlave;
            residualvector[index]   = coeff*vMaster[0];
            residualvector[index+1] = coeff*vMaster[1];
            residualvector[index+2] = coeff*vMaster[2];
        }
        
        //********************************************
        //END OF MASTER: CalculateAndAdd_PressureForce
        //********************************************
        
        //*********************************************
        //BEGIN OF SLAVE: CalculateAndAdd_PressureForce
        //*********************************************
        
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();

        for( unsigned int i=0; i<SlaveNN; i++ )
        {
            int index = MasterNN*dim+dim*i;
            double coeff = normalStress * NSlave[i] * SlaveIntegrationWeight * dASlave;
            residualvector[index]   = coeff*vMaster[0];
            residualvector[index+1] = coeff*vMaster[1];
            residualvector[index+2] = coeff*vMaster[2];
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
        
        if( GetValue( CONTACT_LINK_SLAVE )->GetProperties()[FRICTION_COEFFICIENT] > 0.0 )
        {

	Vector norm_T(2);

	norm_T(0)= sqrt(T(0,0)*T(0,0)+T(0,1)*T(0,1)+T(0,2)*T(0,2));
	norm_T(1)= sqrt(T(1,0)*T(1,0)+T(1,1)*T(1,1)+T(1,2)*T(1,2));

            for( unsigned int i=0; i<MasterNN; i++ )
            {
                for( unsigned int j=0 ; j<dim ; j++)
                {
                Vector Xi = ZeroVector(2);

		Xi[0]=-NMaster[i]*T(0,j)/norm_T(0);

		Xi[1]=-NMaster[i]*T(1,j)/norm_T(1);

                residualvector[ i*dim+j] -= (tangentialStresses[0]*Xi[0]+
                        tangentialStresses[1] * Xi[1]) * SlaveIntegrationWeight* dASlave;
                }
            }
               
        //**********************************************
        //END OF MASTER: CalculateAndAdd_PressureForce
        //**********************************************
        
        //**********************************************
        //BEGIN OF SLAVE: CalculateAndAdd_PressureForce
        //**********************************************
        
            
            for( unsigned int i=0; i<SlaveNN; i++ )
            {
                for( unsigned int j=0;j<dim;j++)
                {
                Vector Xi = ZeroVector(2);

		Xi[0]= NSlave[i]*T(0,j)/norm_T(0);

		Xi[1]= NSlave[i]*T(1,j)/norm_T(1);
                
                    residualvector[MasterNN*dim +i*dim+j] -= (tangentialStresses[0] * Xi[0]+
                            tangentialStresses[1] * Xi[1])* SlaveIntegrationWeight * dASlave;
                }
            }
            
            
        //**********************************************
        //END OF SLAVE: CalculateAndAdd_PressureForce
        //**********************************************
        
             
        }
        
    }
    
    /**
     * This function calculates updates the local and global coordinates 
     * of the master contact partner in order to follow the movement of
     * the slave surface along the master surface
     */
    void ContactLink3D::UpdateMasterLocalPoint( )
    {
        double Xi1 = GetValue( MASTER_CONTACT_LOCAL_POINT )[0];
        double Xi2 = GetValue( MASTER_CONTACT_LOCAL_POINT )[1];
        double deltaXi1 = 0.0;
        double deltaXi2 = 0.0;
        for( int k=0; k<1000; k++ )
        {
	// KRATOS_WATCH( GetValue( MASTER_CONTACT_GLOBAL_POINT ) );
            //setting up tangential vectors
            Vector t1 = ZeroVector(3);//first tangential vector
            Vector t2 = ZeroVector(3);//second tangential vector
                    //derivatives of tangential vectors
            Vector dt11 = ZeroVector(3);
            Vector dt12 = ZeroVector(3);
            Vector dt21 = ZeroVector(3);
            Vector dt22 = ZeroVector(3);
            
                    //retrieving first order derivatives in current solution point 
            Matrix DN = ZeroMatrix(GetValue( CONTACT_LINK_MASTER )->GetGeometry().PointsNumber(),2);
            GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionsLocalGradients( DN, GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT ) );
                    //retrieving second order derivatives in current solution point 
            GeometryType::ShapeFunctionsSecondDerivativesType D2N;
            GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionsSecondDerivatives( D2N, GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT ) );
            for( unsigned  int n=0; n<GetValue( CONTACT_LINK_MASTER )->GetGeometry().PointsNumber(); n++ )
            {
                        //contribution to tangential vectors
                t1[0] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).X0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
			*DN(n,0);
                t1[1] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Y0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
			*DN(n,0);
                t1[2] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Z0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
			*DN(n,0);
                t2[0] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).X0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
			*DN(n,1);
                t2[1] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Y0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
			*DN(n,1);
                t2[2] +=( GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Z0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
			*DN(n,1);
                        //contribution to derivatives of tangential vectors
                dt11[0] +=  (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).X0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
			*D2N[n](0,0);
                dt11[1] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Y0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
			*D2N[n](0,0);
                dt11[2] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Z0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
			*D2N[n](0,0);
                dt12[0] +=  (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).X0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
			*D2N[n](0,1);
                dt12[1] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Y0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
			*D2N[n](0,1);
                dt12[2] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Z0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
			*D2N[n](0,1);
                dt21[0] +=  (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).X0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
			*D2N[n](1,0);
                dt21[1] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Y0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
			*D2N[n](1,0);
                dt21[2] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Z0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
			*D2N[n](1,0);
                dt22[0] +=  (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).X0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
			*D2N[n](1,1);
                dt22[1] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Y0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
			*D2N[n](1,1);
                dt22[2] += (GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).Z0()
			+GetValue( CONTACT_LINK_MASTER )->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
			*D2N[n](1,1);
            }
                    //defining auxiliary terms
            double A1 = ((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0])*t1[0])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1])*t1[1])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2])*t1[2]);
            double A2 = ((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0])*t2[0])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1])*t2[1])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2])*t2[2]);
            double B11 = (-t1[0]*t1[0]-t1[1]*t1[1]-t1[2]*t1[2])
                        + ((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0])*dt11[0])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1])*dt11[1])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2])*dt11[2]);
            double B12 = (-t2[0]*t1[0]-t2[1]*t1[1]-t2[2]*t1[2])
                        + ((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0])*dt12[0])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1])*dt12[1])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2])*dt12[2]);
            double B21 = (-t1[0]*t2[0]-t1[1]*t2[1]-t1[2]*t2[2])
                        + ((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0])*dt21[0])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1])*dt21[1])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2])*dt21[2]);
            double B22 = (-t2[0]*t2[0]-t2[1]*t2[1]-t2[2]*t2[2])
                        + ((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[0]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[0])*dt22[0])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[1]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[1])*dt22[1])
                        +((GetValue( SLAVE_CONTACT_GLOBAL_POINT )[2]-GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT )[2])*dt22[2]);
                    //calculating update for Xi
            deltaXi1 = -A1*B22/(B11*B22-B12*B21)+A2*B12/(B11*B22-B12*B21);
            deltaXi2 =  A2*B21/(B11*B22-B12*B21)-A2*B11/(B11*B22-B12*B21);
                    //updating Xi
            Xi1 += deltaXi1;
            Xi2 += deltaXi2;
                    //updating LocalPoint
            GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT )[0] = Xi1;
            GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT )[1] = Xi2;
                    //updating rResult
            GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT ) = ZeroVector( 3 );
            GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT ) = GlobalCoordinates(GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_CURRENT_GLOBAL_POINT ), GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT ) );

            if( fabs(deltaXi1) < 1e-7 && fabs(deltaXi2) < 1e-7 )
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
    void ContactLink3D::EquationIdVector( EquationIdVectorType& rResult, 
                                          ProcessInfo& CurrentProcessInfo)
    {
        //determining size of DOF list
        //dimension of space
        unsigned int dim = 3;
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        
        unsigned int index;
        rResult.resize((MasterNN+SlaveNN)*dim,false);
        for(  unsigned int i=0; i<MasterNN; i++ )
        {
            index = i*dim;
            rResult[index]   = (GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].GetDof(DISPLACEMENT_X)).EquationId();
            rResult[index+1] = (GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());
            rResult[index+2] = (GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId());
        }
        for( unsigned  int i=0; i<SlaveNN; i++ )
        {
            index = MasterNN*dim+i*dim;
            rResult[index]   = (GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].GetDof(DISPLACEMENT_X)).EquationId();
            rResult[index+1] = (GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());
            rResult[index+2] = (GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId());
        }
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * Setting up the DOF list for the current partners.
     * All conditions are assumed to be defined in 3D space with 3 DOFs per Node.
     * All DOF are given Master first, Slave second
     */
    void ContactLink3D::GetDofList( DofsVectorType& ConditionalDofList,
                                    ProcessInfo& CurrentProcessInfo)
    {
        //determining size of DOF list
        //dimension of space
        unsigned int dim = 3;
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        
        ConditionalDofList.resize((MasterNN+SlaveNN)*dim);
        unsigned int index;
        //setting up master DOFs
        for( unsigned int i=0; i<MasterNN; i++ )
        {
            index = i*dim;
            ConditionalDofList[index]   = (GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ConditionalDofList[index+1] = (GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            ConditionalDofList[index+2] = (GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
        //setting up slave DOFs
        for( unsigned int i=0; i<SlaveNN; i++ )
        {
            index = MasterNN*dim+i*dim;
            ConditionalDofList[index]   = (GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ConditionalDofList[index+1] = (GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            ConditionalDofList[index+2] = (GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
    }

//new functions includes

     Point<3>& ContactLink3D::GlobalCoordinates(Condition::Pointer Surface, Point<3>& rResult, Point<3> const& LocalCoordinates)
    {
	noalias(rResult)= ZeroVector(3);

	for(IndexType i = 0 ; i < Surface->GetGeometry().size() ; i++)
	{
		double shape_func= Surface->GetGeometry().ShapeFunctionValue(i,LocalCoordinates);

		rResult(0) += shape_func* 
				((Surface->GetGeometry()[i]).X0()
				+(Surface->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_X));

		rResult(1) += shape_func* 
				((Surface->GetGeometry()[i]).Y0()
				+(Surface->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Y));

		rResult(2) += shape_func* 
				((Surface->GetGeometry()[i]).Z0()
				+(Surface->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Z));
	}

	return rResult;
     }

	/**
    * returns the relative tangential velocity between the quadrature point on the slave surface and its 
    * closest point projection
    * @param T Matrix of the Tanegntial Vectors on the Master Surface in the current configuration
    * @return tangential velocity
    */
	Vector ContactLink3D::GetRelativTangentialVelocity(Matrix& T)
	{
		Vector result(2);

		Vector slave_velo(3);
		Vector master_velo(3);

		noalias(slave_velo)= ZeroVector(3);
		noalias(master_velo)= ZeroVector(3);

		for(IndexType i = 0 ; i < GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size() ; i++)
		{
			double shape_func= GetValue( CONTACT_LINK_SLAVE )->GetGeometry().ShapeFunctionValue(i,GetValue( SLAVE_CONTACT_LOCAL_POINT ));
			slave_velo+= ((GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_DT))*
				shape_func;
		}
		for(IndexType i = 0 ; i < GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() ; i++)
		{
			double shape_func= GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionValue(i,GetValue( MASTER_CONTACT_LOCAL_POINT ));
			master_velo+= ((GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_DT))*
				shape_func;
		}
		Vector norm_T(2);

		norm_T(0)= sqrt(T(0,0)*T(0,0)+T(0,1)*T(0,1)+T(0,2)*T(0,2));

		norm_T(1)= sqrt(T(1,0)*T(1,0)+T(1,1)*T(1,1)+T(1,2)*T(1,2));

		result(0)= ((slave_velo(0)-master_velo(0))*T(0,0)+(slave_velo(1)-master_velo(1))*T(0,1)
			+(slave_velo(2)-master_velo(2))*T(0,2))/norm_T(0);

		result(1)= ((slave_velo(0)-master_velo(0))*T(1,0)+(slave_velo(1)-master_velo(1))*T(1,1)
			+(slave_velo(2)-master_velo(2))*T(1,2))/norm_T(1);

		return result;
	}
	/**
    * returns the relative velocity between the quadrature point on the slave surface and its 
    * closest point projection
    * @return relative velocity
    */
	Vector ContactLink3D::GetRelativVelocity()
	{
		Vector result(3);

		Vector slave_velo(3);
		Vector master_velo(3);

		noalias(slave_velo)= ZeroVector(3);
		noalias(master_velo)= ZeroVector(3);
		for(IndexType i = 0 ; i < GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size() ; i++)
		{
			double shape_func= GetValue( CONTACT_LINK_SLAVE )->GetGeometry().ShapeFunctionValue(i,GetValue( SLAVE_CONTACT_LOCAL_POINT ));
			slave_velo+= ((GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_DT))*
				shape_func;
		}
		for(IndexType i = 0 ; i < GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() ; i++)
		{
			double shape_func= GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionValue(i,GetValue( MASTER_CONTACT_LOCAL_POINT ));
			master_velo+= ((GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_DT))*
				shape_func;
		}
		result(0)= (slave_velo(0)-master_velo(0));

		result(1)= (slave_velo(1)-master_velo(1));

		result(2)= (slave_velo(2)-master_velo(2));

		return result;
	}
} // Namespace Kratos
