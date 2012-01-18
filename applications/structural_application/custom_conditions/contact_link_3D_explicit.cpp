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
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2008-10-23 12:26:09 $
//   Revision:            $Revision: 1.7 $
//
//
// System includes 

// External includes 

// Project includes 
#include "includes/define.h"
#include "custom_conditions/contact_link_3D_explicit.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "geometries/plane.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    ContactLink3DExplicit::ContactLink3DExplicit( IndexType NewId, 
                                  GeometryType::Pointer pGeometry)
    : Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    ContactLink3DExplicit::ContactLink3DExplicit( IndexType NewId, GeometryType::Pointer pGeometry,  
                                  PropertiesType::Pointer pProperties
                                )
    : Condition( NewId, pGeometry, pProperties )
    {
    }
    
    ContactLink3DExplicit::ContactLink3DExplicit( IndexType NewId, GeometryType::Pointer pGeometry,  
                                  PropertiesType::Pointer pProperties,
                                  Condition::Pointer Master, 
                                  Condition::Pointer Slave,
                                  Point<3>& MasterContactLocalPoint,
                                  Point<3>& SlaveContactLocalPoint,
                                  int SlaveIntegrationPointIndex
                                )
    : Condition( NewId, pGeometry, pProperties )
    {
      
        
        GetValue( CONTACT_LINK_MASTER )                     = Master;
        GetValue( CONTACT_LINK_SLAVE )                      = Slave;
        GetValue( MASTER_CONTACT_LOCAL_POINT )              = MasterContactLocalPoint;
        GetValue( MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT ) = MasterContactLocalPoint;
        GetValue( MASTER_CONTACT_CURRENT_LOCAL_POINT )      = MasterContactLocalPoint;
        GetValue( SLAVE_CONTACT_LOCAL_POINT )               = SlaveContactLocalPoint;
	GetValue( CONTACT_LINK_MASTER ) = Master;
	GetValue( CONTACT_LINK_SLAVE )  = Slave;

	const int& size =  GetGeometry().IntegrationPoints().size();
	GetValue( LAMBDAS ).resize(size, false);
	GetValue( DELTA_LAMBDAS ).resize(size, false);
	noalias(GetValue( LAMBDAS ))        = ZeroVector(size);
	noalias(GetValue( DELTA_LAMBDAS ))  = ZeroVector(size); //GetGeometry().IntegrationPoints().size());


	GetValue(NORMAL_STRESS).resize(GetGeometry().IntegrationPoints().size(),false);
	GetValue(NORMAL_STRESS)=ZeroVector( GetGeometry().IntegrationPoints().size());
	GetValue(TANGENTIAL_STRESS).resize(GetGeometry().IntegrationPoints().size(),false);
	GetValue(TANGENTIAL_STRESS)=ZeroVector( GetGeometry().IntegrationPoints().size());
	
	
        //Test for calculating coordinates at time step midpoint
        //GetValue( MASTER_CONTACT_GLOBAL_POINT ) = GlobalCoordinates(GetValue( CONTACT_LINK_MASTER ), GetValue( MASTER_CONTACT_GLOBAL_POINT ), GetValue( MASTER_CONTACT_LOCAL_POINT ) );
	
	//Test for calculating coordinates at time step midpoint
	//GetValue( SLAVE_CONTACT_GLOBAL_POINT )  = GlobalCoordinates(GetValue( CONTACT_LINK_SLAVE ), GetValue( SLAVE_CONTACT_GLOBAL_POINT ), GetValue( SLAVE_CONTACT_LOCAL_POINT ) );
	
        //GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX ) = SlaveIntegrationPointIndex;

        //GetValue(CONTACT_LINK_M ).resize(2,2,false );
        //noalias(GetValue(CONTACT_LINK_M )) = ZeroMatrix(2,2 );
	
    }
    
    //********************************************************
    //**** Operations ****************************************
    //********************************************************
            
    
    Condition::Pointer ContactLink3DExplicit::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new ContactLink3DExplicit(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    /**
     * Destructor. Never to be called manually
     */
    ContactLink3DExplicit::~ContactLink3DExplicit()
    {
    }
    
    //************************************************************************************
    //************************************************************************************
    
    /**
     * returns the tangential vectors of the current surface in an arbitrary point due to currnt configuration
     */
    Matrix ContactLink3DExplicit::TangentialVectors( Condition::Pointer Surface, 
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

    //************************************************************************************
    //************************************************************************************
    
    /**
     * returns the tangential vectors of the current surface in an arbitrary point due to original configuration
     */
    Matrix ContactLink3DExplicit::TangentialVectors_inOrigin( Condition::Pointer Surface, 
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
     * returns the tangential Vectors of the current surface in an arbitrary point due to the Derivates in global coordinates
     */
    
    Matrix ContactLink3DExplicit::TangentialVectorsGlobal( Condition::Pointer Surface, 
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
    //************************************************************************************
    //************************************************************************************
    
    /**
     * returns normal vector in arbitrary point 
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector 
     */
    Vector ContactLink3DExplicit::NormalVector( Condition::Pointer Surface,
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
    void ContactLink3DExplicit::CalculateRightHandSide( VectorType& rRightHandSideVector, 
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
    void ContactLink3DExplicit::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                              VectorType& rRightHandSideVector, 
                                              ProcessInfo& rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * This function calculates all system contributions due to the contact problem
     * with regard to the current master and slave partners.
     * All Conditions are assumed to be defined in 3D space and havin 3 DOFs per node 
     */
    void ContactLink3DExplicit::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                      VectorType& rRightHandSideVector, 
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        
         //**********************************************
        //setting up the dimensions of the contributions
        //**********************************************
	unsigned int dim = 3;
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN  = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        unsigned int MatSize  = (MasterNN+SlaveNN)*dim;
//         //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
	    Vector Constraint;
	    rRightHandSideVector.resize(MatSize,false);
	    noalias(rRightHandSideVector) = ZeroVector(MatSize); 
	    Vector& lamdas = GetValue(LAMBDAS);
 	    Calculate(CONSTRAINT_VECTOR, Constraint, rCurrentProcessInfo);
 	    noalias(rRightHandSideVector) = -lamdas[0] * Constraint; 
	    
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
    void ContactLink3DExplicit::CalculateAndAdd_RHS( Vector& residualvector,
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
        return;
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * Setting up the EquationIdVector for the current partners.
     * All conditions are assumed to be defined in 3D space with 3 DOFs per node.
     * All Equation IDs are given Master first, Slave second
     */
    void ContactLink3DExplicit::EquationIdVector( EquationIdVectorType& rResult, 
                                          ProcessInfo& CurrentProcessInfo)
    {
        //determining size of DOF list
        //dimension of space
        unsigned int dim = 3;
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        
        unsigned int index;
        rResult.resize((MasterNN+SlaveNN)*dim);
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
    void ContactLink3DExplicit::GetDofList( DofsVectorType& ConditionalDofList,
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

     Point<3>& ContactLink3DExplicit::GlobalCoordinates(Condition::Pointer Surface, Point<3>& rResult, Point<3> const& LocalCoordinates)
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


	Vector ContactLink3DExplicit::GetRelativTangentialVelocity(Matrix& T)
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

	Vector ContactLink3DExplicit::GetRelativVelocity()
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
	
	void ContactLink3DExplicit::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
         {
            Condition::GeometryType& geom = this->GetGeometry();
            const unsigned int dimension  = geom.WorkingSpaceDimension();
            const unsigned int dim2       = geom.size()*dimension;
            rMassMatrix                   = ZeroMatrix(dim2, dim2); 
            for(unsigned int i = 0; i<geom.size(); i++){
               for(unsigned int j = 0; j<dimension; j++ ){
	           rMassMatrix(3*i + j, 3*i + j) = geom[i].FastGetSolutionStepValue(NODAL_MASS);
                  }
               }  
          }
          
          
    void ContactLink3DExplicit::Calculate( const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY 
      if(rVariable==CONSTRAINT_VECTOR)
      {
	array_1d<double, 3>               Points0;
	array_1d<double, 3>               Points1;
	array_1d<double, 3>               Points2;
	array_1d<double, 3>               Points3;

	Output.resize(12, false);
	GeometryType::CoordinatesArrayType rPoint; 

     
	///calculating shape function values for current master element
	Plane rPlane;
	Vector MasterShapeFunctionValues(GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );
	
	Condition::GeometryType& geom_1    = GetValue( CONTACT_LINK_SLAVE )->GetGeometry();	
	Condition::GeometryType& geom_2    = GetValue( CONTACT_LINK_MASTER )->GetGeometry();	
	noalias(MasterShapeFunctionValues) = ZeroVector( GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );
	
	Points0[0] = geom_1[0].X(); 
	Points0[1] = geom_1[0].Y();
	Points0[2] = geom_1[0].Z();
	
	Points1[0] = geom_2[0].X(); 
	Points1[1] = geom_2[0].Y();
	Points1[2] = geom_2[0].Z();
	Points2[0] = geom_2[1].X();
	Points2[1] = geom_2[1].Y(); 
	Points2[2] = geom_2[1].Z(); 
	Points3[0] = geom_2[2].X();
	Points3[1] = geom_2[2].Y(); 
	Points3[2] = geom_2[2].Z(); 
			
	
//	double compare_distance = rPlane.DistPoint3Triangle3(Points0, Points1, Points2, Points3);
	rPoint = rPlane.mClosestPoint;
	Vector Normal = NormalVector(GetValue(CONTACT_LINK_MASTER),rPoint);
	MasterShapeFunctionValues[0] = rPlane.mTriangleBary[0];
	MasterShapeFunctionValues[1] = rPlane.mTriangleBary[1];
	MasterShapeFunctionValues[2] = rPlane.mTriangleBary[2];
	
	
        Output[0]  =  Normal[0];
	Output[1]  =  Normal[1];
	Output[2]  =  Normal[2];
	Output[3]  = -MasterShapeFunctionValues[0]*Normal[0];
        Output[4]  = -MasterShapeFunctionValues[0]*Normal[1];
	Output[5]  = -MasterShapeFunctionValues[0]*Normal[2];
	Output[6]  = -MasterShapeFunctionValues[1]*Normal[0];
	Output[7]  = -MasterShapeFunctionValues[1]*Normal[1];
	Output[8]  = -MasterShapeFunctionValues[1]*Normal[2];
        Output[9]  = -MasterShapeFunctionValues[2]*Normal[0];
	Output[10] = -MasterShapeFunctionValues[2]*Normal[1];
	Output[11] = -MasterShapeFunctionValues[2]*Normal[2];
	
      }
      KRATOS_CATCH("")
      
    }
    
    void ContactLink3DExplicit::Calculate( const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY 
      if(rVariable==NORMAL)
            {
		Condition::GeometryType& geom_1 = this->GetGeometry();
		Vector Output;
		Calculate( CONSTRAINT_VECTOR, Output, rCurrentProcessInfo);
	        for(unsigned int i = 0; i<geom_1.size(); i++)
		{
		   geom_1[i].FastGetSolutionStepValue(NORMAL_X) += Output[3*i]   * GetValue(LAMBDAS)[0];  
		   geom_1[i].FastGetSolutionStepValue(NORMAL_Y) += Output[3*i+1] * GetValue(LAMBDAS)[0];
		   geom_1[i].FastGetSolutionStepValue(NORMAL_Z) += Output[3*i+2] * GetValue(LAMBDAS)[0];
		}  
		
	    }
      KRATOS_CATCH("")
    }
    
    
     void ContactLink3DExplicit::CalculateOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector< array_1d<double,3> >& Output, const ProcessInfo& rCurrentProcessInfo)
     { 
      KRATOS_TRY 
      if(rVariable==NORMAL)
            {
//		Condition::GeometryType& geom_1 = this->GetGeometry();
		const int& size =  GetGeometry().IntegrationPoints().size();
	        Output.resize(size);
	        GeometryType::CoordinatesArrayType rPoint; 
		Output[0] = NormalVector(GetValue(CONTACT_LINK_MASTER),rPoint);  /// ONLY Tetrahedra
	    }
      KRATOS_CATCH("")
    }
      
     
     
    
	
} // Namespace Kratos
