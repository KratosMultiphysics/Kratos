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
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2009-03-17 14:35:29 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/point_segment_contact_link.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    PointSegmentContactLink::PointSegmentContactLink( IndexType NewId, 
                                  GeometryType::Pointer pGeometry) : 
            Condition( NewId, pGeometry )
    {
      
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    PointSegmentContactLink::PointSegmentContactLink( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties) : 
            Condition( NewId, pGeometry, pProperties )
    {
    }
    
    PointSegmentContactLink::PointSegmentContactLink( IndexType NewId, NodesArrayType const& ThisNodes)
    {      
    }
    
    Condition::Pointer PointSegmentContactLink::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new PointSegmentContactLink(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    
    PointSegmentContactLink::PointSegmentContactLink(IndexType NewId, NodesArrayType& ThisNodes)
    {
    } 
    
    
    PointSegmentContactLink::PointSegmentContactLink( 
	                          IndexType NewId, 
				  GeometryType::Pointer pGeometry,  
                                  PropertiesType::Pointer pProperties,
                                  Condition::Pointer Master, 
                                  Condition::Pointer Slave
                                  //Point<3>& MasterContactLocalPoint,
                                  //Point<3>& SlaveContactLocalPoint,
                                  //int SlaveIntegrationPointIndex
                                )  : Condition( NewId, pGeometry, pProperties )
                                {
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
				  
				}
    
    
    /**
     * Destructor. Never to be called manually
     */
    PointSegmentContactLink::~PointSegmentContactLink()
    {
    }

    
    /**
     * returns normal vector in arbitrary point 
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector 
     */
     Vector PointSegmentContactLink::NormalVector()
     { 
       /// El primer nodo es el slave
       array_1d<double, 3> e3      =   ZeroVector(3);
       array_1d<double, 3> Result  =   ZeroVector(3);
       e3[0] = 0.00; e3[1] = 0.00; e3[2] = 1.00; 
       Condition::GeometryType& geom =  this->GetGeometry();   
       array_1d<double, 3> t         =  geom.GetPoint(1) - geom.GetPoint(2); /// tener normal positiva 
       const double tl               =  norm_2(t);
       noalias(t)                    =  t * (1.00/tl);
       MathUtils<double>::CrossProduct(Result,e3,t);
       return Result;
     }
    
    
    //************************************************************************************
    //************************************************************************************ 
     
    Vector PointSegmentContactLink::TangentialVector()
    {
      
       Condition::GeometryType& geom = this->GetGeometry();
       array_1d<double, 3> t         =  geom.GetPoint(1) - geom.GetPoint(2);
       const double tl               =  norm_2(t);
       noalias(t)                    =  t * (1.00/tl);
       return t;
    }


  	Vector PointSegmentContactLink::GetRelativeVelocity()
	{
		Vector result(3);
		Vector slave_velo(3);
		Vector master_velo(3);
		noalias(slave_velo) = ZeroVector(3);
		noalias(master_velo)= ZeroVector(3);
		
		for(IndexType i = 0 ; i < GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size() ; i++)
	            slave_velo+= ((GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i]).GetSolutionStepValue(VELOCITY));
		
		
		for(IndexType i = 0 ; i < GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() ; i++)
		    master_velo+= ((GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i]).GetSolutionStepValue(VELOCITY));
		
		result[0]= (slave_velo(0)-master_velo(0));
		result[1]= (slave_velo(1)-master_velo(1));
		result[2]= (slave_velo(2)-master_velo(2));
		
		return result;
	}
  
       Vector PointSegmentContactLink::GetRelativeTangentialVelocity()
       {
	        Vector result(3);
		Vector slave_velo(3);
		Vector master_velo(3);
		noalias(slave_velo)  = ZeroVector(3);
		noalias(master_velo) = ZeroVector(3);
	        Vector Relative_Velocity =  GetRelativeVelocity();
		return result;
	 
       }
  
  
  

    //************************************************************************************
    //************************************************************************************
    
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void PointSegmentContactLink::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
    {
      
              //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag  = true;
	
        MatrixType matrix = Matrix();
        CalculateAll(matrix, rRightHandSideVector, 
                      rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, 
                      CalculateResidualVectorFlag);
    }
    
    //************************************************************************************
    //************************************************************************************
    
    /**
     * calculates this contact element's local contributions
     */
    void PointSegmentContactLink::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                              VectorType& rRightHandSideVector, 
                                              ProcessInfo& rCurrentProcessInfo)
    {      
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * calculates the contact related contributions to the system
     * Does nothing as assembling is to be switched to linking objects
     */
    void PointSegmentContactLink::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                      VectorType& rRightHandSideVector, 
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
      
      
        //**********************************************
        //setting up the dimensions of the contributions
        //**********************************************
	unsigned int dim      = 2;
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN  = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();

        int MatSize=(MasterNN+SlaveNN)*dim;

        
        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
	    Vector Constraint;
	    rRightHandSideVector.resize(MatSize,false);
	    noalias(rRightHandSideVector) = ZeroVector(MatSize); 
	    Vector& lamdas = GetValue(LAMBDAS);
	    Calculate(CONSTRAINT_VECTOR, Constraint, rCurrentProcessInfo);
	    noalias(rRightHandSideVector) = -lamdas[0] * Constraint; 
	    
	    
	    
	    /*
	    double gap        = CalculateGap();
 	    Vector Normal     = NormalVector();
	    Vector Tangential = TangentialVector();
	    
            rRightHandSideVector.resize(MatSize,false);
	    noalias(rRightHandSideVector) = ZeroVector(MatSize); //resetting RHS
	    
	    Condition::GeometryType& geom = this->GetGeometry();
	    array_1d<double, 3> r = geom[1] - geom[0]; 
	    double segmentlength  = GetValue( CONTACT_LINK_MASTER )->GetGeometry().Length();
	    double shi            = inner_prod(r,Tangential)/segmentlength;
	    
	    Point<3> rPoint;    
	    rPoint[0] = shi;
	    rPoint[1] = shi;
	    rPoint[2] = shi;
	   
 	    ///calculating shape function values for current master element
 	    Vector MasterShapeFunctionValues( GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );
 	    noalias(MasterShapeFunctionValues) = ZeroVector( GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );

 	    
 	    for( IndexType PointNumber = 0;
 	    PointNumber < GetValue( CONTACT_LINK_MASTER )->GetGeometry().size(); PointNumber++ )
 	    {
 	      MasterShapeFunctionValues[PointNumber]
 	      = GetValue( CONTACT_LINK_MASTER )->GetGeometry().ShapeFunctionValue(PointNumber, rPoint);
 	    }
	   
	    /// contact position on the target facet
	    Condition::GeometryType& segmentgeom = GetValue( CONTACT_LINK_MASTER )->GetGeometry();
	    Point<3> Xts;    
	    noalias(Xts) = MasterShapeFunctionValues[0]*segmentgeom[0] + MasterShapeFunctionValues[1]*segmentgeom[1];
	   
	    double penalty = 1.00;   //200e9 * 50.0;
	    
	    rRightHandSideVector[0] =   0.00; //Normal[0];
            rRightHandSideVector[1] =   0.00; // Normal[1];
	    rRightHandSideVector[2] =   0.00; //-MasterShapeFunctionValues[0]*Normal[0];
	    rRightHandSideVector[3] =   0.00; //-MasterShapeFunctionValues[0]*Normal[1];
	    rRightHandSideVector[4] =   0.00; //-MasterShapeFunctionValues[1]*Normal[0];
	    rRightHandSideVector[5] =   0.00; //-MasterShapeFunctionValues[1]*Normal[1];
	    rRightHandSideVector *=  (gap * penalty) ;
	    */
	    
        }
    }
    
        
    //************************************************************************************
    //************************************************************************************    
        
    void PointSegmentContactLink::CalculateNormalImpenetrabilityConstraint( Vector& rCn){}
    void PointSegmentContactLink::CalculateTangentialImpenetrabilityConstraint( Vector& rCt){}
    void PointSegmentContactLink::CalculateNormalContactForce(Vector& rNormalForce){}
    void PointSegmentContactLink::CalculateTangentialContactForce(Vector& rTangentialForce){}
        
        
    
       
    //************************************************************************************
    //************************************************************************************
    
   void PointSegmentContactLink::EquationIdVector( EquationIdVectorType& rResult, 
                                         ProcessInfo& CurrentProcessInfo 
                                       )
   {
     
        //determining size of DOF list
        //dimension of space
        unsigned int dim      = 2;
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN  = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        
        unsigned int index;
        rResult.resize((MasterNN+SlaveNN)*dim,false);
        for(  unsigned int i=0; i<MasterNN; i++ )
        {
            index = i*dim;
            rResult[index]   = GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index+1] = GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        }
        for( unsigned  int i=0; i<SlaveNN; i++ )
        {
            index = MasterNN*dim+i*dim;
            rResult[index]   = GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index+1] = GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        }
     
   }
       
    //************************************************************************************
    //************************************************************************************
   void PointSegmentContactLink::GetDofList( DofsVectorType& ConditionalDofList,
                                   ProcessInfo& CurrentProcessInfo)
   {
    
        //determining size of DOF list
        //dimension of space
        unsigned int dim = 2;
        unsigned int MasterNN = GetValue( CONTACT_LINK_MASTER )->GetGeometry().size();
        unsigned int SlaveNN  = GetValue( CONTACT_LINK_SLAVE )->GetGeometry().size();
        
        ConditionalDofList.resize((MasterNN+SlaveNN)*dim);
        unsigned int index;
        //setting up master DOFs
        for( unsigned int i=0; i<MasterNN; i++ )
        {
            index = i*dim;
            ConditionalDofList[index]   = (GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ConditionalDofList[index+1] = (GetValue( CONTACT_LINK_MASTER )->GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        }
        //setting up slave DOFs
        for( unsigned int i=0; i<SlaveNN; i++ )
        {
            index = MasterNN*dim+i*dim;
            ConditionalDofList[index]   = (GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ConditionalDofList[index+1] = (GetValue( CONTACT_LINK_SLAVE )->GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        } 
     
   }
   
   void PointSegmentContactLink::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
   {
      const unsigned int& size =  GetGeometry().IntegrationPoints().size();
      rValues.resize(size);
      if(rVariable==NORMAL)
      {
	 for(unsigned int i = 0; i<size; i++)
	    rValues[i] = NormalVector();
      }
   }
   
   double PointSegmentContactLink::CalculateGap()
   {
     
     Condition::GeometryType& geom = this->GetGeometry();
     array_1d<double, 3> r = geom[1] - geom[0];
     array_1d<double, 3> n = NormalVector();
     return inner_prod(r,n);
    
   }
   
   void PointSegmentContactLink::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
   {
     Condition::GeometryType& geom = this->GetGeometry();
     const unsigned int dimension  = geom.WorkingSpaceDimension();
     const unsigned int dim2       = geom.size()*dimension;
     rMassMatrix = ZeroMatrix(dim2, dim2); 
     for(unsigned int i = 0; i<geom.size(); i++){
        for(unsigned int j = 0; j<dimension; j++ ){
	 rMassMatrix(2*i + j, 2*i + j) = geom[i].FastGetSolutionStepValue(NODAL_MASS);
       }
     }
     
   }
   
   
    void PointSegmentContactLink::Calculate( const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
    {
      if(rVariable==CONSTRAINT_VECTOR)
      {
	
	/// La dimension esta dada desde afuera
	const double toler = 1E-15;
	Output.resize(6, false);
	Vector Normal     = NormalVector();
	Vector Tangential = TangentialVector();

	Condition::GeometryType& geom = this->GetGeometry();
	array_1d<double, 3> r         = geom[1] - geom[0]; 
	double segmentlength          = GetValue( CONTACT_LINK_MASTER )->GetGeometry().Length();
	double shi                    = inner_prod(r,Tangential)/segmentlength;
        //const int& ID                 = ((GetValue(CONTACT_LINK_SLAVE )->GetGeometry())(0))->Id();
	 
//   	 if(ID==936 || ID==937 || ID==876)
//  	{
//  	    KRATOS_WATCH(this->Id()) 
//  	    KRATOS_WATCH(ID)
//  	    KRATOS_WATCH(GetValue(CONTACT_LINK_MASTER)->Id())
//  	    KRATOS_WATCH("---------------------------------------") 
//         }
        
	///calculating shape function values for current master element
	Vector MasterShapeFunctionValues(GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );
	noalias(MasterShapeFunctionValues) = ZeroVector( GetValue( CONTACT_LINK_MASTER )->GetGeometry().size() );

        /// WARNING = No se usa funcion de forma del elemto, pues no coinside con tesis 3.28
        MasterShapeFunctionValues[0] = 1.0 - shi;  
        MasterShapeFunctionValues[1] = shi;
        
	if(std::fabs(Normal[0])<toler) Normal[0] = 0.00; 
	if(std::fabs(Normal[1])<toler) Normal[1] = 0.00;
	
	Output[0] =  Normal[0];
	Output[1] =  Normal[1];
	Output[2] = -MasterShapeFunctionValues[0]*Normal[0];
        Output[3] = -MasterShapeFunctionValues[0]*Normal[1];
	Output[4] = -MasterShapeFunctionValues[1]*Normal[0];
	Output[5] = -MasterShapeFunctionValues[1]*Normal[1];
      }
      
    }
    
    
    void PointSegmentContactLink::Calculate( const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY 
      if(rVariable==NORMAL)
            {
		Condition::GeometryType& geom_1 = this->GetGeometry();
		Vector Output;
		Calculate( CONSTRAINT_VECTOR, Output, rCurrentProcessInfo);
	        for(unsigned int i = 0; i<geom_1.size(); i++)
		{
		   geom_1[i].FastGetSolutionStepValue(NORMAL_X) += Output[2*i]   * GetValue(LAMBDAS)[0];  
		   geom_1[i].FastGetSolutionStepValue(NORMAL_Y) += Output[2*i+1] * GetValue(LAMBDAS)[0];
		}  
		
	    }
      KRATOS_CATCH("")
    }
    
     void PointSegmentContactLink::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
     {
          
        const int& size =  GetGeometry().IntegrationPoints().size();
	rValues.resize(size);
        if(rVariable==NORMAL_CONTACT_STRESS)
            {
	        Vector& lamdas   = GetValue(LAMBDAS);
		rValues[0]       =  lamdas[0];
		//KRATOS_WATCH(rValues[0])
		
            }
     }
   
   
} // Namespace Kratos
