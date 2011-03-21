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
#include "custom_conditions/point_point_contact_link.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    PointPointContactLink::PointPointContactLink( IndexType NewId, 
                                  GeometryType::Pointer pGeometry) : 
            Condition( NewId, pGeometry )
    {
      
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    PointPointContactLink::PointPointContactLink( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties) : 
            Condition( NewId, pGeometry, pProperties )
    {
    }
    
    PointPointContactLink::PointPointContactLink( IndexType NewId, NodesArrayType const& ThisNodes)
    {      
    }
    
    Condition::Pointer PointPointContactLink::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new PointPointContactLink(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    
    PointPointContactLink::PointPointContactLink(IndexType NewId, NodesArrayType& ThisNodes)
    {
    } 
    
    
    PointPointContactLink::PointPointContactLink( 
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
				  
				  //NormalVector();
				  
				  
				}
    
    
    /**
     * Destructor. Never to be called manually
     */
    PointPointContactLink::~PointPointContactLink()
    {
    }

    
    /**
     * returns normal vector in arbitrary point 
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector 
     */
     Vector PointPointContactLink::NormalVector()
     {
       
       std::vector<array_1d<double, 3 > >    normales;
       array_1d<double, 3> e3            =   ZeroVector(3);
       array_1d<double, 3> t             =   ZeroVector(3);
       array_1d<double, 3> Result        =   ZeroVector(3);
       
       e3[0] = 0.00; 
       e3[1] = 0.00; 
       e3[2] = 1.00;
      
       
       Condition::GeometryType& geom_master      = (GetValue(CONTACT_LINK_MASTER))->GetGeometry();
       WeakPointerVector<Condition>& neighb_cond = geom_master[0].GetValue(NEIGHBOUR_CONDITIONS);
       
       for(WeakPointerVector< Condition >::iterator cond  = neighb_cond.begin(); cond!= neighb_cond.end(); cond++){
	    Condition::GeometryType& geom = cond->GetGeometry();
	    t = geom[0] - geom[1];
	    t = (1.00 / std::sqrt(inner_prod(t,t))) * t;   
            MathUtils<double>::CrossProduct(Result,e3,t);
	    normales.push_back(Result);
       }      
       
       Result        =   ZeroVector(3);
       for(unsigned int i = 0; i<normales.size(); i++)
	  noalias(Result) += normales[i];
       
       
       // sacando la normal promedio
       Result = (1.00 / std::sqrt(inner_prod(Result,Result))) * Result;  
       return Result;
       
     }
     
    
    //************************************************************************************
    //************************************************************************************ 
     
    Vector PointPointContactLink::TangentialVector()
    {
       array_1d<double, 3> t      =   ZeroVector(3);
       array_1d<double, 3> e3     =   ZeroVector(3);
       array_1d<double, 3> normal =   ZeroVector(3); 
       
       e3[0] = 0.00; 
       e3[1] = 0.00; 
       e3[2] = 1.00;
       
       
       normal = NormalVector();
       MathUtils<double>::CrossProduct(t,normal,e3); 
       return t; 
    }


    //************************************************************************************
    //************************************************************************************
    
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void PointPointContactLink::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
    void PointPointContactLink::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
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
    void PointPointContactLink::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                      VectorType& rRightHandSideVector, 
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
      
      
        //**********************************************
        //setting up the dimensions of the contributions
        //**********************************************
	unsigned int dim = 2;
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
        }
    }
    
        
    //************************************************************************************
    //************************************************************************************    
        
    void PointPointContactLink::CalculateNormalImpenetrabilityConstraint( Vector& rCn){}
    void PointPointContactLink::CalculateTangentialImpenetrabilityConstraint( Vector& rCt){}
    void PointPointContactLink::CalculateNormalContactForce(Vector& rNormalForce){}
    void PointPointContactLink::CalculateTangentialContactForce(Vector& rTangentialForce){}
        
        
    
       
    //************************************************************************************
    //************************************************************************************
    
   void PointPointContactLink::EquationIdVector( EquationIdVectorType& rResult, 
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
   void PointPointContactLink::GetDofList( DofsVectorType& ConditionalDofList,
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
   
   void PointPointContactLink::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
   {    
   }
   
   double PointPointContactLink::CalculateGap()
   { 
     return 0.00;
   }
   
   void PointPointContactLink::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
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
   
   
    void PointPointContactLink::Calculate( const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
    {
      if(rVariable==CONSTRAINT_VECTOR)
      {
	
	// La dimension esta dada desde afuera
	Output.resize(4, false);
	Vector Normal     = NormalVector();
	//Vector Tangential = TangentialVector();


	Output[0] =   Normal[0]; // slave
	Output[1] =   Normal[1];
	Output[2] =  -Normal[0]; // master
        Output[3] =  -Normal[1];
      }
    }
    
    
     void PointPointContactLink::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
     {
          
        const int& size =  GetGeometry().IntegrationPoints().size();
	rValues.resize(size);
        if(rVariable==NORMAL_CONTACT_STRESS )
            {
	        Vector& lamdas   = GetValue(LAMBDAS);
		rValues[0]       =  lamdas[0];
		KRATOS_WATCH(rValues[0])
		
            }
     }
   
   
} // Namespace Kratos
