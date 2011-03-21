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
#include "custom_conditions/master_contact_face_2d.h"
// #include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    MasterContactFace2D::MasterContactFace2D( IndexType NewId, 
                                  GeometryType::Pointer pGeometry) : 
            Condition( NewId, pGeometry )
    {
        GetValue( IS_CONTACT_SLAVE  )  = 0;
	GetValue( IS_CONTACT_MASTER )  = 1;  
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    MasterContactFace2D::MasterContactFace2D( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties) : 
            Condition( NewId, pGeometry, pProperties )
    {
    }
    
    Condition::Pointer MasterContactFace2D::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new MasterContactFace2D(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    
    /// nodearraytype is equal to PointerVector<TPointType> 
    MasterContactFace2D::MasterContactFace2D(IndexType NewId, NodesArrayType& ThisNodes)
    {
    } 
    
    
    
    /**
     * Destructor. Never to be called manually
     */
    MasterContactFace2D::~MasterContactFace2D()
    {
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * returns condition type info
     */
//     int MasterContactFace2D::IsContactType()
//     {/*
//         return( 2 );*/
//     }
    //************************************************************************************
    //************************************************************************************
    /**
     * returns the tangential vectors of the current surface in an arbitrary point
     * TODO: TO BE REIMPLEMENTED!!!
     */
//     Matrix MasterContactFace2D::TangentialVectors( GeometryType::CoordinatesArrayType& rPoint )
//     {
//     }
    //************************************************************************************
    //************************************************************************************
    
    /**
     * returns normal vector in arbitrary point 
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector 
     */
//     Vector MasterContactFace2D::NormalVector( GeometryType::CoordinatesArrayType& rPoint )
//     {
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
    bool MasterContactFace2D::ClosestPoint( GeometryType::CoordinatesArrayType& rResultGlobal, 
                                            GeometryType::CoordinatesArrayType& rResultLocal,
                                            const GeometryType::CoordinatesArrayType& rSlaveContactGlobalPoint,
                                            const GeometryType::CoordinatesArrayType& rCandidateGlobal
                                          )
    {
        return false;
    }
    //************************************************************************************
    //************************************************************************************
    
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void MasterContactFace2D::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
    {
        unsigned int ndof = GetGeometry().size()*2;
        if( rRightHandSideVector.size() != ndof )
            rRightHandSideVector.resize(ndof,false);
        rRightHandSideVector = ZeroVector(ndof);
    }
    
    //************************************************************************************
    //************************************************************************************
    
    /**
     * calculates this contact element's local contributions
     */
    void MasterContactFace2D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                              VectorType& rRightHandSideVector, 
                                              ProcessInfo& rCurrentProcessInfo)
    {
        unsigned int ndof = GetGeometry().size()*2;
        if( rRightHandSideVector.size() != ndof )
            rRightHandSideVector.resize(ndof,false);
        rRightHandSideVector = ZeroVector(ndof);
        if( rLeftHandSideMatrix.size1() != ndof )
            rLeftHandSideMatrix(ndof,ndof);
        rLeftHandSideMatrix = ZeroMatrix(ndof,ndof);
      
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * calculates the contact related contributions to the system
     * Does nothing as assembling is to be switched to linking objects
     */
    void MasterContactFace2D::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                      VectorType& rRightHandSideVector, 
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
    } 
    
    //***********************************************************************
    //***********************************************************************
    /**
     * System matrix contribution due to contact energy
     * TO BE IMPLEMENTED
     */
    void MasterContactFace2D::CalculateAndAddKc( Matrix& K,
            const Vector& N,
            double weight,
            double dA,
            double penalty,
            Vector v
                                               )
    {
    }
    
    //***********************************************************************
    //***********************************************************************
    /**
     *
     */
    void MasterContactFace2D::CalculateAndAdd_PressureForce( Vector& residualvector,
            const Vector& N,
            Vector& v3,
            double pressure,
            double weight, double dA 
                                                     )
    {
    }
       
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
   void MasterContactFace2D::EquationIdVector( EquationIdVectorType& rResult, 
                                         ProcessInfo& CurrentProcessInfo 
                                       )
   {

        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = number_of_nodes*3;
        if(rResult.size() != dim)
            rResult.resize(dim);

        for (unsigned int i=0;i<number_of_nodes;i++)
        {
            int index = i*2;
            rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        }
        KRATOS_CATCH("")
   }
       
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
   void MasterContactFace2D::GetDofList( DofsVectorType& ConditionalDofList,
                                   ProcessInfo& CurrentProcessInfo)
   {
	ConditionalDofList.resize(0);

	for (unsigned int i=0;i<GetGeometry().size();i++)
	{
	ConditionalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
	ConditionalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
	} 
     
   }
   
   
   void MasterContactFace2D::Calculate( const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo)
   {
      if( rVariable == NORMAL )
       {
	 Output = NormalVector();
       }
       
       return;
   }
   
   
     array_1d<double,3> MasterContactFace2D::NormalVector()
     {
       
       array_1d<double, 3> e3      =   ZeroVector(3);
       array_1d<double, 3> Result  =   ZeroVector(3);
       
       e3[0] = 0.00; e3[1] = 0.00; e3[2] = 1.00; 
       
       /// El primer nodo es el slave
       Condition::GeometryType& geom = this->GetGeometry();
       /// tener normal positiva    
       array_1d<double, 3> t         =  geom[0] - geom[1];
       t = (1.00 / std::sqrt(inner_prod(t,t))) * t;   
       MathUtils<double>::CrossProduct(Result,e3,t);
            
       return Result;
       
     }
   
   
   void MasterContactFace2D::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
   {
       const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
       const GeometryType::ShapeFunctionsGradientsType& sf_gradients = GetGeometry().ShapeFunctionsLocalGradients();
       if( rVariable == NORMAL )
       {
           if(rValues.size() != integration_points.size())
                rValues.resize(integration_points.size());
           
	   Matrix T( 2, 3 );
	   Vector Result(3);
           for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
           {
               //setting up result matrix
               noalias(T) = ZeroMatrix( 2, 3 );
               //shape function gradients
               Matrix DN = sf_gradients[PointNumber];
               //calculating tangential vectors
               for( unsigned int n=0; n<GetGeometry().PointsNumber(); n++ )
               {
                   T(0,0) += (GetGeometry().GetPoint(n).X0()
                           +GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))*DN(n,0);
                   T(0,1) += (GetGeometry().GetPoint(n).Y0()
                           +GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))*DN(n,0);
                   T(0,2) += (GetGeometry().GetPoint(n).Z0()
                           +GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))*DN(n,0);
                   T(1,0) += (GetGeometry().GetPoint(n).X0()
                           +GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))*DN(n,1);
                   T(1,1) += (GetGeometry().GetPoint(n).Y0()
                           +GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))*DN(n,1);
                   T(1,2) += (GetGeometry().GetPoint(n).Z()
                           +GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))*DN(n,1);
               }
               
               //calculating normal vector
               Result[0] = T(0,1)*T(1,2)-T(0,2)*T(1,1);
               Result[1] = T(0,2)*T(1,0)-T(0,0)*T(1,2);
               Result[2] = T(0,0)*T(1,1)-T(0,1)*T(1,0);
               SD_MathUtils<double>::Normalize( Result );
               
	       //KRATOS_WATCH(Result)
               rValues[PointNumber][0] = Result[0];
               rValues[PointNumber][1] = Result[1];
               rValues[PointNumber][2] = Result[2];
           }
       }  
   }
   
   
    /*
     Vector MasterContactFace2D::NormalVector()
     {
       
       array_1d<double, 3> e3      =   ZeroVector(3);
       array_1d<double, 3> Result  =   ZeroVector(3);
       
       e3[0] = 0.00; e3[1] = 0.00; e3[2] = 1.00; 
       
       /// El primer nodo es el slave
       Condition::GeometryType& geom = this->GetGeometry();
       /// tener normal positiva    
       array_1d<double, 3> t         =  geom[0] - geom[1];
       t = (1.00 / std::sqrt(inner_prod(t,t))) * t;   
       MathUtils<double>::CrossProduct(Result,e3,t);
            
       return Result;
       
     }
     */
    
   
   
} // Namespace Kratos
