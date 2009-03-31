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
*   Last Modified by:    $Author: hurga $
*   Date:                $Date: 2008-02-29 15:18:45 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/master_contact_face_3D_newmark.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    MasterContactFace3DNewmark::MasterContactFace3DNewmark( IndexType NewId, 
                                  GeometryType::Pointer pGeometry) : 
            Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    MasterContactFace3DNewmark::MasterContactFace3DNewmark( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties) : 
            Condition( NewId, pGeometry, pProperties )
    {
//         for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
//         {
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_X));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Y));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Z));
//         }
        GetValue( LAMBDAS ).resize(GetGeometry().IntegrationPoints().size(),false);
        noalias(GetValue( LAMBDAS )) = ZeroVector( GetGeometry().IntegrationPoints().size());
        GetValue( LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2 ,false);
        noalias(GetValue( LAMBDAS_T )) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( GAPS ).resize( GetGeometry().IntegrationPoints().size(),false);
        noalias(GetValue( GAPS )) = ZeroVector( GetGeometry().IntegrationPoints().size());
        GetValue( DELTA_LAMBDAS ).resize( GetGeometry().IntegrationPoints().size(),false );
        noalias(GetValue( DELTA_LAMBDAS )) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( DELTA_LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2 ,false);
        noalias(GetValue( DELTA_LAMBDAS_T )) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( PENALTY ).resize( GetGeometry().IntegrationPoints().size(),false );
        noalias(GetValue( PENALTY )) =  ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( PENALTY_T ).resize( GetGeometry().IntegrationPoints().size() ,false);
        noalias(GetValue( PENALTY_T )) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( IS_CONTACT_MASTER ) = 1;
	GetValue( IS_CONTACT_SLAVE ) = 0;
    }
    
    Condition::Pointer MasterContactFace3DNewmark::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new MasterContactFace3DNewmark(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    /**
     * Destructor. Never to be called manually
     */
    MasterContactFace3DNewmark::~MasterContactFace3DNewmark()
    {
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * returns condition type info
     */
//     int MasterContactFace3DNewmark::IsContactType()
//     {/*
//         return( 2 );*/
//     }
    //************************************************************************************
    //************************************************************************************
    /**
     * returns the tangential vectors of the current surface in an arbitrary point
     * TODO: TO BE REIMPLEMENTED!!!
     */
//     Matrix MasterContactFace3DNewmark::TangentialVectors( GeometryType::CoordinatesArrayType& rPoint )
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
//     Vector MasterContactFace3DNewmark::NormalVector( GeometryType::CoordinatesArrayType& rPoint )
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
    bool MasterContactFace3DNewmark::ClosestPoint( GeometryType::CoordinatesArrayType& rResultGlobal, 
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
    void MasterContactFace3DNewmark::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
    {
        rRightHandSideVector.resize(0,false);
    }
    
    //************************************************************************************
    //************************************************************************************
    
    /**
     * calculates this contact element's local contributions
     */
    void MasterContactFace3DNewmark::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                              VectorType& rRightHandSideVector, 
                                              ProcessInfo& rCurrentProcessInfo)
    {
        rRightHandSideVector.resize(0,false);
        rLeftHandSideMatrix(0,0);
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * calculates the contact related contributions to the system
     * Does nothing as assembling is to be switched to linking objects
     */
    void MasterContactFace3DNewmark::CalculateAll( MatrixType& rLeftHandSideMatrix, 
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
    void MasterContactFace3DNewmark::CalculateAndAddKc( Matrix& K,
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
    void MasterContactFace3DNewmark::CalculateAndAdd_PressureForce( Vector& residualvector,
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
   void MasterContactFace3DNewmark::EquationIdVector( EquationIdVectorType& rResult, 
                                         ProcessInfo& CurrentProcessInfo 
                                       )
   {
	   rResult.resize(0);
   }
       
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
   void MasterContactFace3DNewmark::GetDofList( DofsVectorType& ConditionalDofList,
                                   ProcessInfo& CurrentProcessInfo)
   {
   }
} // Namespace Kratos
