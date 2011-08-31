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
#include "custom_conditions/slave_contact_point_3d.h"
// #include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    SlaveContactPoint3D::SlaveContactPoint3D( IndexType NewId, 
                                  GeometryType::Pointer pGeometry) : 
            Condition( NewId, pGeometry )
    {
      
        GetValue( IS_CONTACT_SLAVE  )         = 1;
	GetValue( IS_CONTACT_MASTER )         = 0; 
	Condition::GeometryType& geom         = this->GetGeometry();
	geom[0].GetValue(IS_CONTACT_SLAVE  )  = 1;
	
    }
    
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    SlaveContactPoint3D::SlaveContactPoint3D( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties) : 
            Condition( NewId, pGeometry, pProperties )
    {
    }
    
    SlaveContactPoint3D::SlaveContactPoint3D( IndexType NewId, NodesArrayType const& ThisNodes)
    {      
    }
    
    Condition::Pointer SlaveContactPoint3D::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new SlaveContactPoint3D(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    
    // nodearraytype is equal to PointerVector<TPointType> 
    SlaveContactPoint3D::SlaveContactPoint3D(IndexType NewId, NodesArrayType& ThisNodes)
    {
        //PointerGeometryType  pgeom    = PointerGeometryType(new Line2D2 <Node<3> >( ThisNodes)); 
	//GetGeometry() = *pgeom;
	//KRATOS_WATCH(GetGeometry())
	//GetValue(IS_CONTACT_MASTER) = 1;   
	//for(NodesArrayType::iterator inode = ThisNodes.begin() ; inode!=ThisNodes.end(); inode++)  
    } 
    
    
    
    /**
     * Destructor. Never to be called manually
     */
    SlaveContactPoint3D::~SlaveContactPoint3D()
    {
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * returns condition type info
     */
//     int SlaveContactPoint3D::IsContactType()
//     {/*
//         return( 2 );*/
//     }
    //************************************************************************************
    //************************************************************************************
    /**
     * returns the tangential vectors of the current surface in an arbitrary point
     * TODO: TO BE REIMPLEMENTED!!!
     */
//     Matrix SlaveContactPoint3D::TangentialVectors( GeometryType::CoordinatesArrayType& rPoint )
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
//     Vector SlaveContactPoint3D::NormalVector( GeometryType::CoordinatesArrayType& rPoint )
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
    bool SlaveContactPoint3D::ClosestPoint( GeometryType::CoordinatesArrayType& rResultGlobal, 
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
    void SlaveContactPoint3D::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
    {
        unsigned int ndof = GetGeometry().size()*3;
        if( rRightHandSideVector.size() != ndof )
            rRightHandSideVector.resize(ndof,false);
        rRightHandSideVector = ZeroVector(ndof);
    }
    
    //************************************************************************************
    //************************************************************************************
    
    /**
     * calculates this contact element's local contributions
     */
    void SlaveContactPoint3D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                              VectorType& rRightHandSideVector, 
                                              ProcessInfo& rCurrentProcessInfo)
    {
        unsigned int ndof = GetGeometry().size()*3;
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
    void SlaveContactPoint3D::CalculateAll( MatrixType& rLeftHandSideMatrix, 
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
    void SlaveContactPoint3D::CalculateAndAddKc( Matrix& K,
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
    void SlaveContactPoint3D::CalculateAndAdd_PressureForce( Vector& residualvector,
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
   void SlaveContactPoint3D::EquationIdVector( EquationIdVectorType& rResult, 
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
            int index = i*3;
            rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
	    rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }
        KRATOS_CATCH("")
   }
       
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
   void SlaveContactPoint3D::GetDofList( DofsVectorType& ConditionalDofList,
                                   ProcessInfo& CurrentProcessInfo)
   {
	ConditionalDofList.resize(0);

	for (unsigned int i=0;i<GetGeometry().size();i++)
	{
	ConditionalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
	ConditionalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
	ConditionalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
	} 
     
   }
   
   void SlaveContactPoint3D::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
   {    
   }
} // Namespace Kratos
