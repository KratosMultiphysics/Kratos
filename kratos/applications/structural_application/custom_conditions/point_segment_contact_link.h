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
*   Last Modified by:    $Author: Nelson$
*   Date:                $Date: 2009-03-17 14:35:29 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

#if !defined(KRATOS_POINT_SEGMENT_LINK_2D_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_SEGMENT_LINK_2D_CONDITION_H_INCLUDED  



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/serializer.h"
//#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_conditions/master_contact_face_2d.h"
#include "custom_conditions/slave_contact_point_2d.h"

namespace Kratos
{  
    /**
     * Contact surface element for 3D contact problems.
     * Defines a facet of a 3D-Element as a contact surface for
     * master contact surfaces
     * adapted from face2D originally written by riccardo.
     */
    class PointSegmentContactLink : public Condition
    {
        
        public:  
	typedef Condition BaseType;
	typedef BaseType::EquationIdVectorType EquationIdVectorType;
	typedef BaseType::MatrixType LHS_ContributionType;
	typedef MasterContactFace2D MasterContactType;  

	/**
	* REMOVED
	*/

	// typedef PointerVectorSet<StressConditionType, IndexedObject> StressConditionContainerType;
	// typedef PointerVectorSet< EquationIdVectorType, IndexedObject> EquationIdVectorContainerType;
	// typedef PointerVectorSet< LHS_ContributionType, IndexedObject> LHS_ContainerType;

	typedef Condition::GeometryType::Pointer PointerGeometryType; 

	// typedef BaseType::SecondaryCondition SecondaryConditionType;
	// typedef PointerVectorSet< SecondaryConditionType, IndexedObject> SecondaryConditionContainerType;
	  
            // Counted pointer of PointSegmentContactLink
            KRATOS_CLASS_POINTER_DEFINITION(PointSegmentContactLink);
            
            /** 
             * Default constructor.
             */
            PointSegmentContactLink( IndexType NewId, GeometryType::Pointer pGeometry);
            
	    PointSegmentContactLink( IndexType NewId, NodesArrayType const& ThisNodes);
	    
            PointSegmentContactLink( IndexType NewId, GeometryType::Pointer pGeometry, 
                           PropertiesType::Pointer pProperties);
			   
			   
	    PointSegmentContactLink(IndexType NewId, NodesArrayType& ThisNode);
	    
	    
	    PointSegmentContactLink( 
	                          IndexType NewId, 
				  GeometryType::Pointer pGeometry,  
                                  PropertiesType::Pointer pProperties,
                                  Condition::Pointer Master, 
                                  Condition::Pointer Slave
                                  //Point<3>& MasterContactLocalPoint,
                                  //Point<3>& SlaveContactLocalPoint,
                                  //int SlaveIntegrationPointIndex
                                );
	    
            /**
             * Destructor.
             */
            virtual ~PointSegmentContactLink();
      
            /**
             * Operations.
             */
            Condition::Pointer Create( IndexType NewId, 
                                       NodesArrayType const& ThisNodes,  
                                       PropertiesType::Pointer pProperties) const;

            
           
            
            
            /**
             * Calculates the local system contributions for this contact element
             */
            void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                       VectorType& rRightHandSideVector, 
                                       ProcessInfo& rCurrentProcessInfo);
            
            void CalculateRightHandSide( VectorType& rRightHandSideVector, 
                                         ProcessInfo& rCurrentProcessInfo);
            
            void EquationIdVector( EquationIdVectorType& rResult, 
                                   ProcessInfo& rCurrentProcessInfo);
            
            void GetDofList( DofsVectorType& ConditionalDofList,
                             ProcessInfo& CurrentProcessInfo);
            
            void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo);
            
	    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);
	    
	    void  MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

	    void Calculate( const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo); 
      
        protected:
	  
	MasterContactType::Pointer mpMasterContact; 
        
        
        private:
            void CalculateAll( MatrixType& rLeftHandSideMatrix, 
                               VectorType& rRightHandSideVector,
                               ProcessInfo& rCurrentProcessInfo,
                               bool CalculateStiffnessMatrixFlag,
                               bool CalculateResidualVectorFlag);
            
	   
	     double CalculateGap();
	     
	     Vector NormalVector();
	     Vector TangentialVector();
              
	     void CalculateNormalImpenetrabilityConstraint( Vector& rCn);
             void CalculateTangentialImpenetrabilityConstraint( Vector& rCt);
             void CalculateNormalContactForce(Vector& rNormalForce);
             void CalculateTangentialContactForce(Vector& rTangentialForce);
	     Vector GetRelativeVelocity();
	     Vector GetRelativeTangentialVelocity();
	     void Calculate( const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo);
            	
	    friend class Serializer;

	    // A private default constructor necessary for serialization 
	    PointSegmentContactLink(){}; 

	    virtual void save(Serializer& rSerializer) const
	    {
	    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
	    }

	    virtual void load(Serializer& rSerializer)
	    {
	    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
	    }

	     
	     
	     
            /**
             * REMOVED
             */
//            
    }; // Class PointSegmentContactLink 
}  // namespace Kratos.

#endif // KRATOS_CONTACT_FACE_3D_CONDITION_H_INCLUDED  defined 
