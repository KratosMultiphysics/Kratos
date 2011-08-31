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

#if !defined(KRATOS_SLAVE_CONTACT_POINT_3D_CONDITION_H_INCLUDED )
#define  KRATOS_SLAVE_CONTACT_POINT_3D_CONDITION_H_INCLUDED 



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
//#include "includes/ublas_interface.h"
#include "includes/variables.h"


namespace Kratos
{  
    /**
     * Contact surface element for 3D contact problems.
     * Defines a facet of a 3D-Element as a contact surface for
     * master contact surfaces
     * adapted from face2D originally written by riccardo.
     */
    class SlaveContactPoint3D : public Condition
    {
        
        public:  
	typedef Condition BaseType;
	typedef BaseType::EquationIdVectorType EquationIdVectorType;
	typedef BaseType::MatrixType LHS_ContributionType;
	

	/**
	* REMOVED
	*/

	// typedef PointerVectorSet<StressConditionType, IndexedObject> StressConditionContainerType;
	// typedef PointerVectorSet< EquationIdVectorType, IndexedObject> EquationIdVectorContainerType;
	// typedef PointerVectorSet< LHS_ContributionType, IndexedObject> LHS_ContainerType;

	typedef Condition::GeometryType::Pointer PointerGeometryType; 

	// typedef BaseType::SecondaryCondition SecondaryConditionType;
	// typedef PointerVectorSet< SecondaryConditionType, IndexedObject> SecondaryConditionContainerType;
	  
            // Counted pointer of SlaveContactPoint3D
            KRATOS_CLASS_POINTER_DEFINITION(SlaveContactPoint3D);
            
            /** 
             * Default constructor.
             */
            SlaveContactPoint3D( IndexType NewId, GeometryType::Pointer pGeometry);
            
	    SlaveContactPoint3D( IndexType NewId, NodesArrayType const& ThisNodes);
	    
            SlaveContactPoint3D( IndexType NewId, GeometryType::Pointer pGeometry, 
                           PropertiesType::Pointer pProperties);
			   
			   
	    SlaveContactPoint3D(IndexType NewId, NodesArrayType& ThisNode);
	    
	    
            /**
             * Destructor.
             */
            virtual ~SlaveContactPoint3D();
      
            /**
             * Operations.
             */
            Condition::Pointer Create( IndexType NewId, 
                                       NodesArrayType const& ThisNodes,  
                                       PropertiesType::Pointer pProperties) const;
            /**
                                        * Turns back information on whether it is a contact type condition
             */
//             int IsContactType();
            
//             Matrix TangentialVectors( GeometryType::CoordinatesArrayType& rPoint );
            
//             Vector NormalVector( GeometryType::CoordinatesArrayType& rPoint );
            
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
            bool ClosestPoint( 
                    GeometryType::CoordinatesArrayType& rResultGlobal,
                    GeometryType::CoordinatesArrayType& rResultLocal,
                    const GeometryType::CoordinatesArrayType& rCandidateGlobal,
                    const GeometryType::CoordinatesArrayType& rSlaveContactGlobalPoint
                                                  );
            
             /**
              * applies the contact stress from a dedicated slave condition.
              * @param coords the coordinates of the slave condition's partner point on the 
              * current master condition in local coordinates
              * @param Stress the value of the current contact stress in current point
              * @param Weight the integration weight in slave element's integration point
              * @param dA the differential area in slave element's integration point 
              * @param NormalDirection the normal vector on current master surface 
              */
            void AddContactStress( Vector coords, const double Stress, double Weight, double dA,
                                   const Vector NormalDirection, double Penalty, double Gap,
                                   EquationIdVectorType SlaveEquationId, Vector SlaveNcontainer );
            
            
            /**
             * calculates the cross-diagonal terms of the system stiffness matrix
             * resulting from interaction between contact elements
             * (REMOVED)
             */
//             void CalculateCrossElementarySystemContributions( SecondaryConditionContainerType& SecondaryConditions, EquationIdVectorType& PrimaryEquationId, ProcessInfo& rCurrentProcessInfo );
            
            /**
             * calculates the stiffness matrix contributions resulting from linearization
             * of normal vector
             * (REMOVED)
             */
//             void CalculateNormalLinearizationElementarySystemContributions( SecondaryConditionContainerType& SecondaryConditions, EquationIdVectorType& PrimaryEquationId, ProcessInfo& rCurrentProcessInfo );
           
            
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
            
	   
	    
            /**
             * Turn back information as a string.
             * (DEACTIVATED)
             */
            //std::string Info();
      
            /**
             * Print information about this object.
             * (DEACTIVATED)
             */
            //virtual void PrintInfo(std::ostream& rOStream) const;

            /**
             * Print object's data.
             * (DEACTIVATED)
             */
            //virtual void PrintData(std::ostream& rOStream) const;
      
        protected:
	  
	//BaseType::Pointer mpMasterContact; 
        
        
        private:
            void CalculateAll( MatrixType& rLeftHandSideMatrix, 
                               VectorType& rRightHandSideVector,
                               ProcessInfo& rCurrentProcessInfo,
                               bool CalculateStiffnessMatrixFlag,
                               bool CalculateResidualVectorFlag);
            
            void CalculateAndAddKc( Matrix& K,
                                    const Vector& N,
                                    double pressure,
                                    double weight,
                                    double penalty,
                                    Vector v
                                  );
            
            void CalculateAndAdd_PressureForce( Vector& residualvector,
                    const Vector& N,
                    Vector& v3,
                    double pressure,
                    double weight,
                    double DetJ ); 
            
		    
		///@}
		///@name Serialization
		///@{	
		friend class Serializer;

		// A private default constructor necessary for serialization 
		SlaveContactPoint3D(){}; 

		void save(Serializer& rSerializer) const
		{
		   KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );   
		}

		virtual void load(Serializer& rSerializer)
		{
		  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
		}
		    
            /**
             * Assignment operator.
             * (DEACTIVATED)
             */
            //SlaveContactPoint3D& operator=(const SlaveContactPoint3D& rOther);
            
            /**
             * Copy constructor.
             * (DEACTIVATED)
             */
            //SlaveContactPoint3D(const SlaveContactPoint3D& rOther);
            
            //member variables
            /**
             * REMOVED
             */
//             StressConditionContainerType::Pointer mpStressConditions;
    }; // Class SlaveContactPoint3D 
}  // namespace Kratos.

#endif // KRATOS_CONTACT_FACE_3D_CONDITION_H_INCLUDED  defined 
