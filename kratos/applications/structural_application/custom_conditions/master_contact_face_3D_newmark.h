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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2007-03-13 15:01:51 $
*   Revision:            $Revision: 1.6 $
*
* ***********************************************************/

#if !defined(KRATOS_MASTER_CONTACT_FACE_3D_NEWMARK_CONDITION_H_INCLUDED )
#define  KRATOS_MASTER_CONTACT_FACE_3D_NEWMARK_CONDITION_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


namespace Kratos
{
    /**
     * stress condition container
     * (REMOVED)
     */
//     struct StressConditionType:IndexedObject
//     {
//         Vector coords;
//         Vector NormalDirection;
//         double Stress;
//         double Weight;
//         double dA;
//         double Penalty;
//         double Gap;
//         Element::EquationIdVectorType SlaveEquationId;
//         Vector SlaveNcontainer;
//         
//     };
    typedef Condition BaseType;
    /**
     * REMOVED
     */
//     typedef PointerVectorSet<StressConditionType, IndexedObject> StressConditionContainerType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef PointerVectorSet< EquationIdVectorType, IndexedObject>
            EquationIdVectorContainerType;
    typedef BaseType::MatrixType LHS_ContributionType;
    typedef PointerVectorSet< LHS_ContributionType, IndexedObject> LHS_ContainerType;
    
//     typedef BaseType::SecondaryCondition SecondaryConditionType;
//     typedef PointerVectorSet< SecondaryConditionType, IndexedObject> SecondaryConditionContainerType;
    
        
    
    /**
     * Contact surface element for 3D contact problems.
     * Defines a facet of a 3D-Element as a contact surface for
     * master contact surfaces
     * adapted from face2D originally written by riccardo.
     */
    class MasterContactFace3DNewmark : public Condition
    {
        
        public:
            // Counted pointer of MasterContactFace3DNewmark
            KRATOS_CLASS_POINTER_DEFINITION(MasterContactFace3DNewmark);
            
            /** 
             * Default constructor.
             */
            MasterContactFace3DNewmark( IndexType NewId, GeometryType::Pointer pGeometry);
            MasterContactFace3DNewmark( IndexType NewId, GeometryType::Pointer pGeometry, 
                           PropertiesType::Pointer pProperties);
            /**
             * Destructor.
             */
            virtual ~MasterContactFace3DNewmark();
      
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
            
//             Matrix MasterContactFace3DNewmark::TangentialVectors( GeometryType::CoordinatesArrayType& rPoint );
            
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
	      MasterContactFace3DNewmark(){}; 

	      virtual void save(Serializer& rSerializer) const
	      {
	      rSerializer.save("Name","MasterContactFace3DNewmark");
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
            //MasterContactFace3DNewmark& operator=(const MasterContactFace3DNewmark& rOther);
            
            /**
             * Copy constructor.
             * (DEACTIVATED)
             */
            //MasterContactFace3DNewmark(const MasterContactFace3DNewmark& rOther);
            
            //member variables
            /**
             * REMOVED
             */
//             StressConditionContainerType::Pointer mpStressConditions;
    }; // Class MasterContactFace3DNewmark 
}  // namespace Kratos.

#endif // KRATOS_CONTACT_FACE_3D_CONDITION_H_INCLUDED  defined 
