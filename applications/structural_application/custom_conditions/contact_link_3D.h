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
//   Date:                $Date: 2008-05-06 14:55:41 $
//   Revision:            $Revision: 1.5 $
//
//
#if !defined(KRATOS_CONTACT_LINK_3D_CONDITION_H_INCLUDED )
#define  KRATOS_CONTACT_LINK_3D_CONDITION_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "custom_conditions/master_contact_face_3D.h"
#include "custom_conditions/slave_contact_face_3D.h"


namespace Kratos
{
    /**
     * Contact link element for 3D contact problems.
     * This condition links two contact surfaces (one master and
     * one slave element) to a condition which may be assembled
     * at once by the builder.
     * Within this condition all system contribution with regard to
     * the gap function, the contact stress, the normal vectors
     * and their linearizations are calculated.
     */
    class ContactLink3D : public Condition
    {
        public:
            // Counted pointer of ContactLink3D
            KRATOS_CLASS_POINTER_DEFINITION(ContactLink3D);
            
            /** 
             * Default constructor.
             */
            ContactLink3D( IndexType NewId, GeometryType::Pointer pGeometry);
            
            ContactLink3D( IndexType NewId, GeometryType::Pointer pGeometry,  
                           PropertiesType::Pointer pProperties
                         );
            
            
            ContactLink3D( IndexType NewId, GeometryType::Pointer pGeometry,  
                                          PropertiesType::Pointer pProperties,
                                          Condition::Pointer Master, 
                                          Condition::Pointer Slave,
                                          Point<3>& MasterContactLocalPoint,
                                          Point<3>& SlaveContactLocalPoint,
                                          int SlaveIntegrationPointIndex
                         );
            /**
             * Destructor.
             */
            virtual ~ContactLink3D();
      
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

	    void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
            
            void EquationIdVector( EquationIdVectorType& rResult, 
                                   ProcessInfo& rCurrentProcessInfo);
            
            void GetDofList( DofsVectorType& ConditionalDofList,
                             ProcessInfo& CurrentProcessInfo);
            
      
            /**
             * Print information about this object.
             */
            virtual void PrintInfo(std::ostream& rOStream) const;

            /**
             * Print object's data.
             */
            virtual void PrintData(std::ostream& rOStream) const;
      
        protected:
        
        
        private:
            void CalculateAll( MatrixType& rLeftHandSideMatrix, 
                               VectorType& rRightHandSideVector,
                               ProcessInfo& rCurrentProcessInfo,
                               bool CalculateStiffnessMatrixFlag,
                               bool CalculateResidualVectorFlag);
            
            
            /////////////////////////
            ///// 
            //////////////////////////
            /*void CalculateAndAdd_RHS( 
                    Vector& residualvector,
                    const Vector& NMaster,
                    const Vector& NSlave,
                    const Vector& vMaster,
                    double normalStress,
                    double SlaveIntegrationWeight,
            double dASlave ); */
            
    void CalculateAndAdd_RHS( Vector& residualvector,
            const Vector& NMaster,
            const Vector& NSlave,
            const Vector& vMaster,
            const Matrix& T,
            const Vector& tangentialStresses,
            double Gap,
            double normalStress,
            double SlaveIntegrationWeight,
            double dASlave);
            
            
            
            
            /**
             * This function calculates updates the local and global coordinates 
             * of the master contact partner in order to follow the movement of
             * the slave surface along the master surface
             */
            void UpdateMasterLocalPoint( );

            
            Vector NormalVector( Condition::Pointer Surface, 
                                 const GeometryType::CoordinatesArrayType& LocalPoint );
            
            Matrix TangentialVectors( Condition::Pointer Surface,
                                      const GeometryType::CoordinatesArrayType& LocalPoint );
            
            Matrix TangentialVectorsGlobal( Condition::Pointer Surface,
                                      const GeometryType::CoordinatesArrayType& LocalPoint );
            
    	    Matrix TangentialVectors_inOrigin( Condition::Pointer Surface, 
                                             const GeometryType::CoordinatesArrayType& rPoint );

     	    Point<3>& GlobalCoordinates(Condition::Pointer Surface, Point<3>& rResult, Point<3> const& LocalCoordinates);


	   Vector GetRelativTangentialVelocity(Matrix& T);

	   Vector GetRelativVelocity();
            /**
             * Assignment operator.
             * (DEACTIVATED)
             */
            //ContactLink3D& operator=(const ContactLink3D& rOther);
            
            /**
             * Copy constructor.
             * (DEACTIVATED)
             */
            //ContactLink3D(const ContactLink3D& rOther);
            
            
            /**
             * private members
             */
            
//             Condition::Pointer mpSlave;
//             Condition::Pointer mpMaster;
//             Point<3> mMasterContactLocalPoint;
//             Point<3> mSlaveContactLocalPoint;
//             Point<3> mMasterContactGlobalPoint;
//             Point<3> mSlaveContactGlobalPoint;
    }; // Class ContactLink3D 
}  // namespace Kratos.

#endif // KRATOS_CONTACT_LINK_3D_CONDITION_H_INCLUDED  defined 
