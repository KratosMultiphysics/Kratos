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
//   Date:                $Date: 2009-01-14 09:30:38 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_KINEMATIC_LINEAR_INCLUDED )
#define  KRATOS_KINEMATIC_LINEAR_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
// #include "includes/constitutive_law_new.h"


namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name  Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.
    /** Detail class definition.
     */

    class KinematicLinear
                : public Element
    {

        public:
            ///@name Type Definitions
            ///@{
            typedef GeometryData::IntegrationMethod IntegrationMethod;

            typedef ConstitutiveLaw ConstitutiveLawType;

            typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
            /// Counted pointer of KinematicLinear

            KRATOS_CLASS_POINTER_DEFINITION( KinematicLinear );

            ///@}
            ///@name Life Cycle
            ///@{

            /// Default constructor.
            KinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry );
            KinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

            /// Destructor.
            virtual ~KinematicLinear();


            ///@}
            ///@name Operators
            ///@{


            ///@}
            ///@name Operations
            ///@{
            IntegrationMethod GetIntegrationMethod();

            Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

            void Initialize();

            void ResetConstitutiveLaw();

            void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

            void CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

            void EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo );

            void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );

            void MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo );

            void DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo );

            void FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo );

            void InitializeSolutionStep( ProcessInfo& CurrentProcessInfo );

            void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void SetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                              std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void CalculateOnIntegrationPoints( const Variable<double >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo );

            void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo );

            void CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& Output, const ProcessInfo& rCurrentProcessInfo );

            void GetValuesVector( Vector& values, int Step = 0 );

            void GetFirstDerivativesVector( Vector& values, int Step = 0 );

            void GetSecondDerivativesVector( Vector& values, int Step = 0 );


            ///@}
            ///@name Access
            ///@{


            ///@}
            ///@name Inquiry
            ///@{


            ///@}
            ///@name Input and output
            ///@{

            /// Turn back information as a string.
//      virtual String Info() const;

            /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

            /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


            ///@}
            ///@name Friends
            ///@{


            ///@}

        protected:
            ///@name Protected static Member Variables
            ///@{


            ///@}
            ///@name Protected member Variables
            ///@{


            ///@}
            ///@name Protected Operators
            ///@{


            ///@}
            ///@name Protected Operations
            ///@{

	    ///@}
	    ///@name Serialization
	    ///@{	
	    friend class Serializer; 
	    
	    // A private default constructor necessary for serialization 
	    KinematicLinear(){}

	    virtual void save(Serializer& rSerializer) const
	    {
	       rSerializer.save("Name","KinematicLinear");
	       KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
	    }

	    virtual void load(Serializer& rSerializer)
	    {
	       KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
	    }
	      

            ///@}
            ///@name Protected  Access
            ///@{


            ///@}
            ///@name Protected Inquiry
            ///@{


            ///@}
            ///@name Protected LifeCycle
            ///@{


            ///@}

        private:
            ///@name Static Member Variables
            ///@{

            ///@}
            ///@name Member Variables
            ///@{
            std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

            IntegrationMethod mThisIntegrationMethod;

            double mTotalDomainInitialSize;


            std::vector< Matrix > mInvJ0;
            Vector mDetJ0;

            bool mfirst_time_step;
            
            bool mIsInitialized;

            Matrix mInitialDisp;
            ///@}
            ///@name Private Operators
            ///@{
            /** K += weight*Btrans*D*B */
            void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                               ProcessInfo& rCurrentProcessInfo,
                               bool CalculateStiffnessMatrixFlag,
                               bool CalculateResidualVectorFlag );

            void CalculateBodyForces(
                Vector& BodyForce,
                const ProcessInfo& CurrentProcessInfo
            );

            void InitializeVariables();

            void InitializeMaterial();

            void InitializeMaterial( std::vector<std::vector<Matrix> >& C );

            void CalculateAndAddExtForceContribution(
                const Vector& N,
                const ProcessInfo& CurrentProcessInfo,
                Vector& BodyForce,
                VectorType& mResidualVector,
                double weight );

            //************************************************************************************
            //************************************************************************************
            //************************************************************************************
            //************************************************************************************

            //CALCULATE FORCEVECTORS DISPLACEMENT

            void AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, double Weight, double detJ );

            void AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ );

            void CalculateStiffnesMatrix( Matrix& K, const
                                          Matrix& tan_C, const Matrix& B_Operator, double Weight, double detJ );

            void CalculateStressAndTangentialStiffness( Vector& StressVector, Matrix& tanC_U,
                    Vector& StrainVector, const Matrix& B_Operator, int PointNumber, const ProcessInfo& CurrentProcessInfo );

            void CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector );

            Matrix GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, int PointNumber );

            void CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX );

	      
            ///@}
            ///@name Private Operations
            ///@{


            ///@}
            ///@name Private  Access
            ///@{


            ///@}
            ///@name Private Inquiry
            ///@{


            ///@}
            ///@name Un accessible methods
            ///@{

            /// Assignment operator.
            //KinematicLinear& operator=(const KinematicLinear& rOther);

            /// Copy constructor.
            //KinematicLinear(const KinematicLinear& rOther);


            ///@}

    }; // Class KinematicLinear

    ///@}
    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    /*  inline std::istream& operator >> (std::istream& rIStream,
                       KinematicLinear& rThis);
    */
    /// output stream function
    /*  inline std::ostream& operator << (std::ostream& rOStream,
                       const KinematicLinear& rThis)
               {
                       rThis.PrintInfo(rOStream);
                       rOStream << std::endl;
                       rThis.PrintData(rOStream);

                       return rOStream;
    }*/
    ///@}

}  // namespace Kratos.

#endif // KRATOS_UNSATURATED_SOILS_ELEMENT_INCLUDED defined 


