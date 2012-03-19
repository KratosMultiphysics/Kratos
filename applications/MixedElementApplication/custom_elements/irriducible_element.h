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
//   Last Modified by:    $Author: hurga $
//   Date:                $Date: 2007-10-18 16:23:41 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_IRREDUCIBLE_ELEMENT_H_INCLUDED )
#define  KRATOS_IRREDUCIBLE_ELEMENT_H_INCLUDED



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

    class IrriducibleElement
                : public Element
    {

        public:
            ///@name Type Definitions
            ///@{
            ///Reference type definition for constitutive laws
            typedef ConstitutiveLaw ConstitutiveLawType;
            ///Pointer type for constitutive laws
            typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
            ///Type definition for integration methods
            typedef GeometryData::IntegrationMethod IntegrationMethod;

            /// Counted pointer of IrriducibleElement
            KRATOS_CLASS_POINTER_DEFINITION( IrriducibleElement );

            ///@}
            ///@name Life Cycle
            ///@{

            /// Default constructor.
            IrriducibleElement( IndexType NewId, GeometryType::Pointer pGeometry );
            IrriducibleElement( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

            /// Destructor.
            virtual ~IrriducibleElement();

            ///@}
            ///@name Operators
            ///@{
            ///@}
            ///@name Operations
            ///@{

            Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

            void Initialize();

            void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

            void CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

            //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

            void EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo );

            void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );

            void InitializeSolutionStep( ProcessInfo& CurrentProcessInfo );

            void FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo );

            void MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo );

            void DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo );

            void CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo );

            void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo );

            void CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo );

            void SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void GetValuesVector( Vector& values, int Step = 0 );
            void GetFirstDerivativesVector( Vector& values, int Step = 0 );
            void GetSecondDerivativesVector( Vector& values, int Step = 0 );


            void Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo );

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

            /**
            * Calculates the elemental contributions
            * \f$ K^e = w\,B^T\,D\,B \f$ and
            * \f$ r^e \f$
            */
            virtual void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                       ProcessInfo& rCurrentProcessInfo,
                                       bool CalculateStiffnessMatrixFlag,
                                       bool CalculateResidualVectorFlag );
            ///@}
            ///@name Protected Operations
            ///@{
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
            /**
             * Container for constitutive law instances on each integration point
             */
            std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

            double mArea0;
            Matrix mDN_DX;

            
            ///@}
            ///@name Private Operators
            ///@{
 
            void InitializeVariables();

            virtual void InitializeMaterial();

            void CalculateB( Matrix& B,
                             Matrix& DN_DX,
                             unsigned int StrainSize );

            std::string Info() const;

            ///@}
            ///@name Private Operations
            ///@{

            ///@}
            ///@name Private  Access
            ///@{
            ///@}
	    
	    
	    ///@}
	    ///@name Serialization
	    ///@{	
	    friend class Serializer; 

	    // A private default constructor necessary for serialization 
	    IrriducibleElement(){}

	    virtual void save(Serializer& rSerializer) const
	    {
	       KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
	    }

	    virtual void load(Serializer& rSerializer)
	    {
	       KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
	    }
	    
            ///@name Private Inquiry
            ///@{
            ///@}
            ///@name Un accessible methods
            ///@{
            /// Assignment operator.
            //IrriducibleElement& operator=(const IrriducibleElement& rOther);
            /// Copy constructor.
            //IrriducibleElement(const IrriducibleElement& rOther);
            ///@}

    }; // Class IrriducibleElement

    ///@}
    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    /// input stream function
    /*  inline std::istream& operator >> (std::istream& rIStream,
    IrriducibleElement& rThis);
    */
    /// output stream function
    /*  inline std::ostream& operator << (std::ostream& rOStream,
    const IrriducibleElement& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }*/
    ///@}

}  // namespace Kratos.
#endif // KRATOS_IRREDUCIBLE_ELEMENT_H_INCLUDED  defined
