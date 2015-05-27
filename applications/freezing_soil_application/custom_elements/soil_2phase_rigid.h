/*
==============================================================================
KratosR1StructuralApplication
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
//   Last Modified by:    $Author: mengmeng $
//   Date:                $Date: 2009-02-23 16:02:32 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_SOIL_2PHASE_RIGID_INCLUDED )
#define  KRATOS_SOIL_2PHASE_RIGID_INCLUDED

// System includes

#include "boost/smart_ptr.hpp"

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"



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

    class Soil2PhaseRigid
                : public Element
    {
        public:
///@name Type Definitions
///@{
            typedef GeometryData::IntegrationMethod IntegrationMethod;
/// Counted pointer of Soil2PhaseRigid

            KRATOS_CLASS_POINTER_DEFINITION( Soil2PhaseRigid );

///@}
///@name Life Cycle
///@{

/// Default constructor.
            Soil2PhaseRigid( IndexType NewId, GeometryType::Pointer pGeometry );
            Soil2PhaseRigid( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

/// Destructor.
            virtual ~Soil2PhaseRigid();


///@}
///@name Operators
///@{


///@}
///@name Operations
///@{
            IntegrationMethod GetIntegrationMethod();

            Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

            void Initialize();

            void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

            void CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

            //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

            void EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo );

            void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );

            void FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo );

            void CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo );

            void GetValuesVector( Vector& values, int Step );

            void GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );

            void SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );


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
            Soil2PhaseRigid(){};

            virtual void save( Serializer& rSerializer ) const
            {
                rSerializer.save( "Name", "Soil2PhaseRigid" );
                KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer,  Element );
            }

            virtual void load( Serializer& rSerializer )
            {
                KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer,  Element );
            }

            /**
             * This function provides the place to perform checks on the completeness of the input.
             * It is designed to be called only once (or anyway, not often) typically at the beginning
             * of the calculations, so to verify that nothing is missing from the input
             * or that no common error is found.
             * @param rCurrentProcessInfo
             */
            virtual int Check( const ProcessInfo& rCurrentProcessInfo );


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
            IntegrationMethod mThisIntegrationMethod;

            double mPorosity;
            double mScale, mUnitRatio;

            unsigned int mNodesMin;
            unsigned int mNodesMax;
            unsigned int mDimension;
            unsigned int mNodesNumber;
            unsigned int mMatSize;

            Vector mInitialTemperature, mMaterialParameters;

            std::vector< Matrix > mInvJ0;
            Vector mDetJ0;


///@}
///@name Private Operators
///@{
            void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                               ProcessInfo& rCurrentProcessInfo,
                               bool CalculateStiffnessMatrixFlag,
                               bool CalculateResidualVectorFlag );

            void DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo );

            void AddInternalForcesToRHSP( Vector& Help_R_P, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t );

            void AddInternalForcesToRHST( Vector& Help_R_T, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t );
//
            void CalculateStiffnessMatrixPP( Matrix& Help_K_PP, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t );
            void CalculateStiffnessMatrixPT( Matrix& Help_K_PT, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t );
            void CalculateStiffnessMatrixTP( Matrix& Help_K_TP, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t );
            void CalculateStiffnessMatrixTT( Matrix& Help_K_TT, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t );

            void CalculateDampingMatrixPP( Matrix& Help_D_PP, const Vector& N, double weight, double DetJ, double p, double t );
            void CalculateDampingMatrixPT( Matrix& Help_D_PT, const Vector& N, double weight, double DetJ, double p, double t );
            void CalculateDampingMatrixTP( Matrix& Help_D_TP, const Vector& N, double weight, double DetJ, double p, double t );
            void CalculateDampingMatrixTT( Matrix& Help_D_TT, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t );

            void AssembleTimeSpaceStiffnessFromDampSubMatrices( MatrixType& rLeftHandSideMatrix,
                    const Matrix& D_PP,
                    const Matrix& D_PT,
                    const Matrix& D_TP,
                    const Matrix& D_TT );

            void AssembleTimeSpaceStiffnessFromStiffSubMatrices( MatrixType& rLeftHandSideMatrix,
                    const Matrix& K_PP,
                    const Matrix& K_PT,
                    const Matrix& K_TP,
                    const Matrix& K_TT );

            void AssembleTimeSpaceRHSFromSubVectors( VectorType& rRightHandSideVector,
                    const Vector& R_P,
                    const Vector& R_T );
//
            double Getp( const Vector& N );
            Vector Getgradp( const Matrix& DN_DX );
            double Getdotp( const Vector& N );

            double Gett( const Vector& N );
            Vector Getgradt( const Matrix& DN_DX );
            double Getdott( const Vector& N );

            double GetrhoL( double p, double t );
            double GetDrhoLDp( int n, double p, double t );
            double GetDrhoLDt( int n, double p, double t );
            double GetD2rhoLDpDt( double p, double t );
            double GetD3rhoLDpDt2( double p, double t );
            double GetD3rhoLDp2Dt( double p, double t );

            Vector GetvL( const Matrix& DN_DX, double p, double t );
            Vector GetDvLDrhoL( const Matrix& DN_DX, double p, double t );
            double GetDvLDgradp( double p, double t );
            Vector GetDvLDp( const Matrix& DN_DX, double p, double t );
            Vector GetDvLDt( const Matrix& DN_DX, double p, double t );

            double GetDpsiLDp( const Vector& N, int n, double p, double t );
            double GetDpsiLDt( const Vector& N, int n, double p, double t );

            double GetD2psiLDpDt( const Vector& N, double p, double t );
            double GetD3psiLDp2Dt( const Vector& N, double p, double t );
            double GetD3psiLDpDt2( const Vector& N, double p, double t );


            double GetDeSDt( int n );

            double GetDeLDt( const Vector& N, int n, double p, double t );
            double GetDeLDp( const Vector& N, int n, double p, double t );
            double GetD2eLDpDt( const Vector& N, double p, double t );

            Vector Getq( const Matrix& DN_DX );
            double GetDqDgradt();

            Vector GethatpL( const Matrix& DN_DX, double p, double t );
            double GetDhatpLDvL();



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
//Soil2PhaseRigid& operator=(const Soil2PhaseRigid& rOther);

/// Copy constructor.
//Soil2PhaseRigid(const Soil2PhaseRigid& rOther);


///@}

    }; // Class Soil2PhaseRigid

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    /*  inline std::istream& operator >> (std::istream& rIStream,
      Soil2PhaseRigid& rThis);
    */
/// output stream function
    /*  inline std::ostream& operator << (std::ostream& rOStream,
      const Soil2PhaseRigid& rThis)
     {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_SOIL_2PHASE_RIGID_INCLUDED defined 


