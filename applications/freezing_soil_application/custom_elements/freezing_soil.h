/*
==============================================================================
KratosR1StructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-PhysCcs Finite Element AnalysCs
VersCon 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo RossC, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossC@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-UniversCty Bochum, Institute for Structural Mechanics, Germany


PermissCon is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limCtation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissCble
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permissCon  notice  shall  be
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
//   RevisCon:            $RevisCon: 1.1 $
//
//


#if !defined(KRATOS_FREEZING_SOIL_INCLUDED )
#define  KRATOS_FREEZING_SOIL_INCLUDED

// System includes

#include "boost/smart_ptr.hpp"

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

/// Short class definition.
/** Detail class definition.
*/

class FreezingSoil
            : public Element
{

    public:
///@name Type Definitions
///@{
        typedef GeometryData::IntegrationMethod IntegrationMethod;
/// Counted pointer of

        KRATOS_CLASS_POINTER_DEFINITION ( FreezingSoil );

///@}
///@name Life Cycle
///@{

/// Default constructor.
        FreezingSoil ( IndexType NewId, GeometryType::Pointer pGeometry );
        FreezingSoil ( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

/// Destructor.
        virtual ~FreezingSoil();


///@}
///@name Operators
///@{


///@}
///@name Operations
///@{
        IntegrationMethod GetIntegrationMethod();

        Element::Pointer Create ( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

        void Initialize();

	void ResetConstitutiveLaw();

        void CalculateLocalSystem ( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );
        void CalculateRightHandSide ( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );
        //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

        void EquationIdVector ( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo );
        void GetDofList ( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );
        void FinalizeSolutionStep ( ProcessInfo& CurrentProcessInfo );
        void InitializeSolutionStep ( ProcessInfo& CurrentProcessInfo );


        void CalculateOnIntegrationPoints ( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo );

        void GetValuesVector ( Vector& values, int Step );
        void GetValueOnIntegrationPoints ( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );
        void GetValueOnIntegrationPoints ( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );
        void GetValueOnIntegrationPoints ( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );
        void SetValueOnIntegrationPoints ( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );
        void SetValueOnIntegrationPoints ( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );
//         void GetPorosity ( Matrix arg1, double arg2, double arg3, double arg4 );
//         void GetWaterMass ( Matrix arg1, double arg2, double arg3, double arg4, int arg5 );

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
///@name Protected  Access
///@{
///@}
        ///@name Serialization
        ///@{

        friend class Serializer;

        // A private default constructor necessary for serialization
        FreezingSoil() {};

        virtual void save ( Serializer& rSerializer ) const
        {
            rSerializer.save ( "Name", "FreezingSoil" );
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer,  Element );
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer,  Element );
        }

        /**
         * This function provides the place to perform checks on the completeness of the input.
         * It is designed to be called only once (or anyway, not often) typically at the beginning
         * of the calculations, so to verify that nothing is missing from the input
         * or that no common error is found.
         * @param rCurrentProcessInfo
         */
        virtual int Check ( const ProcessInfo& rCurrentProcessInfo );


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
        Geometry< Node<3> >::Pointer  mThisGeometryOther;
        std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
        IntegrationMethod mThisIntegrationMethod;

        Vector mMaterialParameters, mGravity, mGravityDisp, mT0, mT0e, mp0, mp0e, mStrainVol0;
        double mn0, mkS, mgS, malphaS, mkappa0, mm,  mtstar, mrhoS0, mlambdaS, mcS, mrhoL0, mrhoC0, mSf, mTf, mkL, mkC, malphaL, malphaC, metaL, mlambdaL, mlambdaC, mcL, mcC, mgL, mgC;
        double mScaleU, mScaleP, mUnitRatio;
	double mK, mG, mTol;

        unsigned int mNodesDispMin, mNodesDispMax,  mNodesNumberDisp, mMatSizeU, mNumberU;
        unsigned int mNodesOtherMin, mNodesOtherMax, mNodesNumberOther, mMatSizeO, mNumberO;
        unsigned int mDimension, mMatSize, mAddIndexU;

        std::vector< Matrix > mInvJ0;
        Vector mDetJ0;
	bool mPrestressAssigned;

///@}
///@name Private Operators
///@{
        void CalculateAll ( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                            ProcessInfo& rCurrentProcessInfo,
                            bool CalculateStiffnessMatrixFlag,
                            bool CalculateResCdualVectorFlag );

        void DampMatrix ( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo );

        void AddInternalForcesToRHS1 ( Vector& Help_R_1, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& Nu, const Matrix& DNu_DX, Vector& stressVectorEff );
        void AddInternalForcesToRHS2 ( Vector& Help_R_2, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );
        void AddInternalForcesToRHS3 ( Vector& Help_R_3, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );

        void CalculateStiffnessMatrixUU ( Matrix& Help_K_UU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& Nu, const Matrix& DNu_DX, Matrix& CtanEff );
        void CalculateStiffnessMatrixUP ( Matrix& Help_K_UP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& Nu, const Matrix& DNu_DX, const Vector& No );
        void CalculateStiffnessMatrixUT ( Matrix& Help_K_UT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& Nu, const Matrix& DNu_DX, const Vector& No, Vector& DstressVectorEff_Dt );
        void CalculateStiffnessMatrixPU ( Matrix& Help_K_PU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNu_DX );
        void CalculateStiffnessMatrixPP ( Matrix& Help_K_PP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );
        void CalculateStiffnessMatrixPT ( Matrix& Help_K_PT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );
        void CalculateStiffnessMatrixTU ( Matrix& Help_K_TU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNu_DX );
        void CalculateStiffnessMatrixTP ( Matrix& Help_K_TP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );
        void CalculateStiffnessMatrixTT ( Matrix& Help_K_TT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );

        void CalculateDampingMatrixPU ( Matrix& Help_D_PU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNu_DX );
        void CalculateDampingMatrixPP ( Matrix& Help_D_PP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );
        void CalculateDampingMatrixPT ( Matrix& Help_D_PT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );
        void CalculateDampingMatrixTU ( Matrix& Help_D_TU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNu_DX );
        void CalculateDampingMatrixTP ( Matrix& Help_D_TP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );
        void CalculateDampingMatrixTT ( Matrix& Help_D_TT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX );
 
        /// +++++++++++++++++++++++++++++++++++++++++++
        double KnoneckerDelta ( int i, int j );
        Matrix GetBu ( const Matrix& DNu_DX );
        Vector Getu ( const Vector& Nu );
        Matrix GetGradu ( const Matrix& DNu_DX );
        double GetDivdotu ( const Matrix& DNu_DX );
        double Getp ( const Vector& No );
        double Getdotp ( const Vector& No );
        Vector GetGradp ( const Matrix& DNo_DX ); 
        double Gett ( const Vector& No );
        double Getdott ( const Vector& No );
        Vector GetGradt ( const Matrix& DNo_DX );
        /// +++++++++++++++++++++++++++++++++++++++++++
        //C1
        Matrix GetElasticTangent ( double t, int index );
        Vector GetstrainVector ( Matrix Grad_u );
        double GetstrainVol ( Matrix Grad_u );
        Matrix GetstressEff ( Matrix Grad_u, double p, double t );
        double GetLiquidSaturation ( double t, int index );
        double GetBulkModulus ( double t, int index );
        double GetShearModulus ( double t, int index );
        double GetBiotCoefficient ( double t, int index );
        double GetPorosity ( Matrix Grad_u, double p, double t, double t0, int index );
        double GetWaterDensity ( double p, double t, int index );
        double GetIceDensity ( double p, double t, int index );
        double GetWaterMass ( Matrix Grad_u, double p, double t, double t0, int index );
        double GetIceMass ( Matrix Grad_u, double p, double t, double t0, int index );
        double GetWaterAndIceMass ( Matrix Grad_u, double p, double t, double t0, int index );

        //C2
        Vector GetDarcyFlow ( Vector Grad_p, double p, double t, int index );
        Vector GetMechanicalDissipation ( Vector Grad_p, double p, double t, int index );

        //C3
        double GetSkeletonEntropy ( Matrix Grad_u, double p, double t, int index );
        double GetWaterEntropy ( double t, int index );
        double GetIceEntropy ( double t, int index );
        double GetEntropyChange ( double p, double t, int index );
        Vector GetHeatFlow ( Vector Grad_t, double t, int index );

        /// +++++++++++++++++++++++++++++++++++++++++++
        void Interpolate ( const Variable<double>& rVariable, const ProcessInfo& rCurrentProcessInfo );
        void Interpolate ( const Variable<Kratos::array_1d<double, 3> >& rVariable, const ProcessInfo& rCurrentProcessInfo );
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
///@name Un accessCble methods
///@{

/// AssCgnment operator.
//& operator=(const & rOther);

/// Copy constructor.
//(const & rOther);


///@}

}; // Class

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
  & rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
  const & rThis)
 {
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}*/
///@}

}  // namespace Kratos.

#endif // KRATOS_FREEZING_SOIL_INCLUDED defined 


