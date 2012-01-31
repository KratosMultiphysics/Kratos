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
 *   Date:                $Date: 2008-01-25 08:37:31 $
 *   Revision:            $Revision: 1.10 $
 *
 * ***********************************************************/

#if !defined(KRATOS_CASM_H_INCLUDED )
#define  KRATOS_CASM_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{
/**
 * Defines the interface to CASM soil models
 */

class Casm : public ConstitutiveLaw
{

public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;
    typedef array_1d<double, 81> MaterialTensorType;
    /**
     * Counted pointer of UPCHMCM
     */
    typedef boost::shared_ptr<Casm> Pointer;
    /**
     * Life Cycle
     */
    Casm();

    /**
     * Destructor.
     */
    virtual ~Casm();

    /**
     * Clone function
     * will be called on initialization of the constitutive law
     */
    virtual boost::shared_ptr<BaseType> Clone() const
    {
        boost::shared_ptr<BaseType> p_clone ( new Casm() );
        return p_clone;
    }

    /**
     * Operators
     */
//    virtual template<class TVariableType> typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable) const
//    {
//     if( rThisVariable == DP_EPSILON )
//      return mEpsilon;
//     if( rThisVariable == INSITU_STRESS )
//     {
//      return mInSituStress;
//     }
//    }
//
//    virtual template<class TVariableType> bool Has(const TVariableType& rThisVariable) const
//    {
//     if( rThisVariable == DP_EPSILON )
//      return true;
//     if( rThisVariable == INSITU_STRESS )
//      return true;
//     return false;
//    }

    bool Has ( const Variable<double>& rThisVariable );
    bool Has ( const Variable<Vector>& rThisVariable );
    bool Has ( const Variable<Matrix>& rThisVariable );
    double& GetValue ( const Variable<double>& rThisVariable, double& rValue );
    Vector& GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue );
    Matrix& GetValue ( const Variable<Matrix>& rThisVariable, Matrix& rValue );
    void SetValue ( const Variable<double>& rThisVariable, const double& rValue,
                    const ProcessInfo& rCurrentProcessInfo );
    void SetValue ( const Variable<array_1d<double, 3> >& rThisVariable,
                    const array_1d<double, 3>& rValue, const ProcessInfo& rCurrentProcessInfo );
    void SetValue ( const Variable<Vector>& rThisVariable, const Vector& rValue,
                    const ProcessInfo& rCurrentProcessInfo );
    void SetValue ( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                    const ProcessInfo& rCurrentProcessInfo );


    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial ( const Properties& props,
                              const GeometryType& geom,
                              const Vector& ShapeFunctionsValues );

    void ResetMaterial ( const Properties& props, const GeometryType& geom, const Vector& ShapeFunctionsValues );

    /**
     * Calculates the constitutive matrix for a given strain vector
     * @param StrainVector the current vector of strains the constitutive
     * matrix is to be generated for
     * @param rResult Matrix the result will be stored in
     */
    void CalculateConstitutiveMatrix ( const Vector& StrainVector, Matrix& rResult );

    /**
     * Calculates the stresses for given strain state
     * @param StrainVector the current vector of strains
     * @param rResult the stress vector corresponding to the given strains
     */
    void CalculateStress ( const Vector& StrainVector, Vector& rResult );

    void InitializeSolutionStep ( const Properties& props,
                                  const GeometryType& geom, //this is just to give the array of nodes
                                  const Vector& ShapeFunctionsValues ,
                                  const ProcessInfo& CurrentProcessInfo );

    void FinalizeSolutionStep ( const Properties& props,
                                const GeometryType& geom, //this is just to give the array of nodes
                                const Vector& ShapeFunctionsValues ,
                                const ProcessInfo& CurrentProcessInfo );

    /**
     * Computes the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear strain measure is used)
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param CurrentProcessInfo current ProcessInfo instance
     * @param props the material's Properties object
     * @param geom the element's geometry
     * @param ShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
    void CalculateMaterialResponse ( const Vector& StrainVector,
                                     const Matrix& DeformationGradient,
                                     Vector& StressVector,
                                     Matrix& AlgorithmicTangent,
                                     const ProcessInfo& CurrentProcessInfo,
                                     const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     bool CalculateStresses = true,
                                     int CalculateTangent = true,
                                     bool SaveInternalVariables = true
                                   );

    /**
     * Calculates the cauchy stresses. For a given deformation and stress state
     * the cauchy stress vector is calculated
     * @param Cauchy_StressVector the result vector of cauchy stresses (will be overwritten)
     * @param F the current deformation gradient
     * @param PK2_StressVector the current second Piola-Kirchhoff-Stress vector
     * @param GreenLagrangeStrainVector the current Green-Lagrangian strains
     */
    void CalculateCauchyStresses ( Vector& Cauchy_StressVector,
                                   const Matrix& F,
                                   const Vector& PK2_StressVector,
                                   const Vector& GreenLagrangeStrainVector );


    /**
     * calculates the current stress and the material tangent
     * NOTE: there are two versions of this function: one for a matrix representation
     * of the material tensor and one for a tensorial formulation. Each ConstitutiveLaw
     * HAS TO IMPLEMENT both of them (for convenience, there are conversation functions
     * available in MathUtils for either of them)
     * @param StressTensor the calculated stress tensor
     * @param StrainTensor the given strain tensor
     * @param algorithmicTangent the 4th order algorithmic tangent tensor
     */
    void CalculateStressAndTangentMatrix ( Matrix& StressTensor,
                                           const Matrix& StrainTensor,
                                           MaterialTensorType& algorithmicTangent );

    /**
     * calculates the current stress and the material tangent
     * NOTE: there are two versions of this function: one for a matrix representation
     * of the material tensor and one for a tensorial formulation. Each ConstitutiveLaw
     * HAS TO IMPLEMENT both of them (for convenience, there are conversation functions
     * available in MathUtils for either of them)
     * @param StressVector the calculated stress vector
     * @param StrainVector the given strain vector
     * @param algorithmicTangent the calculated algorithmic tangent matrix
     */
    void CalculateStressAndTangentMatrix ( Vector& StressVector,
                                           const Vector& StrainVector,
                                           Matrix& algorithmicTangent );

    /**
     * converts a strain vector styled variable into its form, which the
     * deviatoric parts are no longer multiplied by 2
     */
//             void Calculate( const Variable<Matrix >& rVariable,
//                             Matrix& rResult, const ProcessInfo& rCurrentProcessInfo );
//
    void Calculate ( const Variable<Vector >& rVariable,
                     Vector& rResult, const ProcessInfo& rCurrentProcessInfo );



    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

protected:

    /**
     * Member Variables
     */

    double Em;
    double Km;
    double Gm;
    double mMth;
    double lnr;
    double pk;
    double qk;
    double pck;
    double dEpspk;
    double dEpsqk;
    double phiCs;
    double theta;

    double dEpsPkold;
    double unit4thSym3D[6][6];
    double unit2nd3D[6];
    double mTangentOperator[36];
    double mElasticTangent[36];
    Matrix mConsistentTangent;
    Matrix mCtangent;
    Vector mOldStrain;
    Vector mCurrentStrainInc;
    Vector mOldStress;
    Vector mCurrentStress;
    Vector mCurrentStrain;
    double mOldHistoryVariables[6];
    double mCurrentHistoryVariables[6];
    double mModelData[17];
    double mCurrentTime;
    Vector mOldPlasticStrains;
    Vector mCurrentPlasticStrains;
    Vector mInSituStress;
    Vector mInSituStressOriginal;
    Vector mStress;
    Vector mStrain;
    Vector mMaterialParameters;
    int miMod;
    int mIsUndr;
    int mNhv;
    int mNpar;
    int mIplastic;
    double mBulkW;
    double mCurrentSuction;
    double mOldSuction;
    double mCurrentWaterContent;
    bool isYielded;

    /**
     * Assignment operator.
     */
    //Casm& operator=(const Casm& rOther);
    /**
     * Copy constructor.
     */
    Matrix InverseC ( Matrix& InvC );
    //Casm(const Casm& rOther);

private:

    friend class Serializer;
    
    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw );
//         rSerializer.save( "mOldHistoryVariables", mOldHistoryVariables );
//         rSerializer.save( "mCurrentHistoryVariables", mCurrentHistoryVariables );
        rSerializer.save( "mOldPlasticStrains", mOldPlasticStrains );
        rSerializer.save( "mCurrentPlasticStrains", mCurrentPlasticStrains );
    }
    
    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw );
//         rSerializer.load( "mOldHistoryVariables", mOldHistoryVariables );
//         rSerializer.load( "mCurrentHistoryVariables", mCurrentHistoryVariables );
        rSerializer.load( "mOldPlasticStrains", mOldPlasticStrains );
        rSerializer.load( "mCurrentPlasticStrains", mCurrentPlasticStrains );
    }


    void returnMapping (double pTr, double qTr);

    void formPlasticPotentialDerivatives (Vector& dG);

    void formResidualDerivativeMatrix(Matrix& dR, Vector& dG);

    void formResidualVector(Vector& R, Vector& dG);

    void formYieldFunctionDerivatives(Vector& dF);

    static void solveByGauss(Matrix& A, Vector& B);

    static void decompose(Matrix& A, int N);

    void CalculateElasticMatrix(Matrix& C, const double E, const double NU);

    void pressureDependentElasticTangent(Matrix& mCtangent);

    void getM(const Vector& devStr, double& Mth);

    void getFactors(double qTr, Vector& factor);

    void CalculateModules(const double& p, double& Km, double& Gm);

    void CalculateUnit4thSym3D();

    void CalculateUnit2nd3D();



}; // Class Casm
}  // namespace Kratos.
#endif // KRATOS_CASM_H_INCLUDED  defined 
