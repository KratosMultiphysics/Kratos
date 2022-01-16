// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Hoang-Giang Bui
//

#if !defined(KRATOS_LINEAR_ELASTIC_H_INCLUDED )
#define  KRATOS_LINEAR_ELASTIC_H_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

template<int TType>
struct LinearElastic_Helper
{
    /**
     * returns the size of the strain vector of the current constitutive law
     */
    static constexpr unsigned int StrainSize()
    {
        #if defined(__clang__)
        throw std::logic_error("Calling unimplemented function");
        return 0;
        #elif defined(__GNUC__) || defined(__GNUG__)
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #else
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #endif
    }

    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */
    static void CalculateElasticMatrix( Matrix& C, const double& E, const double& NU )
    {
        #if defined(__clang__)
        throw std::logic_error("Calling unimplemented function");
        #elif defined(__GNUC__) || defined(__GNUG__)
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #else
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #endif
    }

    /**
     * compute the stress given the strain in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the stress
     */
    static void Apply( const Vector& StrainVector, Vector& StressVector, const double& E, const double& NU )
    {
        #if defined(__clang__)
        throw std::logic_error("Calling unimplemented function");
        #elif defined(__GNUC__) || defined(__GNUG__)
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #else
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #endif
    }

    /**
     * calculates the inversed linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the inversed linear elastic constitutive matrix
     */
    static void CalculateInversedElasticMatrix( Matrix& C, const double& E, const double& NU )
    {
        #if defined(__clang__)
        throw std::logic_error("Calling unimplemented function");
        #elif defined(__GNUC__) || defined(__GNUG__)
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #else
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #endif
    }

    /**
     * compute the strain given the stress in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the stress
     */
    static void ApplyInverse( const Vector& StressVector, Vector& StrainVector, const double& E, const double& NU )
    {
        #if defined(__clang__)
        throw std::logic_error("Calling unimplemented function");
        #elif defined(__GNUC__) || defined(__GNUG__)
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #else
        KRATOS_THROW_ERROR(std::logic_error, "Calling unimplemented function", __FUNCTION__)
        #endif
    }
};

/**
 * Defines the general linear elastic material law, applicable for plane_strain, plane_stress, axisymmetric and 3d
 * Because this constitutive_law is defined very generally, it can be very memory-consuming
 */
template<int TType>
class LinearElastic : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ConstitutiveLaw BaseType;

    /**
     * Counted pointer of LinearElastic
     */
    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic);

    /**
     * Life Cycle
     */
    /**
     * Default constructor.
     */
    LinearElastic();

    /**
     * Destructor.
     */
    virtual ~LinearElastic() {}

    /**
     * Operators
     */

    /**
     * Operations
     */

    BaseType::Pointer Clone() const override
    {
        BaseType::Pointer p_clone( new LinearElastic<TType>() );
        return p_clone;
    }

    ConstitutiveLaw::StrainMeasure GetStrainMeasure() final
    {
        return StrainMeasure_Infinitesimal;
    }

    ConstitutiveLaw::StressMeasure GetStressMeasure() final
    {
        return StressMeasure_Cauchy;
    }

    bool Has( const Variable<int>& rThisVariable ) override;
    bool Has( const Variable<double>& rThisVariable ) override;
    bool Has( const Variable<array_1d<double, 3> >& rThisVariable ) override;
    bool Has( const Variable<Vector>& rThisVariable ) override;
    bool Has( const Variable<Matrix>& rThisVariable ) override;

    int& GetValue( const Variable<int>& rThisVariable, int& rValue ) override;
    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;
    array_1d<double, 3>& GetValue( const Variable<array_1d<double, 3> >& rThisVariable, array_1d<double, 3>& rValue ) override;
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue ) override;
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue ) override;

    double& CalculateValue( Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue ) override;
    array_1d<double, 3>& CalculateValue( Parameters& rParameterValues, const Variable<array_1d<double, 3> >& rThisVariable, array_1d<double, 3>& rValue ) override;
    Vector& CalculateValue( Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue ) override;
    Matrix& CalculateValue( Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue ) override;

    void SetValue( const Variable<int>& rThisVariable, const int& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<array_1d<double, 3 > >& rThisVariable,
                   const array_1d<double, 3 > & rValue, const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues ) override;

    void ResetMaterial();

    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void ResetMaterial( const Properties& props,
                        const GeometryType& geom,
                        const Vector& ShapeFunctionsValues ) override;

    /**
     * to be called at the beginning of each solution step
     * (e.g. from Element::InitializeSolutionStep)
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
     void InitializeSolutionStep( const Properties& props,
                                  const GeometryType& geom, //this is just to give the array of nodes
                                  const Vector& ShapeFunctionsValues,
                                  const ProcessInfo& CurrentProcessInfo ) override;

    /**
     * to be called at the beginning of each iteration
     * (e.g. from Element::InitializeNonLinearIteration)
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void InitializeNonLinearIteration( const Properties& props,
                                       const GeometryType& geom, //this is just to give the array of nodes
                                       const Vector& ShapeFunctionsValues,
                                       const ProcessInfo& CurrentProcessInfo ) override;

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy( Parameters& parameters ) final;

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2( Parameters& rValues ) final;

    /// DEPRECATED interface
    void CalculateMaterialResponse( const Vector& StrainVector,
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
     * to be called at the end of each step iteration
     * (e.g. from Element::FinalizeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeNonLinearIteration( const Properties& props,
                                     const GeometryType& geom, //this is just to give the array of nodes
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo ) override;

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param props the Properties instance of the current element
     * @param geom the geometry of the current element
     * @param ShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeSolutionStep( const Properties& props,
                               const GeometryType& geom, //this is just to give the array of nodes
                               const Vector& ShapeFunctionsValues,
                               const ProcessInfo& CurrentProcessInfo ) override;

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    int Check( const Properties& props,
               const GeometryType& geom,
               const ProcessInfo& CurrentProcessInfo ) override;

    /**
     * returns the size of the strain vector of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    SizeType GetStrainSize() override
    {
        return LinearElastic_Helper<TType>::StrainSize();
    }

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        if (TType == 0)
        {
            return "PlaneStrain";
        }
        else if (TType == 1)
        {
            return "PlaneStress";
        }
        else if (TType == 2)
        {
            return "Axisymmetric";
        }
        else if (TType == 3)
        {
            return "LinearElastic";
        }
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "E: " << mE << ", Nu: " << mNU << ", Prestress: " << mPrestress << ", PrestressFactor: " << mPrestressFactor
                 << ", CurrentStrain: " << mCurrentStrain << ", OldStrain: " << mOldStrain << ", stress_n: " << m_stress_n
                 << ", CurrentStress: " << mCurrentStress;
    }

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw );
        rSerializer.save( "Prestress", mPrestress );
        rSerializer.save( "PrestressFactor", mPrestressFactor );
        rSerializer.save( "CurrentStress", mCurrentStress );
        rSerializer.save( "OldStress", m_stress_n );
        rSerializer.save( "mE", mE );
        rSerializer.save( "mNU", mNU );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw );
        rSerializer.load( "Prestress", mPrestress );
        rSerializer.load( "PrestressFactor", mPrestressFactor );
        rSerializer.load( "CurrentStress", mCurrentStress );
        rSerializer.load( "OldStress", m_stress_n );
        rSerializer.load( "mE", mE );
        rSerializer.load( "mNU", mNU );
    }

    /**
     * Static Member Variables
     */

    /**
     * Calculates the stresses for given strain state
     * @param StrainVector the current vector of strains
     * @param rResult the stress vector corresponding to the given strains
     */
    void CalculateStress( const Vector& StrainVector, Matrix& AlgorithmicTangent, Vector& rResult );

    int mElemId, mGaussId;
    Vector mPrestress;
    double mPrestressFactor;
    Vector mCurrentStress;
    Matrix m_stress_n;
    Vector mCurrentStrain, mOldStrain;
    double mE, mNU;

    /**
     * Un accessible methods
     */
    /**
     * Assignment operator.
     */
    //LinearElastic& operator=(const IsotropicPlaneStressWrinklingNew& rOther);
    /**
     * Copy constructor.
     */
    //LinearElastic(const IsotropicPlaneStressWrinklingNew& rOther);
}; // Class LinearElastic

template<>
struct LinearElastic_Helper<0> // plane strain
{
    static constexpr unsigned int StrainSize()
    {
        return 3;
    }

    static void CalculateElasticMatrix( Matrix& C, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size, false);

        const double c1 = E * ( 1.00 - NU ) / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
        const double c2 = E * NU / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
        const double c3 = 0.5 * E / ( 1 + NU );

        C( 0, 0 ) = c1;
        C( 0, 1 ) = c2;
        C( 0, 2 ) = 0.0;
        C( 1, 0 ) = c2;
        C( 1, 1 ) = c1;
        C( 1, 2 ) = 0.0;
        C( 2, 0 ) = 0.0;
        C( 2, 1 ) = 0.0;
        C( 2, 2 ) = c3;
    }

    static void Apply( const Vector& StrainVector, Vector& StressVector, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(StressVector.size() != strain_size)
            StressVector.resize(strain_size, false);

        const double c1 = E * ( 1.00 - NU ) / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
        const double c2 = E * NU / (( 1.00 + NU ) * ( 1.00 - 2 * NU ) );
        const double c3 = 0.5 * E / ( 1 + NU );

        StressVector(0) = c1 * StrainVector(0) + c2 * StrainVector(1);
        StressVector(1) = c2 * StrainVector(0) + c1 * StrainVector(1);
        StressVector(2) = c3 * StrainVector(2);
    }

    /// REF: http://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_with_solutions.pdf
    static void CalculateInversedElasticMatrix( Matrix& C, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size, false);

        const double c1 = ( 1.00 - NU * NU ) / E;
        const double c2 = -NU * ( 1.00 + NU ) / E;
        const double c3 = 2 * ( 1 + NU ) / E;

        C( 0, 0 ) = c1;
        C( 0, 1 ) = c2;
        C( 0, 2 ) = 0.0;
        C( 1, 0 ) = c2;
        C( 1, 1 ) = c1;
        C( 1, 2 ) = 0.0;
        C( 2, 0 ) = 0.0;
        C( 2, 1 ) = 0.0;
        C( 2, 2 ) = c3;
    }

    static void ApplyInverse( const Vector& StressVector, Vector& StrainVector, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(StrainVector.size() != strain_size)
            StrainVector.resize(strain_size, false);

        const double c1 = ( 1.00 - NU * NU ) / E;
        const double c2 = -NU * ( 1.00 + NU ) / E;
        const double c3 = 2 * ( 1 + NU ) / E;

        StrainVector(0) = c1 * StressVector(0) + c2 * StressVector(1);
        StrainVector(1) = c2 * StressVector(0) + c1 * StressVector(1);
        StrainVector(2) = c3 * StressVector(2);
    }
};

template<>
struct LinearElastic_Helper<1> // plane stress
{
    static constexpr unsigned int StrainSize()
    {
        return 3;
    }

    static void CalculateElasticMatrix( Matrix& C, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size, false);

        const double c1 = E / (1.00 - NU*NU);
        const double c2 = c1 * NU;
        const double c3 = 0.5 * E / (1 + NU);

        C(0,0) = c1;
        C(0,1) = c2;
        C(0,2) = 0.0;
        C(1,0) = c2;
        C(1,1) = c1;
        C(1,2) = 0.0;
        C(2,0) = 0.0;
        C(2,1) = 0.0;
        C(2,2) = c3;
    }

    static void Apply( const Vector& StrainVector, Vector& StressVector, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(StressVector.size() != strain_size)
            StressVector.resize(strain_size, false);

        const double c1 = E / (1.00 - NU*NU);
        const double c2 = c1 * NU;
        const double c3 = 0.5 * E / (1 + NU);

        StressVector(0) = c1 * StrainVector(0) + c2 * StrainVector(1);
        StressVector(1) = c2 * StrainVector(0) + c1 * StrainVector(1);
        StressVector(2) = c3 * StrainVector(2);
    }

    /// REF: http://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_with_solutions.pdf
    static void CalculateInversedElasticMatrix( Matrix& C, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size, false);

        const double c1 = 1.00 / E;
        const double c2 = -NU * c1;
        const double c3 = 2 * ( 1 + NU ) / E;

        C( 0, 0 ) = c1;
        C( 0, 1 ) = c2;
        C( 0, 2 ) = 0.0;
        C( 1, 0 ) = c2;
        C( 1, 1 ) = c1;
        C( 1, 2 ) = 0.0;
        C( 2, 0 ) = 0.0;
        C( 2, 1 ) = 0.0;
        C( 2, 2 ) = c3;
    }

    static void ApplyInverse( const Vector& StressVector, Vector& StrainVector, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(StrainVector.size() != strain_size)
            StrainVector.resize(strain_size, false);

        const double c1 = 1.00 / E;
        const double c2 = -NU * c1;
        const double c3 = 2 * ( 1 + NU ) / E;

        StrainVector(0) = c1 * StressVector(0) + c2 * StressVector(1);
        StrainVector(1) = c2 * StressVector(0) + c1 * StressVector(1);
        StrainVector(2) = c3 * StressVector(2);
    }
};

template<>
struct LinearElastic_Helper<2> // axisymmetric
{
    static constexpr unsigned int StrainSize()
    {
        return 4;
    }

    static void CalculateElasticMatrix( Matrix& C, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size, false);

        const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
        const double c2 = c1 * ( 1 - NU );
        const double c3 = c1 * NU;
        const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

        C( 0, 0 ) = c2;
        C( 0, 1 ) = c3;
        C( 0, 2 ) = 0.0;
        C( 0, 3 ) = 0.0;
        C( 1, 0 ) = c3;
        C( 1, 1 ) = c2;
        C( 1, 2 ) = 0.0;
        C( 1, 3 ) = 0.0;
        C( 2, 0 ) = 0.0;
        C( 2, 1 ) = 0.0;
        C( 2, 2 ) = c4;
        C( 2, 3 ) = 0.0;
        C( 3, 0 ) = 0.0;
        C( 3, 1 ) = 0.0;
        C( 3, 2 ) = 0.0;
        C( 3, 3 ) = c2;
    }
};

template<>
struct LinearElastic_Helper<3> // 3D
{
    static constexpr unsigned int StrainSize()
    {
        return 6;
    }

    static void CalculateElasticMatrix( Matrix& C, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size, false);

        const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
        const double c2 = c1 * ( 1 - NU );
        const double c3 = c1 * NU;
        const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

        C( 0, 0 ) = c2;
        C( 0, 1 ) = c3;
        C( 0, 2 ) = c3;
        C( 0, 3 ) = 0.0;
        C( 0, 4 ) = 0.0;
        C( 0, 5 ) = 0.0;
        C( 1, 0 ) = c3;
        C( 1, 1 ) = c2;
        C( 1, 2 ) = c3;
        C( 1, 3 ) = 0.0;
        C( 1, 4 ) = 0.0;
        C( 1, 5 ) = 0.0;
        C( 2, 0 ) = c3;
        C( 2, 1 ) = c3;
        C( 2, 2 ) = c2;
        C( 2, 3 ) = 0.0;
        C( 2, 4 ) = 0.0;
        C( 2, 5 ) = 0.0;
        C( 3, 0 ) = 0.0;
        C( 3, 1 ) = 0.0;
        C( 3, 2 ) = 0.0;
        C( 3, 3 ) = c4;
        C( 3, 4 ) = 0.0;
        C( 3, 5 ) = 0.0;
        C( 4, 0 ) = 0.0;
        C( 4, 1 ) = 0.0;
        C( 4, 2 ) = 0.0;
        C( 4, 3 ) = 0.0;
        C( 4, 4 ) = c4;
        C( 4, 5 ) = 0.0;
        C( 5, 0 ) = 0.0;
        C( 5, 1 ) = 0.0;
        C( 5, 2 ) = 0.0;
        C( 5, 3 ) = 0.0;
        C( 5, 4 ) = 0.0;
        C( 5, 5 ) = c4;
    }

    static void Apply( const Vector& StrainVector, Vector& StressVector, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(StressVector.size() != strain_size)
            StressVector.resize(strain_size, false);

        const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
        const double c2 = c1 * ( 1 - NU );
        const double c3 = c1 * NU;
        const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

        StressVector(0) = c2 * StrainVector(0) + c3 * StrainVector(1) + c3 * StrainVector(2);
        StressVector(1) = c3 * StrainVector(0) + c2 * StrainVector(1) + c3 * StrainVector(2);
        StressVector(2) = c3 * StrainVector(0) + c3 * StrainVector(1) + c2 * StrainVector(2);
        StressVector(3) = c4 * StrainVector(3);
        StressVector(4) = c4 * StrainVector(4);
        StressVector(5) = c4 * StrainVector(5);
    }

    /// REF: https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Constitutive_relations
    static void CalculateInversedElasticMatrix( Matrix& C, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(C.size1() != strain_size || C.size2() != strain_size)
            C.resize(strain_size, strain_size, false);

        const double c1 = 1.0 / E;
        const double c2 = -NU * c1;
        const double c3 = 2 * (1 + NU) / E;

        C( 0, 0 ) = c1;
        C( 0, 1 ) = c2;
        C( 0, 2 ) = c2;
        C( 0, 3 ) = 0.0;
        C( 0, 4 ) = 0.0;
        C( 0, 5 ) = 0.0;
        C( 1, 0 ) = c2;
        C( 1, 1 ) = c1;
        C( 1, 2 ) = c2;
        C( 1, 3 ) = 0.0;
        C( 1, 4 ) = 0.0;
        C( 1, 5 ) = 0.0;
        C( 2, 0 ) = c2;
        C( 2, 1 ) = c2;
        C( 2, 2 ) = c1;
        C( 2, 3 ) = 0.0;
        C( 2, 4 ) = 0.0;
        C( 2, 5 ) = 0.0;
        C( 3, 0 ) = 0.0;
        C( 3, 1 ) = 0.0;
        C( 3, 2 ) = 0.0;
        C( 3, 3 ) = c3;
        C( 3, 4 ) = 0.0;
        C( 3, 5 ) = 0.0;
        C( 4, 0 ) = 0.0;
        C( 4, 1 ) = 0.0;
        C( 4, 2 ) = 0.0;
        C( 4, 3 ) = 0.0;
        C( 4, 4 ) = c3;
        C( 4, 5 ) = 0.0;
        C( 5, 0 ) = 0.0;
        C( 5, 1 ) = 0.0;
        C( 5, 2 ) = 0.0;
        C( 5, 3 ) = 0.0;
        C( 5, 4 ) = 0.0;
        C( 5, 5 ) = c3;
    }

    static void ApplyInverse( const Vector& StressVector, Vector& StrainVector, const double& E, const double& NU )
    {
        const unsigned int strain_size = StrainSize();
        if(StrainVector.size() != strain_size)
            StrainVector.resize(strain_size, false);

        const double c1 = 1.0 / E;
        const double c2 = -NU * c1;
        const double c3 = 2 * (1 + NU) / E;

        StrainVector(0) = c1 * StressVector(0) + c2 * StressVector(1) + c2 * StressVector(2);
        StrainVector(1) = c2 * StressVector(0) + c1 * StressVector(1) + c2 * StressVector(2);
        StrainVector(2) = c2 * StressVector(0) + c2 * StressVector(1) + c1 * StressVector(2);
        StrainVector(3) = c3 * StressVector(3);
        StrainVector(4) = c3 * StressVector(4);
        StrainVector(5) = c3 * StressVector(5);
    }
};

} // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_H_INCLUDED  defined

