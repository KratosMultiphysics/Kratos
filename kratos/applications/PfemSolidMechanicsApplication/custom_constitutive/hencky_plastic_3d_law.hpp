//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:             JMCarbonell $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_HENCKY_PLASTIC_3D_LAW_H_INCLUDED)
#define KRATOS_HENCKY_PLASTIC_3D_LAW_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"


namespace Kratos
{

/**
 * Defines a hencky-plastic isotropic constitutive law in 3D 
 * The functionality is limited to large and small displacements 
 */


class HenckyElasticPlastic3DLaw : public HyperElasticPlastic3DLaw
{
protected:

    struct MatrixSplit
    {
        Matrix  EigenValues;
        Matrix  EigenVectors;
    };


public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef FlowRule::Pointer                FlowRulePointer;
    typedef YieldCriterion::Pointer    YieldCriterionPointer;
    typedef HardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;

    /**
     * Counted pointer of HyperElasticPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(HenckyElasticPlastic3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HenckyElasticPlastic3DLaw();


    HenckyElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

    /**
     * Copy constructor.
     */
    HenckyElasticPlastic3DLaw (const HenckyElasticPlastic3DLaw& rOther)  ;


    /**
     * Assignment operator.
     */

    //HyperElasticPlastic3DLaw& operator=(const HyperElasticPlastic3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~HenckyElasticPlastic3DLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension()
    {
        return 3;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize()
    {
        return 6;
    };


/*  bool Has( const Variable<double>& rThisVariable );
    bool Has( const Variable<Vector>& rThisVariable );
    bool Has( const Variable<Matrix>& rThisVariable );

    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );


    void SetValue( const Variable<double>& rVariable,
                   const double& Value,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Vector>& rThisVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Matrix>& rThisVariable,
                   const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
*/    

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& rProps,
                             const GeometryType& rGeom,
                             const Vector& rShapeFunctionsValues );



    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues);



    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    //int Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);



    /**
     * Input and output
     */
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

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
  
//    Matrix mElasticLeftCauchyGreen;
    
//    FlowRulePointer mpFlowRule;

//    YieldCriterionPointer mpYieldCriterion;
	
//    HardeningLawPointer   mpHardeningLaw;
	
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


    /** First and secod term of the CONSISTENT ELASTOPLASTIC MATRIX FOR LARGE DEFORMATIONS 
        in a pullback fashion
    */
    
    virtual void SetConstitutiveMatrixToAppropiateDimension(Matrix& rConstitutiveMatrix);

    /** First and secod term of the CONSISTENT ELASTOPLASTIC MATRIX FOR LARGE DEFORMATIONS 
        in the actual configuration
    */

     Vector& GetStressVectorFromMatrix(const Matrix& rStressMatrix, 
				       Vector& rMainStress, 
				       const Matrix& rEigenVectors);

     virtual void CalculateHenckyMainStrain(const Matrix& rCauchyGreeMatrix, 
					    FlowRule::RadialReturnVariables& rReturnMappingVariables, 
					    Vector& rMainStrain);


    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
     //virtual bool CheckParameters(Parameters& rValues);


private:

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticPlastic3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticPlastic3DLaw )
    }




}; // Class HenckyElasticPlastic3DLaw

} //namespace Kratos

#endif  //KRATOS_HENCKY_PLASTIC_3D_LAW_H_INCLUDED

