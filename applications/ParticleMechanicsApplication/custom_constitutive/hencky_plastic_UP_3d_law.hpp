//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		BSD License 
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//

#if !defined (KRATOS_HENCKY_PLASTIC_UP_3D_LAW_H_INCLUDED)
#define KRATOS_HENCKY_PLASTIC_UP_3D_LAW_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/hencky_plastic_3d_law.hpp"
#include "includes/ublas_interface.h"

namespace Kratos
{

/**
 * Defines a hencky-plastic isotropic constitutive law in 3D 
 * The functionality is limited to large and small displacements 
 */


class HenckyElasticPlasticUP3DLaw : public HenckyElasticPlastic3DLaw
{
//protected:

    //struct MatrixSplit
    //{
        //Matrix  EigenValues;
        //Matrix  EigenVectors;
    //};


public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef MPMFlowRule::Pointer                MPMFlowRulePointer;
    typedef YieldCriterion::Pointer    YieldCriterionPointer;
    typedef HardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;

    /**
     * Counted pointer of HenckyElasticPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(HenckyElasticPlasticUP3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HenckyElasticPlasticUP3DLaw();


    HenckyElasticPlasticUP3DLaw(MPMFlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

    /**
     * Copy constructor.
     */
    HenckyElasticPlasticUP3DLaw (const HenckyElasticPlasticUP3DLaw& rOther)  ;


    /**
     * Assignment operator.
     */

    //HenckyElasticPlastic3DLaw& operator=(const HenckyElasticPlastic3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~HenckyElasticPlasticUP3DLaw();

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

	virtual void GetLawFeatures(Features& rFeatures);
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
    //void InitializeMaterial( const Properties& rProps,
                             //const GeometryType& rGeom,
                             //const Vector& rShapeFunctionsValues );



    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    //void CalculateMaterialResponseKirchhoff (Parameters & rValues);



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
    
    //virtual Matrix SetConstitutiveMatrixToAppropiateDimension(Matrix& rConstitutiveMatrix, const Matrix& rElastoPlasticTangentMatrix);
    
    //virtual Vector SetStressMatrixToAppropiateVectorDimension(Vector& rStressVector, const Matrix& rStressMatrix );
    
    virtual void CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables& rElasticVariables);
    
    //virtual void CalculateElastoPlasticTangentMatrix( const FlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen,const double& rAlpha, Matrix& rElastoPlasticMatrix, const MaterialResponseVariables& rElasticVariables);

	
	
    void GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables);


    void CalculateElastoPlasticTangentMatrix( const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen,const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables);

	 /**
     * Calculates the GreenLagrange strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    virtual void CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
            Vector& rStrainVector );


    /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    virtual void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
                                         Vector& rStrainVector );
                                         
    virtual void CalculatePrincipalStressTrial(const MaterialResponseVariables & rElasticVariables,Parameters & rValues, 
						const MPMFlowRule::RadialReturnVariables& rReturnMappingVariables, 
					    Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix);
	//double& TensorComponent(double & rCabcd,
						 //const Matrix& rMA, const Matrix& rMB, 
						 //const unsigned int& a, const unsigned int& b,
						 //const unsigned int& c, const unsigned int& d);
	
	//virtual void MyTensorProduct(const Matrix& rMA, const Matrix& rMB, 
						      //Matrix& rEigenbasesProductMatrix);
						      
    //double& TensorComponent2(double & rCabcd,
						 //const Matrix& rMA, const Matrix& rMB, 
						 //const unsigned int& a, const unsigned int& b,
						 //const unsigned int& c, const unsigned int& d);
	
	//virtual void MyTensorProduct2(const Matrix& rMA, const Matrix& rMB, 
						      //Matrix& rEigenbasesProductMatrix);
						     
	//double& TensorComponent3(double & rCabcd,
						 //const Matrix& rMA, 
						 //const unsigned int& a, const unsigned int& b,
						 //const unsigned int& c, const unsigned int& d);
	
	//virtual void MyTensorProduct3(const Matrix& rMA, 
						      //Matrix& rEigenbasesProductMatrix);
						      
	//double& TensorComponent4(double & rCabcd,
						 //const Matrix& rMA, 
						 //const unsigned int& a, const unsigned int& b,
						 //const unsigned int& c, const unsigned int& d);
	
	//virtual void MyTensorProduct4(const Matrix& rMA, 
						      //Matrix& rEigenbasesProductMatrix);
						      
						      
	//virtual Matrix CalculateEigenbases(const FlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rEigenbasesMatrix);
	
	//virtual void CalculateElastoPlasticTangentMatrix( const FlowRule::RadialReturnVariables & rReturnMappingVariables, 
													  //const Matrix& rTrialElasticLeftCauchyGreen, const Matrix& rStressMatrix,
													  //Matrix& rElastoPlasticTangentMatrix, Matrix& rConsistentMatrix );

    /** First and secod term of the CONSISTENT ELASTOPLASTIC MATRIX FOR LARGE DEFORMATIONS 
        in the actual configuration
    */

     //Vector& GetStressVectorFromMatrix(const Matrix& rStressMatrix, 
				       //Vector& rMainStress, 
				       //const Matrix& rEigenVectors);

     //virtual void CalculateHenckyMainStrain(const Matrix& rCauchyGreeMatrix, 
					    //FlowRule::RadialReturnVariables& rReturnMappingVariables, 
					    //Vector& rMainStrain);


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

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HenckyElasticPlastic3DLaw )
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HenckyElasticPlastic3DLaw )
    }




}; // Class HenckyElasticPlasticUP3DLaw

} //namespace Kratos

#endif  //KRATOS_HENCKY_PLASTIC_3D_LAW_H_INCLUDED

