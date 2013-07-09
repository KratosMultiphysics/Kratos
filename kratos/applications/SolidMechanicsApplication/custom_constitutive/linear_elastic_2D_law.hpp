//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_LINEAR_ELASTIC_2D_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_2D_LAW_H_INCLUDED

// System includes 

// External includes 

// Project includes
#include "includes/constitutive_law.h"


namespace Kratos
{
  /**
   * Defines a hyper elastic isotropic constitutive law in 2D
   * This material law is defined by the parameters:
   * 1) E  (Young's modulus) 
   * 2) NU (Poisson ratio)
   * As there are no further parameters the functionality is limited 
   * to large displacements elasticity.
   */

  class LinearElastic2DLaw : public ConstitutiveLaw
  {
  public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of LinearElastic2DLaw
     */
    
    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic2DLaw);
    
    /**
     * Life Cycle 
     */

    /**
     * Default constructor.
     */
    LinearElastic2DLaw();
			
    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;
    
    /**
     * Copy constructor.
     */
    LinearElastic2DLaw (const LinearElastic2DLaw& rOther);
   

    /**
     * Assignment operator.
     */

    //LinearElastic2DLaw& operator=(const LinearElastic2DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~LinearElastic2DLaw();
			
    /**
     * Operators 
     */
    
    /**
     * Operations needed by the base class:
     */

    SizeType WorkingSpaceDimension() { return 2; };
    SizeType GetStrainSize()         { return 3; };

    bool Has( const Variable<double>& rThisVariable );
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
    /**
     * Material parameters are inizialized
     */ 
    void InitializeMaterial( const Properties& props,
			     const GeometryType& geom,
			     const Vector& ShapeFunctionsValues );

    /**
     * As this constitutive law describes only linear elastic material properties
     * this function is rather useless and in fact does nothing
     */ 		
    void InitializeSolutionStep( const Properties& props,
				 const GeometryType& geom, //this is just to give the array of nodes
				 const Vector& ShapeFunctionsValues ,
				 const ProcessInfo& CurrentProcessInfo);
			
    void FinalizeSolutionStep( const Properties& props,
			       const GeometryType& geom, //this is just to give the array of nodes
			       const Vector& ShapeFunctionsValues ,
			       const ProcessInfo& CurrentProcessInfo);
    

    /**
     * Calculates the cauchy stresses. For a given deformation and stress state
     * the cauchy stress vector is calculated
     * @param Cauchy_StressVector the result vector of cauchy stresses (will be overwritten)
     * @param F the current deformation gradient
     * @param PK2_StressVector the current second Piola-Kirchhoff-Stress vector
     * @param GreenLagrangeStrainVector the current Green-Lagrangian strains
     */
    void CalculateCauchyStresses( Vector& rCauchy_StressVector,
				  const Matrix& rF,
				  const Vector& rPK2_StressVector,
				  const Vector& rGreenLagrangeStrainVector);

					           
    /**
     * Computes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues 
     * @see   Parameters
     */          
    void CalculateMaterialResponsePK1 (Parameters & rValues);

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues 
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters & rValues);

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues 
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues);


    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues 
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters & rValues);
    
    
   /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues 
     * @see   Parameters
     */
    void FinalizeMaterialResponsePK1 (Parameters & rValues);
 
   /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues 
     * @see   Parameters
     */
    void FinalizeMaterialResponsePK2 (Parameters & rValues);

   /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues 
     * @see   Parameters
     */
    void FinalizeMaterialResponseKirchhoff (Parameters & rValues);

   /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues 
     * @see   Parameters
     */
    void FinalizeMaterialResponseCauchy (Parameters & rValues);

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    int Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

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
    /**
     * there are no protected class members
     */
		
  private:

    ///@}
    ///@Member Variables
    ///@{
			
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

    virtual void load(Serializer& rSerializer)
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }


    /**
     * Calculates the stresses for given strain state
     * @param rStrainVector
     * @param rConstitutiveMatrix
     * @param rStressVector the stress vector corresponding to the deformation
     */
    void CalculateStress( const Vector &rStrainVector,
			  const Matrix &rConstitutiveMatrix,
			  Vector& rStressVector);

	       
    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */


    void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix, 
				       const double &rYoungModulus, 
				       const double &rPoissonCoefficient );

   /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized 
     * @param Parameters
     * @return
     */
    bool CheckParameters(Parameters& rValues);
   
    /**
     * Unaccessible methods 
     */


  }; // Class LinearElastic2DLaw 
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_2D_LAW_H_INCLUDED  defined 
