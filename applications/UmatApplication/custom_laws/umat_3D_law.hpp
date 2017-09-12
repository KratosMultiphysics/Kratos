//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:      JMCarbonell  $
//   Last modified by:    $Co-Author:   LlMonforte   $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_UMAT_3D_LAW_H_INCLUDED )
#define  KRATOS_UMAT_3D_LAW_H_INCLUDED

// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

  class KRATOS_API(UMAT_APPLICATION) Umat3DLaw : public ConstitutiveLaw
  {
  protected:

    using VoigtIndexType = const unsigned int(*)[2];

  public:

    ///@name Type Definitions
    ///@{
    typedef ProcessInfo                                           ProcessInfoType;
    typedef ConstitutiveLaw                                              BaseType;

    typedef array_1d<double, 81>                               MaterialTensorType;
    typedef array_1d<double, 3 >                                   PlaneArrayType;
    typedef array_1d<double, 6 >                                   SpaceArrayType;

    /// Pointer definition of Umat3DLaw
    KRATOS_CLASS_POINTER_DEFINITION( Umat3DLaw );
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Umat3DLaw();

    /// Copy constructor.
    Umat3DLaw(const Umat3DLaw& rOther);

    /// Clone.
    ConstitutiveLaw::Pointer Clone() const override;

    /// Assignment operator.
    Umat3DLaw& operator=(const Umat3DLaw& rOther);

    /// Destructor.
    virtual ~Umat3DLaw();

    ///@}
    ///@name Operators
    ///@{

    
    ///@}
    ///@name Operations
    ///@{

    /**
     * Material parameters are inizialized
     */
    virtual void InitializeMaterial ( const Properties& props,
				      const GeometryType& geom,
				      const Vector& ShapeFunctionsValues ) override;

    /**
     * Step Initialize
     */
    virtual void InitializeSolutionStep ( const Properties& props,
					  const GeometryType& geom,
					  const Vector& ShapeFunctionsValues ,
					  const ProcessInfo& CurrentProcessInfo ) override;

    /**
     * Step Finalize
     */
    virtual void FinalizeSolutionStep ( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues ,
					const ProcessInfo& CurrentProcessInfo ) override;


    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseCauchy (Parameters & rValues) override;

    
    /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues
     * @see   Parameters
     */
    virtual void FinalizeMaterialResponseCauchy(Parameters & rValues) override;


    /**
     * This function is designed to be called once to check compatibility with element and the constitutive law
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    
    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
      return 3;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
      return 6;
    };

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Umat3DLaw";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Umat3DLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "Umat3DLaw Data";
    }


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

    //member variables for umat
    double* STRESS;
    double* STATEV;

    double* STRESS_FINALIZED;
    double* STATEV_FINALIZED;

    double* STRAN_FINALIZED; 

    int* NSTATV;

    double* PROPS;
    int* NPROPS;

    Matrix mPreviousDeformationGradient;

    //static member variables for umat
    //double* STRAN;
    //double* DSTRAN;


    //variable for wrapper
    int* MaterialNumber;

    
    ///@}
    ///@name Protected Operators
    ///@{
    
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Get voigt index tensor:
     */
    virtual VoigtIndexType GetVoigtIndexTensor()
    {
      return this->msIndexVoigt3D6C;
    }

    /**
     * Load previous information:
     */
    void LoadPreviousInformation();
    
    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
    virtual bool CheckParameters(Parameters& rValues);


    ///@}
    
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
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save ( Serializer& rSerializer ) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, ConstitutiveLaw )
    }

    virtual void load ( Serializer& rSerializer ) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, ConstitutiveLaw )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
  }; // Class Umat3DLaw
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_UMAT_3D_LAW_H_INCLUDED  defined 

