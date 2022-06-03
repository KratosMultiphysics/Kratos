// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_SMALL_STRAIN_UMAT_3D_LAW_H_INCLUDED )
#define  KRATOS_SMALL_STRAIN_UMAT_3D_LAW_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include "includes/define.h"

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/constitutive_law.h"

// Application includes
#include "geo_mechanics_application_variables.h"


namespace Kratos
{
/*
   - structure of UMATs:
   - Function to get stress, stiffness matrix
     void umat(double* STRESS,  double* STATEV,  double** DDSDDE,  double* SSE,     double* SPD,     double* SCD,
               double* RPL,     double* DDSDDT,  double* DRPLDE,   double* DRPLDT,  double* STRAN,   double* DSTRAN,
               double* TIME,    double* DTIME,   double* TEMP,     double* DTEMP,   double* PREDEF,  double* DPRED,
               char* MATERL,    int* NDI,        int* NSHR,        int* NTENS,      int* NSTATV,     double* PROPS,
               int* NPROPS,     double* COORDS,  double** DROT,    double* PNEWDT,  double* CELENT,  double** DFGRD0,
               double** DFGRD1, int* NOEL,       int* NPT,         double* KSLAY,   double* KSPT,    int* KSTEP,
               int* KINC);

*/

typedef void(*pF_UMATMod) (double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
                          double* rpl, double* ddsddt, double* drplde, double* drpldt, double* stran, double* dstran,
                          double* time, double* dtime, double* temp, double* dtemp, double* predef, double* dpred,
                          char* materl, int* ndi, int* nshr, int* ntens, int* nstatv, const double* props, int* nprops,
                          double* coords, double** drot, double* pnewdt, double* celent, double** dfgrd0,
                          double** dfgrd1, int* noel, int* npt, double* kslay, double* kspt, int* kstep,
                          int* kinc);


   ///@addtogroup ConstitutiveModelsApplication
   ///@{

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

   class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUMAT3DLaw: public ConstitutiveLaw
   {
   public:
      // The process info type definition
      typedef ProcessInfo            ProcessInfoType;

      // The base class ConstitutiveLaw type definition
      typedef ConstitutiveLaw        BaseType;

      /// The size type definition
      typedef std::size_t             SizeType;

      /// Static definition of the dimension
      static constexpr SizeType Dimension = N_DIM_3D;

      /// Static definition of the VoigtSize
      static constexpr SizeType VoigtSize = VOIGT_SIZE_3D;

      /// Pointer definition of SmallStrainUMAT3DLaw
      KRATOS_CLASS_POINTER_DEFINITION( SmallStrainUMAT3DLaw );


      //@}
      //@name Life Cycle
      //@{

      //----------------------------------------------------------------------------------------
      /**
       * @brief Default constructor.
       */
      SmallStrainUMAT3DLaw();

      /**
       * @brief Clone method
       */
      ConstitutiveLaw::Pointer Clone() const override;

      /**
       * Copy constructor.
       */
      SmallStrainUMAT3DLaw(SmallStrainUMAT3DLaw const& rOther);

      /**
       * @brief Destructor.
       */
      virtual ~SmallStrainUMAT3DLaw();


      // Assignment operator:
      SmallStrainUMAT3DLaw& operator=(SmallStrainUMAT3DLaw const& rOther);

      //----------------------------------------------------------------------------------------
      /**
       * @brief This function is designed to be called once to check compatibility with element
       * @param rFeatures The Features of the law
       */
      void GetLawFeatures(Features& rFeatures) override;

      /**
       * @brief Dimension of the law:
       */
      virtual SizeType WorkingSpaceDimension() override
      {
         return Dimension;
      }

      /**
       * @brief Voigt tensor size:
       */
      virtual SizeType GetStrainSize() override
      {
         return VoigtSize;
      }

      /**
       * @brief Returns the expected strain measure of this constitutive law (by default Green-Lagrange)
       * @return the expected strain measure
       */
      virtual StrainMeasure GetStrainMeasure() override
      {
         return StrainMeasure_Infinitesimal;
      }

      /**
       * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
       * @return the expected stress measure
       */
      virtual StressMeasure GetStressMeasure() override
      {
         return StressMeasure_Cauchy;
      }

      /**
       * @brief Computes the material response:
       * @details PK1 stresses and algorithmic ConstitutiveMatrix
       * @param rValues The internal values of the law
       * @see   Parameters
       */
      void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters & rValues) override;

      /**
       * @brief Computes the material response:
       * @details PK2 stresses and algorithmic ConstitutiveMatrix
       * @param rValues The internal values of the law
       * @see   Parameters
       */
      void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters & rValues) override;

      /**
       * @brief Computes the material response:
       * @details Kirchhoff stresses and algorithmic ConstitutiveMatrix
       * @param rValues The internal values of the law
       * @see   Parameters
       */
      void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters & rValues) override;

      /**
       * @brief Computes the material response:
       * @details Cauchy stresses and algorithmic ConstitutiveMatrix
       * @param rValues The internal values of the law
       * @see   Parameters
       */
      void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters & rValues) override;

      /**
       * @brief Updates the material response:
       * @details Cauchy stresses and Internal Variables
       * @param rValues The internal values of the law
       * @see   Parameters
       */
      void FinalizeMaterialResponseCauchy   (ConstitutiveLaw::Parameters & rValues) override;
      void FinalizeMaterialResponsePK1      (ConstitutiveLaw::Parameters & rValues) override;
      void FinalizeMaterialResponsePK2      (ConstitutiveLaw::Parameters & rValues) override;
      void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters & rValues) override;

      /**
       * @brief It calculates the value of a specified variable (double case)
       * @param rParameterValues the needed parameters for the CL calculation
       * @param rThisVariable the variable to be returned
       * @param rValue a reference to the returned value
       * @return rValue output: the value of the specified variable
       */
      virtual double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                     const Variable<double>& rThisVariable,
                                     double& rValue) override;

      /**
       * @brief It calculates the value of a specified variable (Vector case)
       * @param rParameterValues the needed parameters for the CL calculation
       * @param rThisVariable the variable to be returned
       * @param rValue a reference to the returned value
       * @return rValue output: the value of the specified variable
       */
      virtual Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                     const Variable<Vector>& rThisVariable,
                                     Vector& rValue) override;

      /**
       * @brief It calculates the value of a specified variable (Matrix case)
       * @param rParameterValues the needed parameters for the CL calculation
       * @param rThisVariable the variable to be returned
       * @param rValue a reference to the returned value
       * @return rValue output: the value of the specified variable
       */
      virtual Matrix& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                     const Variable<Matrix>& rThisVariable,
                                     Matrix& rValue) override;


      // @brief This function provides the place to perform checks on the completeness of the input.
      // @details It is designed to be called only once (or anyway, not often) typically at the beginning 
      //          of the calculations, so to verify that nothing is missing from the input or that 
      //          no common error is found.
      int Check(const Properties& rMaterialProperties,
                const GeometryType& rElementGeometry,
                const ProcessInfo& rCurrentProcessInfo) override;

      /**
       * This is to be called at the very beginning of the calculation
       * (e.g. from InitializeElement) in order to initialize all relevant
       * attributes of the constitutive law
       * @param rMaterialProperties the Properties instance of the current element
       * @param rElementGeometry the geometry of the current element
       * @param rShapeFunctionsValues the shape functions values in the current integration point
       */
      void InitializeMaterial(const Properties& rMaterialProperties,
                              const GeometryType& rElementGeometry,
                              const Vector& rShapeFunctionsValues) override;


      /**
       * Initialize the material response in terms of Cauchy stresses
       * @see Parameters
       */
      void InitializeMaterialResponseCauchy   (ConstitutiveLaw::Parameters& rValues) override;
      void InitializeMaterialResponsePK1      (ConstitutiveLaw::Parameters& rValues) override;
      void InitializeMaterialResponsePK2      (ConstitutiveLaw::Parameters& rValues) override;
      void InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

      /**
       * @brief It calculates the strain vector
       * @param rValues The internal values of the law
       * @param rStrainVector The strain vector in Voigt notation
       */
      virtual void CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector);

    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
     void ResetMaterial(const Properties& rMaterialProperties,
                        const GeometryType& rElementGeometry,
                        const Vector& rShapeFunctionsValues) override;

   //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    double& GetValue( const Variable<double> &rThisVariable, double &rValue ) override;
    int& GetValue( const Variable<int> &rThisVariable, int &rValue ) override;
    Vector& GetValue( const Variable<Vector> &rThisVariable, Vector &rValue ) override;

    void SetValue( const Variable<double>& rVariable,
                   const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValue( const Variable<Vector>& rVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;

   //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

      ///@}
      ///@name Inquiry
      ///@{


      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const override
      {
         std::stringstream buffer;
         buffer << "SmallStrainUMAT3DLaw";
         return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override
      {
         rOStream << "SmallStrainUMAT3DLaw";
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override
      {
         rOStream << "SmallStrainUMAT3DLaw Data";
      }


      ///@}
      ///@name Friends
      ///@{


      ///@}

   protected:
      ///@name Protected member Variables
      ///@{
      array_1d<double, VOIGT_SIZE_3D> mStressVector;
      array_1d<double, VOIGT_SIZE_3D> mStressVectorFinalized;

      array_1d<double, VOIGT_SIZE_3D> mDeltaStrainVector;
      array_1d<double, VOIGT_SIZE_3D> mStrainVectorFinalized;

      double mMatrixD[VOIGT_SIZE_3D][VOIGT_SIZE_3D];

      ///@}
      ///@name Protected Operators
      ///@{


      ///@}
      ///@name Protected Operations
      ///@{

      ///@}
      ///@name Protected  Access
      ///@{
      virtual void UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters &rValues);
      virtual void UpdateInternalStrainVectorFinalized(ConstitutiveLaw::Parameters &rValues);
      virtual void SetExternalStressVector(Vector& rStressVector);
      virtual void SetInternalStressVector(const Vector& rStressVector);
      virtual void SetInternalStrainVector(const Vector& rStrainVector);
      virtual void CopyConstitutiveMatrix(ConstitutiveLaw::Parameters &rValues, Matrix& rConstitutiveMatrix);


      void CalculateConstitutiveMatrix(ConstitutiveLaw::Parameters &rValues, Matrix& rConstitutiveMatrix);
      void CalculateStress(ConstitutiveLaw::Parameters &rValues, Vector& rStressVector);

      int GetStateVariableIndex(const Variable<double>& rThisVariable);

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
      pF_UMATMod          pUserMod;

      bool mIsModelInitialized;
      bool mIsUMATLoaded;

      std::vector<int> mProjectDirectory;

      Vector mStateVariables;
      Vector mStateVariablesFinalized;


      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{


      ///@}
      ///@name Private  Access
      ///@{

      // to load UMAT and functions
      bool loadUMAT(const Properties& rMaterialProperties);
      bool loadUMATWindows(const Properties& rMaterialProperties);
      bool loadUMATLinux(const Properties& rMaterialProperties);

      // Set number of MaterialParameters
      void CallUMAT( ConstitutiveLaw::Parameters &rValues);

      // Set state variables to the initial values
      void ResetStateVariables(const Properties& rMaterialProperties);


      ///@}
      ///@name Serialization
      ///@{
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const override
      {
         KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
         rSerializer.save("InitializedModel",           mIsModelInitialized);
         rSerializer.save("StressVectorFinalized",      mStressVectorFinalized);
         rSerializer.save("StrainVectorFinalized",      mStrainVectorFinalized);
         rSerializer.save("StateVariablesFinalized",    mStateVariablesFinalized);
      }

      virtual void load(Serializer& rSerializer) override
      {
         KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
         rSerializer.load("InitializedModel",           mIsModelInitialized);
         rSerializer.load("StressVectorFinalized",      mStressVectorFinalized);
         rSerializer.load("StrainVectorFinalized",      mStrainVectorFinalized);
         rSerializer.load("StateVariablesFinalized",    mStateVariablesFinalized);
      }

      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{

      ///@}

   }; // Class SmallStrainUMAT3DLaw

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{

   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SMALL_STRAIN_UMAT_3D_LAW_H_INCLUDED  defined 


