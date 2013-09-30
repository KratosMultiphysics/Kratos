/* *********************************************************
 *
 *   Last Modified by:    $Author:   JMCarbonell$
 *   Date:                $Date:     12-06-2012$
 *   Revision:            $Revision: 1.0$
 *
 * ***********************************************************/

#include "includes/constitutive_law.h"


namespace Kratos
{
    const unsigned int ConstitutiveLaw::msIndexVoigt3D6C [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };
    const unsigned int ConstitutiveLaw::msIndexVoigt2D4C [4][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1} };
    const unsigned int ConstitutiveLaw::msIndexVoigt2D3C [3][2] = { {0, 0}, {1, 1}, {0, 1} };

    /**
     * Flags related to the Parameters of the Contitutive Law
     */
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, COMPUTE_STRAIN,              0 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, COMPUTE_STRESS,              1 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, COMPUTE_CONSTITUTIVE_TENSOR, 2 );
      
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, ISOCHORIC_TENSOR_ONLY,       3 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, VOLUMETRIC_TENSOR_ONLY,      4 );
      
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, TOTAL_TENSOR,                5 );
      
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, INITIAL_CONFIGURATION,       6 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, LAST_KNOWN_CONFIGURATION,    7 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, FINAL_CONFIGURATION,         8 );

    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, FINALIZE_MATERIAL_RESPONSE,  9 );
  

    /**
     * Flags related to the Features of the Contitutive Law
     */
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, FINITE_STRAINS,             10 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, INFINITESIMAL_STRAINS,     11 );

    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, THREE_DIMENSIONAL_LAW,      12 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, PLANE_STRAIN_LAW,           13 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, PLANE_STRESS_LAW,           14 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, AXISYMMETRIC_LAW,           15 );

    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, U_P_LAW,                    16 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, ISOTROPIC,                  17 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, ANISOTROPIC,                18 );


    /**
     * Constructor.
     */
    ConstitutiveLaw::ConstitutiveLaw() : Flags()
    {
    }

 
    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      ConstitutiveLaw::Pointer p_clone(new ConstitutiveLaw());
     *      return p_clone;
     */
    ConstitutiveLaw::Pointer ConstitutiveLaw::Clone() const
    {
        KRATOS_ERROR(std::logic_error, "Called the virtual function for Clone", "");
    }

    /**
     * @return the working space dimension of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    ConstitutiveLaw::SizeType ConstitutiveLaw::WorkingSpaceDimension()
    {
        KRATOS_ERROR(std::logic_error, "Called the virtual function for WorkingSpaceDimension", "");
    }

    /**
     * returns the size of the strain vector of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    ConstitutiveLaw::SizeType ConstitutiveLaw::GetStrainSize()
    {
        KRATOS_ERROR(std::logic_error, "Called the virtual function for GetStrainSize", "");
    }

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
     bool ConstitutiveLaw::Has(const Variable<double>& rThisVariable)
    {
        return false;
    }

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
     bool ConstitutiveLaw::Has(const Variable<Vector>& rThisVariable)
    {
        return false;
    }

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
     bool ConstitutiveLaw::Has(const Variable<Matrix>& rThisVariable)
    {
        return false;
    }

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * NOTE: fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
     */
     bool ConstitutiveLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
        return false;
    }

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * NOTE: fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
     */
     bool ConstitutiveLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
    {
        return false;
    }

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
     double& ConstitutiveLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
        return rValue;
    }

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
     Vector& ConstitutiveLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
        return rValue;
    }

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @return the value of the specified variable
     */
     Matrix& ConstitutiveLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
        return rValue;
    }

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
     array_1d<double, 3 > & ConstitutiveLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                            array_1d<double, 3 > & rValue)
    {
        return rValue;
    }

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
     array_1d<double, 6 > & ConstitutiveLaw::GetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                            array_1d<double, 6 > & rValue)
    {
        return rValue;
    }

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
     void ConstitutiveLaw::SetValue(const Variable<double>& rVariable,
                          const double& rValue,
                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "Called the virtual function for SetValue", "");
    }

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
     void ConstitutiveLaw::SetValue(const Variable<Vector >& rVariable,
                          const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "Called the virtual function for SetValue", "");
    }

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
     void ConstitutiveLaw::SetValue(const Variable<Matrix >& rVariable,
                          const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "Called the virtual function for SetValue", "");
    }

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
     void ConstitutiveLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                          const array_1d<double, 3 > & rValue,
                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "Called the virtual function for SetValue", "");
    }

    /**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
     void ConstitutiveLaw::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                          const array_1d<double, 6 > & rValue,
                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "Called the virtual function for SetValue", "");
    }

    /**
     * Is called to check whether the provided material parameters in the Properties
     * match the requirements of current constitutive model.
     * @param rMaterialProperties the current Properties to be validated against.
     * @return true, if parameters are correct; false, if parameters are insufficient / faulty
     * NOTE: this has to implemented by each constitutive model. Returns false in base class since
     * no valid implementation is contained here.
     */
     bool ConstitutiveLaw::ValidateInput(const Properties& rMaterialProperties)
    {
        return false;
    }

    /**
     * returns the expected strain measure of this constitutive law (by default linear strains)
     * @return the expected strain measure
     */
    ConstitutiveLaw::StrainMeasure ConstitutiveLaw::GetStrainMeasure()
    {
        return StrainMeasure_Infinitesimal;
    }

    /**
     * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
     ConstitutiveLaw::StressMeasure ConstitutiveLaw::GetStressMeasure()
    {
        return StressMeasure_PK1;
    }

    /**
     * returns whether this constitutive model is formulated in incremental strains/stresses
     * NOTE: by default, all constitutive models should be formulated in total strains
     * @return true, if formulated in incremental strains/stresses, false otherwise
     */
     bool ConstitutiveLaw::IsIncremental()
    {
        return false;
    }

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
     void ConstitutiveLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues)
    {
    }

    /**
     * to be called at the beginning of each solution step
     * (e.g. from Element::InitializeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
     void ConstitutiveLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
                                        const GeometryType& rElementGeometry, //this is just to give the array of nodes
                                        const Vector& rShapeFunctionsValues,
                                        const ProcessInfo& rCurrentProcessInfo)
    {
    }

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
     void ConstitutiveLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo)
    {
    }

    /**
     * to be called at the beginning of each step iteration
     * (e.g. from Element::InitializeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void ConstitutiveLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
            const GeometryType& rElementGeometry,
            const Vector& rShapeFunctionsValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "Calling virtual function for InitializeNonLinearIteration", "");
    }



    /**
     * to be called at the end of each step iteration
     * (e.g. from Element::FinalizeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void ConstitutiveLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
					    const GeometryType& rElementGeometry,
					    const Vector& rShapeFunctionsValues,
					    const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "Calling virtual function for FinalizeNonLinearIteration", "");
    }

    /**
     * Computes the material response in terms of stresses and constitutive tensor
     * @see Parameters
     * @see StressMeasures
     */

    void ConstitutiveLaw::CalculateMaterialResponse(Parameters& rValues,const StressMeasure& rStressMeasure)
    {
      switch(rStressMeasure)
	{
	case StressMeasure_PK1:         CalculateMaterialResponsePK1(rValues);
	  break;
      
	case StressMeasure_PK2:         CalculateMaterialResponsePK2(rValues);
	  break;
	  
	case StressMeasure_Kirchhoff: 	CalculateMaterialResponseKirchhoff(rValues);
	  break;

	case StressMeasure_Cauchy:	CalculateMaterialResponseCauchy(rValues);
	  break;
	  
	default:
	  KRATOS_ERROR(std::logic_error, " Stress Measure not Defined ", "");
	  break;

	}
    }


    /**
     * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */

     void ConstitutiveLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
    {
      KRATOS_ERROR(std::logic_error, "Calling virtual function for CalculateMaterialResponsePK1", "");
    }

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */

     void ConstitutiveLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
    {
      KRATOS_ERROR(std::logic_error, "Calling virtual function for CalculateMaterialResponsePK2", "");
    }

    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */

     void ConstitutiveLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
    {
      KRATOS_ERROR(std::logic_error, "Calling virtual function for CalculateMaterialResponseKirchhoff", "");
    }

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */

     void ConstitutiveLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
    {
      KRATOS_ERROR(std::logic_error, "Calling virtual function for CalculateMaterialResponseCauchy", "");
    }


    
    /**
     * Updates the material response,  called by the element in FinalizeSolutionStep.
     * @see Parameters
     * @see StressMeasures
     */

    void ConstitutiveLaw::FinalizeMaterialResponse(Parameters& rValues,const StressMeasure& rStressMeasure)
    {
      switch(rStressMeasure)
	{
	case StressMeasure_PK1:         FinalizeMaterialResponsePK1(rValues);
	  break;
      
	case StressMeasure_PK2:         FinalizeMaterialResponsePK2(rValues);
	  break;
	  
	case StressMeasure_Kirchhoff: 	FinalizeMaterialResponseKirchhoff(rValues);
	  break;

	case StressMeasure_Cauchy:	FinalizeMaterialResponseCauchy(rValues);
	  break;
	  
	default:
	  KRATOS_ERROR(std::logic_error, " Stress Measure not Defined ", "");
	  break;

	}
    }


    /**
     * Updates the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */

     void ConstitutiveLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
    {
      KRATOS_ERROR(std::logic_error, "Calling virtual function for FinalizeMaterialResponsePK1", "");
    }

    /**
     * Updates the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */

     void ConstitutiveLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
    {
      KRATOS_ERROR(std::logic_error, "Calling virtual function for FinalizeMaterialResponsePK2", "");
    }

    /**
     * Updates the material response in terms of Kirchhoff stresses
     * @see Parameters
     */

     void ConstitutiveLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
    {
      KRATOS_ERROR(std::logic_error, "Calling virtual function for FinalizeMaterialResponseKirchhoff", "");
    }

    /**
     * Updates the material response in terms of Cauchy stresses
     * @see Parameters
     */

     void ConstitutiveLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
    {
      KRATOS_ERROR(std::logic_error, "Calling virtual function for FinalizeMaterialResponseCauchy", "");
    }



    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
     void ConstitutiveLaw::ResetMaterial(const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry,
                               const Vector& rShapeFunctionsValues)
    {
        KRATOS_ERROR(std::logic_error, "Calling virtual function for ResetMaterial", "");
    }




    /**
     * Methods to transform strain Vectors:
     * @param rStrainVector the strain tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStrainInitial the measure of stress of the given  rStrainVector
     * @param rStrainFinal the measure of stress of the returned rStrainVector
     */
     Vector& ConstitutiveLaw::TransformStrains (Vector& rStrainVector,
				       const Matrix &rF,
				       StrainMeasure rStrainInitial,
				       StrainMeasure rStrainFinal)
    {

      switch(rStrainInitial)
	{
        case StrainMeasure_GreenLagrange:         

	  switch(rStrainFinal)
	    {
	    case StrainMeasure_GreenLagrange:         
	      break;
	      
	    case StrainMeasure_Almansi:
	      {
	      Matrix StrainMatrix = MathUtils<double>::StrainVectorToTensor( rStrainVector ); //is it defined ¿?
	      
	      CoVariantPushForward (StrainMatrix,rF);  //Almansi

	      rStrainVector = MathUtils<double>::StrainTensorToVector( StrainMatrix , rStrainVector.size() ); //is it defined ¿?
	      }
	      break;

	    case StrainMeasure_Hencky_Material:
    	      KRATOS_ERROR(std::logic_error,"Hencky strain has no transformation coded", "");
	      break;

	    case StrainMeasure_Hencky_Spatial:
	      KRATOS_ERROR(std::logic_error,"Hencky strain has no transformation coded", "");
	      break;
	      	      
	    default:
	      KRATOS_ERROR(std::logic_error,"FINAL STRAIN NOT DEFINED in StrainTransformation", "");
	      break;
	    }
	  
	  break;

	case StrainMeasure_Almansi: 

	  switch(rStrainFinal)
	    {
	    case StrainMeasure_GreenLagrange:        
	      {
	      Matrix StrainMatrix = MathUtils<double>::StrainVectorToTensor( rStrainVector ); //is it defined ¿?
	      
	      CoVariantPullBack (StrainMatrix,rF);  //GreenLagrange

	      rStrainVector = MathUtils<double>::StrainTensorToVector( StrainMatrix , rStrainVector.size() ); //is it defined ¿?
	      }
	      break;
	      
	    case StrainMeasure_Almansi:
	      break;
	      
	    case StrainMeasure_Hencky_Material:
    	      KRATOS_ERROR(std::logic_error,"Hencky strain has no transformation coded", "");
	      break;

	    case StrainMeasure_Hencky_Spatial:
	      KRATOS_ERROR(std::logic_error,"Hencky strain has no transformation coded", "");
	      break;

	    default:
	      KRATOS_ERROR(std::logic_error,"FINAL STRAIN NOT DEFINED in StrainTransformation", "");
	      break;
	    }

	  break;

	case StrainMeasure_Hencky_Material:
	  KRATOS_ERROR(std::logic_error,"Hencky strain has no transformation coded", "");
	  break;

	case StrainMeasure_Hencky_Spatial:
	  KRATOS_ERROR(std::logic_error,"Hencky strain has no transformation coded", "");
	  break;

	default:
	  KRATOS_ERROR(std::logic_error,"Measure of strain NOT DEFINED in Strains Transformation", "");
	  break;
	}


      return rStrainVector;

    }


    /**
     * Methods to transform stress Matrices:
     * @param rStressMatrix the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressInitial the measure of stress of the given  rStressMatrix 
     * @param rStressFinal the measure of stress of the returned rStressMatrix
     */
     Matrix& ConstitutiveLaw::TransformStresses (Matrix& rStressMatrix,
						 const Matrix &rF,
						 const double &rdetF,
						 StressMeasure rStressInitial,
						 StressMeasure rStressFinal)
    {
      Vector StressVector;

      StressVector = MathUtils<double>::StressTensorToVector( rStressMatrix );

      StressVector=TransformStresses(StressVector,rF,rdetF,rStressInitial,rStressFinal);

      rStressMatrix = MathUtils<double>::StressVectorToTensor( StressVector );
      
      return rStressMatrix;
    }


    /**
     * Methods to transform stress Vectors:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressInitial the measure of stress of the given  rStressVector
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
     Vector& ConstitutiveLaw::TransformStresses (Vector& rStressVector,
						 const Matrix &rF,
						 const double &rdetF,
						 StressMeasure rStressInitial,
						 StressMeasure rStressFinal)
    {
      
      switch(rStressInitial)
	{
        case StressMeasure_PK1:         

	  TransformPK1Stresses(rStressVector,rF,rdetF,rStressFinal);
	  
	  break;

	case StressMeasure_PK2: 

	  TransformPK2Stresses(rStressVector,rF,rdetF,rStressFinal);

	  break;
	  
	case StressMeasure_Kirchhoff:

	  TransformKirchhoffStresses(rStressVector,rF,rdetF,rStressFinal);

	  break;

	case StressMeasure_Cauchy:

	  TransformCauchyStresses(rStressVector,rF,rdetF,rStressFinal);

	  break;
	  
	default:
	  KRATOS_ERROR(std::logic_error,"INITIAL STRESS NOT DEFINED in StressTransformation", "");
	  break;
	}


      return rStressVector;

    }


    /**
     * Methods to transform stress Vectors specialized with the initial stress Measure PK1:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
    Vector& ConstitutiveLaw::TransformPK1Stresses (Vector& rStressVector,
						  const Matrix &rF,
						  const double &rdetF,
						  StressMeasure rStressFinal)
     {
       unsigned int size = rF.size1(); //WorkingSpaceDimension();

       switch(rStressFinal)
	 {
	 case StressMeasure_PK1:         
	   break;
	      
	 case StressMeasure_PK2:
	   {
	   Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );
	   Matrix InvF ( size , size );
	   double J;
	   MathUtils<double>::InvertMatrix( rF, InvF , J );
	      
	   StressMatrix = prod( InvF, StressMatrix ); //PK2

	   rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() ); 
	   }
	   break;
	      
	 case StressMeasure_Kirchhoff:
	   {
	   Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );
	   Matrix InvF ( size , size );
	   double J;
	   MathUtils<double>::InvertMatrix( rF, InvF , J );
	      
	   StressMatrix = prod( InvF, StressMatrix ); //PK2

	   ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

	   rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix , rStressVector.size() );
	   }
	   break;
	      
	 case StressMeasure_Cauchy:
	   {
	   Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );
	   Matrix InvF ( size , size );
	   double J;
	   MathUtils<double>::InvertMatrix( rF, InvF , J );
	      
	   StressMatrix = prod( InvF, StressMatrix ); //PK2

	   ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

	   StressMatrix/=J; //Cauchy
	      
	   rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() ); //is it defined ¿?
	   }
	   break;
	      
	 default:
	   KRATOS_ERROR(std::logic_error,"FINAL STRESS NOT DEFINED in StressTransformation", "");
	   break;
	 }


       return rStressVector;

     }

    /**
     * Methods to transform stress Vectors specialized with the initial stress Measure PK2:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
   Vector& ConstitutiveLaw::TransformPK2Stresses (Vector& rStressVector,
				 const Matrix &rF,
				 const double &rdetF,
				 StressMeasure rStressFinal)
     {

       switch(rStressFinal)
	 {
	 case StressMeasure_PK1:        
	   {
	   Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

	   StressMatrix = prod( rF, StressMatrix ); //PK1

	   rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix , rStressVector.size() ); 
	   }
	   break;
	      
	 case StressMeasure_PK2:
	   break;
	      
	 case StressMeasure_Kirchhoff:
	   {
	   Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

	   ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

	   rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix , rStressVector.size() ); 
	   }
	   break;
	      
	 case StressMeasure_Cauchy:
	   {
	     
	   Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );
     
	   ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

	   if(rdetF!=0)
	     StressMatrix/=rdetF; //Cauchy
	      
	   rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix , rStressVector.size() ); 

	   }
	   break;
	      
	 default:
	   KRATOS_ERROR(std::logic_error,"FINAL STRESS NOT DEFINED in StressTransformation", "");
	   break;
	 }

      return rStressVector;

     }

    /**
     * Methods to transform stress Vectors specialized with the initial stress Measure Kirchooff:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
    Vector& ConstitutiveLaw::TransformKirchhoffStresses (Vector& rStressVector,
					const Matrix &rF,
					const double &rdetF,
					StressMeasure rStressFinal)
      {

	switch(rStressFinal)
	  {
	  case StressMeasure_PK1:        
	    {
	    Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

	    ContraVariantPullBack (StressMatrix,rF);  //PK2

	    StressMatrix = prod( rF, StressMatrix ); //PK1

	    rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix , rStressVector.size() ); 
	    }
	    break;
	      
	  case StressMeasure_PK2:
	    {
	    Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

	    ContraVariantPullBack (StressMatrix,rF);  //PK2

	    rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix , rStressVector.size() ); 
	    }
	    break;
	      
	  case StressMeasure_Kirchhoff:
	    break;
	      
	  case StressMeasure_Cauchy:
	    {
	    if(rdetF!=0)
	      rStressVector/=rdetF; //Cauchy
	    }
	    break;
	      
	  default:
	    KRATOS_ERROR(std::logic_error,"FINAL STRESS NOT DEFINED in StressTransformation", "");
	    break;
	  }

      return rStressVector;

      }

    /**
     * Methods to transform stress Vectors specialized with the initial stress Measure Cauchy:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
    Vector& ConstitutiveLaw::TransformCauchyStresses (Vector& rStressVector,
				     const Matrix &rF,
				     const double &rdetF,
				     StressMeasure rStressFinal)
      {

	switch(rStressFinal)
	  {
	  case StressMeasure_PK1:        
	    {
	    rStressVector*=rdetF; //Kirchhoff

	    Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

	    ContraVariantPullBack (StressMatrix,rF);  //PK2

	    StressMatrix = prod( rF, StressMatrix ); //PK1
	      
	    rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix , rStressVector.size() ); 
	    }
	    break;
	      
	  case StressMeasure_PK2:
	    {
	    rStressVector*=rdetF; //Kirchhoff

	    Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

	    ContraVariantPullBack (StressMatrix,rF);  //PK2

	    rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix , rStressVector.size() ); 
	    }
	    break;
	      
	  case StressMeasure_Kirchhoff:

	    rStressVector*=rdetF; //Kirchhoff

	    break;
	      
	  case StressMeasure_Cauchy:
	    break;
	      
	  default:
	    KRATOS_ERROR(std::logic_error,"FINAL STRESS NOT DEFINED in StressTransformation", "");
	    break;
	  }

      return rStressVector;

      }


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void ConstitutiveLaw::GetLawFeatures(Features& rFeatures)
    {

	KRATOS_ERROR(std::logic_error, "Calling virtual function for GetConstitutiveLawFeatures", "");
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
     int ConstitutiveLaw::Check(const Properties& rMaterialProperties,
				const GeometryType& rElementGeometry,
				const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        return 0;
        KRATOS_CATCH("");
    }


   //*** PROTECTED METHODS: ***//

    /**
     * This method performs a contra-variant push-forward between to tensors
     * i.e. 2nd PK stress to Kirchhoff stress
     */

    void ConstitutiveLaw::ContraVariantPushForward( Matrix& rMatrix,
				   const Matrix& rF)  //i.e. 2nd PK stress to Kirchhoff stress
    {
      unsigned int size = rF.size1(); //WorkingSpaceDimension();
      Matrix temp ( size , size );
      
      noalias( temp )     = prod( rF, rMatrix );
      noalias( rMatrix )  = prod( temp, trans( rF ) );
	
    }

    /**
     * This method performs a contra-variant pull-back between to tensors
     * i.e. Kirchhoff stress to 2nd PK stress
     */

    void ConstitutiveLaw::ContraVariantPullBack( Matrix& rMatrix,
				const Matrix& rF)     //i.e. Kirchhoff stress to 2nd PK stress
    {
      unsigned int size = rF.size1(); //WorkingSpaceDimension();
      Matrix InvF ( size , size );
      double J;
      MathUtils<double>::InvertMatrix( rF, InvF , J );
    
      Matrix temp ( size , size );

      noalias( temp )    = prod( InvF, rMatrix );
      noalias( rMatrix ) = prod( temp, trans( InvF ) );     
    }
     
    /**
     * This method performs a co-variant push-forward between to tensors
     * i.e. Green-Lagrange strain to Almansi strain
     */       
      
    void ConstitutiveLaw::CoVariantPushForward( Matrix& rMatrix,
			       const Matrix& rF)      //i.e. Green-Lagrange strain to Almansi strain
    {
      unsigned int size = rF.size1(); //WorkingSpaceDimension();
      Matrix InvF ( size , size );
      double J;
      MathUtils<double>::InvertMatrix( rF, InvF , J );
      
      Matrix temp ( size , size );
      
      noalias( temp )     = prod( trans( InvF ), rMatrix );
      noalias( rMatrix )  = prod( temp, InvF );
    }

    /**
     * This method performs a co-variant pull-back between to tensors
     * i.e. Almansi strain to Green-Lagrange strain
     */

    void ConstitutiveLaw::CoVariantPullBack( Matrix& rMatrix,
			    const Matrix& rF)         //i.e. Almansi strain to Green-Lagrange strain
    {
      
      unsigned int size = rF.size1(); //WorkingSpaceDimension();
      Matrix temp ( size , size );
      
      noalias( temp )     = prod( trans( rF ), rMatrix );
      noalias( rMatrix )  = prod( temp, rF );
      
    }



    //*** OUTDATED METHODS: ***//



    /**
     * Computes the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear strain measure is used)
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param rCurrentProcessInfo current ProcessInfo instance
     * @param rMaterialProperties the material's Properties object
     * @param rElementGeometry the element's geometry
     * @param rShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
     void ConstitutiveLaw::CalculateMaterialResponse(const Vector& StrainVector,
                                           const Matrix& DeformationGradient,
                                           Vector& StressVector,
                                           Matrix& AlgorithmicTangent,
                                           const ProcessInfo& rCurrentProcessInfo,
                                           const Properties& rMaterialProperties,
                                           const GeometryType& rElementGeometry,
                                           const Vector& rShapeFunctionsValues,
                                           bool CalculateStresses,
                                           int CalculateTangent,
                                           bool SaveInternalVariables)
    {
        KRATOS_ERROR(std::logic_error, "Calling virtual function for CalculateMaterialResponse", "");
    }

    /**
     * Computes the volumetric part of the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param rCurrentProcessInfo current ProcessInfo instance
     * @param rMaterialProperties the material's Properties object
     * @param rElementGeometry the element's geometry
     * @param rShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * NOTE: the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
     void ConstitutiveLaw::CalculateVolumetricResponse(const double VolumetricStrain,
					     const Matrix& DeformationGradient,
					     double& VolumetricStress,
					     double& AlgorithmicBulk,
					     const ProcessInfo& rCurrentProcessInfo,
					     const Properties& rMaterialProperties,
					     const GeometryType& rElementGeometry,
					     const Vector& rShapeFunctionsValues,
					     bool CalculateStresses,
					     int CalculateTangent,
					     bool SaveInternalVariables)
    {
        KRATOS_ERROR(std::logic_error, "Calling virtual function for CalculateVolumetricResponse", "");
    }

    /**
     * Computes the deviatoric part of the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param rCurrentProcessInfo current ProcessInfo instance
     * @param rMaterialProperties the material's Properties object
     * TODO: add proper definition for algorithmic tangent
     */
     void ConstitutiveLaw::CalculateDeviatoricResponse(const Vector& StrainVector,
					     const Matrix& DeformationGradient,
					     Vector& StressVector,
					     Matrix& AlgorithmicTangent,
					     const ProcessInfo& rCurrentProcessInfo,
					     const Properties& rMaterialProperties,
					     const GeometryType& rElementGeometry,
					     const Vector& rShapeFunctionsValues,
					     bool CalculateStresses,
					     int CalculateTangent,
					     bool SaveInternalVariables)
    {
        KRATOS_ERROR(std::logic_error, "Calling virtual function for CalculateDeviatoricResponse", "");
    }


    // VM
     void ConstitutiveLaw::CalculateCauchyStresses(Vector& Cauchy_StressVector,
                                         const Matrix& F,
                                         const Vector& PK2_StressVector,
                                         const Vector& GreenLagrangeStrainVector)
    {
    }




} /* namespace Kratos.*/
