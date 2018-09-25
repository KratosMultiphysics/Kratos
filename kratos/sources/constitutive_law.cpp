//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Nelson Maireni Lafontaine
//                   Josep Maria Carbonell
//

#include "includes/constitutive_law.h"


namespace Kratos
{
    const unsigned int ConstitutiveLaw::msIndexVoigt3D6C [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };
    const unsigned int ConstitutiveLaw::msIndexVoigt2D4C [4][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1} };
    const unsigned int ConstitutiveLaw::msIndexVoigt2D3C [3][2] = { {0, 0}, {1, 1}, {0, 1} };


    /**
     * Flags related to the Parameters of the Constitutive Law
     */
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, USE_ELEMENT_PROVIDED_STRAIN,  0 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, COMPUTE_STRESS,               1 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, COMPUTE_CONSTITUTIVE_TENSOR,  2 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, COMPUTE_STRAIN_ENERGY,        3 );

    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, ISOCHORIC_TENSOR_ONLY,        4 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, VOLUMETRIC_TENSOR_ONLY,       5 );

    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, MECHANICAL_RESPONSE_ONLY,     6 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, THERMAL_RESPONSE_ONLY,        7 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, INCREMENTAL_STRAIN_MEASURE,   8 );

    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, INITIALIZE_MATERIAL_RESPONSE, 9 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, FINALIZE_MATERIAL_RESPONSE,  10 );


    /**
     * Flags related to the Features of the Constitutive Law
     */
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, FINITE_STRAINS,              1 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, INFINITESIMAL_STRAINS,       2 );

    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, THREE_DIMENSIONAL_LAW,       3 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, PLANE_STRAIN_LAW,            4 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, PLANE_STRESS_LAW,            5 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, AXISYMMETRIC_LAW,            6 );

    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, U_P_LAW,                     7 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, ISOTROPIC,                   8 );
    KRATOS_CREATE_LOCAL_FLAG( ConstitutiveLaw, ANISOTROPIC,                 9 );

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
    KRATOS_ERROR <<  "Called the virtual function for Clone"<< std::endl;;
}

/**
 * Create function (should be implemented by any derived class)
 * @return a pointer to a new instance of this constitutive law
 */
ConstitutiveLaw::Pointer ConstitutiveLaw::Create(Kratos::Parameters NewParameters) const
{
    const std::string& name = NewParameters["name"].GetString();
    return KratosComponents<ConstitutiveLaw>::Get(name).Clone();
}

/**
 * @return the working space dimension of the current constitutive law
 * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
 */
ConstitutiveLaw::SizeType ConstitutiveLaw::WorkingSpaceDimension()
{
    KRATOS_ERROR <<  "Called the virtual function for WorkingSpaceDimension"<< std::endl;;
}

/**
 * returns the size of the strain vector of the current constitutive law
 * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
 */
ConstitutiveLaw::SizeType ConstitutiveLaw::GetStrainSize()
{
    KRATOS_ERROR <<  "Called the virtual function for GetStrainSize"<< std::endl;;
}

/**
 * returns whether this constitutive Law has specified variable
 * @param rThisVariable the variable to be checked for
 * @return true if the variable is defined in the constitutive law
 */
bool ConstitutiveLaw::Has(const Variable<bool>& rThisVariable)
{
    return false;
}

/**
 * returns whether this constitutive Law has specified variable
 * @param rThisVariable the variable to be checked for
 * @return true if the variable is defined in the constitutive law
 */
bool ConstitutiveLaw::Has(const Variable<int>& rThisVariable)
{
    return false;
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
bool& ConstitutiveLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
int& ConstitutiveLaw::GetValue(const Variable<int>& rThisVariable, int& rValue)
{
    return rValue;
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
array_1d<double, 3 > & ConstitutiveLaw::GetValue(const Variable<array_1d<double, 3 > >& rThisVariable,
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
array_1d<double, 6 > & ConstitutiveLaw::GetValue(const Variable<array_1d<double, 6 > >& rThisVariable,
        array_1d<double, 6 > & rValue)
{
    return rValue;
}

/**
 * @brief Sets the value of a specified variable (bool)
 * @param rThisVariable the variable to be returned
 * @param Value new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
void ConstitutiveLaw::SetValue(const Variable<bool>& rThisVariable,
                               const bool& Value,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;;
}

/**
 * @brief Sets the value of a specified variable (int)
 * @param rThisVariable the variable to be returned
 * @param Value new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
void ConstitutiveLaw::SetValue(const Variable<int>& rThisVariable,
                               const int& Value,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;;
}

/**
 * @brief Sets the value of a specified variable (double)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
void ConstitutiveLaw::SetValue(const Variable<double>& rVariable,
                               const double& rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;;
}

/**
 * @brief Sets the value of a specified variable (Vector)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
void ConstitutiveLaw::SetValue(const Variable<Vector >& rVariable,
                               const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;;
}

/**
 * @brief Sets the value of a specified variable (Matrix)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
void ConstitutiveLaw::SetValue(const Variable<Matrix >& rVariable,
                               const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;;
}

/**
 * @brief Sets the value of a specified variable (array of 3 components)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
void ConstitutiveLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                               const array_1d<double, 3 > & rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;;
}

/**
 * @brief Sets the value of a specified variable (array of 6 components)
 * @param rVariable the variable to be returned
 * @param rValue new value of the specified variable
 * @param rCurrentProcessInfo the process info
 */
void ConstitutiveLaw::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                               const array_1d<double, 6 > & rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;;
}


/**
 * @brief Calculates the value of a specified variable (bool)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
bool& ConstitutiveLaw::CalculateValue(Parameters& rParameterValues, const Variable<bool>& rThisVariable, bool& rValue)
{
    return rValue;
}

/**
 * @brief Calculates the value of a specified variable (int)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
int& ConstitutiveLaw::CalculateValue(Parameters& rParameterValues, const Variable<int>& rThisVariable, int& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable (double)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
double& ConstitutiveLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable (Vector)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
Vector& ConstitutiveLaw::CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable (Matrix)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
Matrix& ConstitutiveLaw::CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
    return rValue;
}

/**
 * returns the value of a specified variable (array of 3 components)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
array_1d<double, 3 > & ConstitutiveLaw::CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue)
{
    return rValue;
}


  /**
 * returns the value of a specified variable (array of 6 components)
 * @param rParameterValues the needed parameters for the CL calculation
 * @param rThisVariable the variable to be returned
 * @param rValue a reference to the returned value
 * @param rValue output: the value of the specified variable
 */
array_1d<double, 6 > & ConstitutiveLaw::CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double, 6 > >& rVariable,
        array_1d<double, 6 > & rValue)
{
    return rValue;
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
 * @brief This is to be called at the very beginning of the calculation (this initializes in an specific integration point)
 * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
 * @param rMaterialProperties the Properties instance of the current element
 * @param rElementGeometry the geometry of the current element
 * @param ThisIntegrationMethod The integration method considered
 * @param IntegrationPointIndex The current integration point index
 */
void ConstitutiveLaw::InitializeMaterialOnIntegrationPoints(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const GeometryData::IntegrationMethod ThisIntegrationMethod,
    const IndexType IntegrationPointIndex
    )
{
    const Matrix& r_shape_functions = rElementGeometry.ShapeFunctionsValues(ThisIntegrationMethod);
    this->InitializeMaterial(rMaterialProperties, rElementGeometry, row(r_shape_functions, IntegrationPointIndex));
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
    case StressMeasure_PK1:
        CalculateMaterialResponsePK1(rValues);
        break;

    case StressMeasure_PK2:
        CalculateMaterialResponsePK2(rValues);
        break;

    case StressMeasure_Kirchhoff:
        CalculateMaterialResponseKirchhoff(rValues);
        break;

    case StressMeasure_Cauchy:
        CalculateMaterialResponseCauchy(rValues);
        break;

    default:
        KRATOS_ERROR <<  " Stress Measure not Defined "<< std::endl;;
        break;

    }
}


/**
 * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
 * @see Parameters
 */

void ConstitutiveLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for CalculateMaterialResponsePK1"<< std::endl;;
}

/**
 * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
 * @see Parameters
 */

void ConstitutiveLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for CalculateMaterialResponsePK2"<< std::endl;;
}

/**
 * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
 * @see Parameters
 */

void ConstitutiveLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for CalculateMaterialResponseKirchhoff"<< std::endl;;
}

/**
 * Computes the material response in terms of Cauchy stresses and constitutive tensor
 * @see Parameters
 */

void ConstitutiveLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for CalculateMaterialResponseCauchy"<< std::endl;;
}


    /**
     * Initialize the material response,  called by the element in InitializeSolutionStep.
     * @see Parameters
     * @see StressMeasures
     */

    void ConstitutiveLaw::InitializeMaterialResponse(Parameters& rValues,const StressMeasure& rStressMeasure)
    {
      switch(rStressMeasure)
	{
	case StressMeasure_PK1:         InitializeMaterialResponsePK1(rValues);
	  break;

	case StressMeasure_PK2:         InitializeMaterialResponsePK2(rValues);
	  break;

	case StressMeasure_Kirchhoff: 	InitializeMaterialResponseKirchhoff(rValues);
	  break;

	case StressMeasure_Cauchy:	InitializeMaterialResponseCauchy(rValues);
	  break;

	default:
	  KRATOS_THROW_ERROR(std::logic_error, " Stress Measure not Defined ", "");
	  break;

	}
    }


    /**
     * Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */

     void ConstitutiveLaw::InitializeMaterialResponsePK1 (Parameters& rValues)
    {
      KRATOS_THROW_ERROR(std::logic_error, "Calling virtual function for InitializeMaterialResponsePK1", "");
    }

    /**
     * Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */

     void ConstitutiveLaw::InitializeMaterialResponsePK2 (Parameters& rValues)
    {
      KRATOS_THROW_ERROR(std::logic_error, "Calling virtual function for InitializeMaterialResponsePK2", "");
    }

    /**
     * Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */

     void ConstitutiveLaw::InitializeMaterialResponseKirchhoff (Parameters& rValues)
    {
      KRATOS_THROW_ERROR(std::logic_error, "Calling virtual function for InitializeMaterialResponseKirchhoff", "");
    }

    /**
     * Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */

     void ConstitutiveLaw::InitializeMaterialResponseCauchy (Parameters& rValues)
    {
      KRATOS_THROW_ERROR(std::logic_error, "Calling virtual function for InitializeMaterialResponseCauchy", "");
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
    case StressMeasure_PK1:
        FinalizeMaterialResponsePK1(rValues);
        break;

    case StressMeasure_PK2:
        FinalizeMaterialResponsePK2(rValues);
        break;

    case StressMeasure_Kirchhoff:
        FinalizeMaterialResponseKirchhoff(rValues);
        break;

    case StressMeasure_Cauchy:
        FinalizeMaterialResponseCauchy(rValues);
        break;

    default:
        KRATOS_ERROR <<  " Stress Measure not Defined "<< std::endl;;
        break;
    }
}


void ConstitutiveLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for FinalizeMaterialResponsePK1"<< std::endl;;
}

/**
 * Updates the material response in terms of 2nd Piola-Kirchhoff stresses
 * @see Parameters
 */

void ConstitutiveLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for FinalizeMaterialResponsePK2"<< std::endl;;
}

/**
 * Updates the material response in terms of Kirchhoff stresses
 * @see Parameters
 */

void ConstitutiveLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for FinalizeMaterialResponseKirchhoff"<< std::endl;;
}

/**
 * Updates the material response in terms of Cauchy stresses
 * @see Parameters
 */

void ConstitutiveLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for FinalizeMaterialResponseCauchy"<< std::endl;;
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
    KRATOS_ERROR <<  "Calling virtual function for ResetMaterial"<< std::endl;;
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
            Matrix StrainMatrix = MathUtils<double>::StrainVectorToTensor( rStrainVector );

            CoVariantPushForward (StrainMatrix,rF);  //Almansi

            rStrainVector = MathUtils<double>::StrainTensorToVector( StrainMatrix, rStrainVector.size() );
        }
        break;

        case StrainMeasure_Hencky_Material:
            KRATOS_ERROR << "Hencky strain has no transformation coded"<< std::endl;;
            break;

        case StrainMeasure_Hencky_Spatial:
            KRATOS_ERROR << "Hencky strain has no transformation coded"<< std::endl;;
            break;

        default:
            KRATOS_ERROR << "FINAL STRAIN NOT DEFINED in StrainTransformation"<< std::endl;;
            break;
        }

        break;

    case StrainMeasure_Almansi:

        switch(rStrainFinal)
        {
        case StrainMeasure_GreenLagrange:
        {
            Matrix StrainMatrix = MathUtils<double>::StrainVectorToTensor( rStrainVector );

            CoVariantPullBack (StrainMatrix,rF);  //GreenLagrange

            rStrainVector = MathUtils<double>::StrainTensorToVector( StrainMatrix, rStrainVector.size() );
        }
        break;

        case StrainMeasure_Almansi:
            break;

        case StrainMeasure_Hencky_Material:
            KRATOS_ERROR << "Hencky strain has no transformation coded"<< std::endl;;
            break;

        case StrainMeasure_Hencky_Spatial:
            KRATOS_ERROR << "Hencky strain has no transformation coded"<< std::endl;;
            break;

        default:
            KRATOS_ERROR << "FINAL STRAIN NOT DEFINED in StrainTransformation"<< std::endl;;
            break;
        }

        break;

    case StrainMeasure_Hencky_Material:
        KRATOS_ERROR << "Hencky strain has no transformation coded"<< std::endl;;
        break;

    case StrainMeasure_Hencky_Spatial:
        KRATOS_ERROR << "Hencky strain has no transformation coded"<< std::endl;;
        break;

    default:
        KRATOS_ERROR << "Measure of strain NOT DEFINED in Strains Transformation"<< std::endl;;
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
        KRATOS_ERROR << "INITIAL STRESS NOT DEFINED in StressTransformation"<< std::endl;;
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
        Matrix InvF ( size, size );
        double J;
        MathUtils<double>::InvertMatrix( rF, InvF, J );

        StressMatrix = prod( InvF, StressMatrix ); //PK2

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_Kirchhoff:
    {
        Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );
        Matrix InvF ( size, size );
        double J;
        MathUtils<double>::InvertMatrix( rF, InvF, J );

        StressMatrix = prod( InvF, StressMatrix ); //PK2

        ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_Cauchy:
    {
        Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );
        Matrix InvF ( size, size );
        double J;
        MathUtils<double>::InvertMatrix( rF, InvF, J );

        StressMatrix = prod( InvF, StressMatrix ); //PK2

        ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

        StressMatrix/=J; //Cauchy

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    default:
        KRATOS_ERROR << "FINAL STRESS NOT DEFINED in StressTransformation"<< std::endl;;
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

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_PK2:
        break;

    case StressMeasure_Kirchhoff:
    {
        Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

        ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_Cauchy:
    {

        Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

        ContraVariantPushForward (StressMatrix,rF); //Kirchhoff

        if(rdetF!=0)
            StressMatrix/=rdetF; //Cauchy

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );

    }
    break;

    default:
        KRATOS_ERROR << "FINAL STRESS NOT DEFINED in StressTransformation"<< std::endl;;
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

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_PK2:
    {
        Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

        ContraVariantPullBack (StressMatrix,rF);  //PK2

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
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
        KRATOS_ERROR << "FINAL STRESS NOT DEFINED in StressTransformation"<< std::endl;;
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

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_PK2:
    {
        rStressVector*=rdetF; //Kirchhoff

        Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector );

        ContraVariantPullBack (StressMatrix,rF);  //PK2

        rStressVector = MathUtils<double>::StressTensorToVector( StressMatrix, rStressVector.size() );
    }
    break;

    case StressMeasure_Kirchhoff:

        rStressVector*=rdetF; //Kirchhoff

        break;

    case StressMeasure_Cauchy:
        break;

    default:
        KRATOS_ERROR << "FINAL STRESS NOT DEFINED in StressTransformation"<< std::endl;;
        break;
    }

    return rStressVector;

}


/**
 * Methods to transform Constitutive Matrices:
 * @param rConstitutiveMatrix the constitutive matrix
 * @param rF the DeformationGradientF matrix between the configurations
 */

/**
 * This method performs a pull-back of the constitutive matrix
 */
void ConstitutiveLaw::PullBackConstitutiveMatrix ( Matrix& rConstitutiveMatrix,
        const Matrix & rF )
{
    Matrix OriginalConstitutiveMatrix = rConstitutiveMatrix;

    rConstitutiveMatrix.clear();

    Matrix InverseF ( 3, 3 );
    double detF = 0;
    MathUtils<double>::InvertMatrix( rF, InverseF, detF);

    ConstitutiveMatrixTransformation( rConstitutiveMatrix, OriginalConstitutiveMatrix, InverseF );
}


/**
 * This method performs a push-forward of the constitutive matrix
 */
void ConstitutiveLaw::PushForwardConstitutiveMatrix ( Matrix& rConstitutiveMatrix,
        const Matrix & rF )
{
    Matrix OriginalConstitutiveMatrix = rConstitutiveMatrix;

    rConstitutiveMatrix.clear();

    ConstitutiveMatrixTransformation( rConstitutiveMatrix, OriginalConstitutiveMatrix, rF );
}


/**
 * This function is designed to be called once to check compatibility with element
 * @param rFeatures
 */
void ConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{

    KRATOS_ERROR <<  "Calling virtual function for GetConstitutiveLawFeatures"<< std::endl;;
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
    Matrix temp ( size, size );

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
    Matrix InvF ( size, size );
    double J;
    MathUtils<double>::InvertMatrix( rF, InvF, J );

    Matrix temp ( size, size );

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
    Matrix InvF ( size, size );
    double J;
    MathUtils<double>::InvertMatrix( rF, InvF, J );

    Matrix temp ( size, size );

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
    Matrix temp ( size, size );

    noalias( temp )     = prod( trans( rF ), rMatrix );
    noalias( rMatrix )  = prod( temp, rF );

}


/**
 * This method performs a pull-back or a push-forward between two constitutive matrices
 */
void ConstitutiveLaw::ConstitutiveMatrixTransformation ( Matrix& rConstitutiveMatrix,
        const Matrix& rOriginalConstitutiveMatrix,
        const Matrix & rF )
{
    unsigned int size = rOriginalConstitutiveMatrix.size1();
    if(  size == 6 )
    {

        for(unsigned int i=0; i<6; i++)
        {
            for(unsigned int j=0; j<6; j++)
            {
                rConstitutiveMatrix( i, j ) = TransformConstitutiveComponent(rConstitutiveMatrix( i, j ), rOriginalConstitutiveMatrix, rF,
                                              this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
            }

        }
    }
    else if( size == 4 )
    {


        for(unsigned int i=0; i<4; i++)
        {
            for(unsigned int j=0; j<4; j++)
            {
                rConstitutiveMatrix( i, j ) = TransformConstitutiveComponent(rConstitutiveMatrix( i, j ), rOriginalConstitutiveMatrix, rF,
                                              this->msIndexVoigt2D4C[i][0], this->msIndexVoigt2D4C[i][1], this->msIndexVoigt2D4C[j][0], this->msIndexVoigt2D4C[j][1]);
            }

        }
    }
    else if( size == 3 )
    {


        for(unsigned int i=0; i<3; i++)
        {
            for(unsigned int j=0; j<3; j++)
            {
                rConstitutiveMatrix( i, j ) = TransformConstitutiveComponent(rConstitutiveMatrix( i, j ), rOriginalConstitutiveMatrix, rF,
                                              this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
            }

        }
    }


}



/**
 * This method performs a pull-back or a push-forward between two constitutive tensor components
 */
double& ConstitutiveLaw::TransformConstitutiveComponent(double & rCabcd,
        const Matrix & rConstitutiveMatrix,
        const Matrix & rF,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)

{

    rCabcd = 0;
    double Cijkl=0;

    unsigned int dimension = rF.size1();

    //Cabcd
    for(unsigned int j=0; j<dimension; j++)
    {
        for(unsigned int l=0; l<dimension; l++)
        {
            for(unsigned int k=0; k<dimension; k++)
            {
                for(unsigned int i=0; i<dimension; i++)
                {
                    //Cijkl
                    rCabcd +=rF(a,i)*rF(b,j)*rF(c,k)*rF(d,l)*GetConstitutiveComponent(Cijkl,rConstitutiveMatrix,i,j,k,l);
                }
            }
        }
    }

    return rCabcd;

}


/**
 * This method gets the constitutive tensor components
 * from a consitutive matrix supplied in voigt notation
 */
double& ConstitutiveLaw::GetConstitutiveComponent(double & rCabcd,
        const Matrix& rConstitutiveMatrix,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{
    // matrix indices
    unsigned int k=0, l= 0;

    unsigned int size = rConstitutiveMatrix.size1();

    if( size == 3 )
    {

        //index k
        for(unsigned int i=0; i<3; i++)
        {
            if( a == b )
            {
                if( this->msIndexVoigt2D3C[i][0] == a && this->msIndexVoigt2D3C[i][1] == b )
                {
                    k = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt2D3C[i][0] == a && this->msIndexVoigt2D3C[i][1] == b) ||
                        (this->msIndexVoigt2D3C[i][1] == a && this->msIndexVoigt2D3C[i][0] == b) )
                {
                    k = i;
                    break;
                }
            }
        }

        //index l
        for(unsigned int i=0; i<3; i++)
        {
            if( c == d )
            {
                if( this->msIndexVoigt2D3C[i][0] == c && this->msIndexVoigt2D3C[i][1] == d )
                {
                    l = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt2D3C[i][0] == c && this->msIndexVoigt2D3C[i][1] == d) ||
                        (this->msIndexVoigt2D3C[i][1] == c && this->msIndexVoigt2D3C[i][0] == d) )
                {
                    l = i;
                    break;
                }
            }
        }


    }
    else if( size == 4 )
    {

        //index k
        for(unsigned int i=0; i<4; i++)
        {
            if( a == b )
            {
                if( this->msIndexVoigt2D4C[i][0] == a && this->msIndexVoigt2D4C[i][1] == b )
                {
                    k = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt2D4C[i][0] == a && this->msIndexVoigt2D4C[i][1] == b) ||
                        (this->msIndexVoigt2D4C[i][1] == a && this->msIndexVoigt2D4C[i][0] == b) )
                {
                    k = i;
                    break;
                }
            }
        }

        //index l
        for(unsigned int i=0; i<4; i++)
        {
            if( c == d )
            {
                if( this->msIndexVoigt2D4C[i][0] == c && this->msIndexVoigt2D4C[i][1] == d )
                {
                    l = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt2D4C[i][0] == c && this->msIndexVoigt2D4C[i][1] == d) ||
                        (this->msIndexVoigt2D4C[i][1] == c && this->msIndexVoigt2D4C[i][0] == d) )
                {
                    l = i;
                    break;
                }
            }
        }

    }
    else if( size == 6 )
    {

        //index k
        for(unsigned int i=0; i<6; i++)
        {
            if( a == b )
            {
                if( this->msIndexVoigt3D6C[i][0] == a && this->msIndexVoigt3D6C[i][1] == b )
                {
                    k = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt3D6C[i][0] == a && this->msIndexVoigt3D6C[i][1] == b) ||
                        (this->msIndexVoigt3D6C[i][1] == a && this->msIndexVoigt3D6C[i][0] == b) )
                {
                    k = i;
                    break;
                }
            }
        }

        //index l
        for(unsigned int i=0; i<6; i++)
        {
            if( c == d )
            {
                if( this->msIndexVoigt3D6C[i][0] == c && this->msIndexVoigt3D6C[i][1] == d )
                {
                    l = i;
                    break;
                }
            }
            else
            {
                if( (this->msIndexVoigt3D6C[i][0] == c && this->msIndexVoigt3D6C[i][1] == d) ||
                        (this->msIndexVoigt3D6C[i][1] == c && this->msIndexVoigt3D6C[i][0] == d) )
                {
                    l = i;
                    break;
                }
            }
        }
    }

    rCabcd = rConstitutiveMatrix(k,l);

    return rCabcd;
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
    KRATOS_ERROR <<  "Calling virtual function for CalculateMaterialResponse"<< std::endl;;
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
    KRATOS_ERROR <<  "Calling virtual function for CalculateVolumetricResponse"<< std::endl;;
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
    KRATOS_ERROR <<  "Calling virtual function for CalculateDeviatoricResponse"<< std::endl;;
}


// VM
void ConstitutiveLaw::CalculateCauchyStresses(Vector& Cauchy_StressVector,
        const Matrix& F,
        const Vector& PK2_StressVector,
        const Vector& GreenLagrangeStrainVector)
{
}




} /* namespace Kratos.*/
