/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
 *   Last Modified by:    $Author:   Massimo Petracca$
 *   Date:                $Date:     30-10-2013$
 *   Revision:            $Revision: 1.0$
 *
 * ***********************************************************/

#if !defined(PLANESTRESS_FROM_3D_CONSTITUTIVE_LAW_ADAPTER )
#define  PLANESTRESS_FROM_3D_CONSTITUTIVE_LAW_ADAPTER

/* System includes */

/* External includes */

/* Project includes */
#include "multiscale_application_variables.h"
#include "constitutive_law_adapter.h"

#define CHECK_STRAIN_CALCULATION \
	if(rValues.GetOptions().Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) { \
		std::cout << "ERROR: Constitutive Law Adapters work only with strains. The option USE_ELEMENT_PROVIDED_STRAIN is not supported" << std::endl; \
		return; \
	}

namespace Kratos
{

/**
 * Base class of constitutive law adapters.
 */
template<class TAdapter>
class PlaneStressFrom3DConstitutiveLawAdapter : public ConstitutiveLawAdapter< TAdapter >
{

public:

	KRATOS_CLASS_POINTER_DEFINITION(PlaneStressFrom3DConstitutiveLawAdapter);

	typedef typename TAdapter::Pointer TAdapterPointer;

	typedef ConstitutiveLawAdapter<TAdapter> MyBase;

	typedef ConstitutiveLaw::Parameters Parameters;

	typedef ConstitutiveLaw::SizeType SizeType;

	typedef ConstitutiveLaw::GeometryType GeometryType;

public:

    /**
     * Constructor.
     */
	PlaneStressFrom3DConstitutiveLawAdapter(const TAdapterPointer& theAdaptee)
		: MyBase(theAdaptee)
		, mInitialized(false)
		, m_init_strain()
	{
	}

    /**
     * Destructor.
     */
    virtual ~PlaneStressFrom3DConstitutiveLawAdapter()
	{
	}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      ConstitutiveLaw::Pointer p_clone(new PlaneStressFrom3DConstitutiveLawAdapter());
     *      return p_clone;
     */
    virtual ConstitutiveLaw::Pointer Clone() const
	{
        return PlaneStressFrom3DConstitutiveLawAdapter::Pointer( new PlaneStressFrom3DConstitutiveLawAdapter( MyBase::mpAdaptee->Clone() ) );
	}

    /**
     * @return the working space dimension of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType WorkingSpaceDimension()
    {
        return 2;
    }

    /**
     * returns the size of the strain vector of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType GetStrainSize()
    {
        return 3;
    }

	/**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		if (rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
		{
			if (m_error_code != 0.0)
				return m_error_code;
			else
				return MyBase::GetValue(rThisVariable, rValue);
		}
		else
			return MyBase::GetValue(rThisVariable, rValue);
	}

	/**
	* sets the value of a specified variable
	* @param rVariable the variable to be returned
	* @param rValue new value of the specified variable
	* @param rCurrentProcessInfo the process info
	*/
	virtual void SetValue(const Variable<Vector >& rVariable,
		const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == INITIAL_STRAIN) {
			if (rValue.size() == m_init_strain.size())
				noalias(m_init_strain) = rValue;
		}
		else
		{
			MyBase::mpAdaptee->SetValue(rVariable, rValue, rCurrentProcessInfo);
		}
	}

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    virtual void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues)
	{
		MyBase::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
		if(!mInitialized)
		{
			mEz.clear();
			mEz_converged.clear();
			m_error_code = 0.0;
			m_init_strain = ZeroVector(this->GetStrainSize());
			mInitialized = true;
		}
	}

    /**
     * to be called at the beginning of each solution step
     * (e.g. from Element::InitializeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void InitializeSolutionStep(const Properties& rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const Vector& rShapeFunctionsValues,
                                        const ProcessInfo& rCurrentProcessInfo)
	{
		MyBase::InitializeSolutionStep(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
		mEz = mEz_converged;
	}

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void FinalizeSolutionStep(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo)
	{
		MyBase::FinalizeSolutionStep(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
		mEz_converged = mEz;
	}

    /**
     * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
	{
		CalculateAdaptedMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_PK1);
	}

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues)
	{
		CalculateAdaptedMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_PK2);
	}

    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
	{
		CalculateAdaptedMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Kirchhoff);
	}

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
	{
		CalculateAdaptedMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy);
	}

    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void ResetMaterial(const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry,
                               const Vector& rShapeFunctionsValues)
	{
		MyBase::ResetMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
		mInitialized = false;
		mEz.clear();
		mEz_converged.clear();
	}

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    virtual void GetLawFeatures(ConstitutiveLaw::Features& rFeatures)
	{
		MyBase::GetLawFeatures(rFeatures);
		rFeatures.mOptions.Set(ConstitutiveLaw::PLANE_STRESS_LAW);
        rFeatures.mSpaceDimension = WorkingSpaceDimension();
		rFeatures.mStrainSize = GetStrainSize();
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
    virtual int Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		MyBase::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

        if(MyBase::mpAdaptee->GetStrainSize() != 6)
			KRATOS_THROW_ERROR( std::logic_error, "PlaneStressFrom3DConstitutiveLawAdapter - the strain size of the Adaptee material should be 6", "");

		return 0;
		KRATOS_CATCH("")
	}

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


private:

    ///@name Static Member Variables
    ///@{
    ///@}

    ///@name Member Variables
    ///@{
	bool mInitialized;
	array_1d<double, 3> mEz;
	array_1d<double, 3> mEz_converged;
	double m_error_code;
	Vector m_init_strain;

    ///@}

    ///@name Private Operators
    ///@{
    ///@}

    ///@name Private Operations
    ///@{

	void CalculateAdaptedMaterialResponse(ConstitutiveLaw::Parameters& rValues, ConstitutiveLaw::StressMeasure rStressMeasure)
	{
		CHECK_STRAIN_CALCULATION;

		m_error_code = 0.0;

		// some parameters
		const int maxiter = 100;
		const double relative_tolerance = 1.0E-6;
		const double always_converged_tolerance = 1.0E-9;

		// get references (for the adapted shell material)
		/*Vector&*/ Vector StrainVector = rValues.GetStrainVector();
		noalias(StrainVector) -= m_init_strain;
		Vector& StressVector = rValues.GetStressVector();
		Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
		if(StressVector.size() != 3)
			StressVector.resize(3, false);
		if(ConstitutiveMatrix.size1() != 3 || ConstitutiveMatrix.size2() != 3)
			ConstitutiveMatrix.resize(3, 3, false);

	    // construct the parameters for the 3D adaptee
		//ConstitutiveLaw::Parameters rValues3D( rValues );
		ConstitutiveLaw::Parameters rValues3D;
		rValues3D.SetProcessInfo(rValues.GetProcessInfo());
		rValues3D.SetMaterialProperties(rValues.GetMaterialProperties());
		rValues3D.SetElementGeometry(rValues.GetElementGeometry());
		rValues3D.SetShapeFunctionsValues(rValues.GetShapeFunctionsValues());

		Vector strain_3d(6);
		Matrix tangent_3d(6,6,0.0);
		Vector stress_3d(6, 0.0);

		strain_3d(0) = StrainVector(0);
		strain_3d(1) = StrainVector(1);
		strain_3d(2) = mEz(0);
		strain_3d(3) = StrainVector(2);
		strain_3d(4) = mEz(1);
		strain_3d(5) = mEz(2);

		rValues3D.SetStrainVector(strain_3d);
		rValues3D.SetStressVector(stress_3d);
		rValues3D.SetConstitutiveMatrix(tangent_3d);

		Flags& options3D = rValues3D.GetOptions();
		options3D.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
		options3D.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

		// begin
		int iter(0);
        double tolerance = relative_tolerance;
		Matrix Czz(3, 3, 0.0);
		Matrix invCzz(3, 3);
		double dummy_det(0.0);
		array_1d<double, 3> Szz;
		for(iter = 0; iter < maxiter; iter++)
		{
			// calculate 3d material response
            MyBase::mpAdaptee->CalculateMaterialResponse(rValues3D, rStressMeasure);

			// copy the condensed components
			Czz(0,0) = tangent_3d(2,2);  Czz(0,1) = tangent_3d(2,4);  Czz(0,2) = tangent_3d(2,5);
			Czz(1,0) = tangent_3d(4,2);  Czz(1,1) = tangent_3d(4,4);  Czz(1,2) = tangent_3d(4,5);
			Czz(2,0) = tangent_3d(5,2);  Czz(2,1) = tangent_3d(5,4);  Czz(2,2) = tangent_3d(5,5);
			Szz(0) = stress_3d(2);
			Szz(1) = stress_3d(4);
			Szz(2) = stress_3d(5);

			// initialize tolerance
			double Szz_norm = norm_2(Szz);
			if(iter == 0) {
				tolerance = Szz_norm * relative_tolerance;
				if(tolerance < always_converged_tolerance)
					tolerance = always_converged_tolerance;
            }

			// check convergence
			if(Szz_norm <= tolerance) {
				break;
			}

			// solve for the condensed strains
			MathUtils<double>::InvertMatrix3(Czz, invCzz, dummy_det);
			noalias(mEz) -= prod(invCzz, Szz);

			// copy the updated condensed strains
			strain_3d(2) = mEz(0);
			strain_3d(4) = mEz(1);
			strain_3d(5) = mEz(2);
		}
		if(iter >= maxiter) {
			//std::cout << "PlaneStress from 3d material adapter - Maximum iteration reached!\n";
			m_error_code = 1.0;
		}

		/*std::stringstream ss;
		ss << "Ez: " << mEz(0) << ", " << mEz(1) << ", " << mEz(2) << std::endl;
		ss << "Sz: " << Szz(0) << ", " << Szz(1) << ", " << Szz(2) << std::endl;
		std::cout << ss.str();*/

		ConstitutiveMatrix(0, 0) = tangent_3d(0, 0);   ConstitutiveMatrix(0, 1) = tangent_3d(0, 1);  ConstitutiveMatrix(0, 2) = tangent_3d(0, 3);
		ConstitutiveMatrix(1, 0) = tangent_3d(1, 0);   ConstitutiveMatrix(1, 1) = tangent_3d(1, 1);  ConstitutiveMatrix(1, 2) = tangent_3d(1, 3);
		ConstitutiveMatrix(2, 0) = tangent_3d(3, 0);   ConstitutiveMatrix(2, 1) = tangent_3d(3, 1);  ConstitutiveMatrix(2, 2) = tangent_3d(3, 3);

		Matrix L(3, 3);
		L(0, 0) = tangent_3d(2, 0);    L(0, 1) = tangent_3d(2, 1);    L(0, 2) = tangent_3d(2, 3);
		L(1, 0) = tangent_3d(4, 0);    L(1, 1) = tangent_3d(4, 1);    L(1, 2) = tangent_3d(4, 3);
		L(2, 0) = tangent_3d(5, 0);    L(2, 1) = tangent_3d(5, 1);    L(2, 2) = tangent_3d(5, 3);

		Matrix LT(3, 3);
		LT(0, 0) = tangent_3d(0, 2);   LT(0, 1) = tangent_3d(0, 4);   LT(0,  2) = tangent_3d(0, 5);
		LT(1, 0) = tangent_3d(1, 2);   LT(1, 1) = tangent_3d(1, 4);   LT(1,  2) = tangent_3d(1, 5);
		LT(2, 0) = tangent_3d(3, 2);   LT(2, 1) = tangent_3d(3, 4);   LT(2,  2) = tangent_3d(3, 5);

		MathUtils<double>::InvertMatrix3(Czz, invCzz, dummy_det);
		Matrix LTinvC(3, 3);
		noalias( LTinvC ) = prod( LT, invCzz );

		noalias(ConstitutiveMatrix) -= prod(LTinvC, L);

		// just in case of non convergence, we modify the residual
		StressVector(0) = stress_3d(0);
		StressVector(1) = stress_3d(1);
		StressVector(2) = stress_3d(3);
		noalias(StressVector) += prod(LTinvC, Szz);
		/*std::stringstream ss;
		ss << "ConstitutiveMatrix: " << ConstitutiveMatrix << std::endl;
		ss << "LT: " << stress_3d << std::endl;
		ss << "invCzz: " << stress_3d << std::endl;
		ss << "stress_3d: " << stress_3d << std::endl;
		ss << "LTinvC: " << LTinvC << std::endl;
		ss << "Szz: " << Szz << std::endl;
		ss << "StressVector: " << StressVector << std::endl;
		std::cout << ss.str();*/
	}

    ///@}

    ///@name Private  Access
    ///@{
    ///@}

    ///@name Serialization
    ///@{

    friend class Serializer;

	PlaneStressFrom3DConstitutiveLawAdapter(){}

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MyBase );
		rSerializer.save("ez", mEz);
		rSerializer.save("ezc", mEz_converged);
		rSerializer.save("ini", mInitialized);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MyBase );
		rSerializer.load("ez", mEz);
		rSerializer.load("ezc", mEz_converged);
		rSerializer.load("ini", mInitialized);
    }

    ///@}

}; /* Class PlaneStressFrom3DConstitutiveLawAdapter */

} /* namespace Kratos.*/
#endif /* PLANESTRESS_FROM_3D_CONSTITUTIVE_LAW_ADAPTER  defined */
