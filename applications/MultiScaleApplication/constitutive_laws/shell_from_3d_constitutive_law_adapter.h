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

#if !defined(KRATOS_SHELL_FROM_3D_CONSTITUTIVE_LAW_ADAPTER )
#define  KRATOS_SHELL_FROM_3D_CONSTITUTIVE_LAW_ADAPTER

/* System includes */

/* External includes */

/* Project includes */
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
class ShellFrom3DConstitutiveLawAdapter : public ConstitutiveLawAdapter< TAdapter >
{

public:

	KRATOS_CLASS_POINTER_DEFINITION(ShellFrom3DConstitutiveLawAdapter);
	
	typedef typename TAdapter::Pointer TAdapterPointer;
	
	typedef ConstitutiveLawAdapter<TAdapter> MyBase;

	typedef ConstitutiveLaw::Parameters Parameters;
	
	typedef ConstitutiveLaw::SizeType SizeType;
	
	typedef ConstitutiveLaw::GeometryType GeometryType;
	
public:

    /**
     * Constructor.
     */
	ShellFrom3DConstitutiveLawAdapter(const TAdapterPointer& theAdaptee)
		: MyBase(theAdaptee)
		, mInitialized(false)
	{
	}
	
    /**
     * Destructor.
     */
    virtual ~ShellFrom3DConstitutiveLawAdapter()
	{
	}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      ConstitutiveLaw::Pointer p_clone(new ShellFrom3DConstitutiveLawAdapter());
     *      return p_clone;
     */
    virtual ConstitutiveLaw::Pointer Clone() const
	{
        return ShellFrom3DConstitutiveLawAdapter::Pointer( new ShellFrom3DConstitutiveLawAdapter( MyBase::mpAdaptee->Clone() ) );
	}

    /**
     * @return the working space dimension of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType WorkingSpaceDimension()
    {
        return 3;
    }

    /**
     * returns the size of the strain vector of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType GetStrainSize()
    {
        return 5;
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
			mEz = 0.0;
			mEz_converged = 0.0;
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
		mEz = 0.0;
		mEz_converged = 0.0;
	}

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    virtual void GetLawFeatures(ConstitutiveLaw::Features& rFeatures)
	{
		MyBase::GetLawFeatures(rFeatures);
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
			KRATOS_THROW_ERROR( std::logic_error, "ShellFrom3DConstitutiveLawAdapter - the strain size of the Adaptee material should be 6", "");
		
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
	double mEz;
	double mEz_converged;
    ///@}
    
    ///@name Private Operators
    ///@{
    ///@}
    
    ///@name Private Operations
    ///@{

	void CalculateAdaptedMaterialResponse(ConstitutiveLaw::Parameters& rValues, ConstitutiveLaw::StressMeasure rStressMeasure)
	{
		CHECK_STRAIN_CALCULATION;

		// get references (for the adapted shell material)
		Vector& StrainVector = rValues.GetStrainVector();
		Vector& StressVector = rValues.GetStressVector();
		Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
		if(StressVector.size() != 5)
			StressVector.resize(5, false);
		if(ConstitutiveMatrix.size1() != 5 || ConstitutiveMatrix.size2() != 5)
			ConstitutiveMatrix.resize(5, 5, false);

	    // construct the parameters for the 3D adaptee
		ConstitutiveLaw::Parameters rValues3D( rValues );

		Vector strain_3d(6);
		Matrix tangent_3d(6,6,0.0);
		Vector stress_3d(6, 0.0);

		double nu = rValues.GetMaterialProperties()[POISSON_RATIO];
		mEz = (StrainVector(0)*nu + StrainVector(1)*nu)/(nu - 1.0);

		strain_3d(0) = StrainVector(0);
		strain_3d(1) = StrainVector(1);
		strain_3d(2) = mEz;
		strain_3d(3) = StrainVector(2);
		strain_3d(4) = StrainVector(3);
		strain_3d(5) = StrainVector(4);
		
		rValues3D.SetStrainVector(strain_3d);
		rValues3D.SetStressVector(stress_3d);
		rValues3D.SetConstitutiveMatrix(tangent_3d);
		
		Flags& options3D = rValues3D.GetOptions();
		options3D.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
		options3D.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

		MyBase::mpAdaptee->CalculateMaterialResponse(rValues3D, rStressMeasure);

		for(int i = 0; i < 2; i++)
			for(int j = 0; j < 2; j++)
				ConstitutiveMatrix(i,j) = tangent_3d(i,j);

		for(int i = 3; i < 6; i++)
			for(int j = 3; j < 6; j++)
				ConstitutiveMatrix(i-1,j-1) = tangent_3d(i,j);

		for(int i = 0; i < 2; i++) {
			for(int j = 3; j < 6; j++) {
				ConstitutiveMatrix(i, j-1) = tangent_3d(i,j);
				ConstitutiveMatrix(j-1, i) = tangent_3d(j,i);
			}
		}
		StressVector(0) = stress_3d(0);
		StressVector(1) = stress_3d(1);
		StressVector(2) = stress_3d(3);
		StressVector(3) = stress_3d(4);
		StressVector(4) = stress_3d(5);
	}

	//void CalculateAdaptedMaterialResponse(ConstitutiveLaw::Parameters& rValues, ConstitutiveLaw::StressMeasure rStressMeasure)
	//{
	//	CHECK_STRAIN_CALCULATION;

	//	// some parameters
	//	const int maxiter = 30;
	//	const double relative_tolerance = 1.0E-5;
	//	const double always_converged_tolerance = 1.0E-5;

	//	// get references (for the adapted shell material)
	//	Vector& StrainVector = rValues.GetStrainVector();
	//	Vector& StressVector = rValues.GetStressVector();
	//	Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	//	if(StressVector.size() != 5)
	//		StressVector.resize(5, false);
	//	if(ConstitutiveMatrix.size1() != 5 || ConstitutiveMatrix.size2() != 5)
	//		ConstitutiveMatrix.resize(5, 5, false);

	//    // construct the parameters for the 3D adaptee
	//	ConstitutiveLaw::Parameters rValues3D( rValues );

	//	Vector strain_3d(6);
	//	Matrix tangent_3d(6,6,0.0);
	//	Vector stress_3d(6, 0.0);

	//	strain_3d(0) = StrainVector(0);
	//	strain_3d(1) = StrainVector(1);
	//	strain_3d(2) = mEz; // from previous step (if any)
	//	strain_3d(3) = StrainVector(2);
	//	strain_3d(4) = StrainVector(3);
	//	strain_3d(5) = StrainVector(4);
	//	
	//	rValues3D.SetStrainVector(strain_3d);
	//	rValues3D.SetStressVector(stress_3d);
	//	rValues3D.SetConstitutiveMatrix(tangent_3d);
	//	
	//	Flags& options3D = rValues3D.GetOptions();
	//	options3D.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
	//	options3D.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

	//	// begin
	//	int iter(0);
	//	double tolerance = relative_tolerance;
	//	bool converged = false;
	//	double Czz(0.0);
	//	double Szz(0.0);
	//	for(iter = 0; iter < maxiter; iter++)
	//	{
	//		// calculate 3d material response
 //           MyBase::mpAdaptee->CalculateMaterialResponse(rValues3D, rStressMeasure);

	//		Czz = tangent_3d(2,2);
	//		Szz = stress_3d(2);

	//		// initialize tolerance
	//		if(iter == 0) {
	//			tolerance = std::abs(Szz) * relative_tolerance;
	//			if(tolerance < always_converged_tolerance)
	//				tolerance = always_converged_tolerance;
	//		}

	//		// check convergence
	//		if(std::abs(Szz) <= tolerance) {
	//			converged = true;
	//			break;
	//		}

	//		mEz -= Szz / Czz;
	//		strain_3d(2) = mEz;
	//	}
	//	/*if(!converged) {
	//		std::cout << "Shell from 3d material adapter - Maximum iteration reached! : " << Szz << "\n";
	//	}*/
	//	if(converged)
	//	{
	//		Szz = 0.0;
	//	}

	//	for(int i = 0; i < 2; i++)
	//		for(int j = 0; j < 2; j++)
	//			ConstitutiveMatrix(i,j) = tangent_3d(i,j);

	//	for(int i = 3; i < 6; i++)
	//		for(int j = 3; j < 6; j++)
	//			ConstitutiveMatrix(i-1,j-1) = tangent_3d(i,j);

	//	for(int i = 0; i < 2; i++) {
	//		for(int j = 3; j < 6; j++) {
	//			ConstitutiveMatrix(i, j-1) = tangent_3d(i,j);
	//			ConstitutiveMatrix(j-1, i) = tangent_3d(j,i);
	//		}
	//	}

	//	double invCzz = 1.0 / Czz;

	//	Vector L(5);
	//	L(0) = tangent_3d(2,0);
	//	L(1) = tangent_3d(2,1);
	//	L(2) = tangent_3d(2,3);
	//	L(3) = tangent_3d(2,4);
	//	L(4) = tangent_3d(2,5);

	//	Vector LTinvC(5);
	//	LTinvC(0) = tangent_3d(0,2) * invCzz;
	//	LTinvC(1) = tangent_3d(1,2) * invCzz;
	//	LTinvC(2) = tangent_3d(3,2) * invCzz;
	//	LTinvC(3) = tangent_3d(4,2) * invCzz;
	//	LTinvC(4) = tangent_3d(5,2) * invCzz;

	//	noalias(ConstitutiveMatrix) -= outer_prod(LTinvC, L);

	//	// just in case of non convergence, we modify the residual
	//	StressVector(0) = stress_3d(0) + LTinvC(0) * Szz;
	//	StressVector(1) = stress_3d(1) + LTinvC(1) * Szz;
	//	StressVector(2) = stress_3d(3) + LTinvC(2) * Szz;
	//	StressVector(3) = stress_3d(4) + LTinvC(3) * Szz;
	//	StressVector(4) = stress_3d(5) + LTinvC(4) * Szz;
	//}



    ///@}
    
    ///@name Private  Access
    ///@{
    ///@}

    ///@name Serialization
    ///@{

    friend class Serializer;

	ShellFrom3DConstitutiveLawAdapter(){}

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

}; /* Class ShellFrom3DConstitutiveLawAdapter */

} /* namespace Kratos.*/
#endif /* KRATOS_SHELL_FROM_3D_CONSTITUTIVE_LAW_ADAPTER  defined */
