/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-11-04 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_ADAPTER_H_INCLUDED)
#define RVE_ADAPTER_H_INCLUDED

#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "includes/model_part.h"
#include "includes/node.h"
#include "geometries/point_3d.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "utilities/timer.h"

#include "multiscale_application.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "rve_macroscale_status.h"
#include "rve_adapter_settings.h"
#include "math_helpers.h"
#include "time_line.h"

// macro used to print useful informations.
//#define RVE_IS_VERBOSE

// macro used to activate tesing modifications to improve the calculation
// of the homogenized tangent operator by perturbation methods
#define RVE_TEST_PERFORMANCE_TANGENT_OPERATOR

// macro used to use central differences rather than forward differences
// to compute the tangent operator
//#define RVE_TANGENT_OPERATOR_PERTURBATION_TYPE__CENTRAL_DIFF

//#define RVE_THERMAL_FIRST_TEST

namespace Kratos
{

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver,
		 class TRveSettings = RveAdapterSettings_3D
         >
class RveAdapter
{

public:

	KRATOS_CLASS_POINTER_DEFINITION(RveAdapter);
	
	typedef SolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > TSolvingStrategy;

	typedef typename TSolvingStrategy::Pointer TSolvingStrategyPointer;

	typedef ConstitutiveLaw::SizeType SizeType;

	typedef ConstitutiveLaw::GeometryType GeometryType;

public:

	/**
	* Creates a new empty RveAdapter
	*/
	RveAdapter()
		: mpModelPart(ModelPart::Pointer())
		, mpStrategy(TSolvingStrategyPointer())
		, mpMacroscaleStatus(RveMacroscaleStatus::Pointer())
		, mRveGenerated(false)
		, mRveGenerationRequested(true)
		, mInitialized(false)
		, mSolutionStepFinalized(true)
		, mIntegrationErrorCode(0.0)
		, mMacroCharacteristicLength(0.0)
		, mMicroCharacteristicLength(0.0)
		, mUpdateConstitutiveTensor(true)
		, mNDofs(0)
#ifdef RVE_THERMAL_FIRST_TEST
		, mTemp(0.0)
#endif // RVE_THERMAL_FIRST_TEST

	{
	}

	/**
	* Destructor
	*/
	virtual ~RveAdapter()
	{
	}
	
public:
	
	/**
	* Returns a pointer to the model part used for the microscale calculation.
	* @return a pointer to the model part.
	*/
	inline const ModelPart::Pointer& GetModelPart()const
	{
		return mpModelPart;
	}

	/**
	* Returns a string representation of this object.
	* @return the string representation of this object.
	*/
	inline virtual std::string GetInfo()const
	{
		std::stringstream ss;
		ss << "==============================================================" << std::endl;
		ss << "RveAdapter base class:" << std::endl;
		ss << "==============================================================" << std::endl;
		ss << "Rve Generated: " << std::boolalpha << mRveGenerated << std::endl;
		ss << "Rve Generation Requested: " << std::boolalpha << mRveGenerationRequested << std::endl;
		ss << "==============================================================" << std::endl;
		ss << "Model Part:" << std::endl;
		ss << "--------------------------------------------------------------" << std::endl;
		if(mpModelPart == NULL)
			ss << "NULL" << std::endl;
		else
			ss << *mpModelPart << std::endl;
		ss << "==============================================================" << std::endl;
		ss << "Solving Strategy:" << std::endl;
		ss << "--------------------------------------------------------------" << std::endl;
		if(mpStrategy == NULL)
			ss << "NULL" << std::endl;
		else
			ss << "Solving Strategy instance" << std::endl;
		ss << "==============================================================" << std::endl;
		ss << "Macroscale Status:" << std::endl;
		ss << "--------------------------------------------------------------" << std::endl;
		if(mpMacroscaleStatus == NULL)
			ss << "NULL" << std::endl;
		else
			ss << *mpMacroscaleStatus << std::endl;
		ss << "==============================================================" << std::endl;
		return ss.str();
	}

	/**
	* Sets all the data necessary for this RveAdapter.
	* This method should be called before any other calculation method.
	* This method is meant to be called by the RveModeler when it finds out that this
	* RveAdaptors requests the the rve generation.
	* @param pNewModelPart the pointer to the new model part 
	*                      (already built by the RveModeler)
	* @param pNewSolvingStrategy the pointer to the new solving strategy
	*                            (already built by the RveModeler)
	* @param pNewMacroscaleStatus the pointer to the new RveMacroscaleStatus object 
	*                             (already assigned to the RveConditions in the input model part)
	*/
	virtual void SetRveData(const ModelPart::Pointer& pNewModelPart,
		                    const TSolvingStrategyPointer& pNewSolvingStrategy,
		                    const RveMacroscaleStatus::Pointer& pNewMacroscaleStatus)
	{
		KRATOS_TRY

		if(pNewModelPart == mpModelPart) return;

		if(pNewModelPart == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapter - The input ModelPart is NULL", "");

		if(pNewSolvingStrategy == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapter - The input SolvingStrategy is NULL", "");

		if(pNewMacroscaleStatus == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapter - The input RveMacroscaleStatus is NULL", "");

		mpModelPart = pNewModelPart;
		mpStrategy = pNewSolvingStrategy;
		mpMacroscaleStatus = pNewMacroscaleStatus;
		mRveGenerated = true;
		mRveGenerationRequested = false;

		if(mInitialized)
		{
			// In this case the initialization has been already performed,
			// so we have the macro-characteristic length,
			// and we can now compute the micro-characteristic length and assign
			// the multiplier
			InitializeRveModelPart();
			CalculateMicroCharacteristicLength();
			AssignCharacteristicLengthMultiplier();
			mNDofs = CalculateTotalNumberOfDofs();
		}

		KRATOS_CATCH("")
	}

	// ========================================================================
	//
	// HERE THE INTERFACE FOR THE CONSTITUTIVE LAW
	//
	// ========================================================================

	/**
    * @return the working space dimension of the current constitutive law
    * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
    */
    virtual SizeType WorkingSpaceDimension()
	{
		return TRveSettings::WorkingSpaceDimension();
	}

    /**
    * returns the size of the strain vector of the current constitutive law
    * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
    */
    virtual SizeType GetStrainSize()
	{
		return TRveSettings::GetStrainSize();
	}

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    */
    virtual bool Has(const Variable<double>& rThisVariable)
	{
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			return true;
#ifdef RVE_THERMAL_FIRST_TEST
		if(rThisVariable == TEMPERATURE)
			return true;
#endif // RVE_THERMAL_FIRST_TEST
		return false;
	}

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    */
    virtual bool Has(const Variable<Vector>& rThisVariable)
	{
		return false;
	}

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    */
    virtual bool Has(const Variable<Matrix>& rThisVariable)
	{
		if(rThisVariable == PK2_STRESS_TENSOR || rThisVariable == CAUCHY_STRESS_TENSOR)
		   return true;
		return false;
	}

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    * NOTE: fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
    */
    virtual bool Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    * NOTE: fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
    */
    virtual bool Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return false;
	}

    /**
    * returns the value of a specified variable
    * @param rThisVariable the variable to be returned
    * @param rValue a reference to the returned value
    * @param rValue output: the value of the specified variable
    */
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		rValue = 0.0;
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			rValue = mIntegrationErrorCode;
#ifdef RVE_THERMAL_FIRST_TEST
		if(rThisVariable == TEMPERATURE)
			rValue = mTemp;
#endif // RVE_THERMAL_FIRST_TEST
		return rValue;
	}

    /**
    * returns the value of a specified variable
    * @param rThisVariable the variable to be returned
    * @param rValue a reference to the returned value
    * @return the value of the specified variable
    */
    virtual Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
	{
		return rValue;
	}

    /**
    * returns the value of a specified variable
    * @param rThisVariable the variable to be returned
    * @return the value of the specified variable
    */
    virtual Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		if(rThisVariable == PK2_STRESS_TENSOR || rThisVariable == CAUCHY_STRESS_TENSOR)
		{
			StressVectorToTensor(mStressVector, rValue);
		}
		return rValue;
	}

    /**
    * returns the value of a specified variable
    * @param rThisVariable the variable to be returned
    * @param rValue a reference to the returned value
    * @return the value of the specified variable
    */
    virtual array_1d<double, 3 > & GetValue(const Variable<array_1d<double, 3 > >& rVariable,
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
    virtual array_1d<double, 6 > & GetValue(const Variable<array_1d<double, 6 > >& rVariable,
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
    virtual void SetValue(const Variable<double>& rVariable,
                          const double& rValue,
                          const ProcessInfo& rCurrentProcessInfo)
	{
#ifdef RVE_THERMAL_FIRST_TEST
		if(rVariable == TEMPERATURE)
			mTemp = rValue;
#endif // RVE_THERMAL_FIRST_TEST
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
	}
 
    /**
    * sets the value of a specified variable
    * @param rVariable the variable to be returned
    * @param rValue new value of the specified variable
    * @param rCurrentProcessInfo the process info
    */
    virtual void SetValue(const Variable<Matrix >& rVariable,
                          const Matrix& rValue, 
                          const ProcessInfo& rCurrentProcessInfo)
	{
	}

    /**
    * sets the value of a specified variable
    * @param rVariable the variable to be returned
    * @param rValue new value of the specified variable
    * @param rCurrentProcessInfo the process info
    */
    virtual void SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                          const array_1d<double, 3 > & rValue,
                          const ProcessInfo& rCurrentProcessInfo)
	{
	}

    /**
    * sets the value of a specified variable
    * @param rVariable the variable to be returned
    * @param rValue new value of the specified variable
    * @param rCurrentProcessInfo the process info
    */
    virtual void SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                          const array_1d<double, 6 > & rValue,
                          const ProcessInfo& rCurrentProcessInfo)
	{
	}

    /**
    * Is called to check whether the provided material parameters in the Properties
    * match the requirements of current constitutive model.
    * @param rMaterialProperties the current Properties to be validated against.
    * @return true, if parameters are correct; false, if parameters are insufficient / faulty
    * NOTE: this has to be implemented by each constitutive model. Returns false in base class since
    * no valid implementation is contained here.
    */
    virtual bool ValidateInput(const Properties& rMaterialProperties)
	{
		// this is just to make this class as close as possible to the constitutive law interface.
		// for the moment this method seems useless, so we always return true.
		return true;
	}

    /**
    * returns the expected strain measure of this constitutive law (by default linear strains)
    * @return the expected strain measure
    */
    virtual ConstitutiveLaw::StrainMeasure GetStrainMeasure()
	{
		//KRATOS_THROW_ERROR(std::logic_error, "RveAdapter - This method should be implemented by any derived class", "");
		return ConstitutiveLaw::StrainMeasure_Infinitesimal; // for testing...
	}

    /**
    * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
    * @return the expected stress measure
    */
    virtual ConstitutiveLaw::StressMeasure GetStressMeasure()
	{
		//KRATOS_THROW_ERROR(std::logic_error, "RveAdapter - This method should be implemented by any derived class", "");
		return ConstitutiveLaw::StressMeasure_Cauchy; // for testing...
	}

    /**
    * returns whether this constitutive model is formulated in incremental strains/stresses
    * NOTE: by default, all constitutive models should be formulated in total strains
    * @return true, if formulated in incremental strains/stresses, false otherwise
    */
    virtual bool IsIncremental()
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
    virtual void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues)
	{
		if(!mInitialized)
		{
			mStressVector = ZeroVector(GetStrainSize());

			mIntegrationErrorCode = 0.0;
			mInitialized = true;

			mMacroCharacteristicLength = rElementGeometry.Length();
			if(mRveGenerated)
			{
				// In this case the rve generation has been already performed,
				// so we can now compute the micro-characteristic length and assign
				// the multiplier
				InitializeRveModelPart();
				CalculateMicroCharacteristicLength();
				AssignCharacteristicLengthMultiplier();
				mNDofs = CalculateTotalNumberOfDofs();
			}

			mUpdateConstitutiveTensor = true;
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
		const ProcessInfo& macroProcessInfo = rCurrentProcessInfo;
		ProcessInfo& microProcessInfo = mpModelPart->GetProcessInfo();

		double time   = macroProcessInfo[TIME];
		double dtime  = macroProcessInfo[DELTA_TIME];
		int    nsteps = macroProcessInfo[TIME_STEPS];

		mpModelPart->CloneTimeStep(time);

		microProcessInfo[TIME]       = time;
		microProcessInfo[DELTA_TIME] = dtime;
		microProcessInfo[TIME_STEPS] = nsteps;

		mIntegrationErrorCode = 0.0;
		mSolutionStepFinalized = false;

		mUpdateConstitutiveTensor = true;
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
		RestoreSolutionVector(mU);

		mpModelPart->GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL] = RVE_STRATEGY_FINALIZE_STEP_LEVEL___FINALIZE_ONLY;
		mpStrategy->Solve();
		mpModelPart->GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL] = RVE_STRATEGY_FINALIZE_STEP_LEVEL___DO_NOTHING;

#ifdef RVE_IS_VERBOSE
		std::stringstream ss;
		ss << " FINALIZE: N_ITER = " << mpModelPart->GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
		std::cout << ss.str();  
#endif // RVE_IS_VERBOSE

		mSolutionStepFinalized = true;
	}
 
    /**
    * to be called at the beginning of each step iteration
    * (e.g. from Element::InitializeNonLinearIteration)
    * @param rMaterialProperties the Properties instance of the current element
    * @param rElementGeometry the geometry of the current element
    * @param rShapeFunctionsValues the shape functions values in the current integration point
    * @param the current ProcessInfo instance
    */
    virtual void InitializeNonLinearIteration(const Properties& rMaterialProperties,
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
    virtual void FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                            const GeometryType& rElementGeometry,
                                            const Vector& rShapeFunctionsValues,
                                            const ProcessInfo& rCurrentProcessInfo)
	{
	}
    
    /**
    * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
    * @see Parameters
    */
    virtual void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
	{
		CalculateRveResponse(rValues);
	}
    
    /**
    * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
    * @see Parameters
    */
    virtual void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues)
	{
		CalculateRveResponse(rValues);
	}
    
    /**
    * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
    * @see Parameters
    */
    virtual void CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
	{
		CalculateRveResponse(rValues);
	}
    
    /**
    * Computes the material response in terms of Cauchy stresses and constitutive tensor
    * @see Parameters
    */
    virtual void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
	{
		CalculateRveResponse(rValues);
	}
    
    /**
    * Updates the material response in terms of 1st Piola-Kirchhoff stresses
    * @see Parameters
    */
    virtual void FinalizeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
	{
	}
    
    /**
    * Updates the material response in terms of 2nd Piola-Kirchhoff stresses
    * @see Parameters
    */
    virtual void FinalizeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues)
	{
	}

    /**
    * Updates the material response in terms of Kirchhoff stresses
    * @see Parameters
    */
    virtual void FinalizeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
	{
	}

    /**
    * Updates the material response in terms of Cauchy stresses
    * @see Parameters
    */
    virtual void FinalizeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
	{
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
		mStressVector = ZeroVector(GetStrainSize());
		mInitialized = false;
	}

    /**
    * This function is designed to be called once to check compatibility with element
    * @param rFeatures
    */
    virtual void GetLawFeatures(ConstitutiveLaw::Features& rFeatures)
	{
		TRveSettings::GetLawFeatures(rFeatures);
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
		if(mRveGenerated)
		{
			if(mpModelPart == NULL || mpStrategy == NULL || mpMacroscaleStatus == NULL)
			{
				KRATOS_THROW_ERROR(std::logic_error, "RveAdapter - RveGenerated flag set to TRUE, but some data is missing", "");
			}
		}
		return 0;
		KRATOS_CATCH("")
	}
	
public:

	/**
	* Gets a value indicating whether this RveManager has the necessary data (model part, rve macroscale status).
	* @return true, if data for this RveAdapter has been generared, false otherwise
	*/
	inline const bool RveGenerated()const { return mRveGenerated; }

	/**
	* Gets a value indicating whether this RveManager needs the generation of new data (model part, rve macroscale status).
	* @return true, if this RveAdapter needs new data, false otherwise.
	*/
	inline const bool RveGenerationRequested()const { return mRveGenerationRequested; }

	bool TestMaterialResponse(const Vector& E, bool compute_constitutive_tensor)
	{
		SizeType strain_size = GetStrainSize();
		SizeType ndim = WorkingSpaceDimension();

		mStressVector = ZeroVector(strain_size);

		std::stringstream ss;

		Vector S;
		Matrix C;
		Vector Eps(E);

		Node<3>::Pointer pNode(new Node<3>(0, 0.0, 0.0, 0.0));
		Geometry< Node<3> >::PointsArrayType nodesArray;
		nodesArray.push_back(pNode);
		GeometryType::Pointer pGeom( new Point3D< Node<3> >( nodesArray ) );
		Properties::Pointer pProp(new Properties(0));
		Matrix F( IdentityMatrix(ndim, ndim) );
		Matrix F0( IdentityMatrix(ndim, ndim) );
		double detF = 1.0;
		double detF0 = 1.0;
		Vector N(1);
		Matrix dN(1, 2, 0.0);
		N(0) = 1.0;

		ProcessInfo pinfo;
		ConstitutiveLaw::Parameters params(*pGeom, *pProp, pinfo);
		params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
		params.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, compute_constitutive_tensor);
		params.SetConstitutiveMatrix(C);
		params.SetStrainVector(Eps);
		params.SetStressVector(S);
		params.SetDeformationGradientF(F);
		params.SetDeformationGradientF0(F0);
		params.SetDeterminantF(detF);
		params.SetDeterminantF0(detF0);
		params.SetShapeFunctionsValues(N);
		params.SetShapeFunctionsDerivatives(dN);

		double tempTime = Timer::GetTime();

		// initialize
		this->InitializeMaterial(params.GetMaterialProperties(), 
			                     params.GetElementGeometry(), 
								 params.GetShapeFunctionsValues());

		// initialize solution step
		this->InitializeSolutionStep(params.GetMaterialProperties(), 
			                         params.GetElementGeometry(),
			                         params.GetShapeFunctionsValues(), 
									 params.GetProcessInfo());

		// solve
		this->CalculateRveResponse( params );

		// finalize solution step
		RestoreSolutionVector(mU);
		mpModelPart->GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL] = RVE_STRATEGY_FINALIZE_STEP_LEVEL___FINALIZE_ONLY;
		mpStrategy->Solve();
		mpModelPart->GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL] = RVE_STRATEGY_FINALIZE_STEP_LEVEL___DO_NOTHING;
		mSolutionStepFinalized = true;

		double elapsed = Timer::GetTime() - tempTime;
		ss << "ELAPSED TIME: " << elapsed << std::endl;

		ss << std::endl;
		ss << "=============================================================================" << std::endl;
		ss << "Strain:" << std::endl;
		ss << MathHelpers::VectorToString(E, 4, std::scientific);
		ss << "=============================================================================" << std::endl;

		ss << "Stress:" << std::endl;
		ss << MathHelpers::VectorToString(S, 4, std::scientific);
		ss << "=============================================================================" << std::endl;

		ss << "Constitutive Tensor:" << std::endl;
		ss << MathHelpers::MatrixToString(C, 4, std::scientific);
		ss << "=============================================================================" << std::endl;
		
		std::cout << ss.str();

		return mIntegrationErrorCode == 0.0;
	}
	
protected:

	/**
	* Calculates the Response of this RveAdapter.
	* @param rValues material parameters
	*/
	virtual void CalculateRveResponse(ConstitutiveLaw::Parameters& rValues)
	{
		Flags& options = rValues.GetOptions();
		bool compute_constitutive_tensor = options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		mIntegrationErrorCode = 0.0;

		if(CheckForOutputStage(rValues)) return;

		if(!PredictorCalculation(rValues))
		{
			// solve the micro problem
#ifdef RVE_IS_VERBOSE
			std::cout << "----------RVE CALCULATION---------------\n";  
#endif // RVE_IS_VERBOSE

			mpModelPart->GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL] = RVE_STRATEGY_FINALIZE_STEP_LEVEL___DO_NOTHING;
			SolveMicroProblem(rValues);

			if(mU.size() != mNDofs) 
				mU.resize(mNDofs,false);
			SaveSolutionVector(mU);

			if(mIntegrationErrorCode == 0.0)
			{
				ComputeMacroStressVector(rValues);
				noalias(mStressVector) = rValues.GetStressVector();

				if(compute_constitutive_tensor)
				{
#ifdef RVE_IS_VERBOSE
					std::cout << "----------RVE TANGENT OPERATOR---------------\n";  
#endif // RVE_IS_VERBOSE
					mpModelPart->GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL] = RVE_STRATEGY_FINALIZE_STEP_LEVEL___DO_NOTHING;
					ComputeConstitutiveTensor(rValues);
					mpModelPart->GetProcessInfo()[STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL] = RVE_STRATEGY_FINALIZE_STEP_LEVEL___DO_NOTHING;
					mUpdateConstitutiveTensor = false;
				}
				else
				{
					SizeType strain_size = GetStrainSize();
					Matrix& tangent = rValues.GetConstitutiveMatrix();
					if(tangent.size1() != strain_size || tangent.size2() != strain_size)
						tangent.resize(strain_size, strain_size, false);
					noalias(tangent) = mConstitutiveTensor;
				}
			}
			else
			{

#ifdef RVE_IS_VERBOSE
				std::stringstream ss;
				ss << "WARNING - RVE SOLUTION DID NOT CONVERGE" << std::endl;
				std::cout << ss.str();  
#endif // RVE_IS_VERBOSE

				SizeType strain_size = GetStrainSize();
				Vector& stress = rValues.GetStressVector();
				Matrix& tangent = rValues.GetConstitutiveMatrix();
				if(stress.size() != strain_size)
					stress.resize(strain_size, false);
				if(tangent.size1() != strain_size || tangent.size2() != strain_size)
					tangent.resize(strain_size, strain_size, false);
				noalias(stress) = ZeroVector(strain_size);
				noalias(tangent) = IdentityMatrix(strain_size, strain_size);
			}	
		}

		this->AbortSolutionStep();
	}

	/**
	* Returns true if this RVE calculation has been called only for
	* output purposes, False otherwise
	*/
	virtual bool CheckForOutputStage(ConstitutiveLaw::Parameters& rValues)
	{
		const ProcessInfo& macroProcessInfo = rValues.GetProcessInfo();

		bool stepSolved = false;
		if(macroProcessInfo.Has(STRATEGY_SOLUTION_STEP_SOLVED))
			stepSolved = (macroProcessInfo[STRATEGY_SOLUTION_STEP_SOLVED] > 0);

		if(stepSolved) 
		{
			Flags& options = rValues.GetOptions();
			bool compute_constitutive_tensor = options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
			bool compute_stress = compute_constitutive_tensor || options.Is(ConstitutiveLaw::COMPUTE_STRESS);

			SizeType strain_size = GetStrainSize();

			Vector& macroStressVector = rValues.GetStressVector();
			if(macroStressVector.size() != strain_size)
				macroStressVector.resize(strain_size, false);
			noalias(macroStressVector) = mStressVector;

			Matrix& macroConstitutiveMatrix = rValues.GetConstitutiveMatrix();
			if(macroConstitutiveMatrix.size1() != strain_size || macroConstitutiveMatrix.size2() != strain_size)
				macroConstitutiveMatrix.resize(strain_size, strain_size, false);
			noalias(macroConstitutiveMatrix) = mConstitutiveTensor;
		}

#ifdef RVE_IS_VERBOSE
		if(stepSolved)
			std::cout << "Check for output request: TRUE\n";
		else
			std::cout << "Check for output request: FALSE\n";  
#endif // RVE_IS_VERBOSE

		return stepSolved;
	}

	/**
	* Performs (if possible) a predictor step without the solution of a micro-model.
	* This test checks whether this RveAdapter can continue the calculationa without a micro-model.
	* The basic implementation of this method always return false.
	* In order to take advantage of this optimization functionality this method should be
	* re-implemented by a derived class that will perform a proper predictor check.
	* @param rValues material parameters
	* @return true, if this RveAdapter succesfully performed the calculations without the micro-model, false otherwise.
	*/
	virtual bool PredictorCalculation(ConstitutiveLaw::Parameters& rValues)
	{
		if(!mRveGenerated)
		{
			SizeType strain_size = GetStrainSize();

			Vector& StressVector = rValues.GetStressVector();
			if(StressVector.size() != strain_size)
				StressVector.resize(strain_size, false);
			noalias(StressVector) = ZeroVector(strain_size);

			Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
			if(ConstitutiveMatrix.size1() != strain_size || ConstitutiveMatrix.size2() != strain_size)
				ConstitutiveMatrix.resize(strain_size, strain_size, false);
			noalias(ConstitutiveMatrix) = ZeroMatrix(strain_size, strain_size);

			std::stringstream ss;
			ss << "RveAdapter: Cannot calculate the material response without a micromodel" << std::endl;
			std::cout << ss.str();

			// we need to solve the micro-problem
			// but we don't have a micro-model.
			// This should NEVER happen.
			// If this happens it means that there's a bug in the RveModeler
			return true;
		}

		// we need to solve the micro-problem
		return false;
	}

	/**
	* Solves the micro-problem
	* @param rValues material parameters
	*/
	virtual void SolveMicroProblem(ConstitutiveLaw::Parameters& rValues)
	{
#ifdef RVE_THERMAL_FIRST_TEST

		Vector StrainVector = rValues.GetStrainVector();

		double alpha(0.0);
		double totalVolume(0.0);
		for(ModelPart::ElementIterator it = mpModelPart->ElementsBegin(); it != mpModelPart->ElementsEnd(); ++it)
		{
			Element& ielem = *it;
			Element::GeometryType& igeom = ielem.GetGeometry();
			const Properties& props = ielem.GetProperties();
			double ialpha = props[CONDUCTIVITY];
			double ivolume = igeom.DomainSize();
			alpha += ialpha * ivolume;
			totalVolume += ivolume;
		}
		if(totalVolume > 0.0)
			alpha /= totalVolume;
		else
			alpha = 0.0;
		double T0 = 0.0; // TODO: reference temp!!!!
		double Delta_Temp = mTemp - T0;
		double therm_strain = -alpha * Delta_Temp;
		StrainVector[0] += therm_strain;
		StrainVector[1] += therm_strain;
		if(this->GetStrainSize() == 6)
			StrainVector[2] += therm_strain;

		mpMacroscaleStatus->SetStrainVector(StrainVector);

#else

		Vector& StrainVector = rValues.GetStrainVector();
		mpMacroscaleStatus->SetStrainVector(StrainVector);

#endif // RVE_THERMAL_FIRST_TEST

		mpStrategy->Solve();

#ifdef RVE_IS_VERBOSE
		std::stringstream ss;
		ss << " Solve Micro problem: N_ITER = " << mpModelPart->GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
		std::cout << ss.str();  
#endif // RVE_IS_VERBOSE

		if(mpStrategy->IsConverged() == false)
			mIntegrationErrorCode = -1.0;
	}

	/**
	* Computes the homogenized stress vector.
	* A previous call to SolveMicroProblem is required.
	* @param rValues material parameters
	*/
	virtual void ComputeMacroStressVector(ConstitutiveLaw::Parameters& rValues)
	{
		ProcessInfo& processInfo = mpModelPart->GetProcessInfo();

		SizeType strain_size = GetStrainSize();
		SizeType strain_tensor_size = strain_size == 3 ? 2 : 3;

		std::vector< std::pair<SizeType,SizeType> > tmap;
		if(strain_tensor_size == 2) {
			tmap.resize(3);
			tmap[0].first = 0; tmap[0].second = 0;
			tmap[1].first = 1; tmap[1].second = 1;
			tmap[2].first = 0; tmap[2].second = 1;
		}
		else if(strain_tensor_size == 3) {
			tmap.resize(6);
			tmap[0].first = 0; tmap[0].second = 0;
			tmap[1].first = 1; tmap[1].second = 1;
			tmap[2].first = 2; tmap[2].second = 2;
			tmap[3].first = 0; tmap[3].second = 1;
			tmap[4].first = 1; tmap[4].second = 2;
			tmap[5].first = 0; tmap[5].second = 2;
		}

		Vector& macroStressVector = rValues.GetStressVector();
		if(macroStressVector.size() != strain_size)
			macroStressVector.resize(strain_size, false);
		noalias(macroStressVector) = ZeroVector(strain_size);

		double totalVolume(0.0);

		std::vector< Matrix > stressTensors;
		
		for(ModelPart::ElementIterator it = mpModelPart->ElementsBegin(); it != mpModelPart->ElementsEnd(); ++it)
		{
			Element& ielem = *it;
			Element::GeometryType& igeom = ielem.GetGeometry();
			Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();

			ielem.GetValueOnIntegrationPoints(PK2_STRESS_TENSOR, stressTensors, processInfo);

			const GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

			if(stressTensors.size() != ipts.size()) continue;

            for(size_t point_id = 0; point_id < ipts.size(); point_id++)
			{
				double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
				Matrix& igpStressTensor = stressTensors[point_id];
				if(igpStressTensor.size1() == 1 && igpStressTensor.size2() == strain_size)
				{
					for(SizeType j = 0; j < strain_size; j++)
					{
						macroStressVector(j) += igpStressTensor(0, j) * dV;
					}
					totalVolume += dV;
				}
				else if(igpStressTensor.size1() == strain_tensor_size && igpStressTensor.size2() == strain_tensor_size)
				{
					for(SizeType k = 0; k < strain_size; k++)
					{
						SizeType i = tmap[k].first;
						SizeType j = tmap[k].second;
						macroStressVector(k) += igpStressTensor(i, j) * dV;
					}
					totalVolume += dV;
				}
				else // TEMP: to remove...
				{
					for(SizeType k = 0; k < strain_size; k++)
					{
						SizeType i = tmap[k].first;
						SizeType j = tmap[k].second;
						macroStressVector(k) += igpStressTensor(i, j) * dV;
					}
					totalVolume += dV;
				}
			}
		}

		if(totalVolume == 0.0)
			noalias(macroStressVector) = ZeroVector(strain_size);
		else
			macroStressVector /= totalVolume;
	}

#ifdef RVE_TANGENT_OPERATOR_PERTURBATION_TYPE__CENTRAL_DIFF
	/**
	* Computes the algorithmic constitutive tensor.
	* This basic implementation computes the constitutive tensor 
	* numerically by means of perturbation method.
	* @param rValues material parameters
	*/
	virtual void ComputeConstitutiveTensor(ConstitutiveLaw::Parameters& rValues)
	{
		// compute the constitutive tensor using perturbation method

		SizeType strain_size = GetStrainSize();

		Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

		// initialize the output constitutive tensor
		if(ConstitutiveMatrix.size1() != strain_size || ConstitutiveMatrix.size2() != strain_size)
			ConstitutiveMatrix.resize(strain_size, strain_size, false);
		noalias(ConstitutiveMatrix) = ZeroMatrix(strain_size, strain_size);

		// initialize (if not done yet) the stored constitutive tensor
		if(mConstitutiveTensor.size1() != strain_size || mConstitutiveTensor.size2() != strain_size)
			mConstitutiveTensor.resize(strain_size, strain_size);

		// if the update of the tangent operator is not required
		// just copy the stored one an return
		mUpdateConstitutiveTensor = true; // for the moment let's force it!
		if(!mUpdateConstitutiveTensor)
		{
			noalias(ConstitutiveMatrix) = mConstitutiveTensor;
			return;
		}

		Vector StrainVector( rValues.GetStrainVector() ); // a Copy!
		Vector StressVector( rValues.GetStressVector() ); // a Copy!

		// compute the perturbation parameters
		double stress_perturbation_threshold = 1.0E-10;
		double perturbation = 1.0E-10;

		// strain and stress perturbations
		Vector& strainVectorPerturbation = rValues.GetStrainVector(); // a Reference!
		Vector& stressVectorPerturbation = rValues.GetStressVector(); // a Reference!

		Vector S0(strain_size);
		Vector S1(strain_size);

		// compute the stress perturbation for each component of the strain vector
        for(size_t j = 0; j < strain_size; j++)
		{
			// compute the perturbed strain vector
			noalias(strainVectorPerturbation) = StrainVector;
			strainVectorPerturbation(j) = StrainVector(j) - perturbation;

			// solve the perturbed micro-problem
			SolveMicroProblem(rValues);
			ComputeMacroStressVector(rValues);
			noalias( S0 ) = rValues.GetStressVector();

			// compute the perturbed strain vector
			strainVectorPerturbation(j) = StrainVector(j) + perturbation;

			// solve the perturbed micro-problem
			SolveMicroProblem(rValues);
			ComputeMacroStressVector(rValues);
			noalias( S1 ) = rValues.GetStressVector();

			// fill the column of the tangent operator
			// with sigma/perturbation
            for(size_t i = 0; i < strain_size; i++)
			{
				double stress_perturbation_i = S1(i)-S0(i);
				if(std::abs(stress_perturbation_i) > stress_perturbation_threshold)
					ConstitutiveMatrix(i, j) = stress_perturbation_i /( 2.0*perturbation );
				else
					ConstitutiveMatrix(i, j) = 0.0;
			}
		}

		noalias(mConstitutiveTensor) = ConstitutiveMatrix;

		// reset original values
		noalias( rValues.GetStrainVector() ) = StrainVector;
		noalias( rValues.GetStressVector() ) = StressVector;

		// and reset the original (not perturbed) macro scale data
		// it can be used later on if a finalize solution step is required
		mpMacroscaleStatus->SetStrainVector(StrainVector);
	}
#else
	/**
	* Computes the algorithmic constitutive tensor.
	* This basic implementation computes the constitutive tensor 
	* numerically by means of perturbation method.
	* @param rValues material parameters
	*/
	virtual void ComputeConstitutiveTensor(ConstitutiveLaw::Parameters& rValues)
	{
		// compute the constitutive tensor using perturbation method

		SizeType strain_size = GetStrainSize();

		Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

		// initialize the output constitutive tensor
		if(ConstitutiveMatrix.size1() != strain_size || ConstitutiveMatrix.size2() != strain_size)
			ConstitutiveMatrix.resize(strain_size, strain_size, false);
		noalias(ConstitutiveMatrix) = ZeroMatrix(strain_size, strain_size);

		// initialize (if not done yet) the stored constitutive tensor
		if(mConstitutiveTensor.size1() != strain_size || mConstitutiveTensor.size2() != strain_size)
			mConstitutiveTensor.resize(strain_size, strain_size);

		// if the update of the tangent operator is not required
		// just copy the stored one an return
		mUpdateConstitutiveTensor = true; // for the moment let's force it!
		if(!mUpdateConstitutiveTensor)
		{
			noalias(ConstitutiveMatrix) = mConstitutiveTensor;
			return;
		}

		Vector StrainVector( rValues.GetStrainVector() ); // a Copy!
		Vector StressVector( rValues.GetStressVector() ); // a Copy!

		// compute the perturbation parameters
		double stress_perturbation_threshold = 1.0E-10;
		double perturbation = 1.0E-10;

		// strain and stress perturbations
		Vector& strainVectorPerturbation = rValues.GetStrainVector(); // a Reference!
		Vector& stressVectorPerturbation = rValues.GetStressVector(); // a Reference!

		// compute the stress perturbation for each component of the strain vector
        for(size_t j = 0; j < strain_size; j++)
		{
			// compute the perturbed strain vector
			noalias(strainVectorPerturbation) = StrainVector;
			strainVectorPerturbation(j) += perturbation;

			// solve the perturbed micro-problem
			SolveMicroProblem(rValues);

			// compute the macro-stress vector
			ComputeMacroStressVector(rValues);

			// subtract the converged stress vector to
			// obtain the stress variation
			noalias( stressVectorPerturbation ) -= StressVector;

			// fill the column of the tangent operator
			// with sigma/perturbation
			for(size_t i = 0; i < strain_size; i++)
			{
				double stress_perturbation_i = stressVectorPerturbation(i);
				if(std::abs(stress_perturbation_i) > stress_perturbation_threshold)
					ConstitutiveMatrix(i, j) = stress_perturbation_i / perturbation;
				else
					ConstitutiveMatrix(i, j) = 0.0;
			}
		}

		noalias(mConstitutiveTensor) = ConstitutiveMatrix;

		// reset original values
		noalias( rValues.GetStrainVector() ) = StrainVector;
		noalias( rValues.GetStressVector() ) = StressVector;

		// and reset the original (not perturbed) macro scale data
		// it can be used later on if a finalize solution step is required
		mpMacroscaleStatus->SetStrainVector(StrainVector);
	}
#endif // RVE_TANGENT_OPERATOR_PERTURBATION_TYPE__CENTRAL_DIFF

	void InitializeRveModelPart()
	{
		for(ModelPart::ElementIterator it = mpModelPart->ElementsBegin(); it != mpModelPart->ElementsEnd(); ++it)
		{
			Element& ielem = *it;
			ielem.Initialize();
		}
	}

	void AbortSolutionStep()
	{
		bool move_mesh_back = mpStrategy->MoveMeshFlag();

		for (ModelPart::NodeIterator node_iter = mpModelPart->NodesBegin(); 
			 node_iter != mpModelPart->NodesEnd(); 
			 ++node_iter)
		{
			if(move_mesh_back)
			{
				noalias((node_iter)->Coordinates()) = (node_iter)->GetInitialPosition() + 
					                                  (node_iter)->FastGetSolutionStepValue(DISPLACEMENT, 1);
			}

			ModelPart::NodeType& iNode = *node_iter;

			for(ModelPart::NodeType::DofsContainerType::iterator dof_iter = iNode.GetDofs().begin(); 
				dof_iter != iNode.GetDofs().end(); 
				++dof_iter)
			{
				ModelPart::DofType& iDof = *dof_iter;
				iDof.GetSolutionStepValue() = iDof.GetSolutionStepValue(1);
			}
		}
	}

	void SaveSolutionVector(Vector& U)
	{
		for (ModelPart::NodeIterator node_iter = mpModelPart->NodesBegin(); 
			 node_iter != mpModelPart->NodesEnd(); 
			 ++node_iter)
		{
			ModelPart::NodeType& iNode = *node_iter;
			for(ModelPart::NodeType::DofsContainerType::iterator dof_iter = iNode.GetDofs().begin(); 
				dof_iter != iNode.GetDofs().end(); 
				++dof_iter)
			{
				ModelPart::DofType& iDof = *dof_iter;
				U(iDof.EquationId()) = iDof.GetSolutionStepValue();
			}
		}
	}

	void RestoreSolutionVector(const Vector& U)
	{
		for (ModelPart::NodeIterator node_iter = mpModelPart->NodesBegin(); 
			 node_iter != mpModelPart->NodesEnd(); 
			 ++node_iter)
		{
			ModelPart::NodeType& iNode = *node_iter;
			for(ModelPart::NodeType::DofsContainerType::iterator dof_iter = iNode.GetDofs().begin(); 
				dof_iter != iNode.GetDofs().end(); 
				++dof_iter)
			{
				ModelPart::DofType& iDof = *dof_iter;
				iDof.GetSolutionStepValue() = U(iDof.EquationId());
			}
		}
	}

	void StressVectorToTensor(const Vector& a, Matrix& b)
	{
		SizeType vsize = a.size();
		if (vsize==3)
        {
			if(b.size1() != 2 || b.size2() != 2) b.resize(2,2,false);
			b(0,0) = a[0];
			b(0,1) = a[2];
			b(1,0) = a[2];
			b(1,1) = a[1];
        }
		else if (vsize==4)
        {
			if(b.size1() != 3 || b.size2() != 3) b.resize(3,3,false);
			b(0,0) = a[0];
			b(0,1) = a[3];
			b(0,2) = 0.0;
			b(1,0) = a[3];
			b(1,1) = a[1];
			b(1,2) = 0.0;
			b(2,0) = 0.0;
			b(2,1) = 0.0;
			b(2,2) = a[2];
		}
		else if (vsize==6)
        {
			if(b.size1() != 3 || b.size2() != 3) b.resize(3,3,false);
			b(0,0) = a[0];
			b(0,1) = a[3];
			b(0,2) = a[5];
			b(1,0) = a[3];
			b(1,1) = a[1];
			b(1,2) = a[4];
			b(2,0) = a[5];
			b(2,1) = a[4];
			b(2,2) = a[2];
        }
	}

	size_t CalculateTotalNumberOfDofs()
	{
		size_t n(0);
		for (ModelPart::NodeIterator node_iter = mpModelPart->NodesBegin(); 
			 node_iter != mpModelPart->NodesEnd(); 
			 ++node_iter)
		{
			ModelPart::NodeType& iNode = *node_iter;
			n += iNode.GetDofs().size();
		}
		return n;
	}

	void CalculateMicroCharacteristicLength()
	{
		mMicroCharacteristicLength = 0.0;
		
		for(ModelPart::ElementIterator it = mpModelPart->ElementsBegin(); it != mpModelPart->ElementsEnd(); ++it)
		{
			Element& ielem = *it;
			Element::GeometryType& igeom = ielem.GetGeometry();
			Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();

			const GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

			for(size_t point_id = 0; point_id < ipts.size(); point_id++)
			{
				double dV = igeom.DeterminantOfJacobian(point_id, intmethod) * ipts[point_id].Weight();
				mMicroCharacteristicLength += dV;
			}
		}

		SizeType ndim = WorkingSpaceDimension();
		if(ndim == 2)
		{
			mMicroCharacteristicLength = std::sqrt(mMicroCharacteristicLength);
		}
		else // if(ndim == 3)
		{
			mMicroCharacteristicLength = std::pow(mMicroCharacteristicLength, 1.0/3.0);
		}

		//mMicroCharacteristicLength = 62.0;
	}

	void AssignCharacteristicLengthMultiplier()
	{
		double chlen_mult = mMacroCharacteristicLength / mMicroCharacteristicLength;
		std::vector<double> values;

		const ProcessInfo& pinfo = mpModelPart->GetProcessInfo();
		for(ModelPart::ElementIterator it = mpModelPart->ElementsBegin(); it != mpModelPart->ElementsEnd(); ++it)
		{
			Element& ielem = *it;
			Element::GeometryType& igeom = ielem.GetGeometry();
			Element::IntegrationMethod intmethod = ielem.GetIntegrationMethod();

			const GeometryType::IntegrationPointsArrayType& ipts = igeom.IntegrationPoints(intmethod);

			if(values.size() != ipts.size())
				values.resize(ipts.size());

			std::fill(values.begin(), values.end(), chlen_mult);

			ielem.SetValueOnIntegrationPoints(CHARACTERISTIC_LENGTH_MULTIPLIER, values, pinfo);
		}
	}

protected:

	ModelPart::Pointer           mpModelPart;
	TSolvingStrategyPointer      mpStrategy;
	RveMacroscaleStatus::Pointer mpMacroscaleStatus;

	bool   mRveGenerated;
	bool   mRveGenerationRequested;

	bool   mInitialized;
	bool   mSolutionStepFinalized;

	double mIntegrationErrorCode;

	double mMacroCharacteristicLength;
	double mMicroCharacteristicLength;

	bool   mUpdateConstitutiveTensor;

	Vector mStressVector;
	Matrix mConstitutiveTensor;

	Vector mU;
	size_t mNDofs;

#ifdef RVE_THERMAL_FIRST_TEST
	double mTemp;
#endif // RVE_THERMAL_FIRST_TEST


private:

	/**
	* private copy constructor
	*/
	RveAdapter(const RveAdapter& other);
	
	/**
	* private assignment operator
	*/
	RveAdapter& operator = (const RveAdapter& other);

private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
		rSerializer.save("ModelPart", mpModelPart);
		rSerializer.save("Status", mpMacroscaleStatus);
		rSerializer.save("RveGen", mRveGenerated);
		rSerializer.save("RveGenReq", mRveGenerationRequested);
    }

    virtual void load(Serializer& rSerializer)
    {
		rSerializer.load("ModelPart", mpModelPart);
		rSerializer.load("Status", mpMacroscaleStatus);
		rSerializer.load("RveGen", mRveGenerated);
		rSerializer.load("RveGenReq", mRveGenerationRequested);
		if(mpModelPart == NULL)
		{
			mRveGenerationRequested = true;
			mRveGenerated = false;
		}
		else
		{
			mpStrategy = TSolvingStrategyPointer(new TSolvingStrategy(*mpModelPart));
			mRveGenerationRequested = false;
			mRveGenerated = true;
		}
    }
	
};

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver,
		 class TRveSettings
         >
inline std::ostream & operator << (std::ostream& rOStream, const RveAdapter<TSparseSpace, TDenseSpace, TLinearSolver, TRveSettings> & rThis)
{
	return rOStream << rThis.GetInfo();
}

} // namespace Kratos

#endif // RVE_ADAPTER_H_INCLUDED
