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

#if !defined(RVE_ADAPTER_V2_H_INCLUDED)
#define RVE_ADAPTER_V2_H_INCLUDED

#include <limits>
#include <string>
#include <fstream>
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
#include "rve_adapter_settings.h"
#include "math_helpers.h"
#include "time_line.h"

#include "rve_macroscale_data.h"
#include "rve_macroscale_temperature_data.h"
#include "rve_geometry_descriptor.h"
#include "rve_constraint_handler.h"
#include "rve_linear_system_of_equations.h"
#include "rve_homogenizer.h"
#include "rve_utilities.h"

#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//#include "simplexnoise.h"
#include "imperfection_utilities.h"

// macro used to print useful informations.
//#define RVE_V2_IS_VERBOSE
//#define RVE_TIMER_ON

// flag to use or not the linesearch
//#define RVE_V2_USE_LINE_SEARCH

// flag to activate the modified newton with a krylov subspace accelerator
//#define RVE_V2_USES_KRYLOV_NEWTON

//#define RVE_V2_CHECK_HOMG_STRAIN

//#define RVE_V2_INCREMENTAL_REGULARIZATION

#ifdef RVE_V2_USES_KRYLOV_NEWTON
#include "krylov_subspace_accelerator_eigenlib_utilities.h"
#endif // RVE_V2_USES_KRYLOV_NEWTON


namespace Kratos
{


	inline double rve_det_acoustic_tensor_2d(const Matrix& C, double a)
	{
		double sinA = std::sin(a);
		double cosA = std::cos(a);
		double a1 = sinA*sinA;
		double a2 = cosA*cosA;
		double a3 = C(1, 2);
		double a4 = C(2, 0);
		double a5 = C(2, 1);
		double a6 = C(0, 2);
		double a7 = C(1, 0);
		double a8 = C(0, 1);
		double a9 = C(1, 1);
		double a10 = C(0, 0);
		double a11 = C(2, 2);
		double a12 = a2*a2;
		double a13 = a1*a1;
		double detQ = 
				(a1*a4*a8*cosA + a2*a3*a8*cosA - a1*a3*a10*cosA + 
				a1*a6*a7*cosA + a2*a5*a7*cosA - a2*a4*a9*cosA - 
				a1*a5*a10*cosA - a2*a6*a9*cosA)*sinA + 
				a9*a11*a12 - a4*a6*a13 - a3*a5*a12 + a10*a11*a13 + 
				a1*a2*a3*a4 + a1*a2*a5*a6 - a1*a2*a7*a8 - 
				a1*a2*a7*a11 - a1*a2*a8*a11 + a1*a2*a9*a10;
		return detQ;
	}

	inline double rve_min_det_acoustic_tensor_2d(const Matrix& C, double& angle)
	{
		// this returns the absolute minimum
		/*double dq = std::numeric_limits<double>::max();
		angle = 0.0;
		for(int i = 0; i <= 360.0; i++)
		{
			double a = double(i)/180.0*KRATOS_M_PI;
			double dq_trial = rve_det_acoustic_tensor_2d(C,a);
			if(dq_trial < dq) {
				dq = dq_trial;
				angle = a;
			}
		}
		return dq;*/

		// this returns the local maximum of the negatives
		double dq = -std::numeric_limits<double>::max();
		bool found=false;
		angle = 0.0;
		for(int i = 0; i <= 360.0; i++)
		{
			double a = double(i)/180.0*KRATOS_M_PI;
			double dq_trial = rve_det_acoustic_tensor_2d(C,a);
			if(dq_trial < 0.0) {
				if(dq_trial > dq) {
					dq = dq_trial;
					angle = a;
					found=true;
				}
			}
		}
		if(!found)
			dq = 1.0;
		return dq;
	}


template<class TSparseSpace,
         class TDenseSpace,
		 class TRveSettings = RveAdapterSettings_3D
         >
class RveAdapterV2
{

public:

	KRATOS_CLASS_POINTER_DEFINITION(RveAdapterV2);

	typedef ConstitutiveLaw::SizeType SizeType;
	typedef ConstitutiveLaw::GeometryType GeometryType;

	typedef typename TSparseSpace::MatrixType SparseMatrixType;
	typedef typename TSparseSpace::VectorType VectorType;
	typedef typename TDenseSpace::MatrixType  DenseMatrixType;

	typedef ModelPart::DofsArrayType DofsArrayType;

	typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
	typedef typename SchemeType::Pointer SchemePointerType;
	typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > ConvergenceCriteriaType;
	typedef typename ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointerType;

	typedef RveConstraintHandler< TSparseSpace, TDenseSpace > RveConstraintHandlerType;
	typedef typename RveConstraintHandlerType::Pointer RveConstraintHandlerPointerType;

	typedef RveLinearSystemOfEquations< TSparseSpace, TDenseSpace > RveLinearSystemOfEquationsType;
	typedef typename RveLinearSystemOfEquationsType::Pointer RveLinearSystemOfEquationsPointerType;

	typedef RveHomogenizer< TSparseSpace, TDenseSpace > RveHomogenizerType;
	typedef typename RveHomogenizerType::Pointer RveHomogenizerPointerType;

public:

	/**
	* Creates a new empty RveAdapterV2
	*/
	RveAdapterV2()
		: mpModelPart(ModelPart::Pointer())
		, mpMacroscaleData(RveMacroscaleData::Pointer())
		, mpGeometryDescriptor(RveGeometryDescriptor::Pointer())
		, mpConstraintHandler(RveConstraintHandlerPointerType())
		, mpLinearSOE(RveLinearSystemOfEquationsPointerType())
		, mpHomogenizer(RveHomogenizerPointerType())
		, mpScheme(SchemePointerType())
		, mpConvergenceCriteria(ConvergenceCriteriaPointerType())
		, mRveGenerated(false)
		, mRveGenerationRequested(true)
		, mMacroCharacteristicLength(0.0)
		, mMicroCharacteristicLength(0.0)
		, mInitialized(false)
		, mSolutionStepFinalized(true)
		, mIntegrationErrorCode(0.0)
		, mMoveMesh(false)
		, mMaxIterationNumber(10)
		, mU(Vector())
		, mNDofs(0)
		, m_macro_imperfection_factor(0.0)
	{
	}

	/**
	* Destructor
	*/
	virtual ~RveAdapterV2()
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
		ss << "RveAdapterV2 base class:" << std::endl;
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
		ss << "Macroscale Status:" << std::endl;
		ss << "--------------------------------------------------------------" << std::endl;
		if(mpMacroscaleData == NULL)
			ss << "NULL" << std::endl;
		else
			ss << *mpMacroscaleData << std::endl;
		ss << "==============================================================" << std::endl;
		return ss.str();
	}

	/**
	* Sets all the data necessary for this RveAdapterV2.
	* This method should be called before any other calculation method.
	* This method is meant to be called by the RveModeler when it finds out that this
	* RveAdaptors requests the the rve generation.
	*/
	virtual void SetRveData(const ModelPart::Pointer& pNewModelPart,
		                    const RveMacroscaleData::Pointer& pNewMacroScaleData,
		                    const RveGeometryDescriptor::Pointer& pNewGeometryDescriptor,
							const RveConstraintHandlerPointerType& pNewConstraintHandler,
							const RveLinearSystemOfEquationsPointerType& pNewLinearSOE,
							const RveHomogenizerPointerType& pNewHomogenizer,
							const SchemePointerType& pNewScheme,
							const ConvergenceCriteriaPointerType& pNewConvergenceCriteria
							)
	{
		KRATOS_TRY

#ifdef RVE_V2_IS_VERBOSE
		std::cout << "SET RVE DATA - BEGIN\n";
#endif // RVE_V2_IS_VERBOSE

		if(pNewModelPart == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - The input ModelPart is NULL", "");

		if(pNewMacroScaleData == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - The input MacroScaleData is NULL", "");

		if(pNewGeometryDescriptor == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - The input GeometryDescriptor is NULL", "");

		if(pNewConstraintHandler == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - The input ConstraintHandler is NULL", "");

		if(pNewLinearSOE == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - The input LinearSystemOfEquations is NULL", "");

		if(pNewHomogenizer == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - The input Homogenizer is NULL", "");

		if(pNewScheme == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - The input Scheme is NULL", "");

		if(pNewConvergenceCriteria == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - The input Convergence criteria is NULL", "");

		mpModelPart           = pNewModelPart;
		mpMacroscaleData      = pNewMacroScaleData;
		mpGeometryDescriptor  = pNewGeometryDescriptor;
		mpConstraintHandler   = pNewConstraintHandler;
		mpLinearSOE           = pNewLinearSOE;
		mpHomogenizer         = pNewHomogenizer;
		mpScheme              = pNewScheme;
		mpConvergenceCriteria = pNewConvergenceCriteria;

		mRveGenerated = true;
		mRveGenerationRequested = false;

		if(mInitialized)
		{
			// In this case the initialization has been already performed,
			// so we have the macro-characteristic length,
			// and we can now compute the micro-characteristic length and assign
			// the multiplier
			InitializeRve();
			CalculateMicroCharacteristicLength();
			AssignCharacteristicLengthMultiplier();
			mNDofs = RveUtilities::CalculateTotalNumberOfDofs(*mpModelPart);
			if(mU.size() != mNDofs) mU.resize(mNDofs,false);
			mU.clear();
		}

#ifdef RVE_V2_IS_VERBOSE
		std::cout << "SET RVE DATA - END\n";
#endif // RVE_V2_IS_VERBOSE

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
		return true;
	}

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    */
    virtual bool Has(const Variable<Vector>& rThisVariable)
	{
		return true;
	}

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    */
    virtual bool Has(const Variable<Matrix>& rThisVariable)
	{
		return true;
	}

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    * NOTE: fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
    */
    virtual bool Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return true;
	}

    /**
    * returns whether this constitutive Law has specified variable
    * @param rThisVariable the variable to be checked for
    * @return true if the variable is defined in the constitutive law
    * NOTE: fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
    */
    virtual bool Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return true;
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
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE) {
			rValue = mIntegrationErrorCode;
		}
		else {
			mpHomogenizer->HomogenizeVariable(*mpModelPart, *mpGeometryDescriptor, rThisVariable, rValue);
		}
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
		if(rThisVariable == DISCONTINUITY_DIRECTION) {
			rValue = ZeroVector(3);
			double disc_len = 1.0;
			double angle = 0.0;
			double det_Qn = rve_min_det_acoustic_tensor_2d(mConstitutiveTensor,angle);
			if(det_Qn < 0.0) {
				rValue(0) = disc_len*std::cos(angle);
				rValue(1) = disc_len*std::sin(angle);
				m_localized = true;
				m_localization_angle = angle;
			}
			else {
				if(m_localized) {
					rValue(0) = disc_len*std::cos(m_localization_angle);
					rValue(1) = disc_len*std::sin(m_localization_angle);
				}
			}
		}
		else {
			mpHomogenizer->HomogenizeVariable(*mpModelPart, *mpGeometryDescriptor, rThisVariable, rValue);
		}
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
			StressVectorToTensor(mStressVector, rValue);
		else if(rThisVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
			StrainVectorToTensor(mStrainVectorOld, rValue);
		else
			mpHomogenizer->HomogenizeVariable(*mpModelPart, *mpGeometryDescriptor, rThisVariable, rValue);
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
		mpHomogenizer->HomogenizeVariable(*mpModelPart, *mpGeometryDescriptor, rVariable, rValue);
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
		mpHomogenizer->HomogenizeVariable(*mpModelPart, *mpGeometryDescriptor, rVariable, rValue);
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
		//KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - This method should be implemented by any derived class", "");
		return ConstitutiveLaw::StrainMeasure_Infinitesimal; // for testing...
	}

    /**
    * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
    * @return the expected stress measure
    */
    virtual ConstitutiveLaw::StressMeasure GetStressMeasure()
	{
		//KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - This method should be implemented by any derived class", "");
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
#ifdef RVE_V2_IS_VERBOSE
		std::cout << "INITIALIZE MATERIAL - BEGIN\n";
#endif // RVE_V2_IS_VERBOSE
		if(!mInitialized)
		{
			SizeType strain_size = GetStrainSize();
			mFirstStepDone = false;
			mStressVector = ZeroVector(strain_size);
			mStrainVectorOld = ZeroVector(strain_size);
			mStrainVectorOld_trial = ZeroVector(strain_size);
			mConstitutiveTensor = ZeroMatrix(strain_size,strain_size);

			mIntegrationErrorCode = 0.0;
			mInitialized = true;

			mMacroCharacteristicLength = rElementGeometry.Length();

			double impf_scale_factor = ImperfectionUtilties::CalculateRandomImperfectionScaleFactor(rElementGeometry, rShapeFunctionsValues);
			m_macro_imperfection_factor = 1.0 - impf_scale_factor;

			if(mRveGenerated)
			{
				// In this case the rve generation has been already performed,
				// so we can now compute the micro-characteristic length and assign
				// the multiplier
				InitializeRve();
				CalculateMicroCharacteristicLength();
				AssignCharacteristicLengthMultiplier();
				mNDofs = RveUtilities::CalculateTotalNumberOfDofs(*mpModelPart);
				if(mU.size() != mNDofs) mU.resize(mNDofs,false);
				mU.clear();
			}
		}
#ifdef RVE_V2_IS_VERBOSE
		std::cout << "INITIALIZE MATERIAL - END\n";
#endif // RVE_V2_IS_VERBOSE
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
#ifdef RVE_V2_IS_VERBOSE
		std::cout << "RVE INITIALIZE SOLUTION STEP\n"; 
#endif // RVE_V2_IS_VERBOSE
		const ProcessInfo& macroProcessInfo = rCurrentProcessInfo;
		ProcessInfo& microProcessInfo = mpModelPart->GetProcessInfo();

		mStrainVectorOld_trial = mStrainVectorOld;

		double time   = macroProcessInfo[TIME];
		double dtime  = macroProcessInfo[DELTA_TIME];
		int    nsteps = macroProcessInfo[TIME_STEPS];

		if(mSolutionStepFinalized /*&& mSequentialID == 0*/)
			mpModelPart->CloneTimeStep(time);

		microProcessInfo[TIME]       = time;
		microProcessInfo[DELTA_TIME] = dtime;
		microProcessInfo[TIME_STEPS] = nsteps;

		ModelPart& model = *mpModelPart;
		mpScheme->InitializeSolutionStep(model, mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B());
		ModelPart::DofsArrayType& dofSet = mpLinearSOE->DofSet();
		mpConvergenceCriteria->InitializeSolutionStep(model, dofSet, mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B());

		mIntegrationErrorCode = 0.0;
		mSolutionStepFinalized = false;

#ifdef RVE_V2_INCREMENTAL_REGULARIZATION
		if(mRveGenerated)
		{
			CalculateMicroCharacteristicLength();
			AssignCharacteristicLengthMultiplier();
		}
#endif // RVE_V2_INCREMENTAL_REGULARIZATION

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
#ifdef RVE_V2_IS_VERBOSE
		std::cout << "RVE FINALIZE SOLUTION STEP\n"; 
#endif // RVE_V2_IS_VERBOSE
		
		mpConstraintHandler->FinalizeSolutionStep(*mpModelPart, *mpGeometryDescriptor, *mpMacroscaleData);

		mFirstStepDone = true;
		mStrainVectorOld = mStrainVectorOld_trial;

		mpScheme->FinalizeSolutionStep(*mpModelPart, mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B());
		mpScheme->Clean();
		mSolutionStepFinalized = true;

#ifdef RVE_V2_IS_VERBOSE
		std::cout << "RVE FINALIZE SOLUTION STEP - END\n"; 
#endif // RVE_V2_IS_VERBOSE
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
			if( mpModelPart == NULL || 
				mpMacroscaleData == NULL || 
				mpGeometryDescriptor == NULL ||
				mpConstraintHandler == NULL ||
				mpLinearSOE == NULL ||
				mpHomogenizer == NULL ||
				mpScheme == NULL ||
				mpConvergenceCriteria == NULL)
			{
				KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - RveGenerated flag set to TRUE, but some data is missing", "");
			}
		}
		return 0;
		KRATOS_CATCH("")
	}
	
public:

	/**
	* Gets a value indicating whether this RveManager has the necessary data (model part, rve macroscale status).
	* @return true, if data for this RveAdapterV2 has been generared, false otherwise
	*/
	inline const bool RveGenerated()const { return mRveGenerated; }

	/**
	* Gets a value indicating whether this RveManager needs the generation of new data (model part, rve macroscale status).
	* @return true, if this RveAdapterV2 needs new data, false otherwise.
	*/
	inline const bool RveGenerationRequested()const { return mRveGenerationRequested; }

	bool TestMaterialResponse(const Vector& E, bool compute_constitutive_tensor, const SizeType& load_Id)
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
		double detF = 1.0;
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
		params.SetDeterminantF(detF);
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
		RveUtilities::RestoreSolutionVector(*mpModelPart, mU);
		mpLinearSOE->BuildRHS(*mpModelPart, mpScheme);
		mpScheme->FinalizeSolutionStep(*mpModelPart, mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B());
		mpScheme->Clean();
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
	* Initializes the RVE model part and all the components.
	* Called from InitializeMaterial and SetRveData
	*/
	void InitializeRve()
	{
		ModelPart& model = *mpModelPart;

		// bifurcation analysis
		m_localized = false;
		m_localization_angle = 0.0;

		// random imperfections

		/*ModelPart::NodeType& n00 = model.GetNode( mpGeometryDescriptor->CornerNodesIDs()[0] );
		ModelPart::NodeType& n10 = model.GetNode( mpGeometryDescriptor->CornerNodesIDs()[1] );
		ModelPart::NodeType& n01 = model.GetNode( mpGeometryDescriptor->CornerNodesIDs()[3] );
		ModelPart::NodeType& n11 = model.GetNode( mpGeometryDescriptor->CornerNodesIDs()[2] );
		double xc=(n00.X0() + n10.X0() + n11.X0() + n01.X0())/4.0;
		double yc=(n00.Y0() + n10.Y0() + n11.Y0() + n01.Y0())/4.0;
		double lx=n10.X0()-n00.X0();
		double ly=n01.Y0()-n00.Y0();
		double radius = std::min(lx,ly)/2.0;
		for(ModelPart::NodeIterator node_it = mpModelPart->NodesBegin(); node_it != mpModelPart->NodesEnd(); ++node_it) {
			ModelPart::NodeType& inode = *node_it;
			double dx = inode.X0()-xc;
			double dy = inode.Y0()-yc;
			double L  = std::sqrt(dx*dx+dy*dy);
			double AA = 1.0;
			double BB = 0.0;
			double CC = 10.0*radius;
			double fpz_scale_factor = AA*std::exp(-(L-BB)*(L-BB)/(2.0*CC));
			fpz_scale_factor=1.0;
			double noiseval = scaled_raw_noise_2d(-0.1,0.1,inode.X0(),inode.Y0());
			inode.FastGetSolutionStepValue(RANDOM_IMPERFECTION_FACTOR) = fpz_scale_factor*noiseval;
		}*/

		/*for(ModelPart::NodeIterator node_it = model.NodesBegin(); node_it != model.NodesEnd(); ++node_it)
		{
			ModelPart::NodeType& inode = *node_it;
			double noiseval = scaled_raw_noise_2d(-0.1,0.1,inode.X0(),inode.Y0());
			inode.FastGetSolutionStepValue(RANDOM_IMPERFECTION_FACTOR) = noiseval;
		}*/

		if(m_macro_imperfection_factor != 0.0)
		{
			for(ModelPart::NodeIterator node_it = model.NodesBegin(); node_it != model.NodesEnd(); ++node_it)
			{
				ModelPart::NodeType& inode = *node_it;
				inode.FastGetSolutionStepValue(RANDOM_IMPERFECTION_FACTOR) = m_macro_imperfection_factor;
			}
		}

		mpConstraintHandler->AddConditions(model, *mpGeometryDescriptor);

		if(!mpScheme->SchemeIsInitialized())
			mpScheme->Initialize(model);

		if(!mpScheme->ElementsAreInitialized())
			mpScheme->InitializeElements(model);

		if(!mpScheme->ConditionsAreInitialized())
			mpScheme->InitializeConditions(model);

		if(!mpConvergenceCriteria->mConvergenceCriteriaIsInitialized)
			mpConvergenceCriteria->Initialize(model);
	}

	/**
	* Performs: U(n+1,trial) = U(n)
	* Called from RveCalculateMaterialResponse, before the call to Equilibrate
	*/
	void RevertToLastStep()
	{
		for (ModelPart::NodeIterator node_iter = mpModelPart->NodesBegin(); node_iter != mpModelPart->NodesEnd(); ++node_iter)
		{
			ModelPart::NodeType& node = *node_iter;

			if(mMoveMesh)
				noalias(node.GetInitialPosition()) = node.GetInitialPosition() + node.FastGetSolutionStepValue(DISPLACEMENT,1);

			for(ModelPart::NodeType::DofsContainerType::iterator dof_iter = node.GetDofs().begin(); dof_iter != node.GetDofs().end(); ++dof_iter)
			{
				ModelPart::DofType& iDof = *dof_iter;
				iDof.GetSolutionStepValue() = iDof.GetSolutionStepValue(1);
			}
		}
	}

	/**
	* Updates the nodal coordinates if necessary
	*/
	void MoveMesh()
	{
		for (ModelPart::NodeIterator i = mpModelPart->NodesBegin(); i != mpModelPart->NodesEnd(); ++i)
		{
			ModelPart::NodeType& node = *i;
			noalias(node.GetInitialPosition()) = node.GetInitialPosition() + node.FastGetSolutionStepValue(DISPLACEMENT);
		}
	}

	/**
	* Calculates the Response of this RveAdapterV2.
	* @param rValues material parameters
	*/
	virtual void CalculateRveResponse(ConstitutiveLaw::Parameters& rValues)
	{
#ifdef RVE_TIMER_ON
		RveUtilities::RveTimer timer_total;
		RveUtilities::RveTimer timer;
		timer_total.start();  
#endif // RVE_TIMER_ON

		Vector& stress_vector = rValues.GetStressVector();
		Matrix& tangent_matrix = rValues.GetConstitutiveMatrix();
		Flags& options = rValues.GetOptions();
		bool compute_constitutive_tensor = options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		mIntegrationErrorCode = 0.0;

		if(this->CheckForOutputStage(rValues)) {
			mStrainVectorOld_trial = rValues.GetStrainVector();
			return;
		}

		/*this->RevertToLastStep();*/
		/*
		 * optimization: reverting to the last converged state (of the macro-scale analysis)
		 * keeps the number of iterations at the RVE level fixed (more or less).
		 * without reverting to the last converged state, we start the equilibrium-iteration process
		 * with the last displacement as a first guess.
		 * this progressivley diminishes the numer of iterations at the RVE level as the macro-scale
		 * approaches a converged state.
		 * Just for the first 2 iterations it seems better to start from the last converged state
		 * since the prediction of the first iterations might be misleading...
		 */
		if(rValues.GetProcessInfo()[NL_ITERATION_NUMBER] < 3)
			this->RevertToLastStep();

		if(this->PredictorCalculation(rValues)) {
			mStrainVectorOld_trial = rValues.GetStrainVector();
			return;
		}

		ModelPart& model = *mpModelPart;
		RveGeometryDescriptor& geomdescriptor = *mpGeometryDescriptor;
		RveMacroscaleData& macrodata = *mpMacroscaleData;

#ifdef RVE_TIMER_ON
		timer.start();
#endif // RVE_TIMER_ON

		mpMacroscaleData->SetData(rValues,model,geomdescriptor);
		// MAZ_01+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// for perturbation ...
		if(mpMacroscaleData->StrainVectorOld().size() != mStrainVectorOld.size())
			mpMacroscaleData->StrainVectorOld().resize(mStrainVectorOld.size(), false);
		noalias(mpMacroscaleData->StrainVectorOld()) = mStrainVectorOld;
		// MAZ_01+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		mpConstraintHandler->ApplyMacroScaleData(model, geomdescriptor, macrodata);

		mpLinearSOE->Begin(model, geomdescriptor, mpScheme, mpConstraintHandler);

#ifdef RVE_TIMER_ON
		timer.stop();
		double t_01 = timer.value();
		timer.start();  
#endif // RVE_TIMER_ON

		bool converged = false;
		converged = this->Equilibrate();

#ifdef RVE_TIMER_ON
		timer.stop();
		double t_02 = timer.value();
		int nln = mpModelPart->GetProcessInfo()[NL_ITERATION_NUMBER];
#endif // RVE_TIMER_ON

		if(!converged)
		{
			SizeType strain_size = GetStrainSize();
			if(stress_vector.size() != strain_size)
				stress_vector.resize(strain_size, false);
			if(tangent_matrix.size1() != strain_size || tangent_matrix.size2() != strain_size)
				tangent_matrix.resize(strain_size, strain_size, false);
			noalias(stress_vector) = ZeroVector(strain_size);
			noalias(tangent_matrix) = IdentityMatrix(strain_size, strain_size);
			mIntegrationErrorCode = -1.0;
			mStrainVectorOld_trial = rValues.GetStrainVector();
			return;
		}
		
#ifdef RVE_TIMER_ON
		timer.start();  
#endif // RVE_TIMER_ON

		mpHomogenizer->HomogenizeStressTensor(
			model, 
			geomdescriptor, mpLinearSOE, 
			mpConstraintHandler, macrodata, stress_vector);
		noalias(mStressVector) = stress_vector;

#ifdef RVE_TIMER_ON
		timer.stop();
		double t_03 = timer.value();
		timer.start();  
#endif // RVE_TIMER_ON

		if(compute_constitutive_tensor)
		{
			RveUtilities::SaveSolutionVector(model, mU);
			mpHomogenizer->HomogenizeTengentConstitutiveTensor(
				model, geomdescriptor, 
				mpLinearSOE, 
				mpConstraintHandler,
				macrodata, mpScheme,
				stress_vector, tangent_matrix,
				mU,
				mMoveMesh);
			noalias(mConstitutiveTensor) = tangent_matrix;
			RveUtilities::RestoreSolutionVector(model, mU);
			mpLinearSOE->BuildRHS(model, mpScheme);
		}

#ifdef RVE_TIMER_ON
		timer.stop();
		double t_04 = timer.value();  
#endif // RVE_TIMER_ON

		mpLinearSOE->End();

		mStrainVectorOld_trial = rValues.GetStrainVector();

#ifdef RVE_TIMER_ON
		timer_total.stop();
		double t_total = timer_total.value();

		std::stringstream ss;
		ss << "RVE Timings:\n";
		ss << "Begin:       " << t_01 << "; % = " << t_01/t_total << std::endl;
		ss << "Equlibrate:  " << t_02 << "; % = " << t_02/t_total << " - with NLIN: " << nln << std::endl;
		ss << "Hom.Stress:  " << t_03 << "; % = " << t_03/t_total << std::endl;
		ss << "Hom.Tangent: " << t_04 << "; % = " << t_04/t_total << std::endl;
		std::cout << ss.str();  
#endif // RVE_TIMER_ON

#ifdef RVE_V2_CHECK_HOMG_STRAIN
		size_t hs_size = this->GetStrainSize() == 3 ? 2 : 3;
		Matrix homog_strain_tensor(hs_size, hs_size, 0.0);
		mpHomogenizer->HomogenizeVariable(*mpModelPart, *mpGeometryDescriptor, GREEN_LAGRANGE_STRAIN_TENSOR, homog_strain_tensor);
		Vector homog_strain_vector = MathUtils<double>::StrainTensorToVector(homog_strain_tensor);
		std::cout << MathHelpers::VectorToString(homog_strain_vector,4,std::scientific);
		double skew_grad_u = 0.0;
		mpHomogenizer->HomogenizeVariable(*mpModelPart, *mpGeometryDescriptor, TEMPERATURE, skew_grad_u);
		std::cout << "Skwe Grad U = " << skew_grad_u << std::endl;
#endif // RVE_V2_CHECK_HOMG_STRAIN

	}

	/**
	* Returns true if this RVE calculation has been called only for
	* output purposes, False otherwise
	*/
	virtual bool CheckForOutputStage(ConstitutiveLaw::Parameters& rValues)
	{
		Vector strain_increment = rValues.GetStrainVector() - mStrainVectorOld_trial;
		double min_tol = 1.0E-10;
		double rel_tol = 1.0E-6;
		bool zero_increment = norm_2(strain_increment) < std::max(min_tol, rel_tol*norm_2(rValues.GetStrainVector()));
		
		if(mSolutionStepFinalized || zero_increment) 
		{
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

#ifdef RVE_V2_IS_VERBOSE
		if(mSolutionStepFinalized)
			std::cout << "Check for output request: TRUE\n";
		else
			std::cout << "Check for output request: FALSE\n";  
#endif // RVE_V2_IS_VERBOSE

		return mSolutionStepFinalized;
	}

	/**
	* Performs (if possible) a predictor step without the solution of a micro-model.
	* This test checks whether this RveAdapterV2 can continue the calculationa without a micro-model.
	* The basic implementation of this method always return false.
	* In order to take advantage of this optimization functionality this method should be
	* re-implemented by a derived class that will perform a proper predictor check.
	* @param rValues material parameters
	* @return true, if this RveAdapterV2 succesfully performed the calculations without the micro-model, false otherwise.
	*/
	virtual bool PredictorCalculation(ConstitutiveLaw::Parameters& rValues)
	{
		return false;
		//if(!mRveGenerated)
		//{
		//	SizeType strain_size = GetStrainSize();

		//	Vector& StressVector = rValues.GetStressVector();
		//	if(StressVector.size() != strain_size)
		//		StressVector.resize(strain_size, false);
		//	noalias(StressVector) = ZeroVector(strain_size);

		//	Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
		//	if(ConstitutiveMatrix.size1() != strain_size || ConstitutiveMatrix.size2() != strain_size)
		//		ConstitutiveMatrix.resize(strain_size, strain_size, false);
		//	noalias(ConstitutiveMatrix) = ZeroMatrix(strain_size, strain_size);

		//	std::stringstream ss;
		//	ss << "RveAdapterV2: Cannot calculate the material response without a micromodel" << std::endl;
		//	std::cout << ss.str();

		//	// we need to solve the micro-problem
		//	// but we don't have a micro-model.
		//	// This should NEVER happen.
		//	// If this happens it means that there's a bug in the RveModeler
		//	return true;
		//}

		//// we need to solve the micro-problem
		//return false;
	}

	/**
	* Solves the equilibrium of the micro-structure.
	* @return True if convergence is achieved, False otherwise
	*/
	virtual bool Equilibrate()
	{
		if(mpLinearSOE->EquationSystemSize() == 0)
			return true;

		ModelPart& model = *mpModelPart;

		bool converged = false;
		unsigned int iteration_number = 0;
		model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

#ifdef RVE_V2_USES_KRYLOV_NEWTON
		// create the accelerator
		KrylovSubspaceAcceleratorUtilties::KrylovAcceleratorEigenlib krylov_acc;
		int new_numEqns = mpLinearSOE->EquationSystemSize();
		krylov_acc.BeginSolutionStep(new_numEqns);
#endif // RVE_V2_USES_KRYLOV_NEWTON

		while (true)
		{
			/* ----------------------------------------------------------
			* Perfomed operations:
			* - Initialize the nonlinear iteration
			*-----------------------------------------------------------*/

			iteration_number++;
			model.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

			mpScheme->InitializeNonLinIteration(model, mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B());

			/* ----------------------------------------------------------
			* Perfomed operations:
			* - Build the system and solve
			*-----------------------------------------------------------*/

#ifdef RVE_V2_USES_KRYLOV_NEWTON
			if(iteration_number == 1)
			{
				// build the LHS only once...
				mpLinearSOE->Build(model, mpScheme);
			}
			else
			{
				// ... or if the krylov method requires it
				if(krylov_acc.BuildSystemMatrix())
					mpLinearSOE->Build(model, mpScheme);
				else
					mpLinearSOE->BuildRHS(model, mpScheme);
			}
#else
			mpLinearSOE->Build(model, mpScheme);
#endif // RVE_V2_USES_KRYLOV_NEWTON

			mpLinearSOE->Solve();
			
#ifdef RVE_V2_USE_LINE_SEARCH
			PerformLineSearch();
#endif // RVE_V2_USE_LINE_SEARCH

			/* ----------------------------------------------------------
			* Perfomed operations:
			* - Abort iteration process in case of integration errors
			*-----------------------------------------------------------*/

			if(this->CheckIntegrationErrors() != 0.0)
			{
				converged = false;
				break;
			}

#ifdef RVE_V2_USES_KRYLOV_NEWTON
			if(krylov_acc.LeastSquares( mpLinearSOE->X() ) < 0)
			{
				converged = false;
				iteration_number = mMaxIterationNumber;
				break;
			}
			krylov_acc.AccelerateSolution( mpLinearSOE->X() );
#endif // RVE_V2_USES_KRYLOV_NEWTON

			/* ----------------------------------------------------------
			* Perfomed operations:
			* - Finalize the nonlinear iteration
			*-----------------------------------------------------------*/

			// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++
			//mpScheme->Update(model, mpLinearSOE->DofSet(), mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B());
			mpConstraintHandler->Update(model,*mpGeometryDescriptor,*mpMacroscaleData, 
				                        mpLinearSOE->TransformedEquationIds(), *mpScheme,
				                        mpLinearSOE->DofSet(), mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B(),
										mpLinearSOE->R(), mpLinearSOE->EquationSystemSize());
			// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++

			if(mMoveMesh) this->MoveMesh();

			mpScheme->FinalizeNonLinIteration(model, mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B());
			
			/* ----------------------------------------------------------
			* Perfomed operations:
			* - Convergence test
			*-----------------------------------------------------------*/

			converged = mpConvergenceCriteria->PostCriteria(model, mpLinearSOE->DofSet(), mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B());

			if(converged || (iteration_number > mMaxIterationNumber) ) break;
		}

		/*if(!converged)
		{
			std::cout << "WARNING: RVE MAX ITER REACHED!\n";
		}*/

		// MAZ_01 +++++++++++++++++++++++++++++++++++++++++++++++
		if(converged)
		{
			mpLinearSOE->BuildRHS(model, mpScheme); // build reduced
			mpConstraintHandler->PostUpdate(model,*mpGeometryDescriptor,*mpMacroscaleData, 
				                        mpLinearSOE->TransformedEquationIds(), *mpScheme,
				                        mpLinearSOE->DofSet(), mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B(),
										mpLinearSOE->R(), mpLinearSOE->EquationSystemSize());

		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++

		return converged;
	}

	virtual double PerformLineSearch()
	{
		ModelPart& model = *mpModelPart;

		VectorType ReferenceDx = mpLinearSOE->X();
		VectorType B0          = mpLinearSOE->B();

		double R0 = inner_prod(mpLinearSOE->X(),mpLinearSOE->B());
		mpConstraintHandler->Update(model,*mpGeometryDescriptor,*mpMacroscaleData, 
									mpLinearSOE->TransformedEquationIds(), *mpScheme,
									mpLinearSOE->DofSet(), mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B(),
									mpLinearSOE->R(), mpLinearSOE->EquationSystemSize());
		mpLinearSOE->BuildRHS(model, mpScheme);
		double R1= inner_prod(ReferenceDx, mpLinearSOE->B());
		mpLinearSOE->X() *= (-1.0);
		mpConstraintHandler->Update(model,*mpGeometryDescriptor,*mpMacroscaleData, 
									mpLinearSOE->TransformedEquationIds(), *mpScheme,
									mpLinearSOE->DofSet(), mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B(),
									mpLinearSOE->R(), mpLinearSOE->EquationSystemSize());

		double rCurrentAlpha = 1.0;
		if(R0*R1<0.0)
		{
			double R2 = R1;
			if(std::abs(R1) < std::abs(R0)) R2=R0;
			double R0start = R0;
			double alpha = 0.0;
			double nabla = 0.0;
			double delta = 1.0;
			double CurrentAlpha  = 1.0;
			int iterations=0;
			int max_iterations = 10;
			while(std::abs(R2/R0start)>0.3 && iterations<max_iterations && 
				 (R1*R0)<0.0 && std::abs(R1)>1.0e-7 && std::abs(R0)>1.0e-7)
			{
				alpha = 0.5*(nabla+delta);
				CurrentAlpha  = alpha;
				noalias(mpLinearSOE->X()) = ReferenceDx * CurrentAlpha;
				mpConstraintHandler->Update(model,*mpGeometryDescriptor,*mpMacroscaleData, 
									mpLinearSOE->TransformedEquationIds(), *mpScheme,
									mpLinearSOE->DofSet(), mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B(),
									mpLinearSOE->R(), mpLinearSOE->EquationSystemSize());
				mpLinearSOE->BuildRHS(model, mpScheme);
				R2 = inner_prod(ReferenceDx, mpLinearSOE->B());
				mpLinearSOE->X() *= (-1.0);
				mpConstraintHandler->Update(model,*mpGeometryDescriptor,*mpMacroscaleData, 
									mpLinearSOE->TransformedEquationIds(), *mpScheme,
									mpLinearSOE->DofSet(), mpLinearSOE->A(), mpLinearSOE->X(), mpLinearSOE->B(),
									mpLinearSOE->R(), mpLinearSOE->EquationSystemSize());
				if(R2*R1<0.0){
					nabla = alpha;
					R0 = R2;
				}
				else if(R2*R0<0.0){
					delta = alpha;
					R1 = R2;
				}
				else{
					break;
				}  
				iterations++;
			}

			rCurrentAlpha  = CurrentAlpha;
		}
		if(rCurrentAlpha>1.0 || rCurrentAlpha<=0.0) 
			rCurrentAlpha=1.0;
		
		noalias(mpLinearSOE->X()) = ReferenceDx * rCurrentAlpha;
		noalias(mpLinearSOE->B()) = B0;

		return rCurrentAlpha;
	}

	virtual double CheckIntegrationErrors()
	{
		ModelPart & model = *mpModelPart;
		std::vector<double> error_codes;
		for(ModelPart::ElementIterator elem_iter = model.ElementsBegin(); elem_iter != model.ElementsEnd(); elem_iter++)
		{
			Element& elem = *elem_iter;
			elem.GetValueOnIntegrationPoints(CONSTITUTIVE_INTEGRATION_ERROR_CODE, error_codes, model.GetProcessInfo());
			for(size_t i = 0; i < error_codes.size(); i++)
				if(error_codes[i] != 0.0)
					return error_codes[i];
		}
		return 0.0;
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

	void StrainVectorToTensor(const Vector& a, Matrix& b)
	{
		SizeType vsize = a.size();
		if (vsize==3)
        {
			if(b.size1() != 2 || b.size2() != 2) b.resize(2,2,false);
			b(0,0) = a[0];
			b(0,1) = a[2]*0.5;
			b(1,0) = a[2]*0.5;
			b(1,1) = a[1];
        }
		else if (vsize==4)
        {
			if(b.size1() != 3 || b.size2() != 3) b.resize(3,3,false);
			b(0,0) = a[0];
			b(0,1) = a[3]*0.5;
			b(0,2) = 0.0;
			b(1,0) = a[3]*0.5;
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
			b(0,1) = a[3]*0.5;
			b(0,2) = a[5]*0.5;
			b(1,0) = a[3]*0.5;
			b(1,1) = a[1];
			b(1,2) = a[4]*0.5;
			b(2,0) = a[5]*0.5;
			b(2,1) = a[4]*0.5;
			b(2,2) = a[2];
        }
	}

	// RVE Custom operations (coming soon...)

	void CalculateMicroCharacteristicLength()
	{
		// the FE lch at the micro_scale "lch_micro" should be modified as follows:
		// lch_micro_bar = lch_micro * lch_multiplier
		// where
		// lch_multiplier = lch_macro / lch_rve
		SizeType ndim = WorkingSpaceDimension();
		if(ndim == 2)
		{
			mMicroCharacteristicLength = std::sqrt(mpGeometryDescriptor->DomainSize());
		}
		else // if(ndim == 3)
		{
			mMicroCharacteristicLength = std::pow(mpGeometryDescriptor->DomainSize(), 1.0/3.0);
		}
	}

	void AssignCharacteristicLengthMultiplier()
	{
		double chlen_mult = mMacroCharacteristicLength / mMicroCharacteristicLength;
		/*chlen_mult = 1.0;*/
#ifdef RVE_V2_INCREMENTAL_REGULARIZATION
		if(!m_localized)
			chlen_mult=1.0;  
#endif // RVE_V2_INCREMENTAL_REGULARIZATION

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

	// RVE components

	ModelPart::Pointer                    mpModelPart;
	RveMacroscaleData::Pointer            mpMacroscaleData;
	RveGeometryDescriptor::Pointer        mpGeometryDescriptor;
	RveConstraintHandlerPointerType       mpConstraintHandler;
	RveLinearSystemOfEquationsPointerType mpLinearSOE;
	RveHomogenizerPointerType             mpHomogenizer;
	SchemePointerType                     mpScheme;
	ConvergenceCriteriaPointerType        mpConvergenceCriteria;

	// RVE Generator handler (sarà un componente appena possibile...)

	bool   mRveGenerated;
	bool   mRveGenerationRequested;

	// RVE Custom operations (coming soon...)

	double mMacroCharacteristicLength;
	double mMicroCharacteristicLength;

	// RVE management data

	bool   mInitialized;
	bool   mSolutionStepFinalized;
	double mIntegrationErrorCode;
	bool   mMoveMesh;
	size_t mMaxIterationNumber;

	// RVE stored data

	bool   mFirstStepDone;
	Vector mStrainVectorOld;
	Vector mStrainVectorOld_trial;
	Vector mStressVector;
	Matrix mConstitutiveTensor;
	Vector mU;
	size_t mNDofs;

	// Bifurcation analysis

	bool m_localized;
	double m_localization_angle;

	// random imperfections

	double m_macro_imperfection_factor;

private:

	/**
	* private copy constructor
	*/
	RveAdapterV2(const RveAdapterV2& other);
	
	/**
	* private assignment operator
	*/
	RveAdapterV2& operator = (const RveAdapterV2& other);

private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
		//rSerializer.save("ModelPart", mpModelPart);
		//rSerializer.save("Status", mpMacroscaleStatus);
		//rSerializer.save("RveGen", mRveGenerated);
		//rSerializer.save("RveGenReq", mRveGenerationRequested);
    }

    virtual void load(Serializer& rSerializer)
    {
		//rSerializer.load("ModelPart", mpModelPart);
		//rSerializer.load("Status", mpMacroscaleStatus);
		//rSerializer.load("RveGen", mRveGenerated);
		//rSerializer.load("RveGenReq", mRveGenerationRequested);
		//if(mpModelPart == NULL)
		//{
		//	mRveGenerationRequested = true;
		//	mRveGenerated = false;
		//}
		//else
		//{
		//	mpStrategy = TSolvingStrategyPointer(new TSolvingStrategy(*mpModelPart));
		//	mRveGenerationRequested = false;
		//	mRveGenerated = true;
		//}
    }
	
};

template<class TSparseSpace,
         class TDenseSpace,
		 class TRveSettings
         >
inline std::ostream & operator << (std::ostream& rOStream, const RveAdapterV2<TSparseSpace, TDenseSpace, TRveSettings> & rThis)
{
	return rOStream << rThis.GetInfo();
}

} // namespace Kratos

#endif // RVE_ADAPTER_V2_H_INCLUDED
