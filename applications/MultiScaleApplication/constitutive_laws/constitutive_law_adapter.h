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

#if !defined(KRATOS_CONSTITUTIVE_LAW_ADAPTER )
#define  KRATOS_CONSTITUTIVE_LAW_ADAPTER

/* System includes */

/* External includes */

/* Project includes */
#include "includes/constitutive_law.h"

namespace Kratos
{

/**
 * Base class of constitutive law adapters.
 */
template<class TAdapter>
class ConstitutiveLawAdapter : public ConstitutiveLaw
{

public:

	KRATOS_CLASS_POINTER_DEFINITION(ConstitutiveLawAdapter);
	
	typedef typename TAdapter::Pointer TAdapterPointer;
	
	typedef ConstitutiveLaw::Parameters Parameters;
	
	typedef ConstitutiveLaw::SizeType SizeType;
	
	typedef ConstitutiveLaw::GeometryType GeometryType;
	
public:

    /**
     * Constructor.
     */
	ConstitutiveLawAdapter(const TAdapterPointer& theAdaptee)
		: ConstitutiveLaw()
		, mpAdaptee(theAdaptee)
	{
	}
	
    /**
     * Destructor.
     */
    virtual ~ConstitutiveLawAdapter()
	{
	}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      ConstitutiveLaw::Pointer p_clone(new ConstitutiveLawAdapter());
     *      return p_clone;
     */
    virtual ConstitutiveLaw::Pointer Clone() const
	{
		KRATOS_THROW_ERROR(std::logic_error, "ConstitutiveLawAdapter cannot be cloned. This method should be implemented by any derived class","");
	}

    /**
     * @return the working space dimension of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType WorkingSpaceDimension()
    {
        return mpAdaptee->WorkingSpaceDimension();
    }

    /**
     * returns the size of the strain vector of the current constitutive law
     * NOTE: this function HAS TO BE IMPLEMENTED by any derived class
     */
    virtual SizeType GetStrainSize()
    {
        return mpAdaptee->GetStrainSize();
    }

	/**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    virtual bool Has(const Variable<int>& rThisVariable)
	{
		return mpAdaptee->Has(rThisVariable);
	}

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    virtual bool Has(const Variable<double>& rThisVariable)
	{
		return mpAdaptee->Has(rThisVariable);
	}

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    virtual bool Has(const Variable<Vector>& rThisVariable)
	{
		return mpAdaptee->Has(rThisVariable);
	}

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    virtual bool Has(const Variable<Matrix>& rThisVariable)
	{
		return mpAdaptee->Has(rThisVariable);
	}

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * NOTE: fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
     */
    virtual bool Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return mpAdaptee->Has(rThisVariable);
	}

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * NOTE: fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
     */
    virtual bool Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return mpAdaptee->Has(rThisVariable);
	}

	/**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual int& GetValue(const Variable<int>& rThisVariable, int& rValue)
	{
		return mpAdaptee->GetValue(rThisVariable, rValue);
	}

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		return mpAdaptee->GetValue(rThisVariable, rValue);
	}

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
    virtual Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
	{
		return mpAdaptee->GetValue(rThisVariable, rValue);
	}

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @return the value of the specified variable
     */
    virtual Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		return mpAdaptee->GetValue(rThisVariable, rValue);
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
		return mpAdaptee->GetValue(rVariable, rValue);
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
		return mpAdaptee->GetValue(rVariable, rValue);
	}

	/**
     * sets the value of a specified variable
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(const Variable<int>& rVariable,
                          const int& rValue,
                          const ProcessInfo& rCurrentProcessInfo)
	{
		mpAdaptee->SetValue(rVariable, rValue, rCurrentProcessInfo);
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
		mpAdaptee->SetValue(rVariable, rValue, rCurrentProcessInfo);
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
		mpAdaptee->SetValue(rVariable, rValue, rCurrentProcessInfo);
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
		mpAdaptee->SetValue(rVariable, rValue, rCurrentProcessInfo);
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
		mpAdaptee->SetValue(rVariable, rValue, rCurrentProcessInfo);
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
		mpAdaptee->SetValue(rVariable, rValue, rCurrentProcessInfo);
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
		return mpAdaptee->ValidateInput(rMaterialProperties);
	}

    /**
     * returns the expected strain measure of this constitutive law (by default linear strains)
     * @return the expected strain measure
     */
    virtual StrainMeasure GetStrainMeasure()
	{
		return mpAdaptee->GetStrainMeasure();
	}

    /**
     * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    virtual StressMeasure GetStressMeasure()
	{
		return mpAdaptee->GetStressMeasure();
	}

    /**
     * returns whether this constitutive model is formulated in incremental strains/stresses
     * NOTE: by default, all constitutive models should be formulated in total strains
     * @return true, if formulated in incremental strains/stresses, false otherwise
     */
    virtual bool IsIncremental()
	{
		return mpAdaptee->IsIncremental();
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
		mpAdaptee->InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
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
		mpAdaptee->InitializeSolutionStep(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
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
		mpAdaptee->FinalizeSolutionStep(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
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
		mpAdaptee->InitializeNonLinearIteration(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
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
		mpAdaptee->FinalizeNonLinearIteration(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
	}
    
    /**
     * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
	{
		mpAdaptee->CalculateMaterialResponsePK1(rValues);
	}
    
    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues)
	{
		mpAdaptee->CalculateMaterialResponsePK2(rValues);
	}
    
    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
	{
		mpAdaptee->CalculateMaterialResponseKirchhoff(rValues);
	}
    
    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
	{
		mpAdaptee->CalculateMaterialResponseCauchy(rValues);
	}
    
    /**
     * Updates the material response in terms of 1st Piola-Kirchhoff stresses
     * @see ConstitutiveLaw::Parameters
     */
    virtual void FinalizeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
	{
		mpAdaptee->FinalizeMaterialResponsePK1(rValues);
	}
    
    /**
     * Updates the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see ConstitutiveLaw::Parameters
     */
    virtual void FinalizeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues)
	{
		mpAdaptee->FinalizeMaterialResponsePK2(rValues);
	}

    /**
     * Updates the material response in terms of Kirchhoff stresses
     * @see ConstitutiveLaw::Parameters
     */
    virtual void FinalizeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
	{
		mpAdaptee->FinalizeMaterialResponseKirchhoff(rValues);
	}

    /**
     * Updates the material response in terms of Cauchy stresses
     * @see ConstitutiveLaw::Parameters
     */
    virtual void FinalizeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
	{
		mpAdaptee->FinalizeMaterialResponseCauchy(rValues);
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
		mpAdaptee->ResetMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
	}

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    virtual void GetLawFeatures(Features& rFeatures)
	{
		mpAdaptee->GetLawFeatures(rFeatures);
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
		if(mpAdaptee == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "ConstitutiveLawAdapter - missing the Adaptee", "");
		int retval = mpAdaptee->Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
		return retval;
		KRATOS_CATCH("")
	}

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    
    ///@name Protected member Variables
    ///@{
	
	TAdapterPointer mpAdaptee;
	
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

	ConstitutiveLawAdapter(){}

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw );
		rSerializer.save("Adaptee", mpAdaptee);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw );
		rSerializer.load("Adaptee", mpAdaptee);
    }

    ///@}

}; /* Class ConstitutiveLawAdapter */

} /* namespace Kratos.*/
#endif /* KRATOS_CONSTITUTIVE_LAW_ADAPTER  defined */
