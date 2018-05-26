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

#if !defined(KRATOS_RVE_CONSTITUTIVE_LAW )
#define  KRATOS_RVE_CONSTITUTIVE_LAW

/* System includes */

/* External includes */

/* Project includes */
#include "includes/constitutive_law.h"
#include "constitutive_law_adapter.h"

namespace Kratos
{

/**
 * Base class of constitutive laws.
 */
template<class TRveAdapter>
class RveConstitutiveLaw : public ConstitutiveLawAdapter<TRveAdapter>
{

public:

	KRATOS_CLASS_POINTER_DEFINITION(RveConstitutiveLaw);
	
	typedef ConstitutiveLawAdapter<TRveAdapter> MyBase;

    typedef typename MyBase::GeometryType GeometryType;

	typedef typename MyBase::SizeType SizeType;
	
public:

    /**
     * Constructor.
     */
    RveConstitutiveLaw(const typename MyBase::TAdapterPointer& theAdapter)
		: MyBase(theAdapter)
	{
	}
	
    /**
     * Destructor.
     */
    virtual ~RveConstitutiveLaw()
	{
	}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      ConstitutiveLaw::Pointer p_clone(new RveConstitutiveLaw());
     *      return p_clone;
     */
    virtual ConstitutiveLaw::Pointer Clone() const
	{
		KRATOS_THROW_ERROR(std::logic_error, "RveConstitutiveLaw cannot be cloned. It should be created by the modeler and assigned to Elements","");
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
    virtual void CalculateVolumetricResponse(const double VolumetricStrain,
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
		KRATOS_THROW_ERROR(std::logic_error, "RveConstitutiveLaw - Method not implemented","");
	}

    /**
     * Computes the deviatoric part of the material response in terms of stresses and algorithmic tangent
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


     * TODO: add proper definition for algorithmic tangent
     */
    virtual void CalculateDeviatoricResponse(const Vector& StrainVector,
                                             const Matrix& DeformationGradient,
                                             Vector& StressVector,
                                             Matrix& AlgorithmicTangent,
                                             const ProcessInfo& rCurrentProcessInfo,
                                             const Properties& rMaterialProperties,
                                             const GeometryType& rElementGeometry,
                                             const Vector& rShapeFunctionsValues,
                                             bool CalculateStresses = true,
                                             int CalculateTangent = true,
                                             bool SaveInternalVariables = true)
	{
		KRATOS_THROW_ERROR(std::logic_error, "RveConstitutiveLaw - Method not implemented","");
	}


    // VM

    virtual void CalculateCauchyStresses(Vector& Cauchy_StressVector,
                                         const Matrix& F,
                                         const Vector& PK2_StressVector,
                                         const Vector& GreenLagrangeStrainVector)
	{
		KRATOS_THROW_ERROR(std::logic_error, "RveConstitutiveLaw - Method not implemented","");
	}

	// New methods of this class

	inline const ModelPart* GetModelPart()const
	{
        return MyBase::mpAdaptee->GetModelPart();
	}

	bool TestMaterialResponse(const Vector& strains, bool compute_constitutive_tensor, const SizeType& load_id)
	{
		return MyBase::mpAdaptee->TestMaterialResponse(strains, compute_constitutive_tensor, load_id);
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

	RveConstitutiveLaw(){}

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MyBase );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MyBase );
    }

    ///@}

}; /* Class RveConstitutiveLaw */

} /* namespace Kratos.*/
#endif /* KRATOS_RVE_CONSTITUTIVE_LAW  defined */
