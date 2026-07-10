// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#include "includes/element.h"
#include "sph_application_variables.h"
#include "custom_utilities/custom_kernels/kernel_factory.h"
#include "custom_utilities/compute_kernel_correction_utilities.h"
#include "custom_utilities/sph_element_utilities.h"

namespace Kratos
{
    
using SizeType = std::size_t;

template<class TKernelType>
class KRATOS_API(SPH_APPLICATION) TotalLagrangianMixedvFParticle: public Element
{

protected:

public:

    using BaseType = Element;

    typedef GeometryData::IntegrationMethod IntegrationMethod; 

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFNITION(TotalLagrangianMixedvFParticle);

    // Constructor void 
    TotalLagrangianMixedvFParticle()
    {
    }

    // Constructor using an array of nodes 
    TotalLagrangianMixedvFParticle(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    // Constructor using an array of nodes with properties 
    TotalLagrangianMixedvFParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod =  GeometryData::IntegrationMethod::GI_GAUSS_1;  ///
    }

    // Copy constructor
    TotalLagrangianMixedvFParticle(TotalLagrangianMixedvFParticle const& rOther)
        : BaseType(rOther),
        mThisConstitutiveLaw(rOther.mThisConstitutiveLaw),
        mThisIntegrationMethod(rOther.mThisIntegrationMethod)  ////
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<TotalLagrangianMixedvFParticle>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<TotalLagrangianMixedvFParticle>(NewId, pGeom, pProperties);
    }

    /**
     * @brief It creates a new element pointer and clones the previous element data
     */
    Element::Pointer Clone( IndexType NewId, NodesArrayType const& rThisNodes) const override;

    /**
     * @brief Called to initialize the element
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    IntegrationMethod GetIntegrationMethod() const override
    {
        return mThisIntegrationMethod;   ////
    }

    /**
     * @brief This function tells the position of the particle in the list of neighbours
     */
    int GetNeighbourPosition(const std::vector<Element::Pointer>& rNeighbours) const
    {
        int i = 0; 

        while (i<rNeighbours.size() && this->Id() != rNeighbours[i]->Id()) i++;

        return i;
    }

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation IDs
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;
    
    /**
     *  @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;
    
    /**
     * @brief Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(VectorType& rValues, int step = 0) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

protected:

    ConstitutiveLaw::Pointer mThisConstitutiveLaw;
    IntegrationMethod mThisIntegrationMethod;

    /**
     * @brief This function sets the used constitutive laws
     */
    void SetConstitutiveLaw(const ConstitutiveLaw::Pointer rThisConstitutiveLaw)
    {
        mThisConstitutiveLaw = rThisConstitutiveLaw;
    }

    /**
     * @brief Sets the used integration method
     * @details In SPH the inetgration points are the particles themselves,
     * this method is implement to maintain compatibility and display the results on the integration points 
    */
    void SetIntegrationMethod(const IntegrationMethod ThisIntegrationMethod)
    {
        mThisIntegrationMethod = ThisIntegrationMethod; ////
    }

    /**
     * @brief It initializes the material 
     */
    void InitializeMaterial();

private:



};

}