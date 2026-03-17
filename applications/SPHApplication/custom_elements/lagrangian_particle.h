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
class KRATOS_API(SPH_APPLICATION) LagrangianParticle : public Element
{
public: 

    using BaseType = Element;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LagrangianParticle);

    // Constructor void 
    LagrangianParticle()
    {
    }

    // Constructor using an array of nodes 
    LagrangianParticle(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    // Constructor using an array of nodes with properties 
    LagrangianParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    { 
    }

    // Copy constructor
    LagrangianParticle(LagrangianParticle const& rOther)
        : BaseType(rOther),
        mThisConstitutiveLaw(rOther.mThisConstitutiveLaw)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LagrangianParticle>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LagrangianParticle>(NewId, pGeom, pProperties);
    }

    /**
     * @brief Indicates the amount of DoFs per node (v_x, v_y, J)
     */
    virtual IndexType GetDoFsPerNode() const 
    {
        return 2; // Understand if the particle stays 2D or ca be both 2D and 3D, at monent only 2D
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
     * @brief Returns a vector that includes the values of the DoFs
     */
    virtual void GetNodalValuesVector(VectorType& rNodalValue) const;

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     */
    void EquationIdVector(
        EquationIdVectorType& rElementalDofList,
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

    /**
     * @brief Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(VectorType& rValues, int step = 0) const override;

    /**
     * @brief
     */
    void GlobalSizeVector(
        VectorType& rLocalVector,
        VectorType& rNodalValue,
        IndexType i
    )
    {
        IndexType dofs_per_node = GetDoFsPerNode();
        rLocalVector.clear();
        rLocalVector[0] = rNodalValue[dofs_per_node * i];
        rLocalVector[1] = rNodalValue[dofs_per_node * i + 1];
    };
    
    /**
     * @brief
     */
    virtual void GetShapeFunctionDerivatives(
        MatrixType& rB,
        VectorType& rCGCK,
        const double volume
    );

    /**
     * @brief
     */
    virtual void CalculateStrainVector(
        VectorType& rStrainVector,
        VectorType& rCGCK,
        const VectorType& rNodalValues,
        const double volume
    );
    
    /**
     * @brief This function reads and computes the value of body forces from a neighbouring node
     * @param rElement neighbouring particle 
     */
    virtual array_1d<double, 2> GetLocalBodyForces() const; 

    /**
     * @brief This is called during the assembling process in order to calculate the local system
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS matrix
     * @param rRightHandSideVector The RHS vector
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        );

    /**
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix The elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental damping matrix
      * @param rDampingMatrix The elemental damping matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief These functions calculates the values of variables in the integrations points.
     * In SPH case coincide with the neighbouring particles 
     * @details These functions expect a std::vector of values for the specified variable type
     * @param rVariable This parameter selects the output 
     * @param SPH_KERNEL The function computes the kernel values in the neighbours  
     * @param SPH_KERNEL_GRADIENT The function computes the kernel gradient values in the neighbours  
     */

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

protected:

    ConstitutiveLaw::Pointer mThisConstitutiveLaw;

    /**
     * @brief This function sets the used constitutive laws
     */
    void SetConstitutiveLaw(const ConstitutiveLaw::Pointer rThisConstitutiveLaw)
    {
        mThisConstitutiveLaw = rThisConstitutiveLaw;
    }

    /**
     * @brief It initializes the material
     */
    void InitializeMaterial();

private:

};

}