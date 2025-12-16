// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/shellq4_coordinate_transformation.hpp"
#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class MITC4AndesShellThickElement3D4N
 * @ingroup StructuralMechanicsApplication
 * @brief This shell element combines the MITC4 formulation for bending and shear (Reissner-Mindlin theory) and the ANDES membrane formulation
 * @details The element is defined by 4 nodes in 3D space with 6 degrees of freedom per node (3 translations and 3 rotations).
 * Reference papers:
 * MITC4: "Short communication: A four-node plate bending element based on Mindlin/Reissner plate theory and a mixed interpolation", Bathe and Dvorkin, 1985. Int. Journal for Numerical Methods in Engineering, 21:367-383.
 * ANDES: "Buckling and stability problems for thin shell structures using high performance finite elements", Haugen PhD Thesis, 1991.
 * @author Alejandro Cornejo
 */
template <bool IS_COROTATIONAL>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) MITC4AndesShellThickElement3D4N
    : public Element
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = Element;
    using array_3 = array_1d<double, 3>;
    using bounded_24_vector = array_1d<double, 24>;
    using bounded_3_matrix = BoundedMatrix<double, 3, 3>; // rotation matrix
    using bounded_24_matrix = BoundedMatrix<double, 24, 24>; // stiffness matrix

    static constexpr bool is_corotational = IS_COROTATIONAL;
    using CoordinateTransformationType = std::conditional_t<
        IS_COROTATIONAL,
        ShellQ4_CorotationalCoordinateTransformation,
        ShellQ4_CoordinateTransformation>;

    using CoordinateTransformationPointerType = Kratos::unique_ptr<CoordinateTransformationType>;

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MITC4AndesShellThickElement3D4N);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    MITC4AndesShellThickElement3D4N()
    {
    }

    // Constructor using an array of nodes
    MITC4AndesShellThickElement3D4N(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry),
    mpCoordinateTransformation(Kratos::make_unique<CoordinateTransformationType>(pGeometry))
    {
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
    }

    // Constructor using an array of nodes with properties
    MITC4AndesShellThickElement3D4N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId,pGeometry,pProperties),
        mpCoordinateTransformation(Kratos::make_unique<CoordinateTransformationType>(pGeometry))

    {
        // This is needed to prevent uninitialised integration method in inactive elements
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
    }

    // Copy constructor
    MITC4AndesShellThickElement3D4N(MITC4AndesShellThickElement3D4N const& rOther)
        : BaseType(rOther),
        mThisIntegrationMethod(rOther.mThisIntegrationMethod),
        mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<MITC4AndesShellThickElement3D4N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override
    {
        return Kratos::make_intrusive<MITC4AndesShellThickElement3D4N>(NewId, pGeom, pProperties);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Indicates the amount of DoFs per node (u, v, w, theta_x, theta_y, theta_z)
     */
    IndexType GetDoFsPerNode() const
    {
        return 6;
    }

    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This method initializes the constitutive law vector and the individual constitutive laws too
     * @warning Must be called before any calculation is done
     */
    void InitializeMaterial();

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList The vector containing the dof of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Returns the used integration method
     * @return default integration method of the used Geometry
     */
    IntegrationMethod GetIntegrationMethod() const override
    {
        return mThisIntegrationMethod;
    }

    /**
    * element can be integrated using the GP provided by the geometry or custom ones
    * by default, the base element will use the standard integration provided by the geom
    * @return bool to select if use/not use GPs given by the geometry
    */
    bool UseGeometryIntegrationMethod() const
    {
        return true;
    }

    /**
     * @brief Returns the set of integration points
     */
    const GeometryType::IntegrationPointsArrayType IntegrationPoints() const 
    {
        return GetGeometry().IntegrationPoints();
    }

    /**
     * @brief Returns the set of integration points
     */
    const GeometryType::IntegrationPointsArrayType IntegrationPoints(IntegrationMethod ThisMethod) const
    {
        return GetGeometry().IntegrationPoints(ThisMethod);
    }

    /**
     * @brief This method computes the Strain-Displacement matrix B, used to relate nodal displacements to strains for a quadrilateral
     * It assumes that the coordinates are already in the local coordinate system
     * @details The B matrix includes the bending and shear parts. Size of 8x24 since we have 8 generalized strains and 24 dofs (4 nodes with 6 dofs each)
     */
    void CalculateShearBendingB(
        MatrixType& rB,
        const double Area,
        const array_3& r_coord_1, 
        const array_3& r_coord_2, 
        const array_3& r_coord_3,
        const array_3& r_coord_4
    );

    /**
     * @brief This method computes the Strain-Displacement matrix B for the membrane part.
     * It assumes that the coordinates are already in the local coordinate system
     * @details The B matrix includes the membrane based on the Haugen and Felippa's ANDES quad membrane element
     */
    void CalculateMembraneB(
        MatrixType& rB,
        const double Area,
        const array_3& r_coord_1, 
        const array_3& r_coord_2, 
        const array_3& r_coord_3,
        const array_3& r_coord_4
    );

    /**
     * @brief This method computes the are of the triangle defined by the three given coordinates
     */
    double CalculateArea(
        const array_3 &r_coord_1,
        const array_3 &r_coord_2,
        const array_3 &r_coord_3,
        const array_3 &r_coord_4
    ) const;

    /**
     * @brief This method computes the rotation matrix from global to local coordinates
     * The ortonormal basis vectors are stored in the rows of the rotation matrix
     * T = [e1
     *      e2
     *      e3]
     */
    void CalculateRotationMatrixGlobalToLocal(
        bounded_3_matrix& rRotationMatrix,
        const bool UseInitialConfiguration
    ) const;

    /**
     * @brief This method computes rotates the LHS from local to global coordinates
     */
    void RotateLHSToGlobal(
        MatrixType& rLeftHandSideMatrix,
        const bounded_3_matrix& rRotationMatrix
    ) const;

    /**
     * @brief This method computes rotates the LHS from local to global coordinates
     */
    void RotateRHSToGlobal(
        VectorType& rRHS,
        const bounded_3_matrix& rRotationMatrix
    ) const;

    /**
     * @brief This method computes rotates the LHS from local to global coordinates
     */
    void RotateRHSToLocal(
        VectorType& rRHS,
        const bounded_3_matrix& rRotationMatrix
    ) const;

    /**
     * @brief This method builds the nodal values vector in the order: [u1, v1, w1, theta_x1, theta_y1, theta_z1, ...]
     * It returns in local axes in linear case, and only the deformational movements in corotational case
     */
    void GetNodalValuesVector(VectorType &rNodalValues, const bounded_3_matrix &rT) const;

    /**
     * @brief This function provides a more general interface to the element.
     * @details It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
     * @param rRightHandSideVector container for the desired RHS output
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector the elemental right hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo the current process info instance
     */
    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

    /**
     * @brief Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief this is called for non-linear analysis at the end of the iteration process
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;


    /**
    * @brief This function returns the size of the strain vector
    **/
    IndexType GetStrainSize() const
    {
        // We can call this method after we perform the Initialize
        return mConstitutiveLawVector[0]->GetStrainSize();
    }

    /**
     * @brief Returns a custom 3-point Gauss quadrature in area coordinates for a triangle.
     */
    // GeometryType::IntegrationPointsArrayType CustomTriangleAreaCoordinatesQuadrature(const double Area) const
    // {
    //     GeometryType::IntegrationPointsArrayType integration_points(3);

    //     const double w = Area / 3.0;

    //     integration_points[0] = IntegrationPoint<3>(0.5, 0.5, 0.0, w);
    //     integration_points[1] = IntegrationPoint<3>(0.0, 0.5, 0.5, w);
    //     integration_points[2] = IntegrationPoint<3>(0.5, 0.0, 0.5, w);

    //     return integration_points;
    // }

    /**
     * @brief This method add the body forces contribution to the RHS
     */
    void AddBodyForces(const double Area, VectorType &rRightHandSideVector);

    /**
     * @brief This method computes the mass matrix
     */
    void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief This method computes the damping matrix
     */
    void CalculateDampingMatrix(
        MatrixType &rDampingMatrix,
        const ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief This method returns a material property (e.g. Poisson ratio) without assuming that this property is
     * in the main property. It looks into the subproperties to find the property.
     */
    template <class TDataType> 
    TDataType GetMaterialProperty(const Variable<TDataType>& rVariable, const Properties& rProps)
    {
        if (rProps.Has(rVariable)) {
            return rProps.GetValue(rVariable);
        } else {
            const IndexType number_subprops = rProps.NumberOfSubproperties();
            const auto &r_sub_props_list = rProps.GetSubProperties();
            for (auto& r_subprop : r_sub_props_list) {
                if (r_subprop.Has(rVariable)) {
                    return r_subprop.GetValue(rVariable);
                }
            }
        }
        KRATOS_WARNING("MITC4AndesShellThickElement3D4N") << "The variable requested is not present in ANY subproperty..." << std::endl;
        return 0.0;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CS-DSG3 3N triangle shell Element #" << Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    IntegrationMethod mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2; /// Currently selected integration methods

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws

    Kratos::unique_ptr<CoordinateTransformationType> mpCoordinateTransformation = nullptr;


    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Sets the used integration method
     * @param ThisIntegrationMethod Integration method used
     */
    void SetIntegrationMethod(const IntegrationMethod& rThisIntegrationMethod)
    {
        mThisIntegrationMethod = rThisIntegrationMethod;
    }

    /**
     * @brief Sets the used constitutive laws
     * @param ThisConstitutiveLawVector Constitutive laws used
     */
    void SetConstitutiveLawVector(const std::vector<ConstitutiveLaw::Pointer>& rThisConstitutiveLawVector)
    {
        mConstitutiveLawVector = rThisConstitutiveLawVector;
    }


    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer &rSerializer) const override;

    void load(Serializer &rSerializer) override;

}; // class MITC4AndesShellThickElement3D4N.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
