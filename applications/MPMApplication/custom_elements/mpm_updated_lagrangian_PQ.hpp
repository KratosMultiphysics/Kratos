//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//  Pull request:    https://github.com/KratosMultiphysics/Kratos/pull/6927



#if !defined(KRATOS_UPDATED_LAGRANGIAN_PQ_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_PQ_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/mpm_updated_lagrangian.hpp"

namespace Kratos
{
/// Partitioned Quadrature Large Displacement Lagrangian Element for 3D and 2D geometries. (base class)

/**
 * Implements a Partitioned Quadrature Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */

class MPMUpdatedLagrangianPQ
    : public MPMUpdatedLagrangian
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    /// Counted pointer of LargeDisplacementElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMUpdatedLagrangianPQ );
    ///@}

public:

    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    MPMUpdatedLagrangianPQ();


    /// Default constructors
    MPMUpdatedLagrangianPQ(IndexType NewId, GeometryType::Pointer pGeometry);

    MPMUpdatedLagrangianPQ(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    MPMUpdatedLagrangianPQ(MPMUpdatedLagrangianPQ const& rOther);

    /// Destructor.
    ~MPMUpdatedLagrangianPQ() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MPMUpdatedLagrangianPQ& operator=(MPMUpdatedLagrangianPQ const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    //************* STARTING - ENDING  METHODS

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;


    ///@}
    ///@name Access
    ///@{

    void CalculateOnIntegrationPoints(const Variable<int>& rVariable,
        std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

protected:

     /// Calculation of the External Forces Vector. Fe = N * t + N * b
    void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
            GeneralVariables& rVariables,
            Vector& rVolumeForce,
            const double& rIntegrationWeight) override;

     /// Initialize Material Properties on the Constitutive Law
    void InitializeMaterial(const ProcessInfo& rCurrentProcessInfo) override;

private:

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class MPMUpdatedLagrangianPQ
} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_PQ_H_INCLUDED  defined
