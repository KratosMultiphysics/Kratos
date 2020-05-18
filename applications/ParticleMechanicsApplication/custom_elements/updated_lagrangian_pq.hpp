//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_PQ_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_PQ_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian.hpp"

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

/// Large Displacement Lagrangian Element for 3D and 2D geometries. (base class)

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */

class UpdatedLagrangianPQ
    : public UpdatedLagrangian
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
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UpdatedLagrangianPQ );
    ///@}

public:

    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    UpdatedLagrangianPQ();


    /// Default constructors
    UpdatedLagrangianPQ(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianPQ(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianPQ(UpdatedLagrangianPQ const& rOther);

    /// Destructor.
    ~UpdatedLagrangianPQ() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangianPQ& operator=(UpdatedLagrangianPQ const& rOther);

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

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MPM Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MPM Element #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        GetGeometry().PrintData(rOStream);
    }

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    //MaterialPointVariables mMP;

    /**
     * Container for historical total elastic deformation measure F0 = dx/dX
     */
    //Matrix mDeformationGradientF0;
    /**
     * Container for the total deformation gradient determinants
     */
    //double mDeterminantF0;

    /**
     * Container for constitutive law instances on each integration point
     */
    //ConstitutiveLaw::Pointer mConstitutiveLawVector;


    /**
     * Finalize and Initialize label
     */
    //bool mFinalizedStep;


    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
            GeneralVariables& rVariables,
            Vector& rVolumeForce,
            const double& rIntegrationWeight) override;

    /// Calculation of the Explicit Stresses from velocity gradient.
    void CalculateExplicitStresses(const ProcessInfo& rCurrentProcessInfo,
        GeneralVariables& rVariables) override;

    /**
     * Initialize Material Properties on the Constitutive Law
     */
    void InitializeMaterial() override;

    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    SizeType mMPSubPoints;


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

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class UpdatedLagrangianPQ

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_PQ_H_INCLUDED  defined
