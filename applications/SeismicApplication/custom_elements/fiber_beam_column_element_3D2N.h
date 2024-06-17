//  KRATOS  ____       _               _
//         / ___|  ___(_)___ _ __ ___ (_) ___
//         \___ \ / _ \ / __| '_ ` _ \| |/ __|
//          ___) |  __/ \__ \ | | | | | | (__
//         |____/ \___|_|___/_| |_| |_|_|\___|
//
//  License:     BSD License
//  license:     structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//    Co authors: Long Chen
//

#if !defined(KRATOS_FIBER_BEAM_COLUMN_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_FIBER_BEAM_COLUMN_ELEMENT_3D2N_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_classes/fiber_beam_column_section.h"

namespace Kratos
{

///@}
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
///@name  Kratos Classes
///@{

/**
 * @class FiberBeamColumnElement3D2N
 *
 * @brief A 3D-2node fiber beam-column element for reinforced concrete modeling
 * @details A non-linear flexibility-based fiber element for dynamic loading
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(SEISMIC_APPLICATION) FiberBeamColumnElement3D2N : public Element
{

public:

    ///@}
    ///@name Type Definitions
    ///@{

    typedef Element                                     BaseType;
    typedef BaseType::GeometryType                  GeometryType;
    typedef BaseType::NodesArrayType              NodesArrayType;
    typedef BaseType::PropertiesType              PropertiesType;
    typedef BaseType::IndexType                        IndexType;
    typedef BaseType::SizeType                          SizeType;
    typedef BaseType::MatrixType                      MatrixType;
    typedef BaseType::VectorType                      VectorType;
    typedef BaseType::EquationIdVectorType  EquationIdVectorType;
    typedef BaseType::DofsVectorType              DofsVectorType;

    ///@}
    ///@name Pointer Definitions
    ///@brief Pointer definition of FiberBeamColumnElement3D2N
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FiberBeamColumnElement3D2N);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor
     */
    FiberBeamColumnElement3D2N(IndexType NewId = 0);

    /**
     * Constructor using an array of nodes
     */
    FiberBeamColumnElement3D2N(IndexType NewId, const NodesArrayType& rThisNodes);

    /**
     * Constructor using Geometry
     */
    FiberBeamColumnElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);

    /**
     * Constructor using Properties
     */
    FiberBeamColumnElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
     * Copy Constructor
     */
    FiberBeamColumnElement3D2N(FiberBeamColumnElement3D2N const& rOther);

    /**
     * Destructor
     */
    ~FiberBeamColumnElement3D2N() override;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment Operator
     */
    FiberBeamColumnElement3D2N & operator=(FiberBeamColumnElement3D2N const& rOther);

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

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    /**
     * constructs the sections and the fibers from the properties,
     * initializes the sections and the fibers, calculates the
     * transformation matrix, and initializes the local stiffness matrix.
     */
    void Initialize() override;

    /**
     * called after the non linear iteration so that the equilibrium loop is performed
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * called when the non linear iterations have converged
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Gets the vector of displacements and rotations
     * @param rValues: the vector of the displacements and rotations
     * @param Step: solution step index
     */
    void GetValuesVector(Vector& rValues, int Step = 0) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}


protected:

    ///@name Protected static Member Variables
    ///@{

    static constexpr unsigned int msNumberOfNodes = 2;
    static constexpr unsigned int msDimension = 3;
    static constexpr unsigned int msLocalSize = 5;
    static constexpr unsigned int msGlobalSize = 12;

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
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}


private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::vector<FiberBeamColumnSection> mSections;  // vector of the sections

    Matrix mTransformationMatrix = ZeroMatrix(msGlobalSize, msLocalSize);
    Matrix mLocalStiffnessMatrix = ZeroMatrix(msLocalSize, msLocalSize);

    Vector mDeformationI   = ZeroVector(msLocalSize);  // deformations of the current iteration
    Vector mDeformationIM1 = ZeroVector(msLocalSize);  // deformations of the previous iteration
    Vector mInternalForces = ZeroVector(msLocalSize);  // internal forces from the sections
    Vector mDeformationResiduals = ZeroVector(msLocalSize);

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief returns local stiffness matrix 5x5
     * @param rCurrentProcessInfo Current process info
     */
    void CalculateElementLocalStiffnessMatrix();

    /**
     * @brief returns local internal forces 5x1
     * @param rCurrentProcessInfo Current process info
     */
    Vector CreateElementLocalResistingForces(ProcessInfo&);

    /**
     * @brief returns gauss lobatto integration method based on number of sections
     */
    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    /**
     * @brief This function calculates the transformation matrix to globalize vectors and/or matrices
     * return 12x5 matrix
     */
    void CalculateTransformationMatrix();

    /**
     * @return matrix of the unit vector for the beam local coord sys
     */
    Matrix CreateInitialLocalCoordSys() const;

    /**
     * @return reference length of the element
     * @see StructuralMechanicsUtilities::CalculateReferenceLength2D2N
     */
    double CalculateReferenceLength() const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

};  // class FiberBeamColumnElement3D2N

}  // namespace Kratos


#endif
