// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

#if !defined(KRATOS_FIBER_BEAM_COLUMN_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_FIBER_BEAM_COLUMN_ELEMENT_3D2N_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_conditions/fiber_beam_column_section.hpp"
#include "custom_conditions/fiber_beam_column_uniaxial_fiber.hpp"

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
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) FiberBeamColumnElement3D2N : public Element
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

    /// Constructor
    FiberBeamColumnElement3D2N(IndexType NewId = 0);

    /// Constructor using an array of nodes
    FiberBeamColumnElement3D2N(IndexType NewId, const NodesArrayType& rThisNodes);

    /// Constructor using Geometry
    FiberBeamColumnElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constructor using Properties
    FiberBeamColumnElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Copy Constructor
    FiberBeamColumnElement3D2N(FiberBeamColumnElement3D2N const& rOther);

    /// Destructor
    ~FiberBeamColumnElement3D2N() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment Operator
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

    void Initialize() override;

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

    /// FIXME: Should the derivatives be implemented?

    // /**
    //  * @brief This function calculates the total stiffness matrix for the element
    //  */
    // virtual BoundedMatrix<double,msLocalSize,msLocalSize>
    //     CreateElementStiffnessMatrix(ProcessInfo& rCurrentProcessInfo);

    // void CalculateOnIntegrationPoints(
    //     const Variable<double>& rVariable,
    //     std::vector<double>& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // void GetValueOnIntegrationPoints(
    //     const Variable<double>& rVariable,
    //     std::vector<double>& rValues,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // void GetValueOnIntegrationPoints(
    //     const Variable<array_1d<double, 3 > >& rVariable,
    //     std::vector< array_1d<double, 3 > >& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // /**
    //  * @brief This function updates the internal normal force w.r.t. the current deformations
    //  * @param rinternalForces The current updated internal forces
    //  */
    // virtual void UpdateInternalForces(BoundedVector<double,msLocalSize>& rInternalForces);

    // void CalculateOnIntegrationPoints(
    //     const Variable<Vector>& rVariable,
    //     std::vector<Vector>& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // void GetValueOnIntegrationPoints(
    //     const Variable<Vector>& rVariable,
    //     std::vector<Vector>& rValues,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateOnIntegrationPoints(
    //     const Variable<array_1d<double, 3 > >& rVariable,
    //     std::vector< array_1d<double, 3 > >& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateMassMatrix(
    //     MatrixType& rMassMatrix,
    //     ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateDampingMatrix(
    //     MatrixType& rDampingMatrix,
    //     ProcessInfo& rCurrentProcessInfo) override;


    // /**
    //  * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable rDestinationVariable (double version)
    //  * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
    //  * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
    //  * @param rRHSVector input variable containing the RHS vector to be assembled
    //  * @param rRHSVariable variable describing the type of the RHS vector to be assembled
    //  * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
    //  * @param rCurrentProcessInfo the current process info instance
    //  */
    // void AddExplicitContribution(
    //     const VectorType& rRHSVector,
    //     const Variable<VectorType>& rRHSVariable,
    //     Variable<double >& rDestinationVariable,
    //     const ProcessInfo& rCurrentProcessInfo
    //     ) override;

    // /**
    //  * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable (array_1d<double, 3>) version rDestinationVariable.
    //  * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
    //  * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
    //  * @param rRHSVector input variable containing the RHS vector to be assembled
    //  * @param rRHSVariable variable describing the type of the RHS vector to be assembled
    //  * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
    //  * @param rCurrentProcessInfo the current process info instance
    //  */
    // void AddExplicitContribution(const VectorType& rRHSVector,
    //     const Variable<VectorType>& rRHSVariable,
    //     Variable<array_1d<double, 3> >& rDestinationVariable,
    //     const ProcessInfo& rCurrentProcessInfo
    //     ) override;


    // void GetValuesVector(
    //     Vector& rValues,
    //     int Step = 0) override;

    // void GetSecondDerivativesVector(
    //     Vector& rValues,
    //     int Step = 0) override;

    // void GetFirstDerivativesVector(
    //     Vector& rValues,
    //     int Step = 0) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    // int Check(const ProcessInfo& rCurrentProcessInfo) override;

    // /**
    //  * @brief This function calculates the current Green-Lagrange strain
    //  */
    // double CalculateGreenLagrangeStrain();

    // /**
    //  * @brief This function calculates self-weight forces
    //  */
    // BoundedVector<double,msLocalSize> CalculateBodyForces();

    // /**
    //  * @brief This function assembles the geometric stiffness part of the total stiffness matrix
    //  * @param rGeometricStiffnessMatrix The geometric stiffness matrix
    //  * @param rCurrentProcessInfo The current process information
    //  */
    // void CalculateGeometricStiffnessMatrix(BoundedMatrix<double,msLocalSize,msLocalSize>& rGeometricStiffnessMatrix,
    //     ProcessInfo& rCurrentProcessInfo);

    // /**
    //  * @brief This function assembles the elastic stiffness part of the total stiffness matrix
    //  * @param rElasticStiffnessMatrix The elastic stiffness matrix
    //  * @param rCurrentProcessInfo The current process information
    //  */
    // void CalculateElasticStiffnessMatrix(BoundedMatrix<double,msLocalSize,msLocalSize>& rElasticStiffnessMatrix,
    //     ProcessInfo& rCurrentProcessInfo);

    // /**
    //  * @brief This function calculates the current nodal postion for the transformation matrix
    //  * @param rReferenceCoordinates The current coordinates
    //  */
    // virtual void WriteTransformationCoordinates(
    //     BoundedVector<double,msLocalSize>& rReferenceCoordinates);



    // void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    // void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    // void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    // /**
    //  * @brief This function checks if self weight is present
    //  */
    // bool HasSelfWeight() const;

    // /**
    //  * @brief This function calls the constitutive law to get stresses
    //  * @param rCurrentProcessInfo Current process info
    //  * @param rSaveInternalVariables Boolean to save internal constit. law variables
    //  */
    // virtual BoundedVector<double,msLocalSize> GetConstitutiveLawTrialResponse(
    //     const ProcessInfo& rCurrentProcessInfo);

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

    static constexpr int msNumberOfNodes = 2;
    static constexpr int msDimension = 3;
    static constexpr unsigned int msLocalSize = 5;
    static constexpr unsigned int msGlobalSize = 12;

    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector<FiberBeamColumnSection> mSections;
    /// FIXME: we need the constitutive law on the fiber level ...
    // ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    BoundedMatrix<double,msLocalSize,msLocalSize> CreateElementLocalStiffnessMatrix(ProcessInfo&) const;
    BoundedVector<double,msLocalSize> CreateElementLocalResistingForces(ProcessInfo&) const;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    /**
     * @brief This function calculates the transformation matrix to globalize vectors and/or matrices
     */
    BoundedMatrix<double,msGlobalSize,msLocalSize> CreateTransformationMatrix() const;

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
