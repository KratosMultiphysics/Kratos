// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter,
//                   Vahid Galavi
//

#if !defined(KRATOS_GEO_TRUSS_ELEMENT_BASE_H_INCLUDED)
#define KRATOS_GEO_TRUSS_ELEMENT_BASE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"

namespace Kratos
{
/**
 * @class GeoTrussElementBase
 *
 * @brief This is a 2D-2node truss element with 2 translational dofs per node
 *
 * @author Vahid Galavi
 */
template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoTrussElementBase : public Element
{
protected:
    ConstitutiveLaw::Pointer  mpConstitutiveLaw = nullptr;
    static constexpr SizeType NDof              = TDim * TNumNodes;
    static constexpr SizeType DIM               = TDim;
    static constexpr SizeType NUM_NODES         = TNumNodes;

public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoTrussElementBase);

    using BaseType             = Element;
    using GeometryType         = BaseType::GeometryType;
    using NodesArrayType       = BaseType::NodesArrayType;
    using PropertiesType       = BaseType::PropertiesType;
    using IndexType            = BaseType::IndexType;
    using SizeType             = BaseType::SizeType;
    using MatrixType           = BaseType::MatrixType;
    using VectorType           = BaseType::VectorType;
    using EquationIdVectorType = BaseType::EquationIdVectorType;
    using DofsVectorType       = BaseType::DofsVectorType;
    using FullDofMatrixType    = BoundedMatrix<double, TDim * TNumNodes, TDim * TNumNodes>;
    using FullDofVectorType    = BoundedVector<double, TDim * TNumNodes>;

    GeoTrussElementBase(){};
    GeoTrussElementBase(IndexType NewId, GeometryType::Pointer pGeometry);
    GeoTrussElementBase(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~GeoTrussElementBase() override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function calculates the total stiffness matrix for the element
     */
    virtual void CreateElementStiffnessMatrix(MatrixType&        rLocalStiffnessMatrix,
                                              const ProcessInfo& rCurrentProcessInfo);

    void Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    /**
     * @brief This function updates the internal normal force w.r.t. the current deformations
     * @param rinternalForces The current updated internal forces
     */
    virtual void UpdateInternalForces(FullDofVectorType& rInternalForces, const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief This function calculates the transformation matrix to globalize vectors and/or matrices
     * @param rRotationMatrix The transformation matrix
     */
    void CreateTransformationMatrix(BoundedMatrix<double, TDim * TNumNodes, TDim * TNumNodes>& rRotationMatrix);

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateConsistentMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) const;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable rDestinationVariable (double version)
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const VectorType&           rRHSVector,
                                 const Variable<VectorType>& rRHSVariable,
                                 const Variable<double>&     rDestinationVariable,
                                 const ProcessInfo&          rCurrentProcessInfo) override;

    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable (array_1d<double, 3>) version rDestinationVariable.
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const VectorType&                    rRHSVector,
                                 const Variable<VectorType>&          rRHSVariable,
                                 const Variable<array_1d<double, 3>>& rDestinationVariable,
                                 const ProcessInfo&                   rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief This function calculates the current Green-Lagrange strain
     */
    double CalculateGreenLagrangeStrain() const;

    /**
     * @brief This function calculates self-weight forces
     */
    void CalculateBodyForces(FullDofVectorType& rGlobalBodyForces);

    /**
     * @brief This function assembles the geometric stiffness part of the total stiffness matrix
     * @param rGeometricStiffnessMatrix The geometric stiffness matrix
     * @param rCurrentProcessInfo The current process information
     */
    void CalculateGeometricStiffnessMatrix(BoundedMatrix<double, TDim * TNumNodes, TDim * TNumNodes>& rGeometricStiffnessMatrix,
                                           const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief This function assembles the elastic stiffness part of the total stiffness matrix
     * @param rElasticStiffnessMatrix The elastic stiffness matrix
     * @param rCurrentProcessInfo The current process information
     */
    void CalculateElasticStiffnessMatrix(MatrixType& rElasticStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief This function calculates the current nodal postion for the transformation matrix
     * @param rReferenceCoordinates The current coordinates
     */
    virtual void WriteTransformationCoordinates(FullDofVectorType& rReferenceCoordinates);

    double ReturnTangentModulus1D(const ProcessInfo& rCurrentProcessInfo);

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

private:
    /**
     * @brief This method computes directly the lumped mass vector
     * @param rMassVector The lumped mass vector
     */
    void CalculateLumpedMassVector(VectorType& rMassVector, const ProcessInfo& rCurrentProcessInfo) const override;

    [[nodiscard]] Element::DofsVectorType GetDofs() const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos

#endif
