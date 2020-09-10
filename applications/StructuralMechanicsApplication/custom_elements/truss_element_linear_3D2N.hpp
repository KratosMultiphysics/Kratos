// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//

#if !defined(KRATOS_TRUSS_ELEMENT_LINEAR_3D2N_H_INCLUDED )
#define  KRATOS_TRUSS_ELEMENT_LINEAR_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/truss_element_3D2N.hpp"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos {
/**
 * @class TrussElementLinear3D2N
 *
 * @brief This is a linear 3D-2node truss element with 3 translational dofs per node inheriting from TrussElement3D2N
 *
 * @author Klaus B Sautter
 */

class TrussElementLinear3D2N : public TrussElement3D2N
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TrussElementLinear3D2N);

    TrussElementLinear3D2N() {};
    TrussElementLinear3D2N(IndexType NewId,
                           GeometryType::Pointer pGeometry);
    TrussElementLinear3D2N(IndexType NewId,
                           GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties);


    ~TrussElementLinear3D2N() override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function adds forces from prestressing to the force vector
     * @param rRightHandSideVector The right hand side of the problem
     */
    void AddPrestressLinear(VectorType& rRightHandSideVector);

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector< array_1d<double, 3 > >& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable, std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function calculates the total stiffness matrix for the element
     */
    BoundedMatrix<double,msLocalSize,msLocalSize>
    CreateElementStiffnessMatrix(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function calculates the original nodal postion for the transformation matrix
     * @param rReferenceCoordinates The original coordinates
     */
    void WriteTransformationCoordinates(
        BoundedVector<double,msLocalSize>& rReferenceCoordinates) override;

    /**
     * @brief This function calculates the current linear-Lagrange strain
     */
    double CalculateLinearStrain();


    /**
     * @brief This function updates the internal forces
     * @param rinternalForces The internal forces
     */

    void UpdateInternalForces(
        BoundedVector<double,msLocalSize>& rInternalForces) override;

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

private:

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}


#endif
