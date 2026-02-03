// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

#pragma once


// System includes

// External includes

// Project includes
#include "adjoint_finite_difference_base_element.h"

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

/** \brief AdjointFiniteDifferenceMembraneElement
 *
 * This element is inherited from AdjointFiniteDifferencingBaseElement.
 * It overwrites some functions necessary to do proper finite differencing with membrane elements
 */
template <typename TPrimalElement>
class AdjointFiniteDifferenceMembraneElement
    : public AdjointFiniteDifferencingBaseElement<TPrimalElement>
{
public:

    ///@name Type Definitions
    ///@{

    // redefine the typedefs because of templated base class
    typedef AdjointFiniteDifferencingBaseElement<TPrimalElement> BaseType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;
    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;
    typedef typename BaseType::DofsVectorType DofsVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::IntegrationMethod IntegrationMethod;
    typedef typename BaseType::GeometryDataType GeometryDataType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointFiniteDifferenceMembraneElement);
    ///@}

    ///@name Classes
    ///@{
    ///@}

    ///@name Life Cycle
    ///@{


    AdjointFiniteDifferenceMembraneElement(IndexType NewId = 0)
    : BaseType(NewId, false)
    {
    }

    AdjointFiniteDifferenceMembraneElement(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry, false)
    {
    }

    AdjointFiniteDifferenceMembraneElement(IndexType NewId,
                        typename GeometryType::Pointer pGeometry,
                        typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties, false)
    {
    }


    ///@}

    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceMembraneElement<TPrimalElement>>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                              typename GeometryType::Pointer pGeometry,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceMembraneElement<TPrimalElement>>(
            NewId, pGeometry, pProperties);
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    // void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
    //                                    const ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable,
    //                                 Matrix& rOutput,
    //                                 const ProcessInfo& rCurrentProcessInfo) override;

    // Sensitivity functions

    /**
     * Calculates the pseudo-load contribution of the element w.r.t. all properties which are available at the element.
     * This is done by finite differencing of the RHS of the primal element when perturbing a property value.
     * This operation is thread-save, because the property pointer is exchanged by a local one before pertubation.
     */
    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculates the pseudo-load of the design variable SHAPE_SENSITIVITY (coordinates of nodes) contribution of the element.
     * This is done by finite differencing of the RHS of the primal element when perturbing a nodal coordinate.
     * This operation is currently NOT thread-save!
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

    ///@name Public specialized Access - Temporary
    ///@{
    ///@}

protected:

    ///@name Protected Lyfe Cycle
    ///@{

    ///@}

private:

    ///@name Private Classes
    ///@{
    ///@}

    ///@name Private Operations
    ///@{

    void CheckDofs() const;
    void CheckProperties(const ProcessInfo& rCurrentProcessInfo) const;

    double GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) const override;

    ///@}
    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}

    ///@name Member Variables
    ///@{

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

};

}
