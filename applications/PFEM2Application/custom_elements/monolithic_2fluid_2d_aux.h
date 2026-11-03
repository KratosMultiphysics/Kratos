//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2026-10-27 $
//   Revision:            $Revision: 1.0 $
//
//

#ifndef KRATOS_MONOLITHIC_PFEM2_2D_AUX_ELEM_H_INCLUDED
#define KRATOS_MONOLITHIC_PFEM2_2D_AUX_ELEM_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class MonolithicPFEM22DAux
 * @brief Auxiliary 2D element that calculates the standard FEM inner product (w, f) = \int w * f d\Omega
 * @details Evaluates the consistent mass/weight matrix M_ij = \int N_i N_j d\Omega and the RHS
 *          inner product vector (w, f)_i = \int N_i f d\Omega for the scalar TEMPERATURE DOF,
 *          where f is an interpolated nodal source term (e.g. HEAT_FLUX, HEAT_SOURCE, or BODY_FORCE_X).
 */
class MonolithicPFEM22DAux : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MonolithicPFEM22DAux
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MonolithicPFEM22DAux);

    typedef Element BaseType;
    typedef Properties PropertiesType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef std::vector<std::size_t> EquationIdVectorType;
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;
    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MonolithicPFEM22DAux(IndexType NewId = 0) : Element(NewId) {}

    MonolithicPFEM22DAux(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes) {}

    MonolithicPFEM22DAux(IndexType NewId, GeometryType::Pointer pGeometry);

    MonolithicPFEM22DAux(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~MonolithicPFEM22DAux() override;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    void InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return "MonolithicPFEM22DAux";
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void GetNodalSource(array_1d<double, 3>& rSource) const;

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Unaccessible methods
    ///@{

    MonolithicPFEM22DAux& operator=(MonolithicPFEM22DAux const& rOther);

    MonolithicPFEM22DAux(MonolithicPFEM22DAux const& rOther);

    ///@}
}; // Class MonolithicPFEM22DAux

///@}

} // namespace Kratos

#endif // KRATOS_MONOLITHIC_PFEM2_2D_AUX_ELEM_H_INCLUDED
