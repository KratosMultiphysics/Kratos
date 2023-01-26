//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Mohammad R. Hashemi
//
//


#if !defined(KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_ALGEBRAIC_STABILIZATION )
#define  KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_ALGEBRAIC_STABILIZATION

// System includes

// External includes

// Project includes
#include "elements/levelset_convection_element_simplex.h"

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

/// Formulation is based on the works by Kuzmin et al. (especifically see Comput. Methods Appl. Mech. Engrg. 322 (2017) 23–41)
/// Dirichlet boundary condition for rUnknownVar at velocity inlets/outles is essential to be set for this solver since it is based on the flux
template<unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(KRATOS_CORE) LevelSetConvectionElementSimplexAlgebraicStabilization
    : public LevelSetConvectionElementSimplex<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LevelSetConvectionElementSimplexAlgebraicStabilization);

    typedef Element::GeometryType GeometryType;

    typedef Element::PropertiesType PropertiesType;

    typedef Element::NodesArrayType NodesArrayType;

    typedef Element::MatrixType MatrixType;

    typedef Element::VectorType VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LevelSetConvectionElementSimplexAlgebraicStabilization();

    LevelSetConvectionElementSimplexAlgebraicStabilization(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    LevelSetConvectionElementSimplexAlgebraicStabilization(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    ~LevelSetConvectionElementSimplexAlgebraicStabilization() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

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

    std::string Info() const override
    {
        return "LevelSetConvectionElementSimplexAlgebraicStabilization #";
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << this->Id();
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


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

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

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
    ///@name Un accessible methods
    ///@{


    ///@}

};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


} // namespace Kratos.

#endif // KRATOS_LEVELSET_CONVECTION_ELEMENT_SIMPLEX_ALGEBRAIC_STABILIZATION  defined


