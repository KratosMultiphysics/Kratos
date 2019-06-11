// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Jordi Cotela
//

#if !defined(KRATOS_ADJOINT_HEAT_DIFFUSION_ELEMENT_H_INCLUDED )
#define  KRATOS_ADJOINT_HEAT_DIFFUSION_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{

///@addtogroup ConvectionDiffusionApplication
///@{

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

/// Basic element for the ajdoint thermal diffusion problem.
template< class PrimalElement >
class AdjointHeatDiffusionElement: public PrimalElement
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AdjointHeatDiffusionElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointHeatDiffusionElement);

    using IndexType = typename PrimalElement::IndexType;
    using GeometryType = typename PrimalElement::GeometryType;
    using NodesArrayType = typename PrimalElement::NodesArrayType;
    using MatrixType = typename PrimalElement::MatrixType;
    using VectorType = typename PrimalElement::VectorType;
    using EquationIdVectorType = typename PrimalElement::EquationIdVectorType;
    using DofsVectorType = typename PrimalElement::DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    AdjointHeatDiffusionElement(IndexType NewId, typename GeometryType::Pointer pGeometry,  Properties::Pointer pProperties);

    /// Destructor.
    ~AdjointHeatDiffusionElement() override;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  Properties::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, typename GeometryType::Pointer pGeom,  Properties::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    /*void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;*/

    void GetValuesVector(Vector& rValues, int Step = 0) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

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

    // A private default constructor necessary for serialization
    AdjointHeatDiffusionElement() : PrimalElement()
    {
    }

    void save(Serializer& rSerializer) const override
    {
        using BaseType = PrimalElement;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) override
    {
        using BaseType = PrimalElement;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    AdjointHeatDiffusionElement& operator=(const AdjointHeatDiffusionElement& rOther) = delete;

    /// Copy constructor.
    AdjointHeatDiffusionElement(const AdjointHeatDiffusionElement& rOther) = delete;

    ///@}

}; // Class AdjointHeatDiffusionElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class PrimalElement >
inline std::istream& operator >> (std::istream& rIStream, AdjointHeatDiffusionElement<PrimalElement>& rThis)
{
    return rIStream;
}

/// output stream function
template< class PrimalElement >
inline std::ostream& operator << (std::ostream& rOStream, const AdjointHeatDiffusionElement<PrimalElement>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@}

}  // namespace Kratos.

#endif // KRATOS_ADJOINT_HEAT_DIFFUSION_ELEMENT_H_INCLUDED  defined


