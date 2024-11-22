// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "mixed_laplacian_element.h"

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

template<std::size_t TDim>
class MixedLaplacianShiftedBoundaryElement : public MixedLaplacianElement<TDim, TDim+1>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MixedLaplacianShiftedBoundaryElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MixedLaplacianShiftedBoundaryElement);

    static constexpr std::size_t NumNodes = TDim + 1;

    static constexpr std::size_t BlockSize = TDim + 1;

    static constexpr std::size_t LocalSize = NumNodes*BlockSize;

    typedef MixedLaplacianElement<TDim, NumNodes> BaseType;

    typedef typename BaseType::GeometryType GeometryType;

    typedef typename BaseType::PropertiesType PropertiesType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::MatrixType MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with geometry
    MixedLaplacianShiftedBoundaryElement(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry);

    /// Constructor with geometry and properties
    MixedLaplacianShiftedBoundaryElement(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry,
        typename PropertiesType::Pointer pProperties);

    /// Destructor
    virtual ~MixedLaplacianShiftedBoundaryElement() = default;

    /// Copy constructor
    MixedLaplacianShiftedBoundaryElement(const MixedLaplacianShiftedBoundaryElement& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MixedLaplacianShiftedBoundaryElement& operator=(const MixedLaplacianShiftedBoundaryElement& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        typename BaseType::NodesArrayType const& ThisNodes,
        typename PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        typename GeometryType::Pointer pGeom,
        typename PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


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

    // Protected default constructor necessary for serialization
    MixedLaplacianShiftedBoundaryElement() : MixedLaplacianElement<TDim, NumNodes>() {}

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    std::vector<std::size_t> GetSurrogateFacesIds();

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
}; // Class MixedLaplacianShiftedBoundaryElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<std::size_t TDim>
inline std::istream& operator >> (
    std::istream& rIStream,
    MixedLaplacianShiftedBoundaryElement<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const MixedLaplacianShiftedBoundaryElement<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
}  // namespace Kratos.
