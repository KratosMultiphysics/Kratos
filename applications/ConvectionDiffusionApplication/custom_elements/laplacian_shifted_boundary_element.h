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
#include "laplacian_element.h"

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
class LaplacianShiftedBoundaryElement : public LaplacianElement
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of LaplacianShiftedBoundaryElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LaplacianShiftedBoundaryElement);

    typedef LaplacianElement BaseType;

    static constexpr std::size_t NumNodes = TDim + 1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with geometry
    LaplacianShiftedBoundaryElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    /// Constructor with geometry and properties
    LaplacianShiftedBoundaryElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Destructor
    virtual ~LaplacianShiftedBoundaryElement() = default;

    /// Copy constructor
    LaplacianShiftedBoundaryElement(const LaplacianShiftedBoundaryElement& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    LaplacianShiftedBoundaryElement& operator=(const LaplacianShiftedBoundaryElement& rOther) = delete;

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
    LaplacianShiftedBoundaryElement() : LaplacianElement() {}

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LaplacianElement);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LaplacianElement);
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
}; // Class LaplacianShiftedBoundaryElement

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
    LaplacianShiftedBoundaryElement<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const LaplacianShiftedBoundaryElement<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
}  // namespace Kratos.

