// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo' Antonelli
//

#if !defined(KRATOS_CONV_DIFF_SHIFTED_BOUNDARY_ELEMENT_H_INCLUDED )
#define  KRATOS_CONV_DIFF_SHIFTED_BOUNDARY_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "includes/serializer.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "eulerian_conv_diff.h"




namespace Kratos
{
template< unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) ConvDiffShiftedBoundaryElement
    : public EulerianConvectionDiffusionElement<2,3>      // Solo qui devo specificare il template
{
public:

    typedef EulerianConvectionDiffusionElement<2,3> BaseType;
    // typedef typename BaseType::GeometryType GeometryType;
    // typedef typename GeometryType::Pointer  GeometryTypePtr;
    // typedef typename BaseType::PropertiesType PropertiesType;
    // typedef typename PropertiesType::Pointer  PropertiesTypePtr;
    // typedef typename BaseType::NodesArrayType NodesArrayType;
    // typedef typename NodesArrayType::Pointer  NodesArrayTypePtr;
    // typedef typename Element::VectorType VectorType;

    static constexpr std::size_t NumNodes = TDim + 1;

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConvDiffShiftedBoundaryElement);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default constructor. (// qui NON devo specificare il template)
    ConvDiffShiftedBoundaryElement() : EulerianConvectionDiffusionElement()   
    {
    }

    ConvDiffShiftedBoundaryElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : EulerianConvectionDiffusionElement(NewId, pGeometry)
    {}

    ConvDiffShiftedBoundaryElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : EulerianConvectionDiffusionElement(NewId, pGeometry, pProperties)
    {}


    /// Destructor.
    virtual ~ConvDiffShiftedBoundaryElement(){
    //destructor body
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // It is fundamental? Since it is defined in the parent
    Element::Pointer Create(
        IndexType NewId, 
        NodesArrayType const& ThisNodes, 
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<ConvDiffShiftedBoundaryElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }
    
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<ConvDiffShiftedBoundaryElement>(NewId, pGeom, pProperties);
    }

    ///@}
    ///@name Operators
    ///@{

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


protected:
    ///@name Protected static Member Variables



private:
    ///@name Static Member Variables
    ///@{
    std::vector<std::size_t> GetSurrogateFacesIds();

}; // Class ConvDiffShiftedBoundaryElement


}  // namespace Kratos.

#endif // KRATOS_CONV_DIFF_SHIFTED_BOUNDARY_ELEMENT_H_INCLUDED  defined