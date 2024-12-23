//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//
//
//

#pragma once
#define KRATOS_VECTORIAL_CONVECTION_FRACTIONAL_ELEMENT

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"
#include "custom_elements/data_containers/fluid_element_data.h"
#include "custom_elements/data_containers/two_fluid_fractional_navier_stokes/vectorial_convection_fractional_element_data.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
/*
    The "VectorialConvectionFractionalElement" solves the convection problem of a vector field, specifically the fractional velocity, and convects this field using its own velocity.
    This element is part of the two-fluid Navier-Stokes fractional element TwoFluidNavierStokesFractional.
    The Navier-Stokes momentum conservation equation is split into two problems: the first one convects the fractional velocity (using this element), and the second one resolves the remaining terms of the Navier-Stokes momentum conservation.
    Combining both problems in the continuous framework results in the original Navier-Stokes momentum equation
 */

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
    namespace Internals
    {
        template <class TElementData, bool TDataKnowsAboutTimeIntegration>
        class FluidElementTimeIntegrationDetail;
}
template <class TElementData>
class VectorialConvectionFractionalElement: public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(VectorialConvectionFractionalElement);

    using ElementData = TElementData;
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;
    static constexpr unsigned int Dim = TElementData::Dim;

    static constexpr unsigned int NumNodes = TElementData::NumNodes;

    static constexpr unsigned int BlockSize = Dim;

    static constexpr unsigned int LocalSize = NumNodes * BlockSize;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    VectorialConvectionFractionalElement(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    VectorialConvectionFractionalElement(IndexType NewId, const NodesArrayType &ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    VectorialConvectionFractionalElement(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    VectorialConvectionFractionalElement(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~VectorialConvectionFractionalElement();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;
    Element::Pointer Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;
    virtual void ComputeGaussPointLHSContribution(TElementData &rData,MatrixType &rLHS);
    virtual void ComputeGaussPointRHSContribution(TElementData &rData,VectorType &rRHS);

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;
    void GetDofList(
        DofsVectorType &rElementalDofList,
        const ProcessInfo &rCurrentProcessInfo) const override;
    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;
    virtual void AddTimeIntegratedSystem(
        TElementData &rData,
        MatrixType &rLHS,
        VectorType &rRHS);
    virtual void UpdateIntegrationPointData(TElementData &rData,
        unsigned int IntegrationPointIndex,
        double Weight,
        const typename TElementData::MatrixRowType &rN,
        const typename TElementData::ShapeDerivativesType &rDN_DX) const;
    virtual void CalculateGeometryData(Vector &rGaussWeights,
                                       Matrix &rNContainer,
                                       ShapeFunctionDerivativesArrayType &rDN_DX) const;
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
        return "VectorialConvectionFractionalElement #";
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;


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

    bool mElementTauNodal; // Flag to indicate if the stabilization tau is evaluated at each Gauss point or interpolated

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;
    //         ASGS2D() : Element()
    //         {
    //         }

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


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                                    Fluid2DASGS& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                                    const Fluid2DASGS& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

} // namespace Kratos.
