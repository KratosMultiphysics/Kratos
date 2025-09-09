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

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"
#include "custom_elements/data_containers/fluid_element_data.h"
#include "custom_elements/data_containers/two_fluid_fractional_navier_stokes/two_fluid_navier_stokes_fractional_convection_data.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
/*
    The "TwoFluidNavierStokesFractionalConvection" solves the convection problem of a vector field, specifically the fractional velocity, and convects this field using its own velocity.
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
template <class TElementData>
class TwoFluidNavierStokesFractionalConvection: public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TwoFluidNavierStokesFractionalConvection);

    using ElementData = TElementData;
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;
    static constexpr unsigned int Dim = TElementData::Dim;

    static constexpr unsigned int NumNodes = TElementData::NumNodes;

    static constexpr unsigned int BlockSize = Dim;

    static constexpr unsigned int LocalSize = NumNodes * BlockSize;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    TwoFluidNavierStokesFractionalConvection(IndexType NewId = 0);
    
    /**
     * @brief Constructor using an array of nodes.
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    TwoFluidNavierStokesFractionalConvection(
        IndexType NewId,
        const NodesArrayType &ThisNodes);

    /**
     * @brief Constructor using a geometry object.
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    TwoFluidNavierStokesFractionalConvection(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    /**
     * @brief Constructor using geometry and properties.
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    TwoFluidNavierStokesFractionalConvection(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Destructor.
    virtual ~TwoFluidNavierStokesFractionalConvection();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    /**
     @brief Returns a pointer to a new TwoFluidNavierStokesFractionalConvection element, created using given input.
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId, 
        NodesArrayType const& ThisNodes, 
        PropertiesType::Pointer pProperties) const override;


    /**
     * @brief Returns a pointer to a new FluidElement element, created using given input.
     * @param NewId the ID of the new element
     * @param pGeom a pointer to the geometry to be used to create the element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Given a distance function, computes the time integrated Left Hand Side (LHS)
     * and Right Hand Side elemental contributions for the two-fluid element.
     * @param rLeftHandSideMatrix elemental stiffness matrix
     * @param rRightHandSideVector elemental residual vector
     * @param rCurrentProcessInfo reference to the current process info
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Computes the LHS Gauss pt. contribution
     * This method computes the contribution to the LHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     */
    virtual void ComputeGaussPointLHSContribution(
        TElementData &rData,
        MatrixType &rLHS);

    /**
     * @brief Computes the RHS Gauss pt. contribution
     * This method computes the contribution to the RHS of a Gauss pt.
     * @param rData Reference to the element data container
     * @param rRHS Reference to the Right Hand Side vector to be filled
     */
    virtual void ComputeGaussPointRHSContribution(
        TElementData &rData,
        VectorType &rRHS);

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult, 
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType &rElementalDofList,
        const ProcessInfo &rCurrentProcessInfo) const override;
    
    /**
     * @brief Given a distance function, computes the time integrated Right Hand Side (RHS)
     * elemental contribution for the two-fluid element.
     * @param rRightHandSideVector elemental residual vector
     * @param rCurrentProcessInfo reference to the current process info
     */
    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief Computes time integrated LHS and RHS arrays
     * This method computes both the Left Hand Side and
     * Right Hand Side time integrated contributions.
     * @param rData Reference to the element data container
     * @param rLHS Reference to the Left Hand Side matrix to be filled
     * @param rRHS Reference to the Right Hand Side vector to be filled
     */
    void AddTimeIntegratedSystem(
        TElementData &rData,
        MatrixType &rLHS,
        VectorType &rRHS);

    /**
     * @brief Set up the element's data and for the current integration point.
     * @param[in/out] rData Container for the current element's data.
     *  @param[in] Weight Integration point weight.
     *  @param[in] rN Values of nodal shape functions at the integration point.
     *  @param[in] rDN_DX Values of nodal shape function gradients at the integration point.
     */
    void UpdateIntegrationPointData(
        TElementData &rData,
        unsigned int IntegrationPointIndex,
        double Weight,
        const typename TElementData::MatrixRowType &rN,
        const typename TElementData::ShapeDerivativesType &rDN_DX) const;
    
    /**
     * @brief Computes shape function data for all the gauss points
     *
     * @param rGaussWeights         Gauss point weights
     * @param rNContainer           Gauss point shape functions (each row corresponds to a specific gauss point)
     * @param rDN_DX                Gauss point shape function gradients (vector of matrices)
     */
    virtual void CalculateGeometryData(
        Vector &rGaussWeights,
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
        return "TwoFluidNavierStokesFractionalConvection #";
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
    
private:
    ///@name Static Member Variables
    ///@{
    static constexpr double stab_c2 = 2.0;
    static constexpr double stab_c1= 4.0;


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


} // namespace Kratos.
