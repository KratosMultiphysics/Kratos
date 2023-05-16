// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_GEO_STEADY_STATE_PW_INTERFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_STEADY_STATE_PW_INTERFACE_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/transient_Pw_interface_element.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) SteadyStatePwInterfaceElement :
    public TransientPwInterfaceElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( SteadyStatePwInterfaceElement );

    using BaseType = TransientPwInterfaceElement<TDim, TNumNodes>;

    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;

    using DofsVectorType = Element::DofsVectorType;
    using EquationIdVectorType = Element::EquationIdVectorType;

    /// The definition of the sizetype
    using SizeType = std::size_t;

    using BaseType::mRetentionLawVector;
    using BaseType::mThisIntegrationMethod;
    using BaseType::CalculateRetentionResponse;

    using InterfaceElementVariables = typename BaseType::InterfaceElementVariables;
    using SFGradAuxVariables = typename BaseType::SFGradAuxVariables;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    SteadyStatePwInterfaceElement(IndexType NewId = 0) : TransientPwInterfaceElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    SteadyStatePwInterfaceElement(IndexType NewId,
                                const NodesArrayType& ThisNodes) : TransientPwInterfaceElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    SteadyStatePwInterfaceElement(IndexType NewId,
                                GeometryType::Pointer pGeometry) : TransientPwInterfaceElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    SteadyStatePwInterfaceElement(IndexType NewId,
                                GeometryType::Pointer pGeometry,
                                PropertiesType::Pointer pProperties)
                                : TransientPwInterfaceElement<TDim,TNumNodes>( NewId, pGeometry, pProperties )
    {}

    /// Destructor
    ~SteadyStatePwInterfaceElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;


protected:

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& CurrentProcessInfo,
                       const bool CalculateStiffnessMatrixFlag,
                       const bool CalculateResidualVectorFlag) override;

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                            InterfaceElementVariables& rVariables) override;

    void CalculateAndAddRHS(VectorType& rRightHandSideVector,
                            InterfaceElementVariables& rVariables,
                            unsigned int GPoint) override;


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    /// Assignment operator.
    SteadyStatePwInterfaceElement & operator=(SteadyStatePwInterfaceElement const& rOther);

    /// Copy constructor.
    SteadyStatePwInterfaceElement(SteadyStatePwInterfaceElement const& rOther);

}; // Class SteadyStatePwInterfaceElement

} // namespace Kratos

#endif // KRATOS_GEO_STEADY_STATE_PW_INTERFACE_ELEMENT_H_INCLUDED  defined
