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

#if !defined(KRATOS_GEO_STEADY_STATE_PW_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_STEADY_STATE_PW_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/comparison_utilities.hpp"
#include "custom_elements/transient_Pw_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) SteadyStatePwElement :
    public TransientPwElement<TDim, TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( SteadyStatePwElement );

    typedef TransientPwElement<TDim, TNumNodes> BaseType;

    typedef std::size_t IndexType;
    typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef Element::EquationIdVectorType EquationIdVectorType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;
    using BaseType::mRetentionLawVector;
    using BaseType::mConstitutiveLawVector;
    using BaseType::mIsInitialised;
    using BaseType::CalculateRetentionResponse;

    typedef typename BaseType::ElementVariables ElementVariables;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    SteadyStatePwElement(IndexType NewId = 0) : BaseType( NewId ) {}

    /// Constructor using an array of nodes
    SteadyStatePwElement(IndexType NewId,
              const NodesArrayType& ThisNodes) : BaseType(NewId, ThisNodes) {}

    /// Constructor using Geometry
    SteadyStatePwElement(IndexType NewId,
              GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry) {}

    /// Constructor using Properties
    SteadyStatePwElement(IndexType NewId,
              GeometryType::Pointer pGeometry,
              PropertiesType::Pointer pProperties) : BaseType( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~SteadyStatePwElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create( IndexType NewId,
                             NodesArrayType const& ThisNodes,
                             PropertiesType::Pointer pProperties ) const override;

    Element::Pointer Create( IndexType NewId,
                             GeometryType::Pointer pGeom,
                             PropertiesType::Pointer pProperties ) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "steady-state Pw flow Element #" << this->Id() << "\nRetention law: " << mRetentionLawVector[0]->Info();
        return buffer.str();
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "steady-state Pw flow Element #" << this->Id() << "\nRetention law: " << mRetentionLawVector[0]->Info();
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void CalculateAll( MatrixType &rLeftHandSideMatrix,
                       VectorType &rRightHandSideVector,
                       const ProcessInfo& CurrentProcessInfo,
                       const bool CalculateStiffnessMatrixFlag,
                       const bool CalculateResidualVectorFlag ) override;

    void CalculateAndAddLHS(MatrixType &rLeftHandSideMatrix, ElementVariables &rVariables) override;

    void CalculateAndAddRHS(VectorType &rRightHandSideVector, ElementVariables &rVariables, unsigned int GPoint) override;

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

    // Assignment operator.
    SteadyStatePwElement & operator=(SteadyStatePwElement const& rOther);

    // Copy constructor.
    SteadyStatePwElement(SteadyStatePwElement const& rOther);

}; // Class SteadyStatePwElement

} // namespace Kratos

#endif // KRATOS_GEO_STEADY_STATE_PW_ELEMENT_H_INCLUDED  defined
