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

#if !defined(KRATOS_GEO_U_PW_UPDATED_LAGRANGIAN_AXISYMMETRIC_FIC_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_U_PW_UPDATED_LAGRANGIAN_AXISYMMETRIC_FIC_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/stress_strain_utilities.hpp"
#include "custom_elements/U_Pw_updated_lagrangian_FIC_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwUpdatedLagrangianAxisymmetricFICElement :
    public UPwUpdatedLagrangianFICElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPwUpdatedLagrangianAxisymmetricFICElement );

    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;
    /// The definition of the sizetype
    using SizeType = std::size_t;
    using UPwBaseElement<TDim,TNumNodes>::mConstitutiveLawVector;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwUpdatedLagrangianAxisymmetricFICElement(IndexType NewId = 0) :
        UPwUpdatedLagrangianFICElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    UPwUpdatedLagrangianAxisymmetricFICElement(IndexType NewId,
                                         const NodesArrayType& ThisNodes) :
        UPwUpdatedLagrangianFICElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPwUpdatedLagrangianAxisymmetricFICElement(IndexType NewId,
                                         GeometryType::Pointer pGeometry) :
        UPwUpdatedLagrangianFICElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPwUpdatedLagrangianAxisymmetricFICElement(IndexType NewId,
                                         GeometryType::Pointer pGeometry,
                                         PropertiesType::Pointer pProperties) :
        UPwUpdatedLagrangianFICElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~UPwUpdatedLagrangianAxisymmetricFICElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create( IndexType NewId,
                             NodesArrayType const& ThisNodes,
                             PropertiesType::Pointer pProperties ) const override;

    Element::Pointer Create( IndexType NewId,
                             GeometryType::Pointer pGeom,
                             PropertiesType::Pointer pProperties ) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Small strain axisymmetric U-Pw Element #"
               << this->Id()
               << "\nConstitutive law: "
               << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Small strain axisymmetric U-Pw Element #"
                 << this->Id()
                 << "\nConstitutive law: "
                 << mConstitutiveLawVector[0]->Info();
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateBMatrix( Matrix &rB,
                           const Matrix &GradNpT,
                           const Vector &Np ) override;

    double CalculateIntegrationCoefficient( const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
                                            unsigned int PointNumber,
                                            double detJ) override;

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
    UPwUpdatedLagrangianAxisymmetricFICElement & operator=(UPwUpdatedLagrangianAxisymmetricFICElement const& rOther);

    // Copy constructor.
    UPwUpdatedLagrangianAxisymmetricFICElement(UPwUpdatedLagrangianAxisymmetricFICElement const& rOther);

    // Private Operations

    template < class TValueType >
    inline void ThreadSafeNodeWrite(NodeType& rNode, const Variable<TValueType> &Var, const TValueType Value)
    {
        rNode.SetLock();
        rNode.FastGetSolutionStepValue(Var) = Value;
        rNode.UnSetLock();
    }


}; // Class UPwUpdatedLagrangianAxisymmetricFICElement

} // namespace Kratos

#endif // KRATOS_GEO_U_PW_UPDATED_LAGRANGIAN_AXISYMMETRIC_FIC_ELEMENT_H_INCLUDED  defined
