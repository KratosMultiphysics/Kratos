//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_U_PW_NORMAL_LIQUID_FLUX_FIC_CONDITION_H_INCLUDED )
#define  KRATOS_U_PW_NORMAL_LIQUID_FLUX_FIC_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_conditions/U_Pw_normal_liquid_flux_condition.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwNormalLiquidFluxFICCondition : public UPwNormalLiquidFluxCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( UPwNormalLiquidFluxFICCondition );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    typedef typename UPwNormalLiquidFluxCondition<TDim,TNumNodes>::NormalLiquidFluxVariables NormalLiquidFluxVariables;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwNormalLiquidFluxFICCondition() : UPwNormalLiquidFluxCondition<TDim,TNumNodes>() {}

    // Constructor 1
    UPwNormalLiquidFluxFICCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwNormalLiquidFluxCondition<TDim,TNumNodes>(NewId, pGeometry) {}

    // Constructor 2
    UPwNormalLiquidFluxFICCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPwNormalLiquidFluxCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod = this->GetIntegrationMethod();
    }

    // Destructor
    ~UPwNormalLiquidFluxFICCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct NormalLiquidFluxFICVariables
    {
        double DtPressureCoefficient;
        double ElementLength;
        double BiotModulusInverse;

        array_1d<double,TNumNodes> DtPressureVector;
        BoundedMatrix<double,TNumNodes,TNumNodes> PMatrix;
    };

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateElementLength(double& rElementLength, const GeometryType& Geom);


    void CalculateAndAddLHSStabilization(MatrixType& rLeftHandSideMatrix, NormalLiquidFluxVariables& rVariables, NormalLiquidFluxFICVariables& rFICVariables);

    void CalculateAndAddBoundaryMassMatrix(MatrixType& rLeftHandSideMatrix, NormalLiquidFluxVariables& rVariables, NormalLiquidFluxFICVariables& rFICVariables);


    void CalculateAndAddRHSStabilization(VectorType& rRightHandSideVector, NormalLiquidFluxVariables& rVariables, NormalLiquidFluxFICVariables& rFICVariables);

    void CalculateAndAddBoundaryMassFlow(VectorType& rRightHandSideVector, NormalLiquidFluxVariables& rVariables, NormalLiquidFluxFICVariables& rFICVariables);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }

}; // class UPwNormalLiquidFluxFICCondition.

} // namespace Kratos.

#endif // KRATOS_U_PW_NORMAL_LIQUID_FLUX_FIC_CONDITION_H_INCLUDED defined
