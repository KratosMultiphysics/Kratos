// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohamed Nabi
//
//
//


#if !defined(KRATOS_GEO_T_NORMAL_FLUX_CONDITION_H_INCLUDED )
#define  KRATOS_GEO_T_NORMAL_FLUX_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/T_condition.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TNormalFluxCondition : public TCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TNormalFluxCondition);
    
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using TCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    TNormalFluxCondition() : TCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    TNormalFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : TCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    TNormalFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : TCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~TNormalFluxCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const;
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct NormalFluxVariables
    {
        double NormalFlux;
        double IntegrationCoefficient;
        array_1d<double,TNumNodes> Np;
        array_1d<double,TNumNodes> TVector;
    };
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                    
    void CalculateRHS(VectorType& rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo);
    
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, NormalFluxVariables& rVariables);

    virtual void CalculateIntegrationCoefficient(double& rIntegrationCoefficient,
        const Matrix& Jacobian,
        const double& Weight);
    
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
    
}; // class TNormalFluxCondition.

} // namespace Kratos.

#endif // KRATOS_GEO_T_NORMAL_FLUX_CONDITION_H_INCLUDED defined 
