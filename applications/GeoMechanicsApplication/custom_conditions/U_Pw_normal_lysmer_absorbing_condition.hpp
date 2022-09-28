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
//  Main authors:    Aron Noordam
//


#if !defined(KRATOS_GEO_U_PW_NORMAL_LYSMER_ABSORBING_CONDITION_H_INCLUDED )
#define  KRATOS_GEO_U_PW_NORMAL_LYSMER_ABSORBING_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_conditions/U_Pw_face_load_condition.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwLysmerAbsorbingCondition : public UPwFaceLoadCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwLysmerAbsorbingCondition);
    
    typedef std::size_t IndexType;
    typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwLysmerAbsorbingCondition() : UPwFaceLoadCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPwLysmerAbsorbingCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPwLysmerAbsorbingCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPwFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPwLysmerAbsorbingCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

   void GetFirstDerivativesVector(Vector& rValues, int Step) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    static constexpr SizeType N_DOF = TNumNodes * TDim;
    static constexpr SizeType CONDITION_SIZE = TNumNodes * TDim + TNumNodes;

    struct NormalLysmerAbsorbingVariables
    {
        double IntegrationCoefficient;
        array_1d<double,TNumNodes> Np;
        BoundedMatrix<double, N_DOF, N_DOF> UMatrix;
        double rho; // density of soil mixture
        double Ec; // p wave modulus
        double G; // shear modulus
        double n; // porosity
        double vp; // p wave velocity
        double vs; // shear wave velocity
        double p_factor; // p wave relaxation factor
        double s_factor; // s wave relaxation factor
        double virtual_thickness;
        vector<double> EcNodes;
        vector<double> GNodes;
        vector<double> SaturationNodes;
        vector<double> rhoNodes;

        BoundedMatrix<double, TDim, TDim> CAbsMatrix; // damping part of absorbing matrix;
        BoundedMatrix<double, TDim, TDim> KAbsMatrix; // stiffness part of absorbing matrix;
    };

    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& CurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLhsMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo) override;
    
    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, NormalLysmerAbsorbingVariables& rVariables);

    void CalculateRotationMatrix(BoundedMatrix<double, TDim, TDim>& rRotationMatrix, const Element::GeometryType& Geom);

    void GetNeighbourElementVariables(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo);

    void GetVariables(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo);

    void CalculateNodalDampingMatrix(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo, const Element::GeometryType& Geom);

    void CalculateNodalStiffnessMatrix(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo, const Element::GeometryType& Geom);

    Matrix CalculateExtrapolationMatrixNeighbour(const Element& NeighbourElement);
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    typedef std::unordered_multimap<DenseVector<int>, std::vector<Condition::Pointer>, KeyHasherRange<DenseVector<int>>, KeyComparorRange<DenseVector<int>> > hashmap;

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
    
}; // class UPwLysmerAbsorbingCondition.

} // namespace Kratos.

#endif // KRATOS_GEO_U_PW_NORMAL_LYSMER_ABSORBING_CONDITION_H_INCLUDED defined 
