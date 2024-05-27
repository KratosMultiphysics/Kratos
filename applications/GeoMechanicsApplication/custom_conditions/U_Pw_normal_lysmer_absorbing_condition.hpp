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
#pragma once

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
    
    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwLysmerAbsorbingCondition() : UPwFaceLoadCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPwLysmerAbsorbingCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPwLysmerAbsorbingCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPwFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    /**
    * @brief Gets displacement vector at the condional nodes
    * @param rValues displacement values
    * @param Step buffer index, 0 is current time step; 1 is previous timestep; 2 is timestep before; etc
    */
    void GetValuesVector(Vector& rValues, int Step) const override;

    /**
    * @brief Gets velocity vector at the condional nodes
    * @param rValues velocity values
    * @param Step buffer index, 0 is current time step; 1 is previous timestep; 2 is timestep before; etc
    */
    void GetFirstDerivativesVector(Vector& rValues, int Step) const override;

    /**
    * @brief Calculates the right hand side
    * @param rRightHandSideVector global right hand side vector
    * @param rCurrentProcessInfo Current process information
    */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
       const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * @brief Calculates LHS Damping part of absorbing boundary
    * @param rDampingMatrix Global damping matrix
    * @param rCurrentProcessInfo Current process information
    */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * @brief Calculates LHS and RHS stiffness part of absorbing boundary
    * @param rLhsMatrix Global left hand side matrix
    * @param rRightHandSideVector Global right hand side vector
    * @param rCurrentProcessInfo Current process information
    */
    void CalculateLocalSystem(MatrixType& rLhsMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    static constexpr SizeType N_DOF = TNumNodes * TDim;
    static constexpr SizeType CONDITION_SIZE = TNumNodes * TDim + TNumNodes;

    using ElementMatrixType = BoundedMatrix<double, N_DOF, N_DOF>;
    using DimensionMatrixType = BoundedMatrix<double, TDim, TDim>;

    struct NormalLysmerAbsorbingVariables
    {
        double rho; // density of soil mixture
        double Ec; // p wave modulus
        double G; // shear modulus
        double n; // porosity
        double vp; // p wave velocity
        double vs; // shear wave velocity
        double p_factor; // p wave relaxation factor
        double s_factor; // s wave relaxation factor
        double virtual_thickness;
        Vector EcNodes;
        Vector GNodes;
        Vector SaturationNodes;
        Vector rhoNodes;

        DimensionMatrixType CAbsMatrix; // damping part of absorbing matrix;
        DimensionMatrixType KAbsMatrix; // stiffness part of absorbing matrix;
    };

    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    /**
    * @brief Adds the condition matrix to the global left hand side
    * @param rLeftHandSideMatrix Global Left hand side
    * @param rUUMatrix LHS displacement matrix of the current condition
    */
    void AddLHS(MatrixType& rLeftHandSideMatrix, const ElementMatrixType& rUUMatrix);

    /**
    * @brief Calculates and adds terms to the RHS
    * @param rRigtHandSideVector Global Right hand side
    * @param rStiffnessMatrix condition stiffness matrix
    */
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, const MatrixType& rStiffnessMatrix);

    /**
    * @brief Calculates the rotation matrix of the current condition
    * @param rRotationMatrix rotation matrix of the current condition
    * @param rGeom geometry of the current condition
    */
    void CalculateRotationMatrix(DimensionMatrixType& rRotationMatrix, const Element::GeometryType& rGeom);

    /**
    * @brief This method gets the average of the variables of all the neighbour elements of the condition. 
    * @param rVariables Condition specific variables struct
    * @param rCurrentProcessInfo Current process information
    */
    void GetNeighbourElementVariables(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    /**
    * @brief Gets condition variables.
    * @param rVariables Condition specific variables struct
    * @param rCurrentProcessInfo Current process information
    */
    void GetVariables(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculates the damping constant in all directions
     * @param rVariables Condition specific variables struct
     * @param rGeom condition geometry
     */
    void CalculateNodalDampingMatrix(NormalLysmerAbsorbingVariables& rVariables, const Element::GeometryType& rGeom);

    /**
    * @brief Calculates the stiffness constant in all directions
    * @param rVariables Condition specific variables struct
    * @param rGeom condition geometry
    */
    void CalculateNodalStiffnessMatrix(NormalLysmerAbsorbingVariables& rVariables, const Element::GeometryType& rGeom);

    /**
    * @brief Calculates the extrapolation matrix for neighbour elements.Values from integration points are extrapolated to the nodes
    * @param rNeighbourElement The neighbouring element of the condition
    */
    Matrix CalculateExtrapolationMatrixNeighbour(const Element& rNeighbourElement);
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    using hashmap = std::unordered_multimap<DenseVector<int>, std::vector<Condition::Pointer>, KeyHasherRange<DenseVector<int>>, KeyComparorRange<DenseVector<int>>>;

    /**
    * @brief Calculates the rotation matrix of the current condition for 2D line conditions
    * @param rRotationMatrix rotation matrix of the current condition
    * @param rGeom geometry of the current condition
    */
    void CalculateRotationMatrix2DLine(DimensionMatrixType& rRotationMatrix, const Element::GeometryType& rGeom);

	/**
     * @brief Calculates the stiffness matrix for the current condition
     * @param rStiffnessMatrix The stiffness matrix to be calculated
     * @param rCurrentProcessInfo Current process information
     */
    void CalculateConditionStiffnessMatrix(ElementMatrixType& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo);

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
