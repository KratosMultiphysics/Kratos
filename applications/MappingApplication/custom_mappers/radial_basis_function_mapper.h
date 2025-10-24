//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//

#pragma once

// System includes

// External includes

// Project includes
#include "utilities/rbf_shape_functions_utility.h"
#include "mappers/mapper_define.h"
#include "interpolative_mapper_base.h"
#include "custom_utilities/mapper_backend.h"


namespace Kratos
{
///@name Kratos Classes
///@{


// RadialBasisFunctionMapper
//
// The mapper always forward maps from the master to the slave.
// Normally:
//      master  =   interface origin
//      slave   =   interface destination
//
// However, this can be reversed by setting 'destination_is_slave' = false.
// This yields:
//      master  =   interface destination
//      slave   =   interface origin

template<class TSparseSpace, class TDenseSpace>
class RadialBasisFunctionMapper : public InterpolativeMapperBase<TSparseSpace, TDenseSpace, MapperBackend<TSparseSpace, TDenseSpace>>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of RadialBasisFunctionMapper
    KRATOS_CLASS_POINTER_DEFINITION(RadialBasisFunctionMapper);

    using BaseType = InterpolativeMapperBase<TSparseSpace, TDenseSpace, MapperBackend<TSparseSpace, TDenseSpace>>;

    using MapperUniquePointerType = typename BaseType::MapperUniquePointerType;
    using TMappingMatrixType = typename BaseType::TMappingMatrixType;
    using MappingMatrixUniquePointerType = Kratos::unique_ptr<TMappingMatrixType>;

    using SparseMatrixType = typename TSparseSpace::MatrixType;
    using DenseMatrixType = typename TDenseSpace::MatrixType;

    using MappingSparseSpaceType = typename MapperDefinitions::SparseSpaceType;
    using DenseSpaceType = typename MapperDefinitions::DenseSpaceType;

    using MappingMatrixUtilitiesType = MappingMatrixUtilities<MappingSparseSpaceType, DenseSpaceType>;

    using MapperInterfaceInfoUniquePointerType = typename BaseType::MapperInterfaceInfoUniquePointerType;

    
    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    RadialBasisFunctionMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                        : BaseType(rModelPartOrigin, rModelPartDestination){}


    RadialBasisFunctionMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters);

    /// Destructor.
    ~RadialBasisFunctionMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters) const override
    {
        KRATOS_TRY;

        return Kratos::make_unique<RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>>(
        rModelPartOrigin,
        rModelPartDestination,
        JsonParameters);

        KRATOS_CATCH("");
    }

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        AssignInterfaceEquationIds();

        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

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
        return "RadialBasisFunctionMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RadialBasisFunctionMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

    // Get values
    // Always returns the true origin/destination regardless of 'destination_is_slave'
    ModelPart& GetInterfaceModelPartOrigin() override
    {
        return *mpCouplingInterfaceMaster;
    }

    ModelPart& GetInterfaceModelPartDestination() override
    {
        return *mpCouplingInterfaceSlave;
    }

protected:
    ///@name Protected static Member Variables
    Parameters mLocalMapperSettings;

private:

    ///@name Private Operations
    ///@{
    ModelPart* mpCouplingMP = nullptr;
    ModelPart* mpCouplingInterfaceMaster = nullptr;
    ModelPart* mpCouplingInterfaceSlave = nullptr;

    Parameters mMapperSettings;

    MappingMatrixUniquePointerType mpOriginInterpolationMatrix;
    MappingMatrixUniquePointerType mpDestinationEvaluationMatrix;

    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceSlave->GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceMaster->GetCommunicator());
    }

    void FillCoordinatesMatrix(const ModelPart& ModelPart, const std::vector<Condition::Pointer>& IntegrationPointsPointerVector, DenseMatrixType& rCoordinatesMatrix, bool IsDomainIGA);

    // Get scaling factor stabilising numerics, maximum distance between spline support points in either X or Y direction
    double CalculateScaleFactor(DenseMatrixType& rOriginCoords);

    // This function calculates the number of polynomial terms from the degree
    IndexType CalculateNumberOfPolynomialTermsFromDegree(IndexType PolyDegree, bool ProjectToAerodynamicPanels);
    
    // Evaluate the polynomial required for the radial basis function interpolation
    std::vector<double> EvaluatePolynomialBasis(const array_1d<double, 3>& rCoords, unsigned int degree, bool ProjectToPanelsPlane) const;

    // Create and invert the coefficient matrix of the spline C. This depends only on the positions of the origin support points 
    void CreateAndInvertOriginRBFMatrix(DenseMatrixType& rInvCMatrix, const DenseMatrixType& rOriginCoords, bool ProjectToAerodynamicPanels,
        IndexType Poly_Degree, const std::string& RBFType, double Factor = 1.0, double eps = 1.0);

    // Compute Aij splining matrix, relating origin and destination interpolation points
    void CreateDestinationRBFMatrix(DenseMatrixType& rAMatrix, const DenseMatrixType& rOriginCoords, const DenseMatrixType& rDestinationCoords,
        bool ProjectToAerodynamicPanels, IndexType Poly_Degree, const std::string& RBFType, bool ReturnAOAMatrix, double rbf_shape_parameter);
    
    // For IGA, compute the mapping matrix mapping from origin control points to destination nodes
    std::unique_ptr<TMappingMatrixType> ComputeMappingMatrixIga(const TMappingMatrixType& rMappingMatrixGP, const std::vector<Condition::Pointer>& rOriginIntegrationPoints,
        const ModelPart& rOriginModelPart) const;
        
    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters(R"({
            "echo_level"                    : 0,
            "radial_basis_function_type" : "thin_plate_spline",
            "additional_polynomial_degree": 0,
            "destination_is_slave"          : true,
            "is_origin_iga"             : false,
            "is_destination_iga"             : false,
            "aerodynamic_panel_solver_settings": {
                "project_origin_nodes_to_destination_domain_panel_solver": false,
                "map_structural_displacements_to_panels_angles_of_attack": false
            }
        })");
    }

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems) override
    {
        KRATOS_ERROR_IF(false) << "The RBF mapper does not use local systems!" << std::endl;
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const noexcept override
    {
        return nullptr;
    }


}; // Class RadialBasisFunctionMapper

///@} addtogroup block
}  // namespace Kratos.