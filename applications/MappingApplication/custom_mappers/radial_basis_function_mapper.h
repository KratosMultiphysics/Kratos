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
#include "custom_utilities/mapper_local_system.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos
{
///@name Kratos Classes
///@{

class RadialBasisFunctionMapperLocalSystem : public MapperLocalSystem
{
public:

    explicit RadialBasisFunctionMapperLocalSystem(){}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        // KRATOS_DEBUG_ERROR_IF_NOT(mpGeom) << "Members are not initialized!" << std::endl;
        // return mpGeom->Center(); // check why not compiling...
        KRATOS_ERROR << "not implemented, needs checking" << std::endl;
    }

    MapperLocalSystemUniquePointer Create(GeometryPointerType pGeometry) const override
    {
        return Kratos::make_unique<RadialBasisFunctionMapperLocalSystem>();
    }

    /// Turn back information as a string.
    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override {KRATOS_ERROR << "Not implemented!"<<std::endl;}

private:
    // GeometryPointerType mpGeom;
    // bool mIsProjection; // Set to true is we are projecting the master onto the slave.
    //                     // Set to false if we are projecting the slave onto the slave.
    // bool mIsDualMortar = false;
    // bool mIsDestinationIsSlave = true;

};


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
class RadialBasisFunctionMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of RadialBasisFunctionMapper
    KRATOS_CLASS_POINTER_DEFINITION(RadialBasisFunctionMapper);

    typedef Mapper<TSparseSpace, TDenseSpace> BaseType;

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;

    typedef InterfaceVectorContainer<TSparseSpace, TDenseSpace> InterfaceVectorContainerType;
    typedef Kratos::unique_ptr<InterfaceVectorContainerType> InterfaceVectorContainerPointerType;

    typedef std::size_t IndexType;

    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType MappingMatrixType;
    typedef Kratos::unique_ptr<MappingMatrixType> MappingMatrixUniquePointerType;

    typedef LinearSolver<TSparseSpace, TDenseSpace> LinearSolverType;
    typedef Kratos::shared_ptr<LinearSolverType> LinearSolverSharedPointerType;

    typedef typename TSparseSpace::VectorType TSystemVectorType;
    typedef Kratos::unique_ptr<TSystemVectorType> TSystemVectorUniquePointerType;

    using SparseMatrixType = typename TSparseSpace::MatrixType;
    using DenseMatrixType = typename TDenseSpace::MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    RadialBasisFunctionMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                        : mrModelPartOrigin(rModelPartOrigin),
                          mrModelPartDestination(rModelPartDestination){}


    RadialBasisFunctionMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters);

    /// Destructor.
    ~RadialBasisFunctionMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        // mpModeler->PrepareGeometryModel();

        // AssignInterfaceEquationIds();

        // KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE);
            MappingOptions.Set(MapperFlags::INTERNAL_USE_TRANSPOSE, true);
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
        else if (MappingOptions.Is(MapperFlags::INTERNAL_USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }
    }

    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE);
            MappingOptions.Set(MapperFlags::INTERNAL_USE_TRANSPOSE, true);
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
        else if (MappingOptions.Is(MapperFlags::INTERNAL_USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }
    }

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

    ///@}
    ///@name Access
    ///@{

    MappingMatrixType& GetMappingMatrix() override
    {
        if (mMapperSettings["precompute_mapping_matrix"].GetBool()) return *(mpMappingMatrix.get());
        else KRATOS_ERROR << "'precompute_mapping_matrix' must be 'true' in your parameters to retrieve the computed mapping matrix!" << std::endl;
    }

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        return Kratos::make_unique<RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);
    }

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
        return *mpCouplingInterfaceOrigin;
    }

    ModelPart& GetInterfaceModelPartDestination() override
    {
        return *mpCouplingInterfaceDestination;
    }

protected:
    ///@name Protected static Member Variables
    Parameters mLocalMapperSettings;

    /**
     * @brief This function origin model part
     * @return The origin model part
     */
    ModelPart& GetOriginModelPart()
    {
        return mrModelPartOrigin;
    }

    /**
     * @brief This function destination model part
     * @return The destination model part
     */
    ModelPart& GetDestinationModelPart()
    {
        return mrModelPartDestination;
    }

    /**
     * @brief This function return the interface vector container of the origin model part
     * @return The origin model part
    */
    InterfaceVectorContainerType* pGetInterfaceVectorContainerOrigin()
    {
        return mpInterfaceVectorContainerOrigin.get();
    }

    /**
     * @brief This function return the interface vector container of the destination model part
     * @return The destination model part
    */
    InterfaceVectorContainerType* pGetInterfaceVectorContainerDestination()
    {
        return mpInterfaceVectorContainerDestination.get();
    }

private:

    ///@name Private Operations
    ///@{
    // typename Modeler::Pointer mpModeler = nullptr;
    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;
    ModelPart* mpCouplingMP = nullptr;
    ModelPart* mpCouplingInterfaceOrigin = nullptr;
    ModelPart* mpCouplingInterfaceDestination = nullptr;

    Parameters mMapperSettings;

    MapperUniquePointerType mpInverseMapper = nullptr;

    MappingMatrixUniquePointerType mpMappingMatrix;
    // MappingMatrixUniquePointerType mpMappingMatrixProjector;
    // MappingMatrixUniquePointerType mpMappingMatrixSlave;
    MappingMatrixUniquePointerType mpOriginInterpolationMatrix;
    MappingMatrixUniquePointerType mpDestinationEvaluationMatrix;

    //TSystemVectorUniquePointerType mpTempVector;

    MapperLocalSystemPointerVector mMapperLocalSystemsOrigin;
    MapperLocalSystemPointerVector mMapperLocalSystemsDestination;

    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;

    LinearSolverSharedPointerType mpLinearSolver = nullptr;
    using MappingSparseSpaceType = typename MapperDefinitions::SparseSpaceType;
    using DenseSpaceType = typename MapperDefinitions::DenseSpaceType;
    using MappingMatrixUtilitiesType = MappingMatrixUtilities<MappingSparseSpaceType, DenseSpaceType>;

    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceDestination->GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceOrigin->GetCommunicator());
    }

    void MapInternal(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<double>& rOriginVariable,
                              const Variable<double>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    void MapInternal(const Variable<array_1d<double, 3>>& rOriginVariable,
                     const Variable<array_1d<double, 3>>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<array_1d<double, 3>>& rOriginVariable,
                              const Variable<array_1d<double, 3>>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    void CreateLinearSolver();

    //void CalculateMappingMatrixWithSolver(MappingMatrixType& rConsistentInterfaceMatrix, MappingMatrixType& rProjectedInterfaceMatrix);
    void FillCoordinatesMatrix(const ModelPart& ModelPart, const std::vector<Condition::Pointer>& IntegrationPointsPointerVector, DenseMatrixType& rCoordinatesMatrix, bool IsDomainIGA);

    // Get scaling factor stabilising numerics, maximum distance between spline support points in either X or Y direction
    double CalculateScaleFactor(DenseMatrixType& rOriginCoords);

    // This function calculates the number of polynomial terms from the degree
    IndexType CalculateNumberOfPolynomialTermsFromDegree(IndexType PolyDegree, bool project_origin_nodes_to_destination_domain);
    
    // Evaluate the polynomial required for the radial basis function interpolation
    std::vector<double> EvaluatePolynomialBasis(const array_1d<double, 3>& rCoords, unsigned int degree, bool project_origin_nodes_to_destination_domain) const;

    // Create and invert the coefficient matrix of the spline C. This depends only on the positions of the origin support points 
    void CreateAndInvertOriginRBFMatrix(DenseMatrixType& rInvCMatrix, const DenseMatrixType& rOriginCoords, bool project_origin_nodes_to_destination_domain,
        IndexType Poly_Degree, RBFShapeFunctionsUtility::RBFType RBF_Type, double Factor = 1.0, double eps = 1.0);

    // Compute Aij splining matrix, relating origin and destination interpolation points
    void CreateDestinationRBFMatrix(DenseMatrixType& rAMatrix, const DenseMatrixType& rOriginCoords, const DenseMatrixType& rDestinationCoords,
        bool project_origin_nodes_to_destination_domain, IndexType Poly_Degree, RBFShapeFunctionsUtility::RBFType RBF_Type, bool map_displacements_to_rotations, double rbf_shape_parameter);
    
    // For IGA, compute the mapping matrix mapping from origin control points to destination nodes
    std::unique_ptr<MappingMatrixType> ComputeMappingMatrixIga(const MappingMatrixType& rMappingMatrixGP, const std::vector<Condition::Pointer>& rOriginIntegrationPoints,
        const ModelPart& rOriginModelPart) const;


    Parameters GetMapperDefaultSettings() const 
    {
        return Parameters(R"({
            "echo_level"                    : 0,
            "radial_basis_function_type" : "thin_plate_spline",
            "additional_polynomial_degree": 0,
            "is_destination_slave"          : true,
            "is_origin_iga"             : false,
            "is_destination_iga"             : false,
            "precompute_mapping_matrix"      : true, 
            "destination_solver_settings": {
                "project_origin_nodes_to_destination_domain": false,
                "map_displacements_to_rotations": false
            }
        })");
    }

    ///@}
    ///@name Private  Access
    ///@{

    MapperUniquePointerType& GetInverseMapper()
    {
        if (!mpInverseMapper) {
            InitializeInverseMapper();
        }
        return mpInverseMapper;
    }

    void InitializeInverseMapper()
    {
        KRATOS_ERROR << "Inverse Mapping is not supported yet!" << std::endl;
        mpInverseMapper = this->Clone(mrModelPartDestination,
                                      mrModelPartOrigin,
                                      mMapperSettings);
    }

    ///@}

}; // Class RadialBasisFunctionMapper

///@} addtogroup block
}  // namespace Kratos.