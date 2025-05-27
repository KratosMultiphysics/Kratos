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

#include "linear_solvers/linear_solver.h"

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

    typedef InterpolativeMapperBase<TSparseSpace, TDenseSpace, MapperBackend<TSparseSpace, TDenseSpace>> BaseType;

    typedef typename BaseType::MapperUniquePointerType MapperUniquePointerType;
    typedef typename BaseType::TMappingMatrixType TMappingMatrixType;
    typedef Kratos::unique_ptr<TMappingMatrixType> MappingMatrixUniquePointerType;

    typedef LinearSolver<TSparseSpace, TDenseSpace> LinearSolverType;
    typedef Kratos::shared_ptr<LinearSolverType> LinearSolverSharedPointerType;
    
    typedef typename TSparseSpace::MatrixType SparseMatrixType;
    typedef typename TDenseSpace::MatrixType DenseMatrixType;

    typedef typename MapperDefinitions::SparseSpaceType MappingSparseSpaceType;
    typedef typename MapperDefinitions::DenseSpaceType  DenseSpaceType;

    typedef MappingMatrixUtilities<MappingSparseSpaceType, DenseSpaceType> MappingMatrixUtilitiesType;

    typedef typename BaseType::MapperInterfaceInfoUniquePointerType MapperInterfaceInfoUniquePointerType;

    
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
    // ModelPart& mrModelPartOrigin;
    // ModelPart& mrModelPartDestination;
    ModelPart* mpCouplingMP = nullptr;
    ModelPart* mpCouplingInterfaceMaster = nullptr;
    ModelPart* mpCouplingInterfaceSlave = nullptr;

    Parameters mMapperSettings;

    MappingMatrixUniquePointerType mpOriginInterpolationMatrix;
    MappingMatrixUniquePointerType mpDestinationEvaluationMatrix;

    LinearSolverSharedPointerType mpLinearSolver = nullptr;


    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceSlave->GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceMaster->GetCommunicator());
    }

    void CreateLinearSolver();

    std::vector<double> EvaluatePolynomialBasis(const array_1d<double, 3>& coords, unsigned int degree) const;

    Parameters GetMapperDefaultSettings() const override
    {
        return Parameters(R"({
            "echo_level"                    : 0,
            "radial_basis_function_type" : "thin_plate_spline",
            "additional_polynomial_degree": 0,
            "destination_is_slave"          : true,
            "is_origin_iga"             : false,
            "is_destination_iga"             : false,
             "linear_solver_settings": {
                "solver_type": "skyline_lu_factorization"
            }
        })");
    }

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems) override
    {
    }

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const override
    {
        return nullptr;
    }


}; // Class RadialBasisFunctionMapper

///@} addtogroup block
}  // namespace Kratos.