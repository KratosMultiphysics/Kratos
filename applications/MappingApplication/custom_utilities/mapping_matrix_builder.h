
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


#if !defined(KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED )
#define  KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "containers/system_vector.h"
#include "containers/distributed_system_vector.h"
#include "containers/csr_matrix.h"
#include "containers/distributed_csr_matrix.h"
#include "custom_utilities/mapper_local_system.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template<bool IsDistributed>
class KRATOS_API(MAPPING_APPLICATION) MappingMatrixBuilder
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MappingMatrixBuilder
    KRATOS_CLASS_POINTER_DEFINITION(MappingMatrixBuilder);

    using GraphType = typename std::conditional<IsDistributed, DistributedSparseGraph<std::size_t>, SparseGraph<>>::type;

    using MappingMatrixType = typename std::conditional<IsDistributed,
        DistributedCsrMatrix<>,
        CsrMatrix<>>::type;
    using MappingMatrixPointerType = Kratos::unique_ptr<MappingMatrixType>;

    using InterfaceVectorType = typename std::conditional<IsDistributed,
        DistributedSystemVector<>,
        SystemVector<>>::type;
    using InterfaceVectorPointerType = Kratos::unique_ptr<InterfaceVectorType>;

    using MapperLocalSystemPointerVector = std::vector<Kratos::unique_ptr<MapperLocalSystem>>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MappingMatrixBuilder(const int EchoLevel)
        : mEchoLevel(EchoLevel) { }

    /// Destructor.
    virtual ~MappingMatrixBuilder() = default;

    /// Copy constructor.
    MappingMatrixBuilder(MappingMatrixBuilder const& rOther) = delete;

    /// Assignment operator.
    MappingMatrixBuilder& operator=(MappingMatrixBuilder const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void BuildMappingMatrix(
        MapperLocalSystemPointerVector& rLocalSystems,
        MappingMatrixPointerType& rpMappingMatrix,
        InterfaceVectorPointerType& rpInterfaceVectorOrigin,
        InterfaceVectorPointerType& rpInterfaceVectorDestination
        );

    ///@}

private:
    ///@name Member Variables
    ///@{

    int mEchoLevel = 0;

    ///@}
    ///@name Private Operations
    ///@{


    ///@}

}; // Class MappingMatrixBuilder

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED  defined
