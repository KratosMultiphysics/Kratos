//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPING_MATRIX_BUILDER_H)
#define  KRATOS_MAPPING_MATRIX_BUILDER_H

// System includes

// External includes

// Project includes
#include "custom_utilities/interface_vector_container.h"
#include "custom_utilities/mapper_local_system.h"
#include "custom_utilities/mapper_flags.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Object for constructing the mapping-system, equivalent to the BuilderAndSolver
/** The mapping-system is constructed from the the MapperLocalSystems, this class
 * and its derived objects take care of this.
*/
template<class TSparseSpace, class TDenseSpace>
class MappingMatrixBuilder
{
    public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MappingMatrixBuilder
    KRATOS_CLASS_POINTER_DEFINITION(MappingMatrixBuilder);

    typedef typename TSparseSpace::MatrixType TMappingMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef Kratos::unique_ptr<TMappingMatrixType> TMappingMatrixUniquePointerType;
    typedef Kratos::unique_ptr<TSystemVectorType> TSystemVectorUniquePointerType;

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;
    typedef Kratos::shared_ptr<MapperLocalSystemPointerVector> MapperLocalSystemPointerVectorPointer;

    typedef InterfaceVectorContainer<TSparseSpace, TDenseSpace> InterfaceVectorContainerType;
    typedef Kratos::unique_ptr<InterfaceVectorContainerType> InterfaceVectorContainerPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MappingMatrixBuilder();

    /// Destructor.
    virtual ~MappingMatrixBuilder() = default;


    ///@}
    ///@name Operations
    ///@{

    void BuildMappingMatrix(
        const InterfaceVectorContainerPointerType& rpVectorContainerOrigin,
        const InterfaceVectorContainerPointerType& rpVectorContainerDestination,
        MapperLocalSystemPointerVector& rMapperLocalSystems);

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the mapping-matrix
     * @return The mapping-matrix
     */
    TMappingMatrixType& GetMappingMatrix()
    {
        return *mpMappingMatrix;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "MappingMatrixBuilder";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

private:
    ///@name Member Variables
    ///@{

    int mEchoLevel = 0;
    TMappingMatrixUniquePointerType mpMappingMatrix;

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    // MappingMatrixBuilder& operator=(MappingMatrixBuilder const& rOther) {}

    /// Copy constructor.
    // MappingMatrixBuilder(MappingMatrixBuilder const& rOther) {}

    ///@}

}; // Class MappingMatrixBuilder

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPING_MATRIX_BUILDER_H  defined
