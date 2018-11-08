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

#if !defined(KRATOS_MAPPER_H_INCLUDED)
#define  KRATOS_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "containers/flags.h"


namespace Kratos
{

///@}
///@name Kratos Classes
///@{

/// Base Class for all Mappers
/** This is the base class for every mapper. This is the equivalent to a Kratos-SolvingStrategy.
 * It provides the basic interface for the mapping operations
*/
template<class TSparseSpace, class TDenseSpace>
class Mapper
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of Mapper
    KRATOS_CLASS_POINTER_DEFINITION(Mapper);

    typedef typename TSparseSpace::MatrixType  TMappingMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    Mapper() = default;

    /// Destructor.
    virtual ~Mapper() = default;

    /// Copy Constructor
    Mapper(const Mapper&) = delete;

    /// Move Constructor
    Mapper(Mapper&&) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Copy Assignment Operator
    Mapper& operator=(const Mapper&) = delete;

    /// Move Assignment Operator
    Mapper& operator=(Mapper&&) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Updates the mapping-system after the geometry/mesh has changed
    * After changes in the topology (e.g. remeshing or sliding interfaces)
    * the relations for the mapping have to be recomputed. This means that
    * the search has to be conducted again and the mapping-system has to be
    * rebuilt, hence this is expensive
    * @param MappingOptions flags used to specify how the update has to be done
    * @param SearchRadius search radius used for the search
    */
    virtual void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) = 0;

    /**
    * @brief Mapping from Origin to Destination, Scalar Variable
    * Data is exchanged on the Interface, from the Origin-Modelpart
    * to the Destination-ModelPart (the modelparts were specified in the
    * construction Phase of the Mapper)
    * @param rOriginVariable Variable on the Origin-ModelPart
    * @param rDestinationVariable Variable on the Destination-ModelPart
    * @param MappingOptions flags used to specify options for the mapping
    * @see InverseMap
    */
    virtual void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) = 0;

    /**
    * @brief Mapping from Origin to Destination, Vector Variable
    * Same as Map, but maps an array3-variable
    * @see Map
    */
    virtual void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) = 0;

    /**
    * @brief Mapping from Destination to Origin, Scalar Variable
    * Data is exchanged on the Interface, from the Destination-Modelpart
    * to the Origin-ModelPart (the modelparts were specified in the
    * construction Phase of the Mapper)
    * It does the opposite of Map
    * @param rOriginVariable Variable on the Origin-ModelPart
    * @param rDestinationVariable Variable on the Destination-ModelPart
    * @param MappingOptions flags used to specify options for the mapping
    * @see Map
    */
    virtual void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) = 0;

    /**
    * @brief Mapping from Destination to Origin, Vector Variable
    * Same as InveseMap, but maps an array3-variable
    * @see InverseMap
    */
    virtual void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) = 0;

    /**
    * @brief Cloning the Mapper
    * returns a clone of the current Mapper
    * pure virtual, has to be implemented in every derived mapper,
    * used in the creation of the Mappers
    * @see MapperFactory
    */
    virtual MapperUniquePointerType Clone(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters) const = 0;

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the mapping-matrix
     * @return The mapping-matrix
     */
    const TMappingMatrixType& GetMappingMatrix() const = 0;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Mapper";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Mapper";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}

}; // Class Mapper

/// output stream function
template<class TSparseSpace, class TDenseSpace>
inline std::ostream & operator << (
    std::ostream& rOStream,
    const Mapper<TSparseSpace, TDenseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_MAPPER_H_INCLUDED  defined