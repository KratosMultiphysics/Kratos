//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors: Aditya Ghantasala
//

#if !defined(KRATOS_GATHER_MODELPART_ON_ALL_RANKS_UTILITY)
#define KRATOS_GATHER_MODELPART_ON_ALL_RANKS_UTILITY

// System includes
#include <unordered_map>
#include <vector>


/* Project includes */
#include "includes/define.h"
#include "mpi/utilities/gather_modelpart_utility.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
class KRATOS_API(CHIMERA_APPLICATION) GatherModelPartOnAllRanksUtility
{
public:
    typedef std::size_t IndexType;
    ///@name Type Definitions
    ///@{

    // Needed structures for the ExtractSurfaceMesh operation
    struct KeyComparator
    {
        bool operator()(const vector<IndexType> &lhs,
                        const vector<IndexType> &rhs) const
        {
            if (lhs.size() != rhs.size())
                return false;
            for (IndexType i = 0; i < lhs.size(); i++)
                if (lhs[i] != rhs[i])
                    return false;
            return true;
        }
    };

    struct KeyHasher
    {
        IndexType operator()(const vector<int> &k) const
        {
            IndexType seed = 0.0;
            std::hash<int> hasher;
            for (IndexType i = 0; i < k.size(); i++)
                seed ^= hasher(k[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    // Some type-definitions
    typedef std::unordered_map<vector<IndexType>, IndexType, KeyHasher,
                               KeyComparator>
        hashmap;
    typedef std::unordered_map<vector<IndexType>, vector<IndexType>, KeyHasher,
                               KeyComparator>
        hashmap_vec;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ChimeraHoleCuttingUtility
    KRATOS_CLASS_POINTER_DEFINITION(GatherModelPartOnAllRanksUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    GatherModelPartOnAllRanksUtility() = default;

    /// Destructor.
    ~GatherModelPartOnAllRanksUtility() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Gathers the given modelpart on all the ranks
     * @param rModelPartToGather The modelpart which is to be gathered on all nodes
     * @param rGatheredModelPart The full gathered modelpart from all the ranks.
     */
    static void GatherModelPartOnAllRanks(ModelPart &rModelPartToGather, ModelPart &rGatheredModelPart);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
};

} // namespace Kratos

#endif // KRATOS_GATHER_MODELPART_ON_ALL_RANKS_UTILITY  defined