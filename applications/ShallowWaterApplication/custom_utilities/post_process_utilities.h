//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_POST_PROCESS_UTILITIES_H_INCLUDED
#define KRATOS_POST_PROCESS_UTILITIES_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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

/** 
 * @ingroup ShallowWaterApplication
 * @class PostProcessUtilities
 * @brief This utilities manage the properties in order to allow GiD to print the wet and the dry domains in different layers
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) PostProcessUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::vector<IndexType> IndexVectorType;

    /// Pointer definition of PostProcessUtilities
    KRATOS_CLASS_POINTER_DEFINITION(PostProcessUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PostProcessUtilities(ModelPart& rThisModelPart) : mrModelPart(rThisModelPart) {}

    /// Destructor.
    ~PostProcessUtilities() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief For each existing property, an auxiliary property is created. The original properties mean wet state and the auxiliary properties mean dry state.
     */
    void DefineAuxiliaryProperties();

    /**
     * @brief This method assigns a property according to the flag. Gid will print the dry domain in another layer with an specific color.
     * @param The flag which means dry state
     */
    void AssignDryWetProperties(Flags WetStateFlag);

    /**
     * @brief Restore the properties to the original state (wet state)
     */
    void RestoreDryWetProperties();

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
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    std::map<IndexType, IndexType> mWetToDryPropertiesMap;
    std::map<IndexType, IndexType> mDryToWetPropertiesMap;

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

    /// Assignment operator.
    PostProcessUtilities& operator=(PostProcessUtilities const& rOther) = delete;

    /// Copy constructor.
    PostProcessUtilities(PostProcessUtilities const& rOther) = delete;


    ///@}

}; // Class PostProcessUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                PostProcessUtilities& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const PostProcessUtilities& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_POST_PROCESS_UTILITIES_H_INCLUDED  defined
