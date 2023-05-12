//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
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
 * @class ParticlesUtilities
 * @brief Collection of utilities to compute statistic of particles
 * @ingroup KratosCore
 * @author Riccardo Rossi
 */
class ParticlesUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParticlesUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ParticlesUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParticlesUtilities() = delete;

    /**
     * @brief: this function takes all the nodes in rParticlesModelPart and finds in which element they fall within the volume identified in rVolumeModelPart. Then they add one to all the nodes of the element in which the particle is inside
     * @param rLocator this is a search structure for rVolumeModelPart
     * @param rVolumeModelPart: model part on the top of which we will find the particles
     * @param rParticlesModelPart: model part in which the "particles" are contained (as Nodes)
     * @param rCounterVariable: this is the variable which will be used to store the count (has to be present in the nodes of rVolumeModelPart)
     * @param SearchTolerance: search tolerance used in the spatial search
     */
    template<unsigned int TDim, bool CounterHasHistory=false >
    static void CountParticlesInNodes(
        BinBasedFastPointLocator<TDim>& rLocator,
        ModelPart& rVolumeModelPart,
        const ModelPart& rParticlesModelPart,
        const Variable<double>& rCounterVariable,
        const double SearchTolerance=1e-5
    )
    {
        //reset the counter
        block_for_each(rVolumeModelPart.Nodes(), [&rCounterVariable](auto& rNode)
        {
            if constexpr (CounterHasHistory)
                rNode.FastGetSolutionStepValue(rCounterVariable) = 0.0;
            else
                rNode.SetValue(rCounterVariable,0.0);
        });


        unsigned int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType TLS(max_results);

        //for every interface node (nodes in cut elements)
        block_for_each(rParticlesModelPart.Nodes(), TLS, [&rLocator, &rCounterVariable, SearchTolerance](const auto& rNode, auto& rTLS)
        {

            Vector shape_functions;
            Element::Pointer p_element;
            const bool is_found = rLocator.FindPointOnMesh(rNode.Coordinates(), shape_functions, p_element, rTLS.begin(), rTLS.size(), SearchTolerance);

            if(is_found)
            {

                auto& r_geom = p_element->GetGeometry();
                for(unsigned int i=0; i<r_geom.size(); ++i)
                {
                    if constexpr (CounterHasHistory)
                    {
                        auto& rcounter = r_geom[i].FastGetSolutionStepValue(rCounterVariable);
                        AtomicAdd(rcounter, 1.0);
                    }
                    else
                    {
                        auto& rcounter = r_geom[i].GetValue(rCounterVariable);
                        AtomicAdd(rcounter, 1.0);
                    }

                }
            }

        });

    }

    /**
     * @brief: this function takes all the nodes in rParticlesModelPart and finds in which element they fall within the volume identified in rVolumeModelPart. Particles are then classified in the "ClassificationVariable" of type vector, adding 1.0 to the entry which corresponds to the "type" of particle, intended as an integer between 0 and NumberOfTypes identifying the particle
     * @param rLocator this is a search structure for rVolumeModelPart
     * @param rVolumeModelPart: model part on the top of which we will find the particles
     * @param rParticlesModelPart: model part in which the "particles" are contained (as Nodes)
     * @param NumberOfTypes number of admissible types of particles (note that the particle type will be interpreted as integer)
     * @param rParticleTypeVariable: this is an integer identifying the "type" of particle. Numbers lesser than 0 or higher than NumberofTypes will be silently ignored
     * @param SearchTolerance: search tolerance used in the spatial search
     */
    template<unsigned int TDim, class TScalarType, bool ParticleTypeVariableHasHistory=false>
    static void ClassifyParticlesInElements(
        BinBasedFastPointLocator<TDim>& rLocator,
        ModelPart& rVolumeModelPart,
        const ModelPart& rParticlesModelPart,
        const int NumberOfTypes,
        const Variable<TScalarType>& rParticleTypeVariable=AUX_INDEX,
        const Variable<Vector>& rClassificationVectorVariable=MARKER_LABELS,
        const double SearchTolerance=1e-5
    )
    {
        //reset the counter
        Vector zero = ZeroVector(NumberOfTypes);
        block_for_each(rVolumeModelPart.Elements(), [&rClassificationVectorVariable, &zero](auto& rElement){
            rElement.SetValue(rClassificationVectorVariable,zero);
        });


        unsigned int max_results = 10000;
        auto TLS = std::make_pair(typename BinBasedFastPointLocator<TDim>::ResultContainerType(max_results), Vector());
        //typename BinBasedFastPointLocator<TDim>::ResultContainerType TLS(max_results);

        //for every interface node (nodes in cut elements)
        block_for_each(rParticlesModelPart.Nodes(),
                TLS,
                [&rLocator, &rParticleTypeVariable, &rClassificationVectorVariable, &NumberOfTypes, SearchTolerance]
                (const auto& rNode, auto& rTLS)
                {
                    auto& results = rTLS.first;
                    Vector& shape_functions = rTLS.second;
                    Element::Pointer p_element;
                    const bool is_found = rLocator.FindPointOnMesh(rNode.Coordinates(), shape_functions, p_element, results.begin(), results.size(), SearchTolerance);

                    if(is_found)
                    {
                        int particle_type;
                        if constexpr (ParticleTypeVariableHasHistory)
                            particle_type = static_cast<int>(rNode.FastGetSolutionStepValue(rParticleTypeVariable));
                        else
                            particle_type = static_cast<int>(rNode.GetValue(rParticleTypeVariable));

                        if(particle_type>=0 && particle_type<NumberOfTypes) //we will ignore particles identified by a marker <0 or >NumberOfTypes
                        {
                            auto& rclassification = p_element->GetValue(rClassificationVectorVariable);
                            AtomicAdd(rclassification[particle_type], 1.0);
                        }
                    }
                });
    }


    /**
     * @brief: this function looks if particle is found in the locator. If it is not, "rVariable" is marked with the value "OutsiderValue"
     * @param rLocator this is a search structure for the volume
     * @param rParticlesModelPart: model part in which the "particles" are contained (as Nodes)
     * @param rVariable: variable whose value will be modified on the nodes in rParticlesModelPart
     * @param SearchTolerance: search tolerance used in the spatial search
     */
    template<unsigned int TDim, class TDataType, bool VariableHasHistory >
    static void MarkOutsiderParticles(
        BinBasedFastPointLocator<TDim>& rLocator,
        ModelPart& rParticlesModelPart,
        const Variable<TDataType>& rVariable,
        const TDataType& OutsiderValue,
        const double SearchTolerance=1e-5
    )
    {
        unsigned int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType TLS(max_results);

        //for every interface node (nodes in cut elements)
        block_for_each(rParticlesModelPart.Nodes(), TLS, [&rLocator, &rVariable, &OutsiderValue, SearchTolerance](auto& rNode, auto& rTLS)
        {

            Vector shape_functions;
            Element::Pointer p_element;
            const bool is_found = rLocator.FindPointOnMesh(rNode.Coordinates(), shape_functions, p_element, rTLS.begin(), rTLS.size(), SearchTolerance);

            if(!is_found)
            {
                if constexpr (VariableHasHistory)
                    rNode.FastGetSolutionStepValue(rVariable) = OutsiderValue;
                else
                    rNode.SetValue(rVariable, OutsiderValue);
            }
        });
    }

    /**
     * @brief This function takes a matrix of coordinates and gives the intrpolated value of the variable rInterpolationVariable
    at those positions. Note that the value is only initialized if the particle is found
     * @param rLocator this is a search structure for the volume
     * @param rCoordinates positions at which we want to retrieve the interpolation values
     * @param rInterpolationVariable: variable whose value that will be interpolated
     * @param SearchTolerance: search tolerance used in the spatial search
     * @return the function returns:
           - is_inside a vector of bools telling if the given coordinate falls within the volume
           - values the interpolated values (only valid if the corresponding value of is_inside is true)
    */
    template<unsigned int TDim, class TDataType, bool InterpolationVariableHasHistory>
    static std::pair< DenseVector<bool>, std::vector<TDataType> > InterpolateValuesAtCoordinates(
        BinBasedFastPointLocator<TDim>& rLocator,
        const Matrix& rCoordinates,
        const Variable<TDataType>& rInterpolationVariable,
        const double SearchTolerance
        )
    {
        unsigned int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType TLS(max_results);

        auto interpolations = std::make_pair(DenseVector<bool>(rCoordinates.size1()), std::vector<TDataType>(rCoordinates.size1()));

        // For every interface node (nodes in cut elements)
        const auto zero = rInterpolationVariable.Zero();
        IndexPartition(rCoordinates.size1()).for_each(TLS, [&rLocator, &rCoordinates, &interpolations, &rInterpolationVariable, &zero, SearchTolerance](const auto& i, auto& rTLS)
        {
            Vector shape_functions;
            Element::Pointer p_element;
            const bool is_found = rLocator.FindPointOnMesh(row(rCoordinates,i), shape_functions, p_element, rTLS.begin(), rTLS.size(), SearchTolerance);

            (interpolations.first)[i] = is_found;
            if(is_found)
            {
                auto& r_geom = p_element->GetGeometry();
                (interpolations.second)[i] = zero;
                for(unsigned int k=0; k<r_geom.size(); ++k)
                {
                    if constexpr (InterpolationVariableHasHistory)
                        (interpolations.second)[i] += shape_functions[k]*r_geom[k].FastGetSolutionStepValue(rInterpolationVariable);
                    else
                        (interpolations.second)[i] += shape_functions[k]*r_geom[k].GetValue(rInterpolationVariable);

                }
            }
        });

        return interpolations;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
    std::string Info() const
    {
        return std::string("ParticlesUtilities");
    };

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const {};

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {};


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
    ParticlesUtilities& operator=(ParticlesUtilities const& rOther) = delete;

    /// Copy constructor.
    ParticlesUtilities(ParticlesUtilities const& rOther) = delete;

    ///@}
}; // Class ParticlesUtilities

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ParticlesUtilities& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ParticlesUtilities& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

