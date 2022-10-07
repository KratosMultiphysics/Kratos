//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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

/// Collection of utilities to compute statistic of particles
/** Detail class definition.
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
    ParticlesUtilities(){};

    /// Destructor.
    virtual ~ParticlesUtilities(){};

    template<unsigned int TDim>
    static void CountParticlesInNodes(
            BinBasedFastPointLocator<TDim>& rLocator,
            ModelPart& rVolumeModelPart,
            const ModelPart& rParticlesModelPart,
            const Variable<double>& rCounterVariable
            )
        {

            //reset the counter variable
            block_for_each(rVolumeModelPart.Nodes(), [&rCounterVariable](auto& rNode){
                rNode.GetValue(rCounterVariable) = 0.0;
            });


            const unsigned int max_results = 10000; 
            typename BinBasedFastPointLocator<TDim>::ResultContainerType TLS(max_results);
            
            //for each node in the mesh, count the number of particles in the neighbouring elements
            block_for_each(rParticlesModelPart.Nodes(), TLS, [&rLocator, &rCounterVariable](const auto& rNode, auto& rTLS){

                Vector shape_functions;
                Element::Pointer p_element;
                const bool is_found = rLocator.FindPointOnMesh(rNode.Coordinates(), shape_functions, p_element, rTLS.begin(), rTLS.size(), 1e-5);

                if(is_found){

                    auto& r_geom = p_element->GetGeometry();
                    for(unsigned int i=0; i<r_geom.size(); ++i){
                        auto& rcounter = r_geom[i].GetValue(rCounterVariable);
                        AtomicAdd(rcounter, 1.0);
                    }
                }

            });

        }

    /**
    this function looks if particle is found in the locator. If it is not, rVariable is marked with the value "outsider_value"
    */
    template<unsigned int TDim, class TScalarType>
    static void MarkOutsiders(
            BinBasedFastPointLocator<TDim>& rLocator,
            ModelPart& rParticlesModelPart,
            const Variable<TScalarType>& rVariable,
            const TScalarType outsider_value
            )
        {
            const unsigned int max_results = 10000; 
            typename BinBasedFastPointLocator<TDim>::ResultContainerType TLS(max_results);
            
            //for every particle search if it is inside or outside the volume mesh
            block_for_each(rParticlesModelPart.Nodes(), TLS, [&rLocator, &rVariable, outsider_value](auto& rNode, auto& rTLS){

                Vector shape_functions;
                Element::Pointer p_element;
                bool is_found = rLocator.FindPointOnMesh(rNode.Coordinates(), shape_functions, p_element, rTLS.begin(), rTLS.size(), 1e-5);

                if(!is_found){
                    rNode.SetValue(rVariable, outsider_value);
                }
            });
        }

    /**@brief provides the value of the interpolated variable "rInterpolationVariable" at all the positions indicated in rCoordinates
    */
    template<unsigned int TDim, class TScalarType>
    static std::pair< std::vector<bool>, std::vector<TScalarType> > InterpolateValuesAtPoints(
            BinBasedFastPointLocator<TDim>& rLocator,
            const Matrix& rCoordinates,
            const Variable<TScalarType>& rInterpolationVariable
            )
        {
            const unsigned int max_results = 10000; 
            typename BinBasedFastPointLocator<TDim>::ResultContainerType TLS(max_results);

            auto interpolations = std::make_pair(std::vector<bool>(rCoordinates.size1()), std::vector<TScalarType>(rCoordinates.size1()));            
            
            //for each particle interpolate the given variable from the volume mesh
            IndexPartition(rCoordinates.size1()).for_each(TLS, [&rLocator, &rCoordinates, &interpolations, &rInterpolationVariable](const auto& i, auto& rTLS){

                Vector shape_functions;
                Element::Pointer p_element;
                bool is_found = rLocator.FindPointOnMesh(row(rCoordinates,i), shape_functions, p_element, rTLS.begin(), rTLS.size(), 1e-5);

                (interpolations.first)[i] = is_found;
                if(is_found){

                    auto& r_geom = p_element->GetGeometry();
                    (interpolations.second)[i] = TScalarType();
                    for(unsigned int k=0; k<r_geom.size(); ++k){
                        (interpolations.second)[i] += shape_functions[k]*r_geom[k].FastGetSolutionStepValue(rInterpolationVariable);
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
    virtual std::string Info() const {return std::string("ParticlesUtilities");};

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {};

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {};


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

    /// Assignment operator.
    ParticlesUtilities& operator=(ParticlesUtilities const& rOther);

    /// Copy constructor.
    ParticlesUtilities(ParticlesUtilities const& rOther);


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
                ParticlesUtilities& rThis);

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

