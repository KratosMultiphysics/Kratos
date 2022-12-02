//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
///@addtogroup CoSimulationApplication
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class PointElement
 * @ingroup CoSimulationApplication
 * @brief Custom Point container to be used by the search
 * @details It stores the pointer of a certain element
 * @author Vicente Mataix Ferrandiz
 */
class PointElement
    : public Point
{
public:

    ///@name Type Definitions
    ///@{

    /// Base class definition
    typedef Point BaseType;

    /// Counted pointer of PointElement
    KRATOS_CLASS_POINTER_DEFINITION( PointElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointElement():
        BaseType()
    {}

    PointElement(const double X, const double Y, const double Z)
        : BaseType(X, Y, Z)
    {}

    PointElement(Element::Pointer pElement):
        mpElement(pElement)
    {
        UpdatePoint();
    }

    ///Copy constructor  (not really required)
    PointElement(const PointElement& rRHS):
        BaseType(rRHS),
        mpElement(rRHS.mpElement)
    {
    }

    /// Destructor.
    ~PointElement() override= default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the element associated to the point
     * @return mpElement The reference to the element associated to the point
     */
    Element::Pointer pGetElement()
    {
        return mpElement;
    }

    /**
     * @brief This function updates the database, using as base for the coordinates the condition center
     */
    void UpdatePoint()
    {
        noalias(this->Coordinates()) = mpElement->GetGeometry().Center().Coordinates();
    }

private:
    ///@}
    ///@name Member Variables
    ///@{

    Element::Pointer mpElement = nullptr; // The element instance

    ///@}

}; // Class PointElement

/**
 * @class DataTransfer3D1DUtilities
 * @ingroup CoSimulationApplication
 * @brief This utility includes auxiliary methods to transfer from 3D domains to 1D domains and viceversa
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CO_SIMULATION_APPLICATION) DataTransfer3D1DUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DataTransfer3D1DUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DataTransfer3D1DUtilities);

    /// Node definition
    typedef Node<3> NodeType;

    /// Geometry definition
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataTransfer3D1DUtilities() = delete;

    /// Assignment operator.
    DataTransfer3D1DUtilities& operator=(DataTransfer3D1DUtilities const& rOther) = delete;

    /// Copy constructor.
    DataTransfer3D1DUtilities(DataTransfer3D1DUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method transfer data from the 3D to the 1D
     * @param rModelPart3D The 3D model part
     * @param rModelPart1D The 1D model part
     * @param ThisParameters The parameters containing the configuration
     */
    static void From3Dto1DDataTransfer(
        ModelPart& rModelPart3D,
        ModelPart& rModelPart1D,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /**
     * @brief This method transfer data from the 1D to the 3D
     * @param rModelPart3D The 3D model part
     * @param rModelPart1D The 1D model part
     * @param ThisParameters The parameters containing the configuration
     */
    static void From1Dto3DDataTransfer(
        ModelPart& rModelPart3D,
        ModelPart& rModelPart1D,
        Parameters ThisParameters = Parameters(R"({})")
        );

    ///@}

private:
    ///@name Private Operations
    ///@{

    /**
     * @brief This method interpolates from the 1D to the 3D
     * @param rModelPart3D The 3D model part
     * @param rModelPart1D The 1D model part
     * @param ThisParameters The parameters containing the configuration
     */
    void InterpolateFrom1Dto3D(
        ModelPart& rModelPart3D,
        ModelPart& rModelPart1D,
        Parameters ThisParameters
        );

    /**
     * @brief This method interpolates from the 3D to the 1D
     * @param rModelPart3D The 3D model part
     * @param rModelPart1D The 1D model part
     * @param ThisParameters The parameters containing the configuration
     */
    void InterpolateFrom3Dto1D(
        ModelPart& rModelPart3D,
        ModelPart& rModelPart1D,
        Parameters ThisParameters
        );

    /**
     * @brief This method extrapolates from the 1D to the 3D
     * @param rModelPart3D The 3D model part
     * @param rModelPart1D The 1D model part
     * @param ThisParameters The parameters containing the configuration
     */
    void ExtrapolateFrom1Dto3D(
        ModelPart& rModelPart3D,
        ModelPart& rModelPart1D,
        Parameters ThisParameters
        );

    /**
     * @brief This method extrapolates from the 3D to the 1D
     * @param rModelPart3D The 3D model part
     * @param rModelPart1D The 1D model part
     * @param ThisParameters The parameters containing the configuration
     */
    void ExtrapolateFrom3Dto1D(
        ModelPart& rModelPart3D,
        ModelPart& rModelPart1D,
        Parameters ThisParameters
        );

    /**
     * @brief This method computes the list of variables to interpolate/extrapolate
     * @param ThisParameters The parameters containing the list of variables
     * @param rOriginListVariables The list of origin variables to interpolate/extrapolate
     * @param rDestinationListVariables The list of destination variables to interpolate/extrapolate
     */
    static void GetVariablesList(
        Parameters ThisParameters,
        std::vector<const Variable<double>*>& rOriginListVariables,
        std::vector<const Variable<double>*>& rDestinationListVariables
        );

    /**
     * @brief This method computes maximum length of the elements
     * @param rModelPart The model part to compute
     * @return The maximum length
     */
    static double GetMaxLength(ModelPart& rModelPart);

    /**
     * @brief This method returns the default parameters
     * @return The default parameters
     */
    static Parameters GetDefaultParameters();

    ///@}
}; // Class DataTransfer3D1DUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.