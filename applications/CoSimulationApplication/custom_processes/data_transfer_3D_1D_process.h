// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:         BSD License
//                   license: CoSimulationApplication/license.txt
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
#include "processes/process.h"

/* The mappers includes */
#include "spaces/ublas_space.h"
#include "mappers/mapper_flags.h"
#include "factories/mapper_factory.h" 

namespace Kratos
{
///@addtogroup CoSimulationApplication
///@{

///@}
///@name Functions
///@{

/**
 * @brief This method determines the 1D model part
 * @param rFirstModelPart The first ModelPart
 * @param rSecondModelPart The second ModelPart
 * @return The 1D model part 
 */
ModelPart& Determine1DModelPart(ModelPart& rFirstModelPart, ModelPart& rSecondModelPart)
{
    // In MPI could be a bit troublesome
    if (rFirstModelPart.IsDistributed()) { 
        KRATOS_ERROR << "Not compatible with MPI" << std::endl;
    } else { // In serial we just check model part has 1D entities (local space)
        // Determine if the modelparts have entities

        // First model part
        if (rFirstModelPart.NumberOfElements() > 0 || rFirstModelPart.NumberOfConditions() > 0 ) {
            if (rFirstModelPart.NumberOfElements() > 0) {
                if (rFirstModelPart.ElementsBegin()->GetGeometry().LocalSpaceDimension() == 1) {
                    return rFirstModelPart;
                }
            } else {
                if (rFirstModelPart.ConditionsBegin()->GetGeometry().LocalSpaceDimension() == 1) {
                    return rFirstModelPart;
                }
            }
        }

        // Second model part
        if (rSecondModelPart.NumberOfElements() > 0 || rSecondModelPart.NumberOfConditions() > 0 ) {
            if (rSecondModelPart.NumberOfElements() > 0) {
                if (rSecondModelPart.ElementsBegin()->GetGeometry().LocalSpaceDimension() == 1) {
                    return rSecondModelPart;
                }
            } else {
                if (rSecondModelPart.ConditionsBegin()->GetGeometry().LocalSpaceDimension() == 1) {
                    return rSecondModelPart;
                }
            }
        }
    }

    KRATOS_ERROR << "Impossible to detect 1D model part" << std::endl;
    return rFirstModelPart;
}

/**
 * @brief This method determines the 3D model part
 * @details The counter part of the previous method
 * @param rFirstModelPart The first ModelPart
 * @param rSecondModelPart The second ModelPart
 * @return The 3D model part 
 */
ModelPart& Determine3DModelPart(ModelPart& rFirstModelPart, ModelPart& rSecondModelPart)
{
    // It is the counter part of the previous method
    ModelPart* p_1d_model_part = &Determine1DModelPart(rFirstModelPart, rSecondModelPart);
    if (p_1d_model_part == &rFirstModelPart) {
        return rSecondModelPart;
    } else {
        return rFirstModelPart;
    }
}

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

#define DEFINE_MAPPER_FACTORY_SERIAL                                                                                             \
using SparseSpace = UblasSpace<double, boost::numeric::ublas::compressed_matrix<double>, boost::numeric::ublas::vector<double>>; \
using DenseSpace = UblasSpace<double, DenseMatrix<double>, DenseVector<double>>;                                                 \
using MapperFactoryType = MapperFactory<SparseSpace, DenseSpace>;                                                                \
using MapperType = Mapper<SparseSpace, DenseSpace>;

/**
 * @class DataTransfer3D1DProcess
 * @ingroup CoSimulationApplication
 * @brief This utility includes auxiliary methods to transfer from 3D domains to 1D domains and viceversa
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CO_SIMULATION_APPLICATION) DataTransfer3D1DProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DataTransfer3D1DProcess
    KRATOS_CLASS_POINTER_DEFINITION(DataTransfer3D1DProcess);

    /// Geometry definition
    using GeometryType = Geometry<Node>;

    // Define mapper factory
    DEFINE_MAPPER_FACTORY_SERIAL

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataTransfer3D1DProcess(
        ModelPart& rFirstModelPart,
        ModelPart& rSecondModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~DataTransfer3D1DProcess() override;

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method executes the process
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;
    
    ///@}
private:
    ///@name Private member Variables
    ///@{

    ModelPart& mr3DModelPart;                                       /// The 3D model part
    ModelPart& mr1DModelPart;                                       /// The 1D model part
    Parameters mThisParameters;                                     /// The parameters containing the configuration
    MapperType::Pointer mpMapper = nullptr;                         /// The mapper pointer
    std::vector<const Variable<double>*> mOriginListVariables;      /// The list of variables to be transferred (origin)
    std::vector<const Variable<double>*> mDestinationListVariables; /// The list of variables to be transferred (destination)

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method interpolates from the 3D to the 1D
     */
    void InterpolateFrom3Dto1D();

    /**
     * @brief This method interpolates from the 1D to the 3D
     */
    void InterpolateFrom1Dto3D();

    /**
     * @brief This method computes the list of variables to interpolate
     * @param rOriginListVariables The list of origin variables to interpolate
     * @param rDestinationListVariables The list of destination variables to interpolate
     */
    void GetVariablesList(
        std::vector<const Variable<double>*>& rOriginListVariables,
        std::vector<const Variable<double>*>& rDestinationListVariables
        );

    /**
     * @brief This method computes maximum length of the elements
     * @param rModelPart The model part to compute
     * @return The maximum length
     */
    static double GetMaxLength(ModelPart& rModelPart);

    ///@}
}; // Class DataTransfer3D1DProcess

///@}

///@} addtogroup block

}  // namespace Kratos.