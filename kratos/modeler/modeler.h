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

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_components.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class Modeler
 * @ingroup KratosCore
 * @brief Abstract base class for geometry and mesh manipulation in Kratos.
 * @details The Modeler provides a unified interface for generating, preparing, and manipulating
 * geometric and mesh entities required for computational analysis. It serves as the 
 * base class from which specific geometry or mesh-generating modelers should be derived.
 *
 * Key responsibilities include:
 * - Importing or generating geometry from external data.
 * - Preparing and updating geometry configurations.
 * - Generating ModelParts suitable for analyses.
 * @author Riccardo Rossi
 */
class Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(Modeler);

    /// Size type used for indices and counters
    using SizeType = std::size_t;

    /// Index type definition
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * @brief Default constructor.
    * @param ModelerParameters Configuration parameters for the modeler.
    */
    explicit Modeler(
        Parameters ModelerParameters = Parameters())
        : mParameters(ModelerParameters)
        , mEchoLevel(
            ModelerParameters.Has("echo_level")
            ? ModelerParameters["echo_level"].GetInt()
            : 0)
    {}

    /**
    * @brief Constructor with reference to a Model instance.
    * @param rModel Reference to a Kratos Model.
    * @param ModelerParameters Configuration parameters for the modeler.
    */
    explicit Modeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : mParameters(ModelerParameters)
        , mEchoLevel(
            ModelerParameters.Has("echo_level")
            ? ModelerParameters["echo_level"].GetInt()
            : 0)
    {}

    /** 
    * @brief Destructor. 
    */
    virtual ~Modeler() = default;

    /**
    * @brief Creates a new Modeler instance.
    * @details This method must be implemented in derived classes to create new instances
    * according to the specific needs of each derived Modeler.
    * @param rModel Reference to the Kratos Model.
    * @param ModelParameters Parameters configuring the new Modeler.
    * @return Pointer to the new Modeler instance.
    */
    virtual Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const
    {
        KRATOS_ERROR << "Trying to Create Modeler. Please check derived class 'Create' definition." << std::endl;
    }

    ///@}
    ///@name Modeler Stages at Initialize
    ///@{

    /** 
    * @brief Set up or import the geometry model from external inputs.
    */
    virtual void SetupGeometryModel()
    {}

    /** 
    * @brief Prepare or update the geometry model for further processing. 
    */
    virtual void PrepareGeometryModel()
    {}

    /** 
    * @brief Configure the geometry into a suitable ModelPart for analysis. 
    */
    virtual void SetupModelPart()
    {}

    /**
    * @brief Provides default parameters to ensure consistency across constructors.
    * @details Derived classes should override this method.
    * @return Default parameters as a Parameters object.
    */
    virtual const Parameters GetDefaultParameters() const
    {
        KRATOS_ERROR << "Calling base Modeler class GetDefaultParameters. "
                    << "Please implement this in your derived class." << std::endl;
        const Parameters default_parameters = Parameters(R"({})");
        return default_parameters;
    }

    ///@}
    ///@name Mesh and Node Generation
    ///@{

    /**
    * @brief Generate a new ModelPart from an existing one.
    * @param rOriginModelPart The source ModelPart.
    * @param rDestinationModelPart The resulting ModelPart after generation.
    * @param rReferenceElement Reference element for mesh creation.
    * @param rReferenceBoundaryCondition Reference condition for boundaries.
    */
    virtual void GenerateModelPart(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart,
                                Element const& rReferenceElement,
                                Condition const& rReferenceBoundaryCondition)
    {
        KRATOS_ERROR << "This modeler CAN NOT be used for mesh generation." << std::endl;
    }

    /**
    * @brief Generate mesh elements within a given ModelPart.
    * @param ThisModelPart Target ModelPart for mesh generation.
    * @param rReferenceElement Reference element for the generated mesh.
    * @param rReferenceBoundaryCondition Reference boundary conditions.
    */
    virtual void GenerateMesh(ModelPart& ThisModelPart,
                            Element const& rReferenceElement,
                            Condition const& rReferenceBoundaryCondition)
    {
        KRATOS_ERROR << "This modeler CAN NOT be used for mesh generation." << std::endl;
    }

    /**
    * @brief Generate nodes within the specified ModelPart.
    * @param ThisModelPart ModelPart in which nodes are generated.
    */
    virtual void GenerateNodes(ModelPart& ThisModelPart)
    {
        KRATOS_ERROR << "This modeler CAN NOT be used for node generation." << std::endl;
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
    * @brief Provides a short description of the Modeler instance.
    * @return Short description as a string.
    */
    virtual std::string Info() const
    {
        return "Modeler";
    }

    /**
    * @brief Prints basic information about the Modeler.
    * @param rOStream Output stream to write information.
    */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /**
    * @brief Prints detailed data about the Modeler.
    * @param rOStream Output stream for detailed data.
    */
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
protected:
    ///@name Protected Members
    ///@{

    Parameters mParameters; /// Configuration parameters for the Modeler.
    SizeType mEchoLevel;    /// Verbosity level (0=silent).

    ///@}

}; // Class Modeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Modeler& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Modeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<Modeler>;

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Modeler const& ThisComponent);

} // namespace Kratos
 