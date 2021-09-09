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

#ifndef KRATOS_MESH_MOVING_MODELER_H_INCLUDED
#define KRATOS_MESH_MOVING_MODELER_H_INCLUDED


// System includes


// External includes


// Project includes
#include "modeler/modeler.h"


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
 * @class MeshMovingModeler
 * @ingroup KratosShallowWaterApplication
 * @brief Tools for lagrangian computations
 * @details Generate a model part on the wet domain
 * @author Miguel Maso Sotomayor
*/
class KRATOS_API(SHALLOW_WATER_APPLICATION) MeshMovingModeler : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MeshMovingModeler
    KRATOS_CLASS_POINTER_DEFINITION(MeshMovingModeler);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG(TO_COPY);

    typedef Node<3> NodeType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    MeshMovingModeler();

    /**
     * @brief Constructor
     */
    MeshMovingModeler(Model& rModel, Parameters ModelerParameters = Parameters());

    /**
     * @brief Destructor
     */
    virtual ~MeshMovingModeler() = default;

    /**
     * @brief Creates the Modeler Pointer
     */
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<MeshMovingModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Import or generate geometry models from external input.
    void SetupGeometryModel() override;

    /// Prepare or update the geometry model_part.
    void PrepareGeometryModel() override;

    /// Convert the geometry model or import analysis suitable models.
    void SetupModelPart() override;

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MeshMovingModeler";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Member Variables
    ///@{

    Model* mpModel;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    const Parameters GetDefaultParameters() const;

    ///@}

}; // Class MeshMovingModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                MeshMovingModeler& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const MeshMovingModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MESH_MOVING_MODELER_H_INCLUDED  defined
