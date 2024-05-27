//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//


#if !defined(KRATOS_REFINEMENT_MODELER_H_INCLUDED )
#define  KRATOS_REFINEMENT_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler/modeler.h"

#include "geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_volume_geometry.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(IGA_APPLICATION) RefinementModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(RefinementModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RefinementModeler()
        : Modeler()
    {
    }

    /// Constructor.
    RefinementModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~RefinementModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<RefinementModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void PrepareGeometryModel() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "RefinementModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream & rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream & rOStream) const override
    {
    }

    ///@}

private:
    ///@name Memory
    ///@{

    Model* mpModel = nullptr;

    ///@}
    ///@name Operations
    ///@{

    void ApplyRefinements(
        const Parameters rParameters) const;

    void ApplyRefinement(
        const Parameters rParameters) const;

    void GetGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    Parameters ReadParamatersFile(
        const std::string& rDataFileName) const;

    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

    ///@}
}; // Class RefinementModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    RefinementModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const RefinementModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_IO_MODELER_H_INCLUDED  defined
