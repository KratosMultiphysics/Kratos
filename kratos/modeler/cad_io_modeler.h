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


#if !defined(KRATOS_IGA_IO_MODELER_H_INCLUDED )
#define  KRATOS_IGA_IO_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/math_utils.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) CadIoModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(CadIoModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CadIoModeler()
        : Modeler()
    {
    }

    /// Constructor.
    CadIoModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~CadIoModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<CadIoModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupGeometryModel() override;

    void SetupModelPart() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CadIoModeler";
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

    Parameters ReadParamatersFile(
        const std::string& rDataFileName) const;    

    void CreateTheSnakeCoordinates(Vector& knot_vector_u, Vector& knot_vector_v, double& knot_step_u, double& knot_step_v, const Parameters refinements_parameters, ModelPart& surrogate_model_part);

    void SnakeStep(std::vector<std::vector<int>> &knot_spans_available, int knot_span_u_1st_point, int knot_span_u_2nd_point, int knot_span_v_1st_point,int knot_span_v_2nd_point, double& x_true_boundary1, double& x_true_boundary2, double& y_true_boundary1, double& y_true_boundary2, double& knot_step_u, double& knot_step_v);

    ///@}

private:
    ///@name Iga functionalities
    ///@{
    using PointType = Node;
    using PointTypePointer = Node::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;
    using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    using PointerType = DynamicBins::PointerType;

    Model* mpModel = nullptr;


    bool isPointInsideSkinBoundary(Point& pointToSearch, DynamicBins& testBins, ModelPart& skin_model_part);

    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

    ///@}
}; // Class CadIoModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    CadIoModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const CadIoModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_IO_MODELER_H_INCLUDED  defined
