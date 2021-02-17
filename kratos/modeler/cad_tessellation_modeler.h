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

#if !defined(KRATOS_CAD_TESSELLATION_MODELER_INCLUDED)
#define KRATOS_CAD_TESSELLATION_MODELER_INCLUDED

extern "C"
{
    #ifdef SINGLE
        #define REAL float
    #else /* not SINGLE */
        #define REAL double
    #endif /* not SINGLE */
    void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
}

// System includes

// External includes
//this is not ideal...
#include "../../external_libraries/triangle/triangle.h"

// Project includes
#include "modeler.h"
#include "utilities/nurbs_curve_tessellation.h"
#include "geometries/brep_curve_on_surface.h"

namespace Kratos {

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class KRATOS_API(KRATOS_CORE) CadTessellationModeler : public Modeler {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(CadTessellationModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;




    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CadTessellationModeler() : Modeler()
    {
    }

    /// Constructor.
    CadTessellationModeler(Model& rModel, Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters), mpModel(&rModel)
    {
    }

    /// Destructor.
    virtual ~CadTessellationModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<CadTessellationModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupModelPart() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CadTessellationModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;

    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method returns the triangulation of a NURBS surface
     */
    void ComputeTriangulation();

    /**
     * @brief This method initializes the necessary data structure for triangle
     * @param tr This is a struct used to pass data into and out of the triangulate() procedure
     */
    void InitTriangulationDataStructure(triangulateio& tr);

    /**
     * @brief This method clears the data structure for triangle
     * @param tr This is a struct used to pass data into and out of the triangulate() procedure
     */
    void CleanTriangulationDataStructure(triangulateio& tr);



}; // Class CadTessellationModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, CadTessellationModeler& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const CadTessellationModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_CAD_TESSELLATION_MODELER_INCLUDED  defined
