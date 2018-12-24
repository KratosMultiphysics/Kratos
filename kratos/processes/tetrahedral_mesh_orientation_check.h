//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:           BSD License
//                          Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//

#ifndef KRATOS_TETRAHEDRAL_MESH_ORIENTATION_CHECK_H
#define KRATOS_TETRAHEDRAL_MESH_ORIENTATION_CHECK_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "geometries/geometry.h"

namespace Kratos
{
/** 
 * @class TetrahedralMeshOrientationCheck
 * @ingroup KratosCore
 * @brief Check a triangular or tetrahedral mesh to ensure that local connectivities follow the expected convention.
 * @details This process checks all elements to verify that their Jacobian has positive determinant and face conditions to ensure that all face normals point outwards.
 * @note Note that, as a side result of the procedure used, nodal normals (not normalized) are computed and stored on olution step data NORMAL.
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 */
class TetrahedralMeshOrientationCheck
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    //DEFINITION OF FLAGS TO CONTROL THE BEHAVIOUR
    KRATOS_DEFINE_LOCAL_FLAG(ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS);
    KRATOS_DEFINE_LOCAL_FLAG(COMPUTE_NODAL_NORMALS);
    KRATOS_DEFINE_LOCAL_FLAG(COMPUTE_CONDITION_NORMALS);
    KRATOS_DEFINE_LOCAL_FLAG(MAKE_VOLUMES_POSITIVE);

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedralMeshOrientationCheck);

    typedef ModelPart::ElementType ElementType;
    typedef ModelPart::ConditionType ConditionType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for TetrahedralMeshOrientationCheck Process
    /**
     * @param rModelPart The model part to check.
     * @param ThrowErrors If true, an error will be thrown if the input model part contains malformed elements or conditions.
     */
    TetrahedralMeshOrientationCheck(ModelPart& rModelPart,
                                    bool ThrowErrors,
                                    Flags options = NOT_COMPUTE_NODAL_NORMALS | NOT_COMPUTE_CONDITION_NORMALS | NOT_ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
                                    ):
        Process(),
        mrModelPart(rModelPart),
        mThrowErrors(ThrowErrors), //to be changed to a flag
        mrOptions(options)

    {
    }

    TetrahedralMeshOrientationCheck(ModelPart& rModelPart,
                                    Flags options = NOT_COMPUTE_NODAL_NORMALS | NOT_COMPUTE_CONDITION_NORMALS | NOT_ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
                                    ):
        Process(),
        mrModelPart(rModelPart),
        mThrowErrors(false),
        mrOptions(options)
    {
    }

    /// Destructor.
    ~TetrahedralMeshOrientationCheck() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{


    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override;

    void SwapAll();

    void SwapNegativeElements();

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
        return "TetrahedralMeshOrientationCheck";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TetrahedralMeshOrientationCheck";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


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
    const bool mThrowErrors;
    Flags mrOptions;

    ///@}
    ///@name Private Operations
    ///@{

    bool Orient(Geometry< Node<3> >& rGeom);

    void FaceNormal2D(array_1d<double,3>& An,
                      Geometry<Node<3> >& rGeometry);


    void FaceNormal3D(array_1d<double,3>& An,
                      Geometry<Node<3> >& rGeometry);

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TetrahedralMeshOrientationCheck& operator=(TetrahedralMeshOrientationCheck const& rOther);

    /// Copy constructor.
    TetrahedralMeshOrientationCheck(TetrahedralMeshOrientationCheck const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TetrahedralMeshOrientationCheck& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TetrahedralMeshOrientationCheck& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_TETRAHEDRAL_MESH_ORIENTATION_PROCESS_H
