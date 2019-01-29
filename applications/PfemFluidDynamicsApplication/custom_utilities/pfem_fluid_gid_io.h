//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta    maceli@cimne.upc.edu
//




#if !defined(PFEM_FLUID_GID_IO_BASE_H_INCLUDED)
#define  PFEM_FLUID_GID_IO_BASE_H_INCLUDED

// External includes
#define USE_CONST
#include "gidpost/source/gidpost.h"

// Project includes
#include "includes/gid_io.h"

namespace Kratos
{

template<class TGaussPointContainer = GidGaussPointsContainer, class TMeshContainer = GidMeshContainer>
class PfemFluidGidIO : public GidIO<TGaussPointContainer, TMeshContainer>
{

typedef Node < 3 > NodeType;
typedef Properties PropertiesType;
typedef Element ElementType;
typedef Condition ConditionType;
typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;

public:
    ///pointer definition of PfemFluidGidIO
    KRATOS_CLASS_POINTER_DEFINITION(PfemFluidGidIO);

    ///typedefs
    typedef GidIO<TGaussPointContainer, TMeshContainer> BaseType;

    ///Flags for mesh writing
//             enum WriteDeformedMeshFlag{WriteDeformed, WriteUndeformed};
//             enum WriteConditionsFlag{WriteConditions, WriteElementsOnly};
//             enum MultiFileFlag{SingleFile, MultipleFiles};

    ///Constructor
    ///single stream IO constructor
    PfemFluidGidIO( const std::string& rDatafilename,
           GiD_PostMode Mode,
           MultiFileFlag use_multiple_files_flag,
           WriteDeformedMeshFlag write_deformed_flag,
           WriteConditionsFlag write_conditions_flag
         ) : BaseType(rDatafilename, Mode, use_multiple_files_flag, write_deformed_flag, write_conditions_flag) {}

    ///Destructor.
    ~PfemFluidGidIO() override
    {}

    void WriteNodeMesh( MeshType& rThisMesh ) override
    {
        KRATOS_TRY

        Timer::Start("Writing Mesh");

        GiD_fBeginMesh(BaseType::mMeshFile,  "Kratos Mesh",GiD_3D,GiD_Point,1);
        GiD_fBeginCoordinates(BaseType::mMeshFile);
        for ( MeshType::NodeIterator node_iterator = rThisMesh.NodesBegin();
                node_iterator != rThisMesh.NodesEnd();
                ++node_iterator)
        {
            if ( BaseType::mWriteDeformed == WriteUndeformed )
                GiD_fWriteCoordinates(BaseType::mMeshFile, node_iterator->Id(), node_iterator->X0(),
                                      node_iterator->Y0(), node_iterator->Z0() );
            else if ( BaseType::mWriteDeformed == WriteDeformed )
                GiD_fWriteCoordinates(BaseType::mMeshFile, node_iterator->Id(), node_iterator->X(),
                                      node_iterator->Y(), node_iterator->Z() );
            else
                KRATOS_ERROR << "Undefined WriteDeformedMeshFlag" << std::endl;
        }
        GiD_fEndCoordinates(BaseType::mMeshFile);
        int nodes_id[2];
        GiD_fBeginElements(BaseType::mMeshFile);

        for ( MeshType::NodeIterator node_iterator = rThisMesh.NodesBegin();
                node_iterator != rThisMesh.NodesEnd();
                ++node_iterator)
        {
            nodes_id[0] = node_iterator->Id();
            nodes_id[1] = node_iterator->Is(ISOLATED);
            GiD_fWriteElementMat(BaseType::mMeshFile,node_iterator->Id(), nodes_id);
        }
        GiD_fEndElements(BaseType::mMeshFile);
        GiD_fEndMesh(BaseType::mMeshFile);

        Timer::Stop("Writing Mesh");

        KRATOS_CATCH("")
    }//WriteNodeMesh

protected:


private:
    PfemFluidGidIO& operator=(PfemFluidGidIO const& rOther);
    PfemFluidGidIO(PfemFluidGidIO const& rOther);

}; // Class PfemFluidGidIO


inline std::ostream& operator << (std::ostream& rOStream, const PfemFluidGidIO<>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.


#endif // PFEM_FLUID_GID_IO_BASE_H_INCLUDED  defined
