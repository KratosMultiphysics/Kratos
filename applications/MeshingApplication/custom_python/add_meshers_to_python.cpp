// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_meshers_to_python.h"

#ifdef USE_TETGEN_NONFREE_TPL
    #include "external_includes/tetgen_pfem_refine.h"
    #include "external_includes/tetgen_pfem_refine_vms.h"
    #include "external_includes/tetgen_pfem_refine_face.h"
    #include "external_includes/tetgen_pfem_contact.h"
    #include "external_includes/tetgen_cdt.h"
#else
    #define REAL double
#endif 

#include "external_includes/trigen_pfem_refine.h"
#include "external_includes/trigen_pfem_refine_vms.h"
#include "external_includes/trigen_pfem_refine_segment.h"
#include "external_includes/trigen_glass_forming.h"
#include "external_includes/trigen_droplet_refine.h"
#include "external_includes/trigen_cdt.h"

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

#ifdef USE_TETGEN_NONFREE_TPL
///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                                            ADAPTIVE 3D MESHER                         //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

//tetgen pfem refine
void TetRegenerateMesh(TetGenPfemModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& rModelPart, NodeEraseProcess& NodeErase,bool RemNodes, bool AddNodes, double AlphaShape, double HFactor)
{
    Mesher.ReGenerateMesh(rModelPart,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),NodeErase,RemNodes, AddNodes,  AlphaShape, HFactor     );
}

//tetgen pfem refine
void GenerateCDT(TetGenCDT& Mesher, ModelPart& rModelPart, char* ElementName, bool AddNodes, unsigned int PropertyId)
{
    Mesher.GenerateCDT(rModelPart,KratosComponents<Element>::Get(ElementName),AddNodes,  PropertyId  );
}


///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                                     ADAPTIVE 3D MESHER ----> Using Face               //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

//tetgen pfem refine
void TetRegenerateMeshWithFace(TetGenPfemRefineFace& Mesher, char* ElementName, char* ConditionName, ModelPart& rModelPart, ModelPart::ElementsContainerType& rElements, NodeEraseProcess& NodeErase, bool RemNodes, bool AddNodes, double AlphaShape, double HFactor )
{
    Mesher.ReGenerateMesh(rModelPart,rElements,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),NodeErase, RemNodes, AddNodes, AlphaShape, HFactor     );
}
///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                                  Contact Mesher                                       //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

//tetgen pfem refine
void TetRegenerateMeshContact(TetGenPfemContact& Mesher, char* ElementName, char* ConditionName, ModelPart& rModelPart )
{
    Mesher.ReGenerateMesh(rModelPart,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName)
                         );
}
///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                                      VMS Mesher                                       //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

//tetgen pfem refine
void TetRegenerateMeshVMS(TetGenPfemModelerVms& Mesher, char* ElementName, char* ConditionName, ModelPart& rModelPart, NodeEraseProcess& NodeErase,bool RemNodes, bool AddNodes, double AlphaShape, double HFactor)
{
    Mesher.ReGenerateMesh(rModelPart,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),NodeErase,RemNodes, AddNodes,  AlphaShape, HFactor     );


}
#endif

///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                                            ADAPTIVE 2D MESHER                         //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

//trigen pfem refine
void TriRegenerateMesh(TriGenPFEMModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& rModelPart,NodeEraseProcess& NodeErase, bool RemNodes, bool AddNodes, double AlphaShape, double HFactor )
{
    Mesher.ReGenerateMesh(rModelPart,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),NodeErase, RemNodes, AddNodes, AlphaShape, HFactor     );
}

///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                    ADAPTIVE 2D MESHER specifically for GLASS FORMING                  //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

void TriRegenerateMeshGLASS(TriGenGLASSModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& rModelPart,NodeEraseProcess& NodeErase, bool RemNodes, bool AddNodes, double AlphaShape, double HFactor )
{
    Mesher.ReGenerateMesh(rModelPart,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),NodeErase, RemNodes, AddNodes, AlphaShape, HFactor     );
}



///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                    ADAPTIVE 2D MESHER specifically for DROPLET FORMING                //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

void TriRegenerateMeshDroplet(TriGenDropletModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& rModelPart,NodeEraseProcess& NodeErase, bool RemNodes, bool AddNodes, double AlphaShape, double HFactor )
{
    Mesher.ReGenerateMeshDroplet(rModelPart,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),NodeErase, RemNodes, AddNodes, AlphaShape, HFactor     );
}

///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                             BOUNDARY CONFORMANT 2D MESHER                             //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////


void TriGenCDTFluid(TriGenCDT& Mesher,ModelPart& rModelPart)
{
    Mesher.GenerateCDT(rModelPart,
                       KratosComponents<Element>::Get("Fluid2D")
                      );
}

///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                          ADAPTIVE 2D MESHER -->USING SEGMENT                          //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

//trigen pfem refine segment
void TriRegenerateMeshWithSegment(TriGenPFEMRefineSegment& Mesher, char* ElementName, char* ConditionName, ModelPart& rModelPart,NodeEraseProcess& NodeErase, bool RemNodes, bool AddNodes, double AlphaShape, double HFactor )
{
    Mesher.ReGenerateMesh(rModelPart,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),NodeErase, RemNodes, AddNodes, AlphaShape, HFactor     );
}

///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//                                   VMS 2D Mesher                                       //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

//trigen pfem vms refine
void TriRegenerateMeshVMS(TriGenPFEMModelerVMS& Mesher, char* ElementName, char* ConditionName, ModelPart& rModelPart,NodeEraseProcess& NodeErase, bool RemNodes, bool AddNodes, double AlphaShape, double HFactor )
{
    Mesher.ReGenerateMesh(rModelPart,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),NodeErase, RemNodes, AddNodes, AlphaShape, HFactor     );
}

void  AddMeshersToPython(pybind11::module& m)
{
    
#ifdef USE_TETGEN_NONFREE_TPL
    // Class that allows 3D adaptive remeshing (inserting and erasing nodes)
    py::class_<TetGenPfemModeler, TetGenPfemModeler::Pointer >(m, "TetGenPfemModeler")
    .def(py::init< >())
    .def("ReGenerateMesh",TetRegenerateMesh)
    .def("ReGenerateMesh",&TetGenPfemModeler::ReGenerateMesh)
    ;

    py::class_<TetGenPfemRefineFace, TetGenPfemRefineFace::Pointer >(m, "TetGenPfemRefineFace")
    .def(py::init< >())
    .def("ReGenerateMesh",TetRegenerateMeshWithFace)
    ;

    py::class_<TetGenPfemContact, TetGenPfemContact::Pointer >(m, "TetGenPfemContact")
    .def(py::init< >())
    .def("ReGenerateMesh",TetRegenerateMeshContact)
    ;
    
    py::class_<TetGenCDT, TetGenCDT::Pointer >(m, "TetGenCDT")
    .def(py::init< >())
    .def("GenerateCDT",GenerateCDT)
    ;

    py::class_<TetGenPfemModelerVms, TetGenPfemModelerVms::Pointer >(m, "TetGenPfemModelerVms")
    .def(py::init< >())
    .def("ReGenerateMesh",&TetGenPfemModelerVms::ReGenerateMesh)
    ;
#endif
    
    // Class that allows 2D adaptive remeshing (inserting and erasing nodes)
    py::class_<TriGenPFEMModeler, TriGenPFEMModeler::Pointer >(m, "TriGenPFEMModeler")
    .def(py::init< >())
    .def("ReGenerateMesh",TriRegenerateMesh)
    .def("ReGenerateMesh",&TriGenPFEMModeler::ReGenerateMesh)
    ;

    // Class that allows 2D adaptive remeshing (inserting and erasing nodes) as well as preserving the topology (avoiding "holes" at the boundaries). made for glass simulation
    py::class_<TriGenGLASSModeler, TriGenGLASSModeler::Pointer >(m, "TriGenGLASSModeler")
    .def(py::init< >())
    .def("ReGenerateMeshGlass",TriRegenerateMeshGLASS)
    .def("ReGenerateMeshGlass",&TriGenGLASSModeler::ReGenerateMesh)
    ;

    // Class that allows 2D adaptive remeshing (inserting and erasing nodes) as well as preserving the topology (avoiding "holes" at the boundaries). made for droplet simulation
    py::class_<TriGenDropletModeler, TriGenDropletModeler::Pointer >(m,"TriGenDropletModeler")
    .def(py::init< >())
    .def("ReGenerateMeshDroplet",TriRegenerateMeshDroplet)
    .def("ReGenerateMeshDroplet",&TriGenDropletModeler::ReGenerateMeshDroplet)
    ;

    py::class_<TriGenPFEMModelerVMS, TriGenPFEMModelerVMS::Pointer>(m, "TriGenPFEMModelerVMS")
    .def(py::init< >())
    .def("ReGenerateMesh",&TriGenPFEMModelerVMS::ReGenerateMesh);

    //segment mesher adaptive
    py::class_<TriGenPFEMRefineSegment, TriGenPFEMRefineSegment::Pointer >(m, "TriGenPFEMSegment")
    .def(py::init< >())
    .def("ReGenerateMesh",TriRegenerateMeshWithSegment)
    ;
}

}  // namespace Python.

} // Namespace Kratos

