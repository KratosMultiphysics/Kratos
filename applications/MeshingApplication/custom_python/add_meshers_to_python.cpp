//
//   Project Name:        Kratos
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-22 16:37:10 $
//   Revision:            $Revision: 1.13 $
//
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

// #include "triangle.h"
#include "external_includes/trigen_pfem_refine.h"
#include "external_includes/trigen_pfem_refine_vms.h"
#include "external_includes/trigen_pfem_refine_segment.h"

#include "external_includes/trigen_glass_forming.h"

#include "external_includes/trigen_droplet_refine.h"


//#include "external_includes/trigen_mesh_suite.h"
#include "external_includes/trigen_cdt.h"
//#include "external_includes/trigen_refine.h"

// #include "external_includes/msuite_pfem_refine.h"


namespace Kratos
{

namespace Python
{
    
using namespace pybind11;

#ifdef USE_TETGEN_NONFREE_TPL
///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				ADAPTIVE 3D MESHER					//
//											//
//////////////////////////////////////////////////////////////////////////////////////////

//tetgen pfem refine
void TetRegenerateMesh(TetGenPfemModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part, NodeEraseProcess& node_erase,bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),node_erase,rem_nodes, add_nodes,  alpha_shape, h_factor	);
}

//tetgen pfem refine
void GenerateCDT(TetGenCDT& Mesher, ModelPart& model_part, char* ElementName, bool add_nodes, unsigned int property_id)
{
    Mesher.GenerateCDT(model_part,KratosComponents<Element>::Get(ElementName),add_nodes,  property_id  );
}


///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				ADAPTIVE 3D MESHER ----> Using Face			//
//											//
//////////////////////////////////////////////////////////////////////////////////////////

//tetgen pfem refine
void TetRegenerateMeshWithFace(TetGenPfemRefineFace& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part, ModelPart::ElementsContainerType& rElements, NodeEraseProcess& node_erase, bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor )
{
    Mesher.ReGenerateMesh(model_part,rElements,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),node_erase, rem_nodes, add_nodes, alpha_shape, h_factor	);
}
///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				Contact Mesher			                        //
//											//
//////////////////////////////////////////////////////////////////////////////////////////

//tetgen pfem refine
void TetRegenerateMeshContact(TetGenPfemContact& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part )
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName)
                         );
}
///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				VMS Mesher			                        //
//											//
//////////////////////////////////////////////////////////////////////////////////////////

//tetgen pfem refine
void TetRegenerateMeshVMS(TetGenPfemModelerVms& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part, NodeEraseProcess& node_erase,bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),node_erase,rem_nodes, add_nodes,  alpha_shape, h_factor	);


}
#endif



///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				ADAPTIVE 2D MESHER					//
//											//
//////////////////////////////////////////////////////////////////////////////////////////

//trigen pfem refine
void TriRegenerateMesh(TriGenPFEMModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part,NodeEraseProcess& node_erase, bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor )
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),node_erase, rem_nodes, add_nodes, alpha_shape, h_factor	);
}

///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				ADAPTIVE 2D MESHER specifically for GLASS FORMING 	//
//											//
//////////////////////////////////////////////////////////////////////////////////////////

void TriRegenerateMeshGLASS(TriGenGLASSModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part,NodeEraseProcess& node_erase, bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor )
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),node_erase, rem_nodes, add_nodes, alpha_shape, h_factor	);
}



///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				ADAPTIVE 2D MESHER specifically for DROPLET FORMING 	//
//											//
//////////////////////////////////////////////////////////////////////////////////////////

void TriRegenerateMeshDROPLET(TriGenDropletModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part,NodeEraseProcess& node_erase, bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor )
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),node_erase, rem_nodes, add_nodes, alpha_shape, h_factor	);
}

///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				BOUNDARY CONFORMANT 2D MESHER				//
//											//
//////////////////////////////////////////////////////////////////////////////////////////


void TriGenCDTFluid(TriGenCDT& Mesher,ModelPart& model_part)
{
    Mesher.GenerateCDT(model_part,
                       KratosComponents<Element>::Get("Fluid2D")
                      );
}

///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				ADAPTIVE 2D MESHER -->USING SEGMENT  		//
//											//
//////////////////////////////////////////////////////////////////////////////////////////

//trigen pfem refine segment
void TriRegenerateMeshWithSegment(TriGenPFEMRefineSegment& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part,NodeEraseProcess& node_erase, bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor )
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),node_erase, rem_nodes, add_nodes, alpha_shape, h_factor	);
}

///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				VMS 2D Mesher			                        //
//											//
//////////////////////////////////////////////////////////////////////////////////////////

//trigen pfem vms refine
void TriRegenerateMeshVMS(TriGenPFEMModelerVMS& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part,NodeEraseProcess& node_erase, bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor )
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),node_erase, rem_nodes, add_nodes, alpha_shape, h_factor	);
}



void  AddMeshersToPython(pybind11::module& m)
{
    
#ifdef USE_TETGEN_NONFREE_TPL
    //class that allows 3D adaptive remeshing (inserting and erasing nodes)
    class_<TetGenPfemModeler, TetGenPfemModeler::Pointer >(m, "TetGenPfemModeler")
    .def(init< >())
    .def("ReGenerateMesh",TetRegenerateMesh)
    .def("ReGenerateMesh",&TetGenPfemModeler::ReGenerateMesh)
    ;

    class_<TetGenPfemRefineFace, TetGenPfemRefineFace::Pointer >(m, "TetGenPfemRefineFace")
    .def(init< >())
    .def("ReGenerateMesh",TetRegenerateMeshWithFace)
    ;

    class_<TetGenPfemContact, TetGenPfemContact::Pointer >(m, "TetGenPfemContact")
    .def(init< >())
    .def("ReGenerateMesh",TetRegenerateMeshContact)
    ;
    
    class_<TetGenCDT, TetGenCDT::Pointer >(m, "TetGenCDT")
    .def(init< >())
    .def("GenerateCDT",GenerateCDT)
    ;

    class_<TetGenPfemModelerVms, TetGenPfemModelerVms::Pointer >(m, "TetGenPfemModelerVms")
    .def(init< >())
    .def("ReGenerateMesh",&TetGenPfemModelerVms::ReGenerateMesh)
    ;
#endif
    
    //class that allows 2D adaptive remeshing (inserting and erasing nodes)
    class_<TriGenPFEMModeler, TriGenPFEMModeler::Pointer >(m, "TriGenPFEMModeler")
    .def(init< >())
    .def("ReGenerateMesh",TriRegenerateMesh)
    .def("ReGenerateMesh",&TriGenPFEMModeler::ReGenerateMesh)
    ;

    //class that allows 2D adaptive remeshing (inserting and erasing nodes) as well as preserving the topology (avoiding "holes" at the boundaries). made for glass simulation
    class_<TriGenGLASSModeler, TriGenGLASSModeler::Pointer >(m, "TriGenGLASSModeler")
    .def(init< >())
    .def("ReGenerateMeshGlass",TriRegenerateMeshGLASS)
    .def("ReGenerateMeshGlass",&TriGenGLASSModeler::ReGenerateMesh)
    ;

    
    //class that allows 2D adaptive remeshing (inserting and erasing nodes) as well as preserving the topology (avoiding "holes" at the boundaries). made for droplet simulation
    class_<TriGenDropletModeler, TriGenDropletModeler::Pointer >(m,"TriGenDropletModeler")
    .def(init< >())
    .def("ReGenerateMeshDROPLET",TriRegenerateMeshDROPLET)
    .def("ReGenerateMeshDROPLET",&TriGenDropletModeler::ReGenerateMesh)
    ;


    class_<TriGenPFEMModelerVMS, TriGenPFEMModelerVMS::Pointer>(m, "TriGenPFEMModelerVMS")
    .def(init< >())
    .def("ReGenerateMesh",&TriGenPFEMModelerVMS::ReGenerateMesh);

    //segment mesher adaptive
    class_<TriGenPFEMRefineSegment, TriGenPFEMRefineSegment::Pointer >(m, "TriGenPFEMSegment")
    .def(init< >())
    .def("ReGenerateMesh",TriRegenerateMeshWithSegment)
    ;



}

}  // namespace Python.

} // Namespace Kratos

