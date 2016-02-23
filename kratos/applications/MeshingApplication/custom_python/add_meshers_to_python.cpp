//
//   Project Name:        Kratos
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-22 16:37:10 $
//   Revision:            $Revision: 1.13 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_meshers_to_python.h"
#include "external_includes/tetgen_pfem_refine.h"
#include "external_includes/tetgen_pfem_refine_vms.h"
#include "external_includes/tetgen_pfem_refine_face.h"
#include "external_includes/trigen_pfem_refine.h"
#include "external_includes/trigen_pfem_refine_vms.h"
#include "external_includes/trigen_pfem_refine_segment.h"
#include "external_includes/tetgen_pfem_contact.h"
#include "external_includes/tetgen_cdt.h"

//#include "external_includes/trigen_mesh_suite.h"
#include "external_includes/trigen_cdt.h"
//#include "external_includes/trigen_refine.h"

// #include "external_includes/msuite_pfem_refine.h"


namespace Kratos
{

namespace Python
{
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
//				BOUNDARY CONFORMANT 2D MESHER				//
//											//
//////////////////////////////////////////////////////////////////////////////////////////


void TriGenCDTFluid(TriGenCDT& Mesher,ModelPart& model_part)
{
    Mesher.GenerateCDT(model_part,
                       KratosComponents<Element>::Get("Fluid2D")
                      );
}
/*
void RefineCDT(TriGenCDTrefine& Mesher,ModelPart& model_part)
{
	Mesher.RefineCDT(model_part, true,
			   KratosComponents<Element>::Get("Fluid2D")
			  );
}

void QualityCDT(TriGenCDTrefine& Mesher,ModelPart& model_part)
{
	Mesher.RefineCDT(model_part, false,
			 KratosComponents<Element>::Get("Fluid2D")
			);
}
*/
///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				ADAPTIVE 2D MESHER -->USING MESH SUITE  		//
//											//
//////////////////////////////////////////////////////////////////////////////////////////

//trigen pfem refine
// void MsuiteRegenerateMesh(MSuitePFEMModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part,NodeEraseProcess& node_erase, bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor )
// {
//     Mesher.ReGenerateMesh(model_part,
//                           KratosComponents<Element>::Get(ElementName),
//                           KratosComponents<Condition>::Get(ConditionName),node_erase, rem_nodes, add_nodes, alpha_shape, h_factor	);
// }


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



void  AddMeshersToPython()
{

    using namespace boost::python;
    //class that allows 3D adaptive remeshing (inserting and erasing nodes)
    class_<TetGenPfemModeler >("TetGenPfemModeler",
                               init< >())
    .def("ReGenerateMesh",TetRegenerateMesh)
    .def("ReGenerateMesh",&TetGenPfemModeler::ReGenerateMesh)
    ;

    class_<TetGenPfemRefineFace >("TetGenPfemRefineFace",
                                  init< >())
    .def("ReGenerateMesh",TetRegenerateMeshWithFace)
    ;

    class_<TetGenPfemContact >("TetGenPfemContact",
                               init< >())
    .def("ReGenerateMesh",TetRegenerateMeshContact)
    ;
    
    class_<TetGenCDT >("TetGenCDT",
                               init< >())
    .def("GenerateCDT",GenerateCDT)
    ;

    class_<TetGenPfemModelerVms >("TetGenPfemModelerVms",
                                  init< >())
    .def("ReGenerateMesh",&TetGenPfemModelerVms::ReGenerateMesh)
    ;
    
    //class that allows 2D adaptive remeshing (inserting and erasing nodes)
    class_<TriGenPFEMModeler >("TriGenPFEMModeler",
                               init< >())
    .def("ReGenerateMesh",TriRegenerateMesh)
    .def("ReGenerateMesh",&TriGenPFEMModeler::ReGenerateMesh)
    ;

    class_<TriGenPFEMModelerVMS>("TriGenPFEMModelerVMS",
                                 init< >())
            .def("ReGenerateMesh",&TriGenPFEMModelerVMS::ReGenerateMesh);
    /*
     class_<TriGenModeler >("TriGenModeler",
    	 init< >())
     .def("ReGenerateMesh",TriRegenerateMesh)
    	;
    */
    /*
      class_<TriGenCDT >("TriGenCDT",
    	init< >())
      .def("GenerateFluidElements",TriGenCDTFluid)
      ;

      class_<TriGenCDTrefine >("TriGenCDTrefine",
    		     init< >())
    		  .def("RefineCDT",RefineCDT)
    		  .def("QualityCDT",QualityCDT)
    		  ;
          */
    //class that allows 2D adaptive remeshing (inserting and erasing nodes)
//     class_<MSuitePFEMModeler >("MSuitePFEMModeler",
//                                init< >())
//     .def("ReGenerateMesh",MsuiteRegenerateMesh);
    //segment mesher adaptive
    class_<TriGenPFEMRefineSegment >("TriGenPFEMSegment",
                                     init< >())
    .def("ReGenerateMesh",TriRegenerateMeshWithSegment)
    ;



}

}  // namespace Python.

} // Namespace Kratos

