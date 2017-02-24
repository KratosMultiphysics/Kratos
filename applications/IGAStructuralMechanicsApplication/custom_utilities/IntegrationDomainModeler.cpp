//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//

// System includes


// External includes 


// Project includes
#include "custom_utilities/IntegrationDomainModeler.h"
#include "iga_structural_mechanics_application.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"


namespace Kratos
{
	void IntegrationDomainModeler::SetUp(const String &GeometryData)
	{
		
	}

IntegrationDomainModeler::IntegrationDomainModeler()
{
	//NurbsPatchGeometry2D< Node<3> > mNurbsPatchGeometry();
}
    // vector of pointers for 3D nodes to mNurbsPatchGeometry
    //NurbsPatchGeometry2D< Node<3> > mNurbsPatchGeometry();

    //std::cout <<"The empty constructor of the class NurbsModeler has been called called"<<std::endl;


IntegrationDomainModeler::~IntegrationDomainModeler()
{}

}  // namespace Kratos.


