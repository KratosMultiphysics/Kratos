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
#include "custom_utilities/NurbsShapeFunctionModeler.h"
#include "iga_structural_mechanics_application.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"


namespace Kratos
{
	/**
	 * @brief NurbsShapeFunctionModeler::GeneratePatchesAndElements
	 * This function reads in a NURBS-Geometry using an existing ModelPart and stores
	 * the NURBS-Object as a member variable
	 *
	 * @param NurbsModelPart: Contains all Control Points
	 * @param FirstControlPoint: The first Control Point of the NURBS-Patch
	 * @param ConnectivitiesEnd: The last Control Point of the NURBS-Patch
	 * @param PolynomialDegreeXi: Polynomial Degree in U/Xi direction
	 * @param PolynomialDegreeEta: Polynomial Degree in V/Eta direction
	 * @param KnotsXi: Knots in U/Xi direction
	 * @param KnotsEta: Knots in V/Eta direction
	 * @param NumberOfCpsU: Number of Control Points in U/Xi direction
	 * @param NumberOfCpsV: Number of Control Points in V/Eta direction
	 * @param Weights: Weights of the Control Points
	 */
	void NurbsShapeFunctionModeler::SetUp(ModelPart &NurbsModelPart,
		int dimension,
		Vector ControlPoints, //integer IDs
		Vector Weights,
		int PolynomialDegreeXi,
		int PolynomialDegreeEta,
		Vector KnotsXi,
		Vector KnotsEta)
	{
		// Define a Pointer vector to 3D nodes named ElementControlPoints
		PointerVector< Node<3> > ElementControlPoints;
		//Besser Resize!!
		// Append the vector PatchControlPoints with the CPs of the current patch
		for (int i = 0; i < ControlPoints.size(); i++)
		{
			ElementControlPoints.push_back(NurbsModelPart.pGetNode(ControlPoints[i]));
		}

		//std::cout << "Jetzt: " << ElementControlPoints.size() << std::endl;
		//for (int i = 0; i < ElementControlPoints.size(); i++)
		//{
		//	std::cout << "CP1: " << ElementControlPoints[i] << std::endl;
		//}
		
		// Define a pointer pNurbsElement
		NurbsPatchGeometry2D< Node<3> >::Pointer pNurbsSurface;


		// Create the NURBS-patch
		if (dimension == 2)
		{
			pNurbsSurface = NurbsPatchGeometry2D< Node<3> >::Pointer(
				new NurbsPatchGeometry2D< Node<3> >(
				ElementControlPoints,
				Weights,
				KnotsXi,
				KnotsEta,
				PolynomialDegreeXi,
				PolynomialDegreeEta,
				PolynomialDegreeXi+1,
				PolynomialDegreeEta+1,
				KnotsXi(PolynomialDegreeXi),
				KnotsXi(PolynomialDegreeXi+1),
				KnotsEta(PolynomialDegreeEta),
				KnotsEta(PolynomialDegreeEta+1)
				));
		}

		else if (dimension == 3)
		{

		}

		else
		{
			std::cout << "Error in NurbsModeler::GeneratePatchesAndElements() : This implementation only support dimension 2 or 3." << std::endl;
			return;
		}
		// mNurbsPatchGeometry is a vector of pointers
		mNurbsPatchGeometry=pNurbsSurface;
	}

	void NurbsShapeFunctionModeler::EvaluateShapeFunction(
		const array_1d<double, 3>& LocalCoordinatesOfQuadraturePoint,
		Vector &ShapeFunctionValues,
		Matrix &ShapeFunctionLocalDerivatives) 
	{
		mNurbsPatchGeometry->ShapeFunctionsValues(ShapeFunctionValues, LocalCoordinatesOfQuadraturePoint);
		mNurbsPatchGeometry->ShapeFunctionsDerivativesValues(ShapeFunctionLocalDerivatives, LocalCoordinatesOfQuadraturePoint);
	}

	void NurbsShapeFunctionModeler::EvaluateShapeFunctionSecondOrder(
		const array_1d<double, 3>& LocalCoordinatesOfQuadraturePoint,
		Vector &ShapeFunctionValues,
		Matrix &ShapeFunctionLocalDerivatives)
	{
		mNurbsPatchGeometry->ShapeFunctionsValues(ShapeFunctionValues, LocalCoordinatesOfQuadraturePoint);
		mNurbsPatchGeometry->ShapeFunctionsSecondDerivativesValues(ShapeFunctionLocalDerivatives, LocalCoordinatesOfQuadraturePoint);
	}

NurbsShapeFunctionModeler::NurbsShapeFunctionModeler()
{
	NurbsPatchGeometry2D< Node<3> > mNurbsPatchGeometry();
}
    // vector of pointers for 3D nodes to mNurbsPatchGeometry
    //NurbsPatchGeometry2D< Node<3> > mNurbsPatchGeometry();

    //std::cout <<"The empty constructor of the class NurbsModeler has been called called"<<std::endl;


NurbsShapeFunctionModeler::~NurbsShapeFunctionModeler()
{}

}  // namespace Kratos.


