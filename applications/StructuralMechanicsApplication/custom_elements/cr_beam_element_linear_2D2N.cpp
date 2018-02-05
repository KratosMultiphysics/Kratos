// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license: 	 structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//                   
//                   
//
// System includes

// External includes

// Project includes
#include "custom_elements/cr_beam_element_linear_2D2N.hpp"
#include "custom_utilities/static_condensation_utility.h"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"


namespace Kratos
{

	CrBeamElementLinear2D2N::CrBeamElementLinear2D2N(IndexType NewId,
		GeometryType::Pointer pGeometry)
		: CrBeamElement2D2N(NewId, pGeometry)
	{}

	CrBeamElementLinear2D2N::CrBeamElementLinear2D2N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties)
		: CrBeamElement2D2N(NewId, pGeometry, pProperties)
	{}

	Element::Pointer CrBeamElementLinear2D2N::Create(IndexType NewId,
		NodesArrayType const& rThisNodes,
		PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return Kratos::make_shared<CrBeamElementLinear2D2N>(
            NewId, rGeom.Create(rThisNodes), pProperties);
	}

	CrBeamElementLinear2D2N::~CrBeamElementLinear2D2N() {}

	

	void CrBeamElementLinear2D2N::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
		{
			KRATOS_TRY;
			// Kt
			this->CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
				
			Vector NodalDeformation = ZeroVector(msElementSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(msElementSize);
			rRightHandSideVector -= prod(rLeftHandSideMatrix, NodalDeformation);
	
			rRightHandSideVector += this->CalculateBodyForces();

			KRATOS_CATCH("")
		}

	void CrBeamElementLinear2D2N::CalculateRightHandSide(VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
		{
			KRATOS_TRY;
			Vector NodalDeformation = ZeroVector(msElementSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(msElementSize);
			rRightHandSideVector -= prod(this->K_master, NodalDeformation);



			rRightHandSideVector += this->CalculateBodyForces();
			KRATOS_CATCH("")
		}

	void CrBeamElementLinear2D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
		ProcessInfo& rCurrentProcessInfo) 
		{
			KRATOS_TRY;
			rLeftHandSideMatrix = this->CreateElementStiffnessMatrix_Total();


			//// start static condensation
			if (this->Has(CONDENSED_DOF_LIST))
			{
				Vector dof_list_input = this->GetValue(CONDENSED_DOF_LIST);
				std::vector<int> dofList(dof_list_input.size());
				for (size_t i=0;i<dof_list_input.size();++i) dofList[i]=dof_list_input[i];
				StaticCondensationUtility::CondenseLeftHandSide(*this,rLeftHandSideMatrix,dofList);
			}
			//// end static condensation



			this->GlobalizeMatrix(rLeftHandSideMatrix);
			this->K_master = rLeftHandSideMatrix;
			KRATOS_CATCH("")
		}

	/////////////////////////////////////////////////
	///////////// CUSTOM FUNCTIONS --->>
	/////////////////////////////////////////////////


	bounded_matrix<double,CrBeamElement2D2N::msElementSize,CrBeamElement2D2N::msElementSize>
	 CrBeamElementLinear2D2N::CreateRotationMatrix()
	{
		KRATOS_TRY;
		const double current_element_angle = this->CalculateInitialElementAngle();
		const double c = std::cos(current_element_angle);
		const double s = std::sin(current_element_angle);

		bounded_matrix<double,msElementSize,msElementSize>
		 RotationMatrix = ZeroMatrix(msElementSize,msElementSize);

		RotationMatrix(0, 0) = c;
		RotationMatrix(0, 1) = -s;
		RotationMatrix(1, 0) = s;
		RotationMatrix(1, 1) = c;
		RotationMatrix(2, 2) = 1.00;

		RotationMatrix(3, 3) = c;
		RotationMatrix(3, 4) = -s;
		RotationMatrix(4, 3) = s;
		RotationMatrix(4, 4) = c;
		RotationMatrix(5, 5) = 1.00;

		return RotationMatrix;
		KRATOS_CATCH("")
	}


	void CrBeamElementLinear2D2N::CalculateOnIntegrationPoints(
		const Variable<array_1d<double, 3 > >& rVariable,
		std::vector< array_1d<double, 3 > >& rOutput,
		const ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
		//element with two nodes can only represent results at one node 
		const unsigned int&  write_points_number = GetGeometry()
			.IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
		if (rOutput.size() != write_points_number) {
			rOutput.resize(write_points_number);
		}

		Matrix LeftHandSideMatrix = this->CreateElementStiffnessMatrix_Total();

		Vector NodalDeformation = ZeroVector(msElementSize);
		this->GetValuesVector(NodalDeformation);

		bounded_matrix<double,msElementSize,msElementSize> TransformationMatrix = this->CreateRotationMatrix();
		// calculate local displacements
		NodalDeformation = prod(Matrix(trans(TransformationMatrix)),NodalDeformation);
		

		//// start static back condensation
		if (this->Has(CONDENSED_DOF_LIST))
		{
			Vector dof_list_input = this->GetValue(CONDENSED_DOF_LIST);
			std::vector<int> dofList(dof_list_input.size());
			for (size_t i=0;i<dof_list_input.size();++i) dofList[i]=dof_list_input[i];
			Vector nodal_deformation_temp = NodalDeformation;
			StaticCondensationUtility::ConvertingCondensation(
				*this,nodal_deformation_temp,NodalDeformation,dofList,LeftHandSideMatrix);
		}
		//// end static back condensation

		Vector Stress = prod(LeftHandSideMatrix, NodalDeformation); 

		KRATOS_WATCH(this->Id());
		KRATOS_WATCH(NodalDeformation);
		KRATOS_WATCH(Stress);
		std::cout << "-------------------------------------------" << std::endl;


		//rOutput[GP 1,2,3][x,y,z]

		if (rVariable == MOMENT)
		{
			rOutput[0][0] = 0.00;
			rOutput[1][0] = 0.00;
			rOutput[2][0] = 0.00;

			rOutput[0][1] = 0.00;
			rOutput[1][1] = 0.00;
			rOutput[2][1] = 0.00;

			rOutput[0][2] = 1.0 *Stress[2] * 0.75 - Stress[5] * 0.25;
			rOutput[1][2] = 1.0 *Stress[2] * 0.50 - Stress[5] * 0.50;
			rOutput[2][2] = 1.0 *Stress[2] * 0.25 - Stress[5] * 0.75;

		}
		if (rVariable == FORCE)
		{
			rOutput[0][0] = -1.0 * Stress[0] * 0.75 + Stress[3] * 0.25;
			rOutput[1][0] = -1.0 * Stress[0] * 0.50 + Stress[3] * 0.50;
			rOutput[2][0] = -1.0 * Stress[0] * 0.25 + Stress[3] * 0.75;

			rOutput[0][1] = -1.0 * Stress[1] * 0.75 + Stress[4] * 0.25;
			rOutput[1][1] = -1.0 *Stress[1] * 0.50 + Stress[4] * 0.50;
			rOutput[2][1] = -1.0 *Stress[1] * 0.25 + Stress[4] * 0.75;

			rOutput[0][2] = 0.00;
			rOutput[1][2] = 0.00;
			rOutput[2][2] = 0.00;

		}

		KRATOS_CATCH("")
	}

	void CrBeamElementLinear2D2N::GetValueOnIntegrationPoints(
		const Variable<array_1d<double, 3 > >& rVariable,
		std::vector< array_1d<double, 3 > >& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		this->CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
		KRATOS_CATCH("")
	}


	double CrBeamElementLinear2D2N::CalculateLength() 
	{
		KRATOS_TRY;
		return this->CalculateReferenceLength();
		KRATOS_CATCH("")
	}




	void CrBeamElementLinear2D2N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		rSerializer.save("MasterStiffnessMatrix", this->K_master);
	}

	void CrBeamElementLinear2D2N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
		rSerializer.load("MasterStiffnessMatrix", this->K_master);
	}

} // namespace Kratos.


