//
//   Project Name:        KratosPFEMFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:                   June 2018 $
//   Revision:            $Revision:                     0.0 $
//
//

#ifndef KRATOS_NODAL_TWO_STEP_V_P_STRATEGY_H
#define KRATOS_NODAL_TWO_STEP_V_P_STRATEGY_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/builders_and_solvers/nodal_residualbased_elimination_builder_and_solver.h"
#include "custom_strategies/builders_and_solvers/nodal_residualbased_elimination_builder_and_solver_continuity.h"
#include "custom_strategies/builders_and_solvers/nodal_residualbased_block_builder_and_solver.h"

#include "custom_utilities/solver_settings.h"

#include "custom_strategies/strategies/gauss_seidel_linear_strategy.h"

#include "pfem_fluid_dynamics_application_variables.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>

namespace Kratos
{

///@addtogroup PFEMFluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template <class TSparseSpace,
		  class TDenseSpace,
		  class TLinearSolver>
class NodalTwoStepVPStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
	///@name Type Definitions
	///@{
	KRATOS_CLASS_POINTER_DEFINITION(NodalTwoStepVPStrategy);

	/// Counted pointer of NodalTwoStepVPStrategy
	//typedef boost::shared_ptr< NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver> > Pointer;

	typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

	typedef typename BaseType::TDataType TDataType;

	/// Node type (default is: Node<3>)
	typedef Node<3> NodeType;

	/// Geometry type (using with given NodeType)
	typedef Geometry<NodeType> GeometryType;

	typedef std::size_t SizeType;

	//typedef typename BaseType::DofSetType DofSetType;

	typedef typename BaseType::DofsArrayType DofsArrayType;

	typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

	typedef typename BaseType::TSystemVectorType TSystemVectorType;

	typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

	typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

	typedef typename BaseType::ElementsArrayType ElementsArrayType;

	typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

	typedef TwoStepVPSolverSettings<TSparseSpace, TDenseSpace, TLinearSolver> SolverSettingsType;

	typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

	typedef GlobalPointersVector<Node<3>> NodeWeakPtrVectorType;
	///@}
	///@name Life Cycle
	///@{

	NodalTwoStepVPStrategy(ModelPart &rModelPart,
						   SolverSettingsType &rSolverConfig) : BaseType(rModelPart)
	{
		InitializeStrategy(rSolverConfig);
	}

	NodalTwoStepVPStrategy(ModelPart &rModelPart,
						   /*SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& rSolverConfig,*/
						   typename TLinearSolver::Pointer pVelocityLinearSolver,
						   typename TLinearSolver::Pointer pPressureLinearSolver,
						   bool ReformDofSet = true,
						   double VelTol = 0.0001,
						   double PresTol = 0.0001,
						   int MaxPressureIterations = 1, // Only for predictor-corrector
						   unsigned int TimeOrder = 2,
						   unsigned int DomainSize = 2) : BaseType(rModelPart), // Move Mesh flag, pass as input?
														  mVelocityTolerance(VelTol),
														  mPressureTolerance(PresTol),
														  mMaxPressureIter(MaxPressureIterations),
														  mDomainSize(DomainSize),
														  mTimeOrder(TimeOrder),
														  mReformDofSet(ReformDofSet)
	{
		KRATOS_TRY;

		BaseType::SetEchoLevel(1);

		// Check that input parameters are reasonable and sufficient.
		this->Check();

		bool CalculateNormDxFlag = true;

		bool ReformDofAtEachIteration = false; // DofSet modifiaction is managed by the fractional step strategy, auxiliary strategies should not modify the DofSet directly.

		// Additional Typedefs
		//typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3 > > > VarComponent;
		typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
		typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

		//initializing fractional velocity solution step
		typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;
		typename SchemeType::Pointer pScheme;

		typename SchemeType::Pointer Temp = typename SchemeType::Pointer(new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>());
		/* typename SchemeType::Pointer Temp = typename SchemeType::Pointer(new IncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ()); */
		pScheme.swap(Temp);

		//CONSTRUCTION OF VELOCITY
		BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new NodalResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pVelocityLinearSolver));
		/* BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (pVelocityLinearSolver)); */

		this->mpMomentumStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pVelocityLinearSolver, vel_build, ReformDofAtEachIteration, CalculateNormDxFlag));

		this->mpMomentumStrategy->SetEchoLevel(BaseType::GetEchoLevel());

		vel_build->SetCalculateReactionsFlag(false);

		/* BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pPressureLinearSolver, PRESSURE)); */
		/* BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver >(pPressureLinearSolver)); */
		BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new NodalResidualBasedEliminationBuilderAndSolverContinuity<TSparseSpace, TDenseSpace, TLinearSolver>(pPressureLinearSolver));
		/* BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new NodalResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver >(pPressureLinearSolver)); */

		this->mpPressureStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pPressureLinearSolver, pressure_build, ReformDofAtEachIteration, CalculateNormDxFlag));

		this->mpPressureStrategy->SetEchoLevel(BaseType::GetEchoLevel());

		pressure_build->SetCalculateReactionsFlag(false);

		KRATOS_CATCH("");
	}

	/// Destructor.
	virtual ~NodalTwoStepVPStrategy() {}

	int Check() override
	{
		KRATOS_TRY;

		// Check elements and conditions in the model part
		int ierr = BaseType::Check();
		if (ierr != 0)
			return ierr;

		if (DELTA_TIME.Key() == 0)
			KRATOS_THROW_ERROR(std::runtime_error, "DELTA_TIME Key is 0. Check that the application was correctly registered.", "");
		if (BDF_COEFFICIENTS.Key() == 0)
			KRATOS_THROW_ERROR(std::runtime_error, "BDF_COEFFICIENTS Key is 0. Check that the application was correctly registered.", "");

		ModelPart &rModelPart = BaseType::GetModelPart();

		if (mTimeOrder == 2 && rModelPart.GetBufferSize() < 3)
			KRATOS_THROW_ERROR(std::invalid_argument, "Buffer size too small for fractional step strategy (BDF2), needed 3, got ", rModelPart.GetBufferSize());
		if (mTimeOrder == 1 && rModelPart.GetBufferSize() < 2)
			KRATOS_THROW_ERROR(std::invalid_argument, "Buffer size too small for fractional step strategy (Backward Euler), needed 2, got ", rModelPart.GetBufferSize());

		const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();

		for (ModelPart::ElementIterator itEl = rModelPart.ElementsBegin(); itEl != rModelPart.ElementsEnd(); ++itEl)
		{
			ierr = itEl->Check(rCurrentProcessInfo);
			if (ierr != 0)
				break;
		}

		for (ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond)
		{
			ierr = itCond->Check(rCurrentProcessInfo);
			if (ierr != 0)
				break;
		}

		return ierr;

		KRATOS_CATCH("");
	}

	double Solve() override
	{
		// Initialize BDF2 coefficients
		ModelPart &rModelPart = BaseType::GetModelPart();
		this->SetTimeCoefficients(rModelPart.GetProcessInfo());
		double NormDp = 0.0;
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		double currentTime = rCurrentProcessInfo[TIME];
		double timeInterval = rCurrentProcessInfo[DELTA_TIME];
		bool timeIntervalChanged = rCurrentProcessInfo[TIME_INTERVAL_CHANGED];

		// bool momentumAlreadyConverged=false;
		// bool continuityAlreadyConverged=false;

		unsigned int maxNonLinearIterations = mMaxPressureIter;
		std::cout << "\n                   Solve with nodally_integrated_two_step_vp strategy at t=" << currentTime << "s" << std::endl;

		if (timeIntervalChanged == true && currentTime > 10 * timeInterval)
		{
			maxNonLinearIterations *= 2;
		}
		if (currentTime < 10 * timeInterval)
		{
			if (BaseType::GetEchoLevel() > 1)
				std::cout << "within the first 10 time steps, I consider the given iteration number x3" << std::endl;
			maxNonLinearIterations *= 3;
		}
		if (currentTime < 20 * timeInterval && currentTime >= 10 * timeInterval)
		{
			if (BaseType::GetEchoLevel() > 1)
				std::cout << "within the second 10 time steps, I consider the given iteration number x2" << std::endl;
			maxNonLinearIterations *= 2;
		}

		bool momentumConverged = true;
		bool continuityConverged = false;
		bool fixedTimeStep = false;
		/* boost::timer solve_step_time; */

		this->UnactiveSliverElements();

		this->InitializeSolutionStep();
		for (unsigned int it = 0; it < maxNonLinearIterations; ++it)
		{
			if (BaseType::GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
				std::cout << "----- > iteration: " << it << std::endl;

			if (it == 0)
			{

				this->ComputeNodalVolume();

				this->InitializeNonLinearIterations();
			}
			this->CalcNodalStrainsAndStresses();

			momentumConverged = this->SolveMomentumIteration(it, maxNonLinearIterations, fixedTimeStep);

			this->UpdateTopology(rModelPart, BaseType::GetEchoLevel());
			this->ComputeNodalVolume();
			this->InitializeNonLinearIterations();
			this->CalcNodalStrains();

			if (fixedTimeStep == false)
			{
				continuityConverged = this->SolveContinuityIteration(it, maxNonLinearIterations);
			}

			// if((momentumConverged==true || it==maxNonLinearIterations-1) && momentumAlreadyConverged==false){
			// 	std::ofstream myfile;
			// 	myfile.open ("momentumConvergedIteration.txt",std::ios::app);
			// 	myfile << currentTime << "\t" << it << "\n";
			// 	myfile.close();
			// 	momentumAlreadyConverged=true;
			// }
			// if((continuityConverged==true || it==maxNonLinearIterations-1) && continuityAlreadyConverged==false){
			// 	std::ofstream myfile;
			// 	myfile.open ("continuityConvergedIteration.txt",std::ios::app);
			// 	myfile << currentTime << "\t" << it << "\n";
			// 	myfile.close();
			// 	continuityAlreadyConverged=true;
			// }

			if (it == maxNonLinearIterations - 1 || ((continuityConverged && momentumConverged) && it > 1))
			{
				//this->ComputeErrorL2NormCaseImposedG();
				//this->ComputeErrorL2NormCasePoiseuille();
				this->CalculateAccelerations();
				// std::ofstream myfile;
				// myfile.open ("maxConvergedIteration.txt",std::ios::app);
				// myfile << currentTime << "\t" << it << "\n";
				// myfile.close();
			}
			if ((continuityConverged && momentumConverged) && it > 1)
			{
				rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, false);
				rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, false);
				std::cout << "nodal V-P strategy converged in " << it + 1 << " iterations." << std::endl;
				break;
			}
		}

		if (!continuityConverged && !momentumConverged && BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
			std::cout << "Convergence tolerance not reached." << std::endl;

		if (mReformDofSet)
			this->Clear();

		/* std::cout << "solve_step_time : " << solve_step_time.elapsed() << std::endl; */

		return NormDp;
	}

	void FinalizeSolutionStep() override
	{
		/* this->UpdateStressStrain(); */
	}

	void Initialize() override
	{

		std::cout << "                                 Initialize in nodal_two_step_v_p_strategy" << std::endl;
		ModelPart &rModelPart = BaseType::GetModelPart();
		const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
		unsigned int sizeStrains = 3 * (dimension - 1);

		// #pragma omp parallel
		// 	{
		ModelPart::NodeIterator NodesBegin;
		ModelPart::NodeIterator NodesEnd;
		OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

		for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
		{
			NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
			unsigned int neighbourNodes = neighb_nodes.size();
			unsigned int sizeSDFNeigh = neighbourNodes * dimension;

			if (itNode->SolutionStepsDataHas(NODAL_CAUCHY_STRESS))
			{
				Vector &rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS);
				if (rNodalStress.size() != sizeStrains)
				{
					rNodalStress.resize(sizeStrains, false);
				}
				noalias(rNodalStress) = ZeroVector(sizeStrains);
			}
			else
			{
				std::cout << "THIS node does not have NODAL_CAUCHY_STRESS... " << itNode->X() << " " << itNode->Y() << std::endl;
			}

			if (itNode->SolutionStepsDataHas(NODAL_DEVIATORIC_CAUCHY_STRESS))
			{
				Vector &rNodalStress = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS);
				if (rNodalStress.size() != sizeStrains)
				{
					rNodalStress.resize(sizeStrains, false);
				}
				noalias(rNodalStress) = ZeroVector(sizeStrains);
			}
			else
			{
				std::cout << "THIS node does not have NODAL_DEVIATORIC_CAUCHY_STRESS... " << itNode->X() << " " << itNode->Y() << std::endl;
			}

			if (itNode->SolutionStepsDataHas(NODAL_VOLUME))
			{
				itNode->FastGetSolutionStepValue(NODAL_VOLUME) = 0;
			}
			else
			{
				std::cout << "THIS node does not have NODAL_VOLUME... " << itNode->X() << " " << itNode->Y() << std::endl;
			}

			if (itNode->SolutionStepsDataHas(NODAL_MEAN_MESH_SIZE))
			{
				itNode->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE) = 0;
			}
			else
			{
				std::cout << "THIS node does not have NODAL_MEAN_MESH_SIZE... " << itNode->X() << " " << itNode->Y() << std::endl;
			}

			if (itNode->SolutionStepsDataHas(NODAL_FREESURFACE_AREA))
			{
				itNode->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) = 0;
			}
			else
			{
				std::cout << "THIS node does not have NODAL_FREESURFACE_AREA... " << itNode->X() << " " << itNode->Y() << std::endl;
			}

			if (itNode->SolutionStepsDataHas(NODAL_SFD_NEIGHBOURS))
			{
				Vector &rNodalSFDneighbours = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
				if (rNodalSFDneighbours.size() != sizeSDFNeigh)
				{
					rNodalSFDneighbours.resize(sizeSDFNeigh, false);
				}
				noalias(rNodalSFDneighbours) = ZeroVector(sizeSDFNeigh);
			}
			else
			{
				std::cout << "THIS node does not have NODAL_SFD_NEIGHBOURS... " << itNode->X() << " " << itNode->Y() << std::endl;
			}

			if (itNode->SolutionStepsDataHas(NODAL_SPATIAL_DEF_RATE))
			{
				Vector &rSpatialDefRate = itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
				if (rSpatialDefRate.size() != sizeStrains)
				{
					rSpatialDefRate.resize(sizeStrains, false);
				}
				noalias(rSpatialDefRate) = ZeroVector(sizeStrains);
			}
			else
			{
				std::cout << "THIS node does not have NODAL_SPATIAL_DEF_RATE... " << itNode->X() << " " << itNode->Y() << std::endl;
			}

			if (itNode->SolutionStepsDataHas(NODAL_DEFORMATION_GRAD))
			{
				Matrix &rFgrad = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
				if (rFgrad.size1() != dimension)
				{
					rFgrad.resize(dimension, dimension, false);
				}
				noalias(rFgrad) = ZeroMatrix(dimension, dimension);
			}
			else
			{
				std::cout << "THIS node does not have NODAL_DEFORMATION_GRAD... " << itNode->X() << " " << itNode->Y() << std::endl;
			}

			if (itNode->SolutionStepsDataHas(NODAL_DEFORMATION_GRAD_VEL))
			{
				Matrix &rFgradVel = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
				if (rFgradVel.size1() != dimension)
				{
					rFgradVel.resize(dimension, dimension, false);
				}
				noalias(rFgradVel) = ZeroMatrix(dimension, dimension);
			}
			else
			{
				std::cout << "THIS node does not have NODAL_DEFORMATION_GRAD_VEL... " << itNode->X() << " " << itNode->Y() << std::endl;
			}

			this->AssignFluidMaterialToEachNode(itNode);
		}

		// }
	}

	void UnactiveSliverElements()
	{
		KRATOS_TRY;

		ModelPart &rModelPart = BaseType::GetModelPart();
		const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
		MesherUtilities MesherUtils;
		double ModelPartVolume = MesherUtils.ComputeModelPartVolume(rModelPart);
		double CriticalVolume = 0.001 * ModelPartVolume / double(rModelPart.Elements().size());
		double ElementalVolume = 0;

#pragma omp parallel
		{
			ModelPart::ElementIterator ElemBegin;
			ModelPart::ElementIterator ElemEnd;
			OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
			for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
			{
				unsigned int numNodes = itElem->GetGeometry().size();
				if (numNodes == (dimension + 1))
				{
					if (dimension == 2)
					{
						ElementalVolume = (itElem)->GetGeometry().Area();
					}
					else if (dimension == 3)
					{
						ElementalVolume = (itElem)->GetGeometry().Volume();
					}

					if (ElementalVolume < CriticalVolume)
					{
						// std::cout << "sliver element: it has Volume: " << ElementalVolume << " vs CriticalVolume(meanVol/1000): " << CriticalVolume<< std::endl;
						(itElem)->Set(ACTIVE, false);
					}
					else
					{
						(itElem)->Set(ACTIVE, true);
					}
				}
			}
		}
		KRATOS_CATCH("");
	}

	void AssignFluidMaterialToEachNode(ModelPart::NodeIterator itNode)
	{

		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		const double timeInterval = rCurrentProcessInfo[DELTA_TIME];

		double deviatoricCoeff = itNode->FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
		double volumetricCoeff = timeInterval * itNode->FastGetSolutionStepValue(BULK_MODULUS);

		double currFirstLame = volumetricCoeff - 2.0 * deviatoricCoeff / 3.0;

		itNode->FastGetSolutionStepValue(VOLUMETRIC_COEFFICIENT) = currFirstLame;
		itNode->FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT) = deviatoricCoeff;
	}

	void ComputeNodalVolume()
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		ElementsArrayType &pElements = rModelPart.Elements();
		const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

#ifdef _OPENMP
		int number_of_threads = omp_get_max_threads();
#else
		int number_of_threads = 1;
#endif

		vector<unsigned int> element_partition;
		OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

		//  #pragma omp parallel
		// 	    {
		int k = OpenMPUtils::ThisThread();
		typename ElementsArrayType::iterator ElemBegin = pElements.begin() + element_partition[k];
		typename ElementsArrayType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

		for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++) //MSI: To be parallelized
		{
			Element::GeometryType &geometry = itElem->GetGeometry();
			double elementalVolume = 0;

			if (dimension == 2)
			{
				elementalVolume = geometry.Area() / 3.0;
			}
			else if (dimension == 3)
			{
				elementalVolume = geometry.Volume() * 0.25;
			}
			// index = 0;
			unsigned int numNodes = geometry.size();
			for (unsigned int i = 0; i < numNodes; i++)
			{
				double &nodalVolume = geometry(i)->FastGetSolutionStepValue(NODAL_VOLUME);
				nodalVolume += elementalVolume;
			}
		}
		// }
	}

	void InitializeSolutionStep() override
	{
		this->FillNodalSFDVector();
	}

	void FillNodalSFDVector()
	{

		ModelPart &rModelPart = BaseType::GetModelPart();

		//  #pragma omp parallel
		//  	{
		// 		ModelPart::NodeIterator NodesBegin;
		// 		ModelPart::NodeIterator NodesEnd;
		// 		OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

		// for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
		// 	{

		for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
		{
			InitializeNodalVariablesForRemeshedDomain(itNode);

			SetNeighboursOrderToNode(itNode); // it assigns neighbours to inner nodes, filling NODAL_SFD_NEIGHBOURS_ORDER
		}
	}

	void SetNeighboursOrderToNode(ModelPart::NodeIterator itNode)
	{
		NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
		unsigned int neighbourNodes = neighb_nodes.size() + 1; // +1 becausealso the node itself must be considered as nieghbor node
		Vector &rNodeOrderedNeighbours = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);

		if (rNodeOrderedNeighbours.size() != neighbourNodes)
			rNodeOrderedNeighbours.resize(neighbourNodes, false);

		noalias(rNodeOrderedNeighbours) = ZeroVector(neighbourNodes);

		rNodeOrderedNeighbours[0] = itNode->Id();

		if (neighbourNodes > 1)
		{
			for (unsigned int k = 0; k < neighbourNodes - 1; k++)
			{
				rNodeOrderedNeighbours[k + 1] = neighb_nodes[k].Id();
			}
		}
	}

	void InitializeNodalVariablesForRemeshedDomain(ModelPart::NodeIterator itNode)
	{

		ModelPart &rModelPart = BaseType::GetModelPart();
		const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
		unsigned int sizeStrains = 3 * (dimension - 1);
		NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
		unsigned int neighbourNodes = neighb_nodes.size() + 1;
		unsigned int sizeSDFNeigh = neighbourNodes * dimension;

		if (itNode->SolutionStepsDataHas(NODAL_CAUCHY_STRESS))
		{
			Vector &rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS);
			if (rNodalStress.size() != sizeStrains)
				rNodalStress.resize(sizeStrains, false);
			noalias(rNodalStress) = ZeroVector(sizeStrains);
		}
		if (itNode->SolutionStepsDataHas(NODAL_DEVIATORIC_CAUCHY_STRESS))
		{
			Vector &rNodalDevStress = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS);
			if (rNodalDevStress.size() != sizeStrains)
				rNodalDevStress.resize(sizeStrains, false);
			noalias(rNodalDevStress) = ZeroVector(sizeStrains);
		}
		if (itNode->SolutionStepsDataHas(NODAL_SFD_NEIGHBOURS_ORDER))
		{
			Vector &rNodalSFDneighboursOrder = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
			if (rNodalSFDneighboursOrder.size() != neighbourNodes)
				rNodalSFDneighboursOrder.resize(neighbourNodes, false);
			noalias(rNodalSFDneighboursOrder) = ZeroVector(neighbourNodes);
		}
		if (itNode->SolutionStepsDataHas(NODAL_SFD_NEIGHBOURS))
		{
			Vector &rNodalSFDneighbours = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
			if (rNodalSFDneighbours.size() != sizeSDFNeigh)
				rNodalSFDneighbours.resize(sizeSDFNeigh, false);
			noalias(rNodalSFDneighbours) = ZeroVector(sizeSDFNeigh);
		}
		if (itNode->SolutionStepsDataHas(NODAL_SPATIAL_DEF_RATE))
		{
			Vector &rSpatialDefRate = itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
			if (rSpatialDefRate.size() != sizeStrains)
				rSpatialDefRate.resize(sizeStrains, false);
			noalias(rSpatialDefRate) = ZeroVector(sizeStrains);
		}
		if (itNode->SolutionStepsDataHas(NODAL_DEFORMATION_GRAD))
		{
			Matrix &rFgrad = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
			if (rFgrad.size1() != dimension)
				rFgrad.resize(dimension, dimension, false);
			noalias(rFgrad) = ZeroMatrix(dimension, dimension);
		}
		if (itNode->SolutionStepsDataHas(NODAL_DEFORMATION_GRAD_VEL))
		{
			Matrix &rFgradVel = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
			if (rFgradVel.size1() != dimension)
				rFgradVel.resize(dimension, dimension, false);
			noalias(rFgradVel) = ZeroMatrix(dimension, dimension);
		}
		if (itNode->SolutionStepsDataHas(NODAL_VOLUME))
		{
			itNode->FastGetSolutionStepValue(NODAL_VOLUME) = 0;
		}
		if (itNode->SolutionStepsDataHas(NODAL_MEAN_MESH_SIZE))
		{
			itNode->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE) = 0;
		}
		if (itNode->SolutionStepsDataHas(NODAL_FREESURFACE_AREA))
		{
			itNode->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) = 0;
		}
		if (itNode->SolutionStepsDataHas(NODAL_VOLUMETRIC_DEF_RATE))
		{
			itNode->FastGetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE) = 0;
		}
		if (itNode->SolutionStepsDataHas(NODAL_EQUIVALENT_STRAIN_RATE))
		{
			itNode->FastGetSolutionStepValue(NODAL_EQUIVALENT_STRAIN_RATE) = 0;
		}
	}

	void InitializeNonLinearIterations()
	{

		ModelPart &rModelPart = BaseType::GetModelPart();
		ElementsArrayType &pElements = rModelPart.Elements();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();

#ifdef _OPENMP
		int number_of_threads = omp_get_max_threads();
#else
		int number_of_threads = 1;
#endif

		vector<unsigned int> element_partition;
		OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

		// #pragma omp parallel
		//       {
		int k = OpenMPUtils::ThisThread();
		typename ElementsArrayType::iterator ElemBegin = pElements.begin() + element_partition[k];
		typename ElementsArrayType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

		for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++) //MSI: To be parallelized
		{
			itElem->InitializeNonLinearIteration(rCurrentProcessInfo);
		}
		//    }
	}

	void CalcNodalStrainsAndStresses()
	{
		ModelPart &rModelPart = BaseType::GetModelPart();

		// #pragma omp parallel
		//   {
		ModelPart::NodeIterator NodesBegin;
		ModelPart::NodeIterator NodesEnd;
		OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

		for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
		{

			double nodalVolume = itNode->FastGetSolutionStepValue(NODAL_VOLUME);

			double theta = 0.5;

			if (nodalVolume > 0)
			{
				this->ComputeAndStoreNodalDeformationGradient(itNode, theta);
				this->CalcNodalStrainsAndStressesForNode(itNode);
			}
			else
			{ // if nodalVolume==0
				InitializeNodalVariablesForRemeshedDomain(itNode);
			}
		}
		//   }
	}

	void CalcNodalStrainsAndStressesForNode(ModelPart::NodeIterator itNode)
	{

		ModelPart &rModelPart = BaseType::GetModelPart();

		const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

		double currFirstLame = itNode->FastGetSolutionStepValue(VOLUMETRIC_COEFFICIENT);
		double deviatoricCoeff = itNode->FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);

		Matrix Fgrad = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
		Matrix FgradVel = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
		double detFgrad = 1.0;
		Matrix InvFgrad = ZeroMatrix(dimension, dimension);
		Matrix SpatialVelocityGrad = ZeroMatrix(dimension, dimension);

		if (dimension == 2)
		{
			MathUtils<double>::InvertMatrix2(Fgrad, InvFgrad, detFgrad);
		}
		else if (dimension == 3)
		{
			MathUtils<double>::InvertMatrix3(Fgrad, InvFgrad, detFgrad);
		}

		//it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
		SpatialVelocityGrad = prod(FgradVel, InvFgrad);

		if (dimension == 2)
		{
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] = SpatialVelocityGrad(0, 0);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] = SpatialVelocityGrad(1, 1);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));

			double yieldShear = itNode->FastGetSolutionStepValue(YIELD_SHEAR);
			if (yieldShear > 0)
			{
				itNode->FastGetSolutionStepValue(NODAL_EQUIVALENT_STRAIN_RATE) = sqrt((2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] +
																					   2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] +
																					   4.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2]));
				double adaptiveExponent = itNode->FastGetSolutionStepValue(ADAPTIVE_EXPONENT);
				double equivalentStrainRate = itNode->FastGetSolutionStepValue(NODAL_EQUIVALENT_STRAIN_RATE);
				double exponent = -adaptiveExponent * equivalentStrainRate;
				if (equivalentStrainRate != 0)
				{
					deviatoricCoeff += (yieldShear / equivalentStrainRate) * (1 - exp(exponent));
				}
				if (equivalentStrainRate < 0.00001 && yieldShear != 0 && adaptiveExponent != 0)
				{
					// for gamma_dot very small the limit of the Papanastasiou viscosity is mu=m*tau_yield
					deviatoricCoeff = adaptiveExponent * yieldShear;
				}
			}

			double DefVol = itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] + itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];

			itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE) = DefVol;

			double nodalSigmaTot_xx = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
			double nodalSigmaTot_yy = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];
			double nodalSigmaTot_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];

			double nodalSigmaDev_xx = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] - DefVol / 3.0);
			double nodalSigmaDev_yy = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] - DefVol / 3.0);
			double nodalSigmaDev_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];

			itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS, 0)[0] = nodalSigmaTot_xx;
			itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS, 0)[1] = nodalSigmaTot_yy;
			itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS, 0)[2] = nodalSigmaTot_xy;

			itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS, 0)[0] = nodalSigmaDev_xx;
			itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS, 0)[1] = nodalSigmaDev_yy;
			itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS, 0)[2] = nodalSigmaDev_xy;
		}
		else if (dimension == 3)
		{
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] = SpatialVelocityGrad(0, 0);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] = SpatialVelocityGrad(1, 1);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] = SpatialVelocityGrad(2, 2);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4] = 0.5 * (SpatialVelocityGrad(2, 0) + SpatialVelocityGrad(0, 2));
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5] = 0.5 * (SpatialVelocityGrad(2, 1) + SpatialVelocityGrad(1, 2));

			double yieldShear = itNode->FastGetSolutionStepValue(YIELD_SHEAR);
			if (yieldShear > 0)
			{
				itNode->FastGetSolutionStepValue(NODAL_EQUIVALENT_STRAIN_RATE) = sqrt(2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] +
																					  2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] +
																					  2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] +
																					  4.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3] +
																					  4.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4] +
																					  4.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5]);
				double adaptiveExponent = itNode->FastGetSolutionStepValue(ADAPTIVE_EXPONENT);
				double equivalentStrainRate = itNode->FastGetSolutionStepValue(NODAL_EQUIVALENT_STRAIN_RATE);
				double exponent = -adaptiveExponent * equivalentStrainRate;
				if (equivalentStrainRate != 0)
				{
					deviatoricCoeff += (yieldShear / equivalentStrainRate) * (1 - exp(exponent));
				}
				if (equivalentStrainRate < 0.00001 && yieldShear != 0 && adaptiveExponent != 0)
				{
					// for gamma_dot very small the limit of the Papanastasiou viscosity is mu=m*tau_yield
					deviatoricCoeff = adaptiveExponent * yieldShear;
				}
			}

			double DefVol = itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] + itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] + itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];

			itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE) = DefVol;

			double nodalSigmaTot_xx = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
			double nodalSigmaTot_yy = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];
			double nodalSigmaTot_zz = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];
			double nodalSigmaTot_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3];
			double nodalSigmaTot_xz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4];
			double nodalSigmaTot_yz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5];

			double nodalSigmaDev_xx = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] - DefVol / 3.0);
			double nodalSigmaDev_yy = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] - DefVol / 3.0);
			double nodalSigmaDev_zz = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] - DefVol / 3.0);
			double nodalSigmaDev_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3];
			double nodalSigmaDev_xz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4];
			double nodalSigmaDev_yz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5];

			itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS, 0)[0] = nodalSigmaTot_xx;
			itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS, 0)[1] = nodalSigmaTot_yy;
			itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS, 0)[2] = nodalSigmaTot_zz;
			itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS, 0)[3] = nodalSigmaTot_xy;
			itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS, 0)[4] = nodalSigmaTot_xz;
			itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS, 0)[5] = nodalSigmaTot_yz;

			itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS, 0)[0] = nodalSigmaDev_xx;
			itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS, 0)[1] = nodalSigmaDev_yy;
			itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS, 0)[2] = nodalSigmaDev_zz;
			itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS, 0)[3] = nodalSigmaDev_xy;
			itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS, 0)[4] = nodalSigmaDev_xz;
			itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS, 0)[5] = nodalSigmaDev_yz;
		}
	}

	void CalcNodalStrainsForNode(ModelPart::NodeIterator itNode)
	{

		/* std::cout << "Calc Nodal Strains  " << std::endl; */
		ModelPart &rModelPart = BaseType::GetModelPart();

		const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

		//   Matrix Fgrad=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
		//   Matrix FgradVel=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
		//   double detFgrad=1.0;
		//   Matrix InvFgrad=ZeroMatrix(dimension,dimension);
		//   Matrix SpatialVelocityGrad=ZeroMatrix(dimension,dimension);

		double detFgrad = 1.0;
		Matrix nodalFgrad = ZeroMatrix(dimension, dimension);
		Matrix FgradVel = ZeroMatrix(dimension, dimension);
		Matrix InvFgrad = ZeroMatrix(dimension, dimension);
		Matrix SpatialVelocityGrad = ZeroMatrix(dimension, dimension);

		nodalFgrad = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
		FgradVel = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);

		//Inverse

		if (dimension == 2)
		{
			MathUtils<double>::InvertMatrix2(nodalFgrad, InvFgrad, detFgrad);
		}
		else if (dimension == 3)
		{
			MathUtils<double>::InvertMatrix3(nodalFgrad, InvFgrad, detFgrad);
		}

		//it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
		SpatialVelocityGrad = prod(FgradVel, InvFgrad);

		if (dimension == 2)
		{

			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] = SpatialVelocityGrad(0, 0);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] = SpatialVelocityGrad(1, 1);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));

			itNode->FastGetSolutionStepValue(NODAL_EQUIVALENT_STRAIN_RATE) = sqrt((2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] +
																				   2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] +
																				   4.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2]));

			double DefX = itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
			double DefY = itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];

			double DefVol = DefX + DefY;

			itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE) = DefVol;
		}
		else if (dimension == 3)
		{

			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] = SpatialVelocityGrad(0, 0);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] = SpatialVelocityGrad(1, 1);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] = SpatialVelocityGrad(2, 2);
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4] = 0.5 * (SpatialVelocityGrad(2, 0) + SpatialVelocityGrad(0, 2));
			itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5] = 0.5 * (SpatialVelocityGrad(2, 1) + SpatialVelocityGrad(1, 2));

			itNode->FastGetSolutionStepValue(NODAL_EQUIVALENT_STRAIN_RATE) = sqrt(2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] +
																				  2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] +
																				  2.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] +
																				  4.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3] +
																				  4.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4] +
																				  4.0 * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5] * itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5]);

			double DefX = itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
			double DefY = itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];
			double DefZ = itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];

			double DefVol = DefX + DefY + DefZ;

			itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE) = DefVol;
		}
	}

	void CalcNodalStrains()
	{

		/* std::cout << "Calc Nodal Strains  " << std::endl; */
		ModelPart &rModelPart = BaseType::GetModelPart();

		// #pragma omp parallel
		//     {
		ModelPart::NodeIterator NodesBegin;
		ModelPart::NodeIterator NodesEnd;
		OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

		for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
		{

			double nodalVolume = itNode->FastGetSolutionStepValue(NODAL_VOLUME);

			double theta = 1.0;

			if (nodalVolume > 0)
			{
				this->ComputeAndStoreNodalDeformationGradient(itNode, theta);
				this->CalcNodalStrainsForNode(itNode);
			}
			else
			{ // if nodalVolume==0
				InitializeNodalVariablesForRemeshedDomain(itNode);
			}
		}
		// }

		/* std::cout << "Calc Nodal Strains And Stresses DONE " << std::endl; */
	}

	void ComputeAndStoreNodalDeformationGradient(ModelPart::NodeIterator itNode, double theta)
	{

		KRATOS_TRY;

		ModelPart &rModelPart = BaseType::GetModelPart();
		const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
		Vector nodalSFDneighboursId = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
		Vector rNodalSFDneigh = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
		/* unsigned int idThisNode=nodalSFDneighboursId[0]; */
		const unsigned int neighSize = nodalSFDneighboursId.size();
		Matrix Fgrad = ZeroMatrix(dimension, dimension);
		Matrix FgradVel = ZeroMatrix(dimension, dimension);
		NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);

		if (dimension == 2)
		{

			double dNdXi = rNodalSFDneigh[0];
			double dNdYi = rNodalSFDneigh[1];

			Fgrad(0, 0) += dNdXi * itNode->X();
			Fgrad(0, 1) += dNdYi * itNode->X();
			Fgrad(1, 0) += dNdXi * itNode->Y();
			Fgrad(1, 1) += dNdYi * itNode->Y();

			double VelocityX = itNode->FastGetSolutionStepValue(VELOCITY_X, 0) * theta + itNode->FastGetSolutionStepValue(VELOCITY_X, 1) * (1 - theta);
			double VelocityY = itNode->FastGetSolutionStepValue(VELOCITY_Y, 0) * theta + itNode->FastGetSolutionStepValue(VELOCITY_Y, 1) * (1 - theta);

			FgradVel(0, 0) += dNdXi * VelocityX;
			FgradVel(0, 1) += dNdYi * VelocityX;
			FgradVel(1, 0) += dNdXi * VelocityY;
			FgradVel(1, 1) += dNdYi * VelocityY;

			unsigned int firstRow = 2;

			if (neighSize > 0)
			{
				for (unsigned int i = 0; i < neighSize - 1; i++) //neigh_nodes has one cell less than nodalSFDneighboursId becuase this has also the considered node ID at the beginning
				{
					dNdXi = rNodalSFDneigh[firstRow];
					dNdYi = rNodalSFDneigh[firstRow + 1];
					unsigned int neigh_nodes_id = neighb_nodes[i].Id();
					unsigned int other_neigh_nodes_id = nodalSFDneighboursId[i + 1];
					if (neigh_nodes_id != other_neigh_nodes_id)
					{
						std::cout << "node (x,y)=(" << itNode->X() << "," << itNode->Y() << ") with neigh_nodes_id " << neigh_nodes_id << " different than  other_neigh_nodes_id " << other_neigh_nodes_id << std::endl;
					}
					Fgrad(0, 0) += dNdXi * neighb_nodes[i].X();
					Fgrad(0, 1) += dNdYi * neighb_nodes[i].X();
					Fgrad(1, 0) += dNdXi * neighb_nodes[i].Y();
					Fgrad(1, 1) += dNdYi * neighb_nodes[i].Y();

					VelocityX = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_X, 0) * theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_X, 1) * (1 - theta);
					VelocityY = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Y, 0) * theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Y, 1) * (1 - theta);

					FgradVel(0, 0) += dNdXi * VelocityX;
					FgradVel(0, 1) += dNdYi * VelocityX;
					FgradVel(1, 0) += dNdXi * VelocityY;
					FgradVel(1, 1) += dNdYi * VelocityY;

					firstRow += 2;
				}
			}
		}
		else
		{

			double dNdXi = rNodalSFDneigh[0];
			double dNdYi = rNodalSFDneigh[1];
			double dNdZi = rNodalSFDneigh[2];

			double VelocityX = itNode->FastGetSolutionStepValue(VELOCITY_X, 0) * theta + itNode->FastGetSolutionStepValue(VELOCITY_X, 1) * (1 - theta);
			double VelocityY = itNode->FastGetSolutionStepValue(VELOCITY_Y, 0) * theta + itNode->FastGetSolutionStepValue(VELOCITY_Y, 1) * (1 - theta);
			double VelocityZ = itNode->FastGetSolutionStepValue(VELOCITY_Z, 0) * theta + itNode->FastGetSolutionStepValue(VELOCITY_Z, 1) * (1 - theta);

			Fgrad(0, 0) += dNdXi * itNode->X();
			Fgrad(0, 1) += dNdYi * itNode->X();
			Fgrad(0, 2) += dNdZi * itNode->X();

			Fgrad(1, 0) += dNdXi * itNode->Y();
			Fgrad(1, 1) += dNdYi * itNode->Y();
			Fgrad(1, 2) += dNdZi * itNode->Y();

			Fgrad(2, 0) += dNdXi * itNode->Z();
			Fgrad(2, 1) += dNdYi * itNode->Z();
			Fgrad(2, 2) += dNdZi * itNode->Z();

			FgradVel(0, 0) += dNdXi * VelocityX;
			FgradVel(0, 1) += dNdYi * VelocityX;
			FgradVel(0, 2) += dNdZi * VelocityX;

			FgradVel(1, 0) += dNdXi * VelocityY;
			FgradVel(1, 1) += dNdYi * VelocityY;
			FgradVel(1, 2) += dNdZi * VelocityY;

			FgradVel(2, 0) += dNdXi * VelocityZ;
			FgradVel(2, 1) += dNdYi * VelocityZ;
			FgradVel(2, 2) += dNdZi * VelocityZ;

			unsigned int firstRow = 3;

			if (neighSize > 0)
			{
				for (unsigned int i = 0; i < neighSize - 1; i++)
				{

					dNdXi = rNodalSFDneigh[firstRow];
					dNdYi = rNodalSFDneigh[firstRow + 1];
					dNdZi = rNodalSFDneigh[firstRow + 2];

					VelocityX = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_X, 0) * theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_X, 1) * (1 - theta);
					VelocityY = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Y, 0) * theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Y, 1) * (1 - theta);
					VelocityZ = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Z, 0) * theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Z, 1) * (1 - theta);

					Fgrad(0, 0) += dNdXi * neighb_nodes[i].X();
					Fgrad(0, 1) += dNdYi * neighb_nodes[i].X();
					Fgrad(0, 2) += dNdZi * neighb_nodes[i].X();

					Fgrad(1, 0) += dNdXi * neighb_nodes[i].Y();
					Fgrad(1, 1) += dNdYi * neighb_nodes[i].Y();
					Fgrad(1, 2) += dNdZi * neighb_nodes[i].Y();

					Fgrad(2, 0) += dNdXi * neighb_nodes[i].Z();
					Fgrad(2, 1) += dNdYi * neighb_nodes[i].Z();
					Fgrad(2, 2) += dNdZi * neighb_nodes[i].Z();

					FgradVel(0, 0) += dNdXi * VelocityX;
					FgradVel(0, 1) += dNdYi * VelocityX;
					FgradVel(0, 2) += dNdZi * VelocityX;

					FgradVel(1, 0) += dNdXi * VelocityY;
					FgradVel(1, 1) += dNdYi * VelocityY;
					FgradVel(1, 2) += dNdZi * VelocityY;

					FgradVel(2, 0) += dNdXi * VelocityZ;
					FgradVel(2, 1) += dNdYi * VelocityZ;
					FgradVel(2, 2) += dNdZi * VelocityZ;

					firstRow += 3;
				}
			}
		}

		itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD) = Fgrad;
		itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL) = FgradVel;
		KRATOS_CATCH("");
	}

	void UpdateTopology(ModelPart &rModelPart, unsigned int echoLevel)
	{
		KRATOS_TRY;
		/* this->CalculateDisplacements(); */
		this->CalculateDisplacementsAndResetNodalVariables();
		BaseType::MoveMesh();
		BoundaryNormalsCalculationUtilities BoundaryComputation;
		BoundaryComputation.CalculateWeightedBoundaryNormals(rModelPart, echoLevel);

		KRATOS_CATCH("");
	}

	void CalculatePressureVelocity()
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
		unsigned int timeStep = rCurrentProcessInfo[STEP];

		for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
			 i != rModelPart.NodesEnd(); ++i)
		{
			if (timeStep == 1)
			{
				(i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0) = 0;
				(i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1) = 0;
			}
			else
			{
				double &CurrentPressure = (i)->FastGetSolutionStepValue(PRESSURE, 0);
				double &PreviousPressure = (i)->FastGetSolutionStepValue(PRESSURE, 1);
				double &CurrentPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
				CurrentPressureVelocity = (CurrentPressure - PreviousPressure) / timeInterval;
			}
		}
	}

	void CalculatePressureAcceleration()
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
		unsigned int timeStep = rCurrentProcessInfo[STEP];

		for (ModelPart::NodeIterator i = rModelPart.NodesBegin(); i != rModelPart.NodesEnd(); ++i)
		{
			if (timeStep == 1)
			{
				(i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0) = 0;
				(i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1) = 0;
			}
			else
			{
				double &CurrentPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
				double &PreviousPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1);
				double &CurrentPressureAcceleration = (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0);
				CurrentPressureAcceleration = (CurrentPressureVelocity - PreviousPressureVelocity) / timeInterval;
			}
		}
	}

	void CalculateAccelerations()
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

		for (ModelPart::NodeIterator i = rModelPart.NodesBegin(); i != rModelPart.NodesEnd(); ++i)
		{

			array_1d<double, 3> &CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 0);
			array_1d<double, 3> &PreviousVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);

			array_1d<double, 3> &CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
			array_1d<double, 3> &PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

			if ((i)->IsNot(ISOLATED) && ((i)->IsNot(RIGID) || (i)->Is(SOLID)))
			{
				UpdateAccelerations(CurrentAcceleration, CurrentVelocity, PreviousAcceleration, PreviousVelocity, BDFcoeffs);
			}
			else if ((i)->Is(RIGID))
			{
				array_1d<double, 3> Zeros(3, 0.0);
				(i)->FastGetSolutionStepValue(ACCELERATION, 0) = Zeros;
				(i)->FastGetSolutionStepValue(ACCELERATION, 1) = Zeros;
			}
			else
			{
				(i)->FastGetSolutionStepValue(NODAL_VOLUME) = 0.0;
				(i)->FastGetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE) = 0.0;
				(i)->FastGetSolutionStepValue(NODAL_EQUIVALENT_STRAIN_RATE) = 0;
				(i)->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE) = 0.0;
				(i)->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) = 0.0;
				(i)->FastGetSolutionStepValue(PRESSURE, 0) = 0.0;
				(i)->FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
				(i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0) = 0.0;
				(i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1) = 0.0;
				(i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0) = 0.0;
				(i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1) = 0.0;
				if ((i)->SolutionStepsDataHas(VOLUME_ACCELERATION))
				{
					array_1d<double, 3> &VolumeAcceleration = (i)->FastGetSolutionStepValue(VOLUME_ACCELERATION);
					(i)->FastGetSolutionStepValue(ACCELERATION, 0) = VolumeAcceleration;
					(i)->FastGetSolutionStepValue(VELOCITY, 0) += VolumeAcceleration * rCurrentProcessInfo[DELTA_TIME];
				}
			}
		}
	}

	inline void UpdateAccelerations(array_1d<double, 3> &CurrentAcceleration,
									const array_1d<double, 3> &CurrentVelocity,
									array_1d<double, 3> &PreviousAcceleration,
									const array_1d<double, 3> &PreviousVelocity,
									Vector &BDFcoeffs)
	{
		/* noalias(PreviousAcceleration)=CurrentAcceleration; */
		noalias(CurrentAcceleration) = -BDFcoeffs[1] * (CurrentVelocity - PreviousVelocity) - PreviousAcceleration;
		// std::cout<<"rBDFCoeffs[0] is "<<rBDFCoeffs[0]<<std::endl;//3/(2*delta_t)
		// std::cout<<"rBDFCoeffs[1] is "<<rBDFCoeffs[1]<<std::endl;//-2/(delta_t)
		// std::cout<<"rBDFCoeffs[2] is "<<rBDFCoeffs[2]<<std::endl;//1/(2*delta_t)
	}

	void CalculateDisplacements()
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		const double TimeStep = rCurrentProcessInfo[DELTA_TIME];

		for (ModelPart::NodeIterator i = rModelPart.NodesBegin(); i != rModelPart.NodesEnd(); ++i)
		{

			array_1d<double, 3> &CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 0);
			array_1d<double, 3> &PreviousVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);

			array_1d<double, 3> &CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 0);
			array_1d<double, 3> &PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);

			/* if( i->IsFixed(DISPLACEMENT_X) == false ) */
			CurrentDisplacement[0] = 0.5 * TimeStep * (CurrentVelocity[0] + PreviousVelocity[0]) + PreviousDisplacement[0];

			/* if( i->IsFixed(DISPLACEMENT_Y) == false ) */
			CurrentDisplacement[1] = 0.5 * TimeStep * (CurrentVelocity[1] + PreviousVelocity[1]) + PreviousDisplacement[1];

			/* if( i->IsFixed(DISPLACEMENT_Z) == false ) */
			CurrentDisplacement[2] = 0.5 * TimeStep * (CurrentVelocity[2] + PreviousVelocity[2]) + PreviousDisplacement[2];
		}
	}

	void CalculateDisplacementsAndResetNodalVariables()
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
		const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
		unsigned int sizeStrains = 3 * (dimension - 1);

		// #pragma omp parallel
		//       {
		ModelPart::NodeIterator NodesBegin;
		ModelPart::NodeIterator NodesEnd;
		OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

		for (ModelPart::NodeIterator i = NodesBegin; i != NodesEnd; ++i)
		{
			array_1d<double, 3> &CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 0);
			array_1d<double, 3> &PreviousVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);

			array_1d<double, 3> &CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 0);
			array_1d<double, 3> &PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);

			CurrentDisplacement[0] = 0.5 * TimeStep * (CurrentVelocity[0] + PreviousVelocity[0]) + PreviousDisplacement[0];
			CurrentDisplacement[1] = 0.5 * TimeStep * (CurrentVelocity[1] + PreviousVelocity[1]) + PreviousDisplacement[1];
			if (dimension == 3)
			{
				CurrentDisplacement[2] = 0.5 * TimeStep * (CurrentVelocity[2] + PreviousVelocity[2]) + PreviousDisplacement[2];
			}

			///// reset Nodal variables //////
			Vector &rNodalSFDneighbours = i->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
			unsigned int sizeSDFNeigh = rNodalSFDneighbours.size();
			// unsigned int neighbourNodes=i->GetValue(NEIGHBOUR_NODES).size()+1;
			// unsigned int sizeSDFNeigh=neighbourNodes*dimension;

			i->FastGetSolutionStepValue(NODAL_VOLUME) = 0;
			i->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE) = 0;
			i->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) = 0;
			i->FastGetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE) = 0;
			i->FastGetSolutionStepValue(NODAL_EQUIVALENT_STRAIN_RATE) = 0;

			noalias(rNodalSFDneighbours) = ZeroVector(sizeSDFNeigh);

			Vector &rSpatialDefRate = i->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
			noalias(rSpatialDefRate) = ZeroVector(sizeStrains);

			Matrix &rFgrad = i->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
			noalias(rFgrad) = ZeroMatrix(dimension, dimension);

			Matrix &rFgradVel = i->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
			noalias(rFgradVel) = ZeroMatrix(dimension, dimension);
		}
		//  }
	}

	void UpdatePressureAccelerations()
	{
		this->CalculateAccelerations();
		this->CalculatePressureVelocity();
		this->CalculatePressureAcceleration();
	}

	void Clear() override
	{
		mpMomentumStrategy->Clear();
		mpPressureStrategy->Clear();
	}

	///@}
	///@name Access
	///@{

	void SetEchoLevel(int Level) override
	{
		BaseType::SetEchoLevel(Level);
		int StrategyLevel = Level > 0 ? Level - 1 : 0;
		mpMomentumStrategy->SetEchoLevel(StrategyLevel);
		mpPressureStrategy->SetEchoLevel(StrategyLevel);
	}

	///@}
	///@name Inquiry
	///@{

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const override
	{
		std::stringstream buffer;
		buffer << "NodalTwoStepVPStrategy";
		return buffer.str();
	}

	/// Print information about this object.
	void PrintInfo(std::ostream &rOStream) const override
	{
		rOStream << "NodalTwoStepVPStrategy";
	}

	/// Print object's data.
	void PrintData(std::ostream &rOStream) const override
	{
	}

	///@}
	///@name Friends
	///@{

	///@}

protected:
	///@name Protected Life Cycle
	///@{

	///@}
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	/// Calculate the coefficients for time iteration.
	/**
     * @param rCurrentProcessInfo ProcessInfo instance from the fluid ModelPart. Must contain DELTA_TIME and BDF_COEFFICIENTS variables.
     */
	void SetTimeCoefficients(ProcessInfo &rCurrentProcessInfo)
	{
		KRATOS_TRY;

		if (mTimeOrder == 2)
		{
			//calculate the BDF coefficients
			double Dt = rCurrentProcessInfo[DELTA_TIME];
			double OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

			double Rho = OldDt / Dt;
			double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

			Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
			BDFcoeffs.resize(3, false);

			BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho);		   //coefficient for step n+1 (3/2Dt if Dt is constant)
			BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
			BDFcoeffs[2] = TimeCoeff;								   //coefficient for step n-1 (1/2Dt if Dt is constant)
		}
		else if (mTimeOrder == 1)
		{
			double Dt = rCurrentProcessInfo[DELTA_TIME];
			double TimeCoeff = 1.0 / Dt;

			Vector &BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
			BDFcoeffs.resize(2, false);

			BDFcoeffs[0] = TimeCoeff;  //coefficient for step n+1 (1/Dt)
			BDFcoeffs[1] = -TimeCoeff; //coefficient for step n (-1/Dt)
		}

		KRATOS_CATCH("");
	}

	bool SolveMomentumIteration(unsigned int it, unsigned int maxIt, bool &fixedTimeStep)
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		int Rank = rModelPart.GetCommunicator().MyPID();
		bool ConvergedMomentum = false;
		double NormDv = 0;
		fixedTimeStep = false;
		// build momentum system and solve for fractional step velocity increment
		rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP, 1);

		if (it == 0)
		{
			mpMomentumStrategy->InitializeSolutionStep();
			/* this->SetNeighboursVelocityId(); */
		}

		NormDv = mpMomentumStrategy->Solve();

		if (BaseType::GetEchoLevel() > 1 && Rank == 0)
			std::cout << "-------------- s o l v e d ! ------------------" << std::endl;

		double DvErrorNorm = 0;
		ConvergedMomentum = this->CheckVelocityConvergence(NormDv, DvErrorNorm);

		unsigned int iterationForCheck = 3;
		KRATOS_INFO("TwoStepVPStrategy") << "iteration(" << it << ") Velocity error: " << DvErrorNorm << " velTol: " << mVelocityTolerance << std::endl;

		// Check convergence
		if (it == maxIt - 1)
		{
			std::cout << "         iteration(" << it << ") Final Velocity error: " << DvErrorNorm << " velTol: " << mVelocityTolerance << std::endl;
			fixedTimeStep = this->FixTimeStepMomentum(DvErrorNorm);
		}
		else if (it > iterationForCheck)
		{
			fixedTimeStep = this->CheckMomentumConvergence(DvErrorNorm);
		}
		// 	    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
		// double currentTime = rCurrentProcessInfo[TIME];
		// double tolerance=0.0000000001;
		// if(currentTime>(0.25-tolerance) && currentTime<(0.25+tolerance)){
		// 	std::ofstream myfile;
		//   myfile.open ("velocityConvergenceAt025s.txt",std::ios::app);
		// 	myfile << it << "\t" << DvErrorNorm << "\n";
		//   myfile.close();
		// }
		// else if(currentTime>(0.5-tolerance) && currentTime<(0.5+tolerance)){
		// 	std::ofstream myfile;
		//   myfile.open ("velocityConvergenceAt05s.txt",std::ios::app);
		// 	myfile << it << "\t" << DvErrorNorm << "\n";
		//   myfile.close();
		// }
		// else if(currentTime>(0.75-tolerance) && currentTime<(0.75+tolerance)){
		// 	std::ofstream myfile;
		//   myfile.open ("velocityConvergenceAt075s.txt",std::ios::app);
		// 	myfile << it << "\t" << DvErrorNorm << "\n";
		//   myfile.close();
		// }
		// else if(currentTime>(1.0-tolerance) && currentTime<(1.0+tolerance)){
		// 	std::ofstream myfile;
		//   myfile.open ("velocityConvergenceAt100s.txt",std::ios::app);
		// 	myfile << it << "\t" << DvErrorNorm << "\n";
		//   myfile.close();
		// }

		if (!ConvergedMomentum && BaseType::GetEchoLevel() > 0 && Rank == 0)
			std::cout << "Momentum equations did not reach the convergence tolerance." << std::endl;

		return ConvergedMomentum;
	}

	bool SolveContinuityIteration(unsigned int it, unsigned int maxIt)
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		int Rank = rModelPart.GetCommunicator().MyPID();
		bool ConvergedContinuity = false;
		double NormDp = 0;
		// 2. Pressure solution
		rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP, 5);
		if (it == 0)
		{
			mpPressureStrategy->InitializeSolutionStep();
		}

		NormDp = mpPressureStrategy->Solve();

		if (BaseType::GetEchoLevel() > 0 && Rank == 0)
			std::cout << "The norm of pressure is: " << NormDp << std::endl;

		double DpErrorNorm = 0;
		ConvergedContinuity = this->CheckPressureConvergence(NormDp, DpErrorNorm);

		// Check convergence
		if (it == maxIt - 1)
		{
			std::cout << "                  iteration(" << it << ") Final Pressure error: " << DpErrorNorm << " presTol: " << mPressureTolerance << std::endl;
			ConvergedContinuity = this->FixTimeStepContinuity(DpErrorNorm);
		}
		else
		{
			std::cout << "                             iteration(" << it << ") Pressure error: " << DpErrorNorm << " presTol: " << mPressureTolerance << std::endl;
		}

		// ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
		// double currentTime = rCurrentProcessInfo[TIME];
		// double tolerance=0.0000000001;
		// if(currentTime>(0.25-tolerance) && currentTime<(0.25+tolerance)){
		// 	std::ofstream myfile;
		//   myfile.open ("pressureConvergenceAt025s.txt",std::ios::app);
		// 	myfile << it << "\t" << DpErrorNorm << "\n";
		//   myfile.close();
		// }
		// else if(currentTime>(0.5-tolerance) && currentTime<(0.5+tolerance)){
		// 	std::ofstream myfile;
		//   myfile.open ("pressureConvergenceAt05s.txt",std::ios::app);
		// 	myfile << it << "\t" << DpErrorNorm << "\n";
		//   myfile.close();
		// }
		// else if(currentTime>(0.75-tolerance) && currentTime<(0.75+tolerance)){
		// 	std::ofstream myfile;
		//   myfile.open ("pressureConvergenceAt075s.txt",std::ios::app);
		// 	myfile << it << "\t" << DpErrorNorm << "\n";
		//   myfile.close();
		// }
		// else if(currentTime>(1.0-tolerance) && currentTime<(1.0+tolerance)){
		// 	std::ofstream myfile;
		//   myfile.open ("pressureConvergenceAt100s.txt",std::ios::app);
		// 	myfile << it << "\t" << DpErrorNorm << "\n";
		//   myfile.close();
		// }

		if (!ConvergedContinuity && BaseType::GetEchoLevel() > 0 && Rank == 0)
			std::cout << "Continuity equation did not reach the convergence tolerance." << std::endl;

		return ConvergedContinuity;
	}

	bool CheckVelocityConvergence(const double NormDv, double &errorNormDv)
	{
		ModelPart &rModelPart = BaseType::GetModelPart();

		double NormV = 0.00;
		errorNormDv = 0;

#pragma omp parallel reduction(+ \
							   : NormV)
		{
			ModelPart::NodeIterator NodeBegin;
			ModelPart::NodeIterator NodeEnd;
			OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
			for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
			{
				const array_1d<double, 3> &Vel = itNode->FastGetSolutionStepValue(VELOCITY);

				double NormVelNode = 0;

				for (unsigned int d = 0; d < 3; ++d)
				{
					NormVelNode += Vel[d] * Vel[d];
					NormV += Vel[d] * Vel[d];
				}
			}
		}
		BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormV);

		NormV = sqrt(NormV);

		if (NormV == 0.0)
			NormV = 1.00;

		errorNormDv = NormDv / NormV;

		if (BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
		{
			std::cout << "The norm of velocity increment is: " << NormDv << std::endl;
			std::cout << "The norm of velocity is: " << NormV << std::endl;
			std::cout << "Velocity error: " << errorNormDv << "mVelocityTolerance: " << mVelocityTolerance << std::endl;
		}
		/* else{ */
		/*   std::cout<<"Velocity error: "<< errorNormDv <<" velTol: " << mVelocityTolerance<< std::endl; */
		/* } */

		if (errorNormDv < mVelocityTolerance)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	void ComputeErrorL2NormCaseImposedG()
	{

		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		const double currentTime = rCurrentProcessInfo[TIME];

		double sumErrorL2Velocity = 0;
		double sumErrorL2VelocityX = 0;
		double sumErrorL2VelocityY = 0;
		double sumErrorL2Pressure = 0;
		double sumErrorL2TauXX = 0;
		double sumErrorL2TauYY = 0;
		double sumErrorL2TauXY = 0;

#pragma omp parallel
		{
			ModelPart::NodeIterator NodeBegin;
			ModelPart::NodeIterator NodeEnd;
			OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
			for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
			{
				const double posX = itNode->X();
				const double posY = itNode->Y();
				const double nodalArea = itNode->FastGetSolutionStepValue(NODAL_VOLUME);
				const double velX = itNode->FastGetSolutionStepValue(VELOCITY_X);
				const double velY = itNode->FastGetSolutionStepValue(VELOCITY_Y);
				const double pressure = itNode->FastGetSolutionStepValue(PRESSURE);
				const double tauXX = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[0];
				const double tauYY = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[1];
				const double tauXY = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[2];

				double expectedVelocityX = pow(posX, 2) * (1.0 - posX) * (1.0 - posX) * (2.0 * posY - 6.0 * pow(posY, 2) + 4.0 * pow(posY, 3));
				double expectedVelocityY = -pow(posY, 2) * (1.0 - posY) * (1.0 - posY) * (2.0 * posX - 6.0 * pow(posX, 2) + 4.0 * pow(posX, 3));
				double expectedPressure = -posX * (1.0 - posX);

				double expectedTauXX = 2.0 * (-4.0 * (1 - posX) * posX * (-1.0 + 2.0 * posX) * posY * (1.0 - 3.0 * posY + 2.0 * pow(posY, 2)));
				double expectedTauYY = 2.0 * (4.0 * posX * (1.0 - 3.0 * posX + 2.0 * pow(posX, 2)) * (1 - posY) * posY * (-1.0 + 2.0 * posY));
				double expectedTauXY = (2.0 * (1.0 - 6.0 * posY + 6.0 * pow(posY, 2)) * (1 - posX) * (1 - posX) * pow(posX, 2) - 2.0 * (1.0 - 6.0 * posX + 6.0 * pow(posX, 2)) * (1 - posY) * (1 - posY) * pow(posY, 2));

				double nodalErrorVelocityX = velX - expectedVelocityX;
				double nodalErrorVelocityY = velY - expectedVelocityY;
				double nodalErrorPressure = pressure - expectedPressure;
				double nodalErrorTauXX = tauXX - expectedTauXX;
				double nodalErrorTauYY = tauYY - expectedTauYY;
				double nodalErrorTauXY = tauXY - expectedTauXY;

				sumErrorL2Velocity += (pow(nodalErrorVelocityX, 2) + pow(nodalErrorVelocityY, 2)) * nodalArea;
				sumErrorL2VelocityX += pow(nodalErrorVelocityX, 2) * nodalArea;
				sumErrorL2VelocityY += pow(nodalErrorVelocityY, 2) * nodalArea;
				sumErrorL2Pressure += pow(nodalErrorPressure, 2) * nodalArea;
				sumErrorL2TauXX += pow(nodalErrorTauXX, 2) * nodalArea;
				sumErrorL2TauYY += pow(nodalErrorTauYY, 2) * nodalArea;
				sumErrorL2TauXY += pow(nodalErrorTauXY, 2) * nodalArea;

				// itNode->FastGetSolutionStepValue(NODAL_ERROR_XX)=nodalErrorTauXX;
			}
		}

		double errorL2Velocity = sqrt(sumErrorL2Velocity);
		double errorL2VelocityX = sqrt(sumErrorL2VelocityX);
		double errorL2VelocityY = sqrt(sumErrorL2VelocityY);
		double errorL2Pressure = sqrt(sumErrorL2Pressure);
		double errorL2TauXX = sqrt(sumErrorL2TauXX);
		double errorL2TauYY = sqrt(sumErrorL2TauYY);
		double errorL2TauXY = sqrt(sumErrorL2TauXY);

		std::ofstream myfileVelocity;
		myfileVelocity.open("errorL2VelocityFile.txt", std::ios::app);
		myfileVelocity << currentTime << "\t" << errorL2Velocity << "\n";
		myfileVelocity.close();

		std::ofstream myfileVelocityX;
		myfileVelocityX.open("errorL2VelocityXFile.txt", std::ios::app);
		myfileVelocityX << currentTime << "\t" << errorL2VelocityX << "\n";
		myfileVelocityX.close();

		std::ofstream myfileVelocityY;
		myfileVelocityY.open("errorL2VelocityYFile.txt", std::ios::app);
		myfileVelocityY << currentTime << "\t" << errorL2VelocityY << "\n";
		myfileVelocityY.close();

		std::ofstream myfilePressure;
		myfilePressure.open("errorL2PressureFile.txt", std::ios::app);
		myfilePressure << currentTime << "\t" << errorL2Pressure << "\n";
		myfilePressure.close();

		std::ofstream myfileTauXX;
		myfileTauXX.open("errorL2TauXXFile.txt", std::ios::app);
		myfileTauXX << currentTime << "\t" << errorL2TauXX << "\n";
		myfileTauXX.close();

		std::ofstream myfileTauYY;
		myfileTauYY.open("errorL2TauYYFile.txt", std::ios::app);
		myfileTauYY << currentTime << "\t" << errorL2TauYY << "\n";
		myfileTauYY.close();

		std::ofstream myfileTauXY;
		myfileTauXY.open("errorL2TauXYFile.txt", std::ios::app);
		myfileTauXY << currentTime << "\t" << errorL2TauXY << "\n";
		myfileTauXY.close();
	}

	void ComputeErrorL2NormCasePoiseuille()
	{

		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		const double currentTime = rCurrentProcessInfo[TIME];

		double sumErrorL2VelocityTheta = 0;
		double sumErrorL2TauTheta = 0;

		double r_in = 0.2;
		double R_out = 0.5;
		double kappa = r_in / R_out;
		double omega = 0.5;
		double viscosity = 100.0;

#pragma omp parallel
		{
			ModelPart::NodeIterator NodeBegin;
			ModelPart::NodeIterator NodeEnd;
			OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
			for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
			{
				const double posX = itNode->X();
				const double posY = itNode->Y();
				const double rPos = sqrt(pow(posX, 2) + pow(posY, 2));
				const double cosalfa = posX / rPos;
				const double sinalfa = posY / rPos;
				const double sin2alfa = 2.0 * cosalfa * sinalfa;
				const double cos2alfa = 1.0 - 2.0 * pow(sinalfa, 2);
				const double nodalArea = itNode->FastGetSolutionStepValue(NODAL_VOLUME);
				const double velX = itNode->FastGetSolutionStepValue(VELOCITY_X);
				const double velY = itNode->FastGetSolutionStepValue(VELOCITY_Y);
				const double tauXX = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[0];
				const double tauYY = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[1];
				const double tauXY = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[2];

				double expectedVelocityTheta = pow(kappa, 2) * omega * R_out / (1.0 - pow(kappa, 2)) * (R_out / rPos - rPos / R_out);
				double computedVelocityTheta = sqrt(pow(velX, 2) + pow(velY, 2));
				double nodalErrorVelocityTheta = computedVelocityTheta - expectedVelocityTheta;

				double expectedTauTheta = (2.0 * viscosity * pow(kappa, 2) * omega * pow(R_out, 2)) / (1.0 - pow(kappa, 2)) / pow(rPos, 2);
				double computedTauTheta = +(tauXX - tauYY) * sin2alfa / 2.0 - tauXY * cos2alfa;
				double nodalErrorTauTheta = computedTauTheta - expectedTauTheta;
				itNode->FastGetSolutionStepValue(NODAL_ERROR_XX) = computedVelocityTheta;
				// if(posY>-0.01 && posY<0.01){
				//  std::cout<<"expectedTauTheta "<<expectedTauTheta<<"   computedTauTheta "<<computedTauTheta <<std::endl;
				//  std::cout<<"tauXX "<<tauXX<<"   tauYY "<<tauYY<<"   tauXY "<<tauXY <<std::endl;
				//  std::cout<<"posX  "<<posX <<"   posY  "<<posY <<std::endl;
				//  std::cout<<"\n  ";
				// }

				// if(posX>-0.01 && posX<0.01){
				//  std::cout<<"expectedTauTheta "<<expectedTauTheta<<"   computedTauTheta "<<computedTauTheta <<std::endl;
				//  std::cout<<"tauXX "<<tauXX<<"   tauYY "<<tauYY<<"   tauXY "<<tauXY <<std::endl;
				//  std::cout<<"posX  "<<posX <<"   posY  "<<posY <<std::endl;
				//  std::cout<<"\n  ";
				// }

				sumErrorL2VelocityTheta += pow(nodalErrorVelocityTheta, 2) * nodalArea;
				sumErrorL2TauTheta += pow(nodalErrorTauTheta, 2) * nodalArea;
			}
		}

		double errorL2VelocityTheta = sqrt(sumErrorL2VelocityTheta);
		double errorL2TauTheta = sqrt(sumErrorL2TauTheta);

		std::ofstream myfileVelocity;
		myfileVelocity.open("errorL2Poiseuille.txt", std::ios::app);
		myfileVelocity << currentTime << "\t" << errorL2VelocityTheta << "\t" << errorL2TauTheta << "\n";
		myfileVelocity.close();
	}

	bool CheckPressureConvergence(const double NormDp, double &errorNormDp)
	{
		ModelPart &rModelPart = BaseType::GetModelPart();

		double NormP = 0.00;
		errorNormDp = 0;

		// #pragma omp parallel reduction(+:NormP)
		//         {
		ModelPart::NodeIterator NodeBegin;
		ModelPart::NodeIterator NodeEnd;
		OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
		for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
		{
			const double Pr = itNode->FastGetSolutionStepValue(PRESSURE);
			NormP += Pr * Pr;
		}
		// }

		BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormP);

		NormP = sqrt(NormP);

		if (NormP == 0.0)
			NormP = 1.00;

		errorNormDp = NormDp / NormP;

		if (BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
		{
			std::cout << "         The norm of pressure increment is: " << NormDp << std::endl;
			std::cout << "         The norm of pressure is: " << NormP << std::endl;
			std::cout << "         Pressure error: " << errorNormDp << std::endl;
		}
		/* else{ */
		/*     std::cout<<"         Pressure error: "<<errorNormDp <<" presTol: "<<mPressureTolerance << std::endl; */
		/* } */

		if (errorNormDp < mPressureTolerance)
		{
			return true;
		}
		else
			return false;
	}

	bool FixTimeStepMomentum(const double DvErrorNorm)
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		double currentTime = rCurrentProcessInfo[TIME];
		double timeInterval = rCurrentProcessInfo[DELTA_TIME];
		double minTolerance = 0.005;
		bool fixedTimeStep = false;
		if (currentTime < 3 * timeInterval)
		{
			minTolerance = 10;
		}

		if ((DvErrorNorm > minTolerance || (DvErrorNorm < 0 && DvErrorNorm > 0) || (DvErrorNorm != DvErrorNorm)) &&
			DvErrorNorm != 0 &&
			(DvErrorNorm != 1 || currentTime > timeInterval))
		{
			rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, true);
			std::cout << "NOT GOOD CONVERGENCE!!! I'll reduce the next time interval" << DvErrorNorm << std::endl;
			minTolerance = 0.05;
			if (DvErrorNorm > minTolerance)
			{
				std::cout << "BAD CONVERGENCE!!! I GO AHEAD WITH THE PREVIOUS VELOCITY AND PRESSURE FIELDS" << DvErrorNorm << std::endl;
				fixedTimeStep = true;
				// #pragma omp parallel
				// 	  {
				ModelPart::NodeIterator NodeBegin;
				ModelPart::NodeIterator NodeEnd;
				OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
				for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
				{
					itNode->FastGetSolutionStepValue(VELOCITY, 0) = itNode->FastGetSolutionStepValue(VELOCITY, 1);
					itNode->FastGetSolutionStepValue(PRESSURE, 0) = itNode->FastGetSolutionStepValue(PRESSURE, 1);
					itNode->FastGetSolutionStepValue(ACCELERATION, 0) = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
				}
				// }
			}
		}
		else
		{
			rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, false);
		}
		return fixedTimeStep;
	}

	bool CheckMomentumConvergence(const double DvErrorNorm)
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		double currentTime = rCurrentProcessInfo[TIME];
		double timeInterval = rCurrentProcessInfo[DELTA_TIME];
		double minTolerance = 0.99999;
		bool fixedTimeStep = false;

		if ((DvErrorNorm > minTolerance || (DvErrorNorm < 0 && DvErrorNorm > 0) || (DvErrorNorm != DvErrorNorm)) &&
			DvErrorNorm != 0 &&
			(DvErrorNorm != 1 || currentTime > timeInterval))
		{
			rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, true);
			std::cout << "           BAD CONVERGENCE DETECTED DURING THE ITERATIVE LOOP!!! error: " << DvErrorNorm << " higher than 0.9999" << std::endl;
			std::cout << "      I GO AHEAD WITH THE PREVIOUS VELOCITY AND PRESSURE FIELDS" << std::endl;
			fixedTimeStep = true;
#pragma omp parallel
			{
				ModelPart::NodeIterator NodeBegin;
				ModelPart::NodeIterator NodeEnd;
				OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodeBegin, NodeEnd);
				for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
				{
					itNode->FastGetSolutionStepValue(VELOCITY, 0) = itNode->FastGetSolutionStepValue(VELOCITY, 1);
					itNode->FastGetSolutionStepValue(PRESSURE, 0) = itNode->FastGetSolutionStepValue(PRESSURE, 1);
					itNode->FastGetSolutionStepValue(ACCELERATION, 0) = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
				}
			}
		}
		else
		{
			rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, false);
		}
		return fixedTimeStep;
	}

	bool FixTimeStepContinuity(const double DvErrorNorm)
	{
		ModelPart &rModelPart = BaseType::GetModelPart();
		ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
		double currentTime = rCurrentProcessInfo[TIME];
		double timeInterval = rCurrentProcessInfo[DELTA_TIME];
		double minTolerance = 0.01;
		bool fixedTimeStep = false;
		if (currentTime < 3 * timeInterval)
		{
			minTolerance = 10;
		}

		if ((DvErrorNorm > minTolerance || (DvErrorNorm < 0 && DvErrorNorm > 0) || (DvErrorNorm != DvErrorNorm)) &&
			DvErrorNorm != 0 &&
			(DvErrorNorm != 1 || currentTime > timeInterval))
		{
			fixedTimeStep = true;
			rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, true);
		}
		else
		{
			rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, false);
		}
		return fixedTimeStep;
	}

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

	// private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	double mVelocityTolerance;

	double mPressureTolerance;

	unsigned int mMaxPressureIter;

	unsigned int mDomainSize;

	unsigned int mTimeOrder;

	bool mReformDofSet;

	// Fractional step index.
	/*  1 : Momentum step (calculate fractional step velocity)
      * 2-3 : Unused (reserved for componentwise calculation of frac step velocity)
      * 4 : Pressure step
      * 5 : Computation of projections
      * 6 : End of step velocity
      */
	//    unsigned int mStepId;

	/// Scheme for the solution of the momentum equation
	StrategyPointerType mpMomentumStrategy;

	/// Scheme for the solution of the mass equation
	StrategyPointerType mpPressureStrategy;

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	void InitializeStrategy(SolverSettingsType &rSolverConfig)
	{
		KRATOS_TRY;

		mTimeOrder = rSolverConfig.GetTimeOrder();

		// Check that input parameters are reasonable and sufficient.
		this->Check();

		//ModelPart& rModelPart = this->GetModelPart();

		mDomainSize = rSolverConfig.GetDomainSize();

		mReformDofSet = rSolverConfig.GetReformDofSet();

		BaseType::SetEchoLevel(rSolverConfig.GetEchoLevel());

		// Initialize strategies for each step
		bool HaveVelStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Velocity, mpMomentumStrategy);

		if (HaveVelStrategy)
		{
			rSolverConfig.FindTolerance(SolverSettingsType::Velocity, mVelocityTolerance);
			/* rSolverConfig.FindMaxIter(SolverSettingsType::Velocity,mMaxVelocityIter); */
		}
		else
		{
			KRATOS_THROW_ERROR(std::runtime_error, "NodalTwoStepVPStrategy error: No Velocity strategy defined in FractionalStepSettings", "");
		}

		bool HavePressStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Pressure, mpPressureStrategy);

		if (HavePressStrategy)
		{
			rSolverConfig.FindTolerance(SolverSettingsType::Pressure, mPressureTolerance);
			rSolverConfig.FindMaxIter(SolverSettingsType::Pressure, mMaxPressureIter);
		}
		else
		{
			KRATOS_THROW_ERROR(std::runtime_error, "NodalTwoStepVPStrategy error: No Pressure strategy defined in FractionalStepSettings", "");
		}

		// Check input parameters
		this->Check();

		KRATOS_CATCH("");
	}

	///@}
	///@name Private  Access
	///@{

	///@}
	///@name Private Inquiry
	///@{

	///@}
	///@name Un accessible methods
	///@{

	/// Assignment operator.
	NodalTwoStepVPStrategy &operator=(NodalTwoStepVPStrategy const &rOther) {}

	/// Copy constructor.
	NodalTwoStepVPStrategy(NodalTwoStepVPStrategy const &rOther) {}

	///@}

}; /// Class NodalTwoStepVPStrategy

///@}
///@name Type Definitions
///@{

///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_NODAL_TWO_STEP_V_P_STRATEGY_H
