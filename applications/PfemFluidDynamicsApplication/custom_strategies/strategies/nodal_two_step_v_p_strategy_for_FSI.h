//
//   Project Name:        KratosPFEMFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:                   June 2018 $
//   Revision:            $Revision:                     0.0 $
//
//

#ifndef KRATOS_NODAL_TWO_STEP_V_P_STRATEGY_FOR_FSI_H
#define KRATOS_NODAL_TWO_STEP_V_P_STRATEGY_FOR_FSI_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/builders_and_solvers/nodal_residualbased_elimination_builder_and_solver_for_FSI.h"
#include "custom_strategies/builders_and_solvers/nodal_residualbased_elimination_builder_and_solver_continuity_for_FSI.h"

#include "custom_utilities/solver_settings.h"

#include "custom_strategies/strategies/gauss_seidel_linear_strategy.h"

#include "pfem_fluid_dynamics_application_variables.h"

#include "nodal_two_step_v_p_strategy.h"

#include <stdio.h>
#include <cmath>
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
	class NodalTwoStepVPStrategyForFSI : public NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
	{
	public:
		///@name Type Definitions
		///@{
		KRATOS_CLASS_POINTER_DEFINITION(NodalTwoStepVPStrategyForFSI);

		/// Counted pointer of NodalTwoStepVPStrategy
		// typedef boost::shared_ptr< NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver> > Pointer;

		typedef NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

		typedef typename BaseType::TDataType TDataType;

		/// Node type (default is: Node)
		typedef Node NodeType;

		/// Geometry type (using with given NodeType)
		typedef Geometry<NodeType> GeometryType;

		typedef std::size_t SizeType;

		// typedef typename BaseType::DofSetType DofSetType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

		typedef typename BaseType::ElementsArrayType ElementsArrayType;

		typedef typename ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

		typedef TwoStepVPSolverSettings<TSparseSpace, TDenseSpace, TLinearSolver> SolverSettingsType;

		using NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mVelocityTolerance;
		using NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mPressureTolerance;
		using NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mMaxPressureIter;
		using NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mDomainSize;
		using NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mReformDofSet;
		using NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mpMomentumStrategy;
		using NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::mpPressureStrategy;

		typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

		typedef GlobalPointersVector<Node> NodeWeakPtrVectorType;
		///@}
		///@name Life Cycle
		///@{

		NodalTwoStepVPStrategyForFSI(ModelPart &rModelPart,
									 SolverSettingsType &rSolverConfig) : BaseType(rModelPart)
		{
			NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::InitializeStrategy(rSolverConfig);
		}

		NodalTwoStepVPStrategyForFSI(ModelPart &rModelPart,
									 /*SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& rSolverConfig,*/
									 typename TLinearSolver::Pointer pVelocityLinearSolver,
									 typename TLinearSolver::Pointer pPressureLinearSolver,
									 bool ReformDofSet = true,
									 double VelTol = 0.0001,
									 double PresTol = 0.0001,
									 int MaxPressureIterations = 1, // Only for predictor-corrector
									 unsigned int TimeOrder = 2,
									 unsigned int DomainSize = 2) : BaseType(rModelPart,
																			 pVelocityLinearSolver,
																			 pPressureLinearSolver,
																			 ReformDofSet,
																			 VelTol,
																			 PresTol,
																			 MaxPressureIterations,
																			 TimeOrder,
																			 DomainSize)
		{

			KRATOS_TRY;

			BaseType::SetEchoLevel(1);

			// Check that input parameters are reasonable and sufficient.
			this->Check();

			bool CalculateNormDxFlag = true;

			bool ReformDofAtEachIteration = false; // DofSet modifiaction is managed by the fractional step strategy, auxiliary strategies should not modify the DofSet directly.

			// Additional Typedefs
			typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
			typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

			// initializing fractional velocity solution step
			typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;
			typename SchemeType::Pointer pScheme;

			typename SchemeType::Pointer Temp = typename SchemeType::Pointer(new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>());
			pScheme.swap(Temp);

			// CONSTRUCTION OF VELOCITY
			BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new NodalResidualBasedEliminationBuilderAndSolverForFSI<TSparseSpace, TDenseSpace, TLinearSolver>(pVelocityLinearSolver));

			this->mpMomentumStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pVelocityLinearSolver, vel_build, ReformDofAtEachIteration, CalculateNormDxFlag));

			this->mpMomentumStrategy->SetEchoLevel(BaseType::GetEchoLevel());

			vel_build->SetCalculateReactionsFlag(false);

			BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new NodalResidualBasedEliminationBuilderAndSolverContinuityForFSI<TSparseSpace, TDenseSpace, TLinearSolver>(pPressureLinearSolver));

			this->mpPressureStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pPressureLinearSolver, pressure_build, ReformDofAtEachIteration, CalculateNormDxFlag));

			this->mpPressureStrategy->SetEchoLevel(BaseType::GetEchoLevel());

			pressure_build->SetCalculateReactionsFlag(false);

			KRATOS_CATCH("");
		}

		/// Destructor.
		virtual ~NodalTwoStepVPStrategyForFSI() {}

		bool SolveSolutionStep() override
		{
			// Initialize BDF2 coefficients
			ModelPart &rModelPart = BaseType::GetModelPart();
			ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
			double currentTime = rCurrentProcessInfo[TIME];
			double timeInterval = rCurrentProcessInfo[DELTA_TIME];
			bool timeIntervalChanged = rCurrentProcessInfo[TIME_INTERVAL_CHANGED];
			bool converged = false;

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
			double pressureNorm = 0;
			double velocityNorm = 0;

			FillNodalSFDVector(); // it fills SOLID_NODAL_SFD_NEIGHBOURS_ORDER for solids and NODAL_SFD_NEIGHBOURS_ORDER for fluids and inner solids
			for (unsigned int it = 0; it < maxNonLinearIterations; ++it)
			{
				if (BaseType::GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
					std::cout << "----- > iteration: " << it << std::endl;

				if (it == 0)
				{
					ComputeNodalVolumeAndAssignFlagToElementType(); // it assings NODAL_VOLUME to fluid and SOLID_NODAL_VOLUME to solid. Interface nodes have both
					this->InitializeNonLinearIterations();			// it fills SOLID_NODAL_SFD_NEIGHBOURS for solids and NODAL_SFD_NEIGHBOURS for fluids
				}
				CalcNodalStrainsAndStresses(); // it computes stresses and strains for fluid and solid nodes

				momentumConverged = this->SolveMomentumIteration(it, maxNonLinearIterations, fixedTimeStep, velocityNorm);

				UpdateTopology(rModelPart, BaseType::GetEchoLevel());
				ComputeNodalVolume();
				this->InitializeNonLinearIterations();
				CalcNodalStrains();

				if (fixedTimeStep == false)
				{
					continuityConverged = this->SolveContinuityIteration(it, maxNonLinearIterations, pressureNorm);
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
					// this->ComputeErrorL2NormCaseImposedG();
					// this->ComputeErrorL2NormCasePoiseuille();
					this->CalculateAccelerations();
					// std::ofstream myfile;
					// myfile.open ("maxConvergedIteration.txt",std::ios::app);
					// myfile << currentTime << "\t" << it << "\n";
					// myfile.close();
				}
				bool hybridMethod = false;
				if (hybridMethod == true)
				{
					if (it == maxNonLinearIterations - 1 || ((continuityConverged && momentumConverged) && it > 0))
					{
						this->UpdateElementalStressStrain();
					}
				}

				if ((continuityConverged && momentumConverged) && it > 1)
				{
					rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE, false);
					rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE, false);
					converged = true;
					std::cout << "nodal V-P strategy converged in " << it + 1 << " iterations." << std::endl;
					break;
				}
				if (fixedTimeStep == true)
				{
					break;
				}
			}

			if (!continuityConverged && !momentumConverged && BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
				std::cout << "Convergence tolerance not reached." << std::endl;

			if (mReformDofSet)
				this->Clear();

			return converged;
		}

		void UpdateElementalStressStrain()
		{
			ModelPart &rModelPart = BaseType::GetModelPart();
			const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();

#pragma omp parallel
			{
				ModelPart::ElementIterator ElemBegin;
				ModelPart::ElementIterator ElemEnd;
				OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
				for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
				{
					itElem->InitializeSolutionStep(rCurrentProcessInfo);
				}
			}
		}

		void Initialize() override
		{

			std::cout << "  \n     Initialize in nodal_two_step_v_p_strategy_FSI" << std::endl;
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

				if (itNode->SolutionStepsDataHas(SOLID_NODAL_CAUCHY_STRESS))
				{
					Vector &rSolidNodalStress = itNode->FastGetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS);
					if (rSolidNodalStress.size() != sizeStrains)
					{
						rSolidNodalStress.resize(sizeStrains, false);
					}
					noalias(rSolidNodalStress) = ZeroVector(sizeStrains);
				}
				else
				{
					std::cout << "THIS node does not have SOLID_NODAL_CAUCHY_STRESS... " << itNode->X() << " " << itNode->Y() << std::endl;
				}

				if (itNode->SolutionStepsDataHas(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS))
				{
					Vector &rSolidNodalStress = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS);
					if (rSolidNodalStress.size() != sizeStrains)
					{
						rSolidNodalStress.resize(sizeStrains, false);
					}
					noalias(rSolidNodalStress) = ZeroVector(sizeStrains);
				}
				else
				{
					std::cout << "THIS node does not have SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS... " << itNode->X() << " " << itNode->Y() << std::endl;
				}

				if (itNode->SolutionStepsDataHas(SOLID_NODAL_VOLUME))
				{
					itNode->FastGetSolutionStepValue(SOLID_NODAL_VOLUME) = 0;
				}
				else
				{
					std::cout << "THIS node does not have SOLID_NODAL_VOLUME... " << itNode->X() << " " << itNode->Y() << std::endl;
				}

				if (itNode->SolutionStepsDataHas(SOLID_NODAL_MEAN_MESH_SIZE))
				{
					itNode->FastGetSolutionStepValue(SOLID_NODAL_MEAN_MESH_SIZE) = 0;
				}
				else
				{
					std::cout << "THIS node does not have SOLID_NODAL_MEAN_MESH_SIZE... " << itNode->X() << " " << itNode->Y() << std::endl;
				}

				if (itNode->SolutionStepsDataHas(SOLID_NODAL_FREESURFACE_AREA))
				{
					itNode->FastGetSolutionStepValue(SOLID_NODAL_FREESURFACE_AREA) = 0;
				}
				else
				{
					std::cout << "THIS node does not have SOLID_NODAL_FREESURFACE_AREA... " << itNode->X() << " " << itNode->Y() << std::endl;
				}

				if (itNode->SolutionStepsDataHas(SOLID_NODAL_SFD_NEIGHBOURS))
				{
					Vector &rSolidNodalSFDneighbours = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS);
					if (rSolidNodalSFDneighbours.size() != sizeSDFNeigh)
					{
						rSolidNodalSFDneighbours.resize(sizeSDFNeigh, false);
					}
					noalias(rSolidNodalSFDneighbours) = ZeroVector(sizeSDFNeigh);
				}
				else
				{
					std::cout << "THIS node does not have SOLID_NODAL_SFD_NEIGHBOURS... " << itNode->X() << " " << itNode->Y() << std::endl;
				}

				if (itNode->SolutionStepsDataHas(SOLID_NODAL_SPATIAL_DEF_RATE))
				{
					Vector &rSolidSpatialDefRate = itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE);
					if (rSolidSpatialDefRate.size() != sizeStrains)
					{
						rSolidSpatialDefRate.resize(sizeStrains, false);
					}
					noalias(rSolidSpatialDefRate) = ZeroVector(sizeStrains);
				}
				else
				{
					std::cout << "THIS node does not have SOLID_NODAL_SPATIAL_DEF_RATE... " << itNode->X() << " " << itNode->Y() << std::endl;
				}

				if (itNode->SolutionStepsDataHas(SOLID_NODAL_DEFORMATION_GRAD))
				{
					Matrix &rSolidFgrad = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD);
					if (rSolidFgrad.size1() != dimension)
					{
						rSolidFgrad.resize(dimension, dimension, false);
					}
					noalias(rSolidFgrad) = ZeroMatrix(dimension, dimension);
				}
				else
				{
					std::cout << "THIS node does not have SOLID_NODAL_DEFORMATION_GRAD... " << itNode->X() << " " << itNode->Y() << std::endl;
				}

				if (itNode->SolutionStepsDataHas(SOLID_NODAL_DEFORMATION_GRAD_VEL))
				{
					Matrix &rSolidFgradVel = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL);
					if (rSolidFgradVel.size1() != dimension)
					{
						rSolidFgradVel.resize(dimension, dimension, false);
					}
					noalias(rSolidFgradVel) = ZeroMatrix(dimension, dimension);
				}
				else
				{
					std::cout << "THIS node does not have SOLID_NODAL_DEFORMATION_GRAD_VEL... " << itNode->X() << " " << itNode->Y() << std::endl;
				}

				InitialAssignMaterialToEachNode(itNode);
			}

			// }
		}

		void InitialAssignMaterialToEachNode(ModelPart::NodeIterator itNode)
		{

			ModelPart &rModelPart = BaseType::GetModelPart();
			const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
			const double timeInterval = rCurrentProcessInfo[DELTA_TIME];

			double deviatoricCoeff = 0;
			double volumetricCoeff = 0;

			if (itNode->Is(SOLID))
			{
				const double youngModulus = itNode->FastGetSolutionStepValue(YOUNG_MODULUS);
				const double poissonRatio = itNode->FastGetSolutionStepValue(POISSON_RATIO);

				deviatoricCoeff = timeInterval * youngModulus / (1.0 + poissonRatio) * 0.5;
				volumetricCoeff = timeInterval * poissonRatio * youngModulus / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)) + 2.0 * deviatoricCoeff / 3.0;
			}
			else if (itNode->Is(FLUID) || itNode->Is(RIGID))
			{
				if (rModelPart.GetNodalSolutionStepVariablesList().Has(DYNAMIC_VISCOSITY))
				{
					deviatoricCoeff = itNode->FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
				}
				volumetricCoeff = timeInterval * itNode->FastGetSolutionStepValue(BULK_MODULUS);
			}

			if ((itNode->Is(SOLID) && itNode->Is(RIGID)))
			{
				itNode->FastGetSolutionStepValue(INTERFACE_NODE) = true;
			}
			else
			{
				itNode->FastGetSolutionStepValue(INTERFACE_NODE) = false;
			}

			const double currFirstLame = volumetricCoeff - 2.0 * deviatoricCoeff / 3.0;

			// currFirstLame=deltaT*firstLame
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

			for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++) // MSI: To be parallelized
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

					if (itElem->Is(SOLID))
					{
						double &solidVolume = geometry(i)->FastGetSolutionStepValue(SOLID_NODAL_VOLUME);
						solidVolume += elementalVolume;
						nodalVolume += -elementalVolume;
					}
				}
			}
			// }
		}

		void ComputeNodalVolumeAndAssignFlagToElementType()
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

			for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++) // MSI: To be parallelized
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
				unsigned int fluidNodes = 0;
				unsigned int solidNodes = 0;
				unsigned int interfaceNodes = 0;

				for (unsigned int i = 0; i < numNodes; i++)
				{

					if ((geometry(i)->Is(FLUID) && geometry(i)->IsNot(SOLID)) || (geometry(i)->Is(FLUID) && geometry(i)->FastGetSolutionStepValue(INTERFACE_NODE) == true))
					{
						fluidNodes += 1;
					}
					if (geometry(i)->Is(SOLID))
					{
						solidNodes += 1;
					}
					if (geometry(i)->FastGetSolutionStepValue(INTERFACE_NODE) == true)
					{
						interfaceNodes += 1;
					}
				}

				if (solidNodes == numNodes)
				{
					itElem->Set(SOLID);
				}
				if (interfaceNodes == numNodes)
				{
					itElem->Set(SOLID);
				}
				if (fluidNodes == numNodes)
				{
					itElem->Set(FLUID);
				}
				if (solidNodes == numNodes && fluidNodes == numNodes)
				{
					itElem->Reset(FLUID);
				}

				for (unsigned int i = 0; i < numNodes; i++)
				{

					double &nodalVolume = geometry(i)->FastGetSolutionStepValue(NODAL_VOLUME);
					nodalVolume += elementalVolume;

					if (itElem->Is(SOLID))
					{
						double &solidVolume = geometry(i)->FastGetSolutionStepValue(SOLID_NODAL_VOLUME);
						solidVolume += elementalVolume;
						nodalVolume += -elementalVolume;
					}
				}
			}
			// }
		}

		void SetNeighboursOrderToSolidNode(ModelPart::NodeIterator itNode)
		{
			NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
			unsigned int neighbourNodes = neighb_nodes.size() + 1; // +1 becausealso the node itself must be considered as nieghbor node
			Vector &rNodeOrderedNeighbours = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);

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

		void SetNeighboursOrderToInterfaceNode(ModelPart::NodeIterator itNode)
		{
			NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
			unsigned int neighbourNodes = neighb_nodes.size() + 1;

			unsigned int fluidCounter = 1;
			unsigned int solidCounter = 1;

			if (neighbourNodes > 1)
			{
				for (unsigned int k = 0; k < neighbourNodes - 1; k++)
				{
					if (neighb_nodes[k].IsNot(SOLID) || neighb_nodes[k].FastGetSolutionStepValue(INTERFACE_NODE) == true)
					{
						fluidCounter += 1;
					}
					if (neighb_nodes[k].Is(SOLID))
					{
						solidCounter += 1;
					}
				}
			}

			Vector &rFluidNodeOrderedNeighbours = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
			Vector &rSolidNodeOrderedNeighbours = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);

			if (rFluidNodeOrderedNeighbours.size() != fluidCounter)
				rFluidNodeOrderedNeighbours.resize(fluidCounter, false);

			if (rSolidNodeOrderedNeighbours.size() != solidCounter)
				rSolidNodeOrderedNeighbours.resize(solidCounter, false);

			noalias(rFluidNodeOrderedNeighbours) = ZeroVector(fluidCounter);
			noalias(rSolidNodeOrderedNeighbours) = ZeroVector(solidCounter);

			rFluidNodeOrderedNeighbours[0] = itNode->Id();
			rSolidNodeOrderedNeighbours[0] = itNode->Id();

			fluidCounter = 0;
			solidCounter = 0;

			if (neighbourNodes > 1)
			{
				for (unsigned int k = 0; k < neighbourNodes - 1; k++)
				{
					if (neighb_nodes[k].IsNot(SOLID) || neighb_nodes[k].FastGetSolutionStepValue(INTERFACE_NODE) == true)
					{
						fluidCounter += 1;
						rFluidNodeOrderedNeighbours[fluidCounter] = neighb_nodes[k].Id();
					}
					if (neighb_nodes[k].Is(SOLID))
					{
						solidCounter += 1;
						rSolidNodeOrderedNeighbours[solidCounter] = neighb_nodes[k].Id();
					}
				}
			}
		}

		void InitializeNodalVariablesForSolidRemeshedDomain(ModelPart::NodeIterator itNode)
		{

			ModelPart &rModelPart = BaseType::GetModelPart();
			const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
			unsigned int sizeStrains = 3 * (dimension - 1);
			NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
			unsigned int neighbourNodes = neighb_nodes.size() + 1;
			unsigned int sizeSDFNeigh = neighbourNodes * dimension;

			if (itNode->SolutionStepsDataHas(SOLID_NODAL_CAUCHY_STRESS))
			{
				Vector &rSolidNodalStress = itNode->FastGetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS);
				if (rSolidNodalStress.size() != sizeStrains)
					rSolidNodalStress.resize(sizeStrains, false);
				noalias(rSolidNodalStress) = ZeroVector(sizeStrains);
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS))
			{
				Vector &rSolidNodalDevStress = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS);
				if (rSolidNodalDevStress.size() != sizeStrains)
					rSolidNodalDevStress.resize(sizeStrains, false);
				noalias(rSolidNodalDevStress) = ZeroVector(sizeStrains);
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_SFD_NEIGHBOURS))
			{
				Vector &rSolidNodalSFDneighbours = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS);
				if (rSolidNodalSFDneighbours.size() != sizeSDFNeigh)
					rSolidNodalSFDneighbours.resize(sizeSDFNeigh, false);
				noalias(rSolidNodalSFDneighbours) = ZeroVector(sizeSDFNeigh);
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_SFD_NEIGHBOURS_ORDER))
			{
				Vector &rSolidNodalSFDneighboursOrder = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
				if (rSolidNodalSFDneighboursOrder.size() != neighbourNodes)
					rSolidNodalSFDneighboursOrder.resize(neighbourNodes, false);
				noalias(rSolidNodalSFDneighboursOrder) = ZeroVector(neighbourNodes);
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_SPATIAL_DEF_RATE))
			{
				Vector &rSolidSpatialDefRate = itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE);
				if (rSolidSpatialDefRate.size() != sizeStrains)
					rSolidSpatialDefRate.resize(sizeStrains, false);
				noalias(rSolidSpatialDefRate) = ZeroVector(sizeStrains);
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_DEFORMATION_GRAD))
			{
				Matrix &rSolidFgrad = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD);
				if (rSolidFgrad.size1() != dimension)
					rSolidFgrad.resize(dimension, dimension, false);
				noalias(rSolidFgrad) = ZeroMatrix(dimension, dimension);
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_DEFORMATION_GRAD_VEL))
			{
				Matrix &rSolidFgradVel = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL);
				if (rSolidFgradVel.size1() != dimension)
					rSolidFgradVel.resize(dimension, dimension, false);
				noalias(rSolidFgradVel) = ZeroMatrix(dimension, dimension);
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_VOLUME))
			{
				itNode->FastGetSolutionStepValue(SOLID_NODAL_VOLUME) = 0;
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_MEAN_MESH_SIZE))
			{
				itNode->FastGetSolutionStepValue(SOLID_NODAL_MEAN_MESH_SIZE) = 0;
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_FREESURFACE_AREA))
			{
				itNode->FastGetSolutionStepValue(SOLID_NODAL_FREESURFACE_AREA) = 0;
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_VOLUMETRIC_DEF_RATE))
			{
				itNode->FastGetSolutionStepValue(SOLID_NODAL_VOLUMETRIC_DEF_RATE) = 0;
			}
			if (itNode->SolutionStepsDataHas(SOLID_NODAL_EQUIVALENT_STRAIN_RATE))
			{
				itNode->FastGetSolutionStepValue(SOLID_NODAL_EQUIVALENT_STRAIN_RATE) = 0;
			}
		}

		void CalcNodalStrainsAndStresses()
		{

			ModelPart &rModelPart = BaseType::GetModelPart();

			const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			// #pragma omp parallel
			//   {
			ModelPart::NodeIterator NodesBegin;
			ModelPart::NodeIterator NodesEnd;
			OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

			for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
			{

				const double nodalVolume = itNode->FastGetSolutionStepValue(NODAL_VOLUME);
				const double solidNodalVolume = itNode->FastGetSolutionStepValue(SOLID_NODAL_VOLUME);

				double theta = 0.5;

				if (itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true)
				{

					if (nodalVolume > 0)
					{
						Vector nodalSFDneighboursId = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
						Vector rNodalSFDneigh = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);

						Matrix &interfaceFgrad = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
						Matrix &interfaceFgradVel = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);

						if (interfaceFgrad.size1() != dimension)
							interfaceFgrad.resize(dimension, dimension, false);

						if (interfaceFgradVel.size1() != dimension)
							interfaceFgradVel.resize(dimension, dimension, false);

						noalias(interfaceFgrad) = ZeroMatrix(dimension, dimension);
						noalias(interfaceFgradVel) = ZeroMatrix(dimension, dimension);

						// I have to compute the stresses and strains two times because one time is for the solid and the other for the fluid
						//  Matrix interfaceFgrad=ZeroMatrix(dimension,dimension);
						//  Matrix interfaceFgradVel=ZeroMatrix(dimension,dimension);
						// the following function is more expensive than the general one because there is one loop more over neighbour nodes. This is why I do it here also for fluid interface nodes.
						ComputeAndStoreNodalDeformationGradientForInterfaceNode(itNode, nodalSFDneighboursId, rNodalSFDneigh, theta, interfaceFgrad, interfaceFgradVel);
						// itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)=interfaceFgrad;
						// itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL)=interfaceFgradVel;
						CalcNodalStrainsAndStressesForInterfaceFluidNode(itNode);
					}

					if (solidNodalVolume > 0)
					{
						Vector solidNodalSFDneighboursId = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
						Vector rSolidNodalSFDneigh = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS);

						Matrix &solidInterfaceFgrad = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD);
						Matrix &solidInterfaceFgradVel = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL);

						if (solidInterfaceFgrad.size1() != dimension)
							solidInterfaceFgrad.resize(dimension, dimension, false);

						if (solidInterfaceFgradVel.size1() != dimension)
							solidInterfaceFgradVel.resize(dimension, dimension, false);

						noalias(solidInterfaceFgrad) = ZeroMatrix(dimension, dimension);
						noalias(solidInterfaceFgradVel) = ZeroMatrix(dimension, dimension);

						theta = 1.0;
						// Matrix solidInterfaceFgrad=ZeroMatrix(dimension,dimension);
						// Matrix solidInterfaceFgradVel=ZeroMatrix(dimension,dimension);
						ComputeAndStoreNodalDeformationGradientForInterfaceNode(itNode, solidNodalSFDneighboursId, rSolidNodalSFDneigh, theta, solidInterfaceFgrad, solidInterfaceFgradVel);
						// itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD)=solidInterfaceFgrad;
						// itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL)=solidInterfaceFgradVel;
						CalcNodalStrainsAndStressesForInterfaceSolidNode(itNode);
					}
				}
				else
				{
					if (itNode->Is(SOLID) && solidNodalVolume > 0)
					{
						theta = 1.0;
						ComputeAndStoreNodalDeformationGradientForSolidNode(itNode, theta);
						CalcNodalStrainsAndStressesForSolidNode(itNode);
					}
					else if (nodalVolume > 0)
					{
						theta = 0.5;
						this->ComputeAndStoreNodalDeformationGradient(itNode, theta);
						this->CalcNodalStrainsAndStressesForNode(itNode);
					}
				}
				if (nodalVolume == 0 && solidNodalVolume == 0)
				{ // if nodalVolume==0
					theta = 0.5;
					this->InitializeNodalVariablesForRemeshedDomain(itNode);
					InitializeNodalVariablesForSolidRemeshedDomain(itNode);
				}
			}
			//   }
		}

		void CopyValuesToSolidNonInterfaceNodes(ModelPart::NodeIterator itNode)
		{
			Vector &solidNodalSFDneighboursId = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
			Vector &solidNodalSFDneigh = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS);
			Matrix &solidInterfaceFgrad = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD);
			Matrix &solidInterfaceFgradVel = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL);
			Vector &solidSpatialDefRate = itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE);
			double &volumetricDefRate = itNode->GetSolutionStepValue(SOLID_NODAL_VOLUMETRIC_DEF_RATE);
			Vector &solidCauchyStress = itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS);
			Vector &solidDeviatoricCauchyStress = itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS);

			Vector nodalSFDneighboursId = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
			unsigned int sizeNodalSFDneighboursId = nodalSFDneighboursId.size();
			solidNodalSFDneighboursId.resize(sizeNodalSFDneighboursId, false);

			Vector nodalSFDneigh = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
			unsigned int sizeNodalSFDneigh = nodalSFDneigh.size();
			solidNodalSFDneigh.resize(sizeNodalSFDneigh, false);

			solidNodalSFDneighboursId = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
			solidNodalSFDneigh = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
			solidInterfaceFgrad = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
			solidInterfaceFgradVel = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
			solidSpatialDefRate = itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
			volumetricDefRate = itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE);
			solidCauchyStress = itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS);
			solidDeviatoricCauchyStress = itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS);
		}

		void CalcNodalStrainsAndStressesForInterfaceFluidNode(ModelPart::NodeIterator itNode)
		{

			ModelPart &rModelPart = BaseType::GetModelPart();
			const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
			const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
			const double timeInterval = rCurrentProcessInfo[DELTA_TIME];

			double deviatoricCoeff = 0;

			const double currFirstLame = timeInterval * itNode->FastGetSolutionStepValue(BULK_MODULUS);

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

			// it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
			SpatialVelocityGrad = prod(FgradVel, InvFgrad);

			if (dimension == 2)
			{
				itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] = SpatialVelocityGrad(0, 0);
				itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] = SpatialVelocityGrad(1, 1);
				itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));

				this->ComputeDeviatoricCoefficientForFluid(itNode, dimension, deviatoricCoeff);

				const double DefVol = itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] + itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];

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

				this->ComputeDeviatoricCoefficientForFluid(itNode, dimension, deviatoricCoeff);

				const double DefVol = itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0] + itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1] + itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];

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

		void CalcNodalStrainsAndStressesForInterfaceSolidNode(ModelPart::NodeIterator itNode)
		{

			ModelPart &rModelPart = BaseType::GetModelPart();
			const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
			const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
			const double timeInterval = rCurrentProcessInfo[DELTA_TIME];

			const double youngModulus = itNode->FastGetSolutionStepValue(YOUNG_MODULUS);
			const double poissonRatio = itNode->FastGetSolutionStepValue(POISSON_RATIO);

			const double currFirstLame = timeInterval * poissonRatio * youngModulus / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio));
			double deviatoricCoeff = timeInterval * youngModulus / (1.0 + poissonRatio) * 0.5;

			Matrix Fgrad = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD);
			Matrix FgradVel = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL);
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

			// it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
			SpatialVelocityGrad = prod(FgradVel, InvFgrad);

			if (dimension == 2)
			{
				auto &r_stain_tensor2D = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE);
				r_stain_tensor2D[0] = SpatialVelocityGrad(0, 0);
				r_stain_tensor2D[1] = SpatialVelocityGrad(1, 1);
				r_stain_tensor2D[2] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));

				const double DefVol = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] + itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1];

				itNode->GetSolutionStepValue(SOLID_NODAL_VOLUMETRIC_DEF_RATE) = DefVol;

				double nodalSigmaTot_xx = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0];
				double nodalSigmaTot_yy = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1];
				double nodalSigmaTot_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2];

				double nodalSigmaDev_xx = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] - DefVol / 3.0);
				double nodalSigmaDev_yy = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] - DefVol / 3.0);
				double nodalSigmaDev_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2];

				if (itNode->Is(SOLID))
				{
					nodalSigmaTot_xx += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[0];
					nodalSigmaTot_yy += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[1];
					nodalSigmaTot_xy += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[2];

					nodalSigmaDev_xx += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[0];
					nodalSigmaDev_yy += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[1];
					nodalSigmaDev_xy += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[2];
				}

				auto &r_stress_tensor2D = itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 0);
				r_stress_tensor2D[0] = nodalSigmaTot_xx;
				r_stress_tensor2D[1] = nodalSigmaTot_yy;
				r_stress_tensor2D[2] = nodalSigmaTot_xy;

				auto &r_dev_stress_tensor2D = itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 0);
				r_dev_stress_tensor2D[0] = nodalSigmaDev_xx;
				r_dev_stress_tensor2D[1] = nodalSigmaDev_yy;
				r_dev_stress_tensor2D[2] = nodalSigmaDev_xy;
			}
			else if (dimension == 3)
			{

				auto &r_stain_tensor3D = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE);
				r_stain_tensor3D[0] = SpatialVelocityGrad(0, 0);
				r_stain_tensor3D[1] = SpatialVelocityGrad(1, 1);
				r_stain_tensor3D[2] = SpatialVelocityGrad(2, 2);
				r_stain_tensor3D[3] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));
				r_stain_tensor3D[4] = 0.5 * (SpatialVelocityGrad(2, 0) + SpatialVelocityGrad(0, 2));
				r_stain_tensor3D[5] = 0.5 * (SpatialVelocityGrad(2, 1) + SpatialVelocityGrad(1, 2));

				const double DefVol = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] + itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] + itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2];

				itNode->GetSolutionStepValue(SOLID_NODAL_VOLUMETRIC_DEF_RATE) = DefVol;

				double nodalSigmaTot_xx = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0];
				double nodalSigmaTot_yy = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1];
				double nodalSigmaTot_zz = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2];
				double nodalSigmaTot_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[3];
				double nodalSigmaTot_xz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[4];
				double nodalSigmaTot_yz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[5];

				double nodalSigmaDev_xx = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] - DefVol / 3.0);
				double nodalSigmaDev_yy = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] - DefVol / 3.0);
				double nodalSigmaDev_zz = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2] - DefVol / 3.0);
				double nodalSigmaDev_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[3];
				double nodalSigmaDev_xz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[4];
				double nodalSigmaDev_yz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[5];

				if (itNode->Is(SOLID))
				{
					nodalSigmaTot_xx += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[0];
					nodalSigmaTot_yy += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[1];
					nodalSigmaTot_zz += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[2];
					nodalSigmaTot_xy += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[3];
					nodalSigmaTot_xz += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[4];
					nodalSigmaTot_yz += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[5];

					nodalSigmaDev_xx += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[0];
					nodalSigmaDev_yy += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[1];
					nodalSigmaDev_zz += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[2];
					nodalSigmaDev_xy += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[3];
					nodalSigmaDev_xz += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[4];
					nodalSigmaDev_yz += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[5];
				}

				auto &r_stress_tensor3D = itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 0);
				r_stress_tensor3D[0] = nodalSigmaTot_xx;
				r_stress_tensor3D[1] = nodalSigmaTot_yy;
				r_stress_tensor3D[2] = nodalSigmaTot_zz;
				r_stress_tensor3D[3] = nodalSigmaTot_xy;
				r_stress_tensor3D[4] = nodalSigmaTot_xz;
				r_stress_tensor3D[5] = nodalSigmaTot_yz;

				auto &r_dev_stress_tensor3D = itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 0);
				r_dev_stress_tensor3D[0] = nodalSigmaDev_xx;
				r_dev_stress_tensor3D[1] = nodalSigmaDev_yy;
				r_dev_stress_tensor3D[2] = nodalSigmaDev_zz;
				r_dev_stress_tensor3D[3] = nodalSigmaDev_xy;
				r_dev_stress_tensor3D[4] = nodalSigmaDev_xz;
				r_dev_stress_tensor3D[5] = nodalSigmaDev_yz;
			}
		}

		void CalcNodalStrainsAndStressesForSolidNode(ModelPart::NodeIterator itNode)
		{

			ModelPart &rModelPart = BaseType::GetModelPart();
			const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
			const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
			const double timeInterval = rCurrentProcessInfo[DELTA_TIME];

			const double youngModulus = itNode->FastGetSolutionStepValue(YOUNG_MODULUS);
			const double poissonRatio = itNode->FastGetSolutionStepValue(POISSON_RATIO);

			const double currFirstLame = timeInterval * poissonRatio * youngModulus / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio));
			double deviatoricCoeff = timeInterval * youngModulus / (1.0 + poissonRatio) * 0.5;

			Matrix Fgrad = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD);
			Matrix FgradVel = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL);
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

			// it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
			SpatialVelocityGrad = prod(FgradVel, InvFgrad);

			if (dimension == 2)
			{
				auto &r_stain_tensor2D = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE);
				r_stain_tensor2D[0] = SpatialVelocityGrad(0, 0);
				r_stain_tensor2D[1] = SpatialVelocityGrad(1, 1);
				r_stain_tensor2D[2] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));

				const double DefVol = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] + itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1];

				itNode->GetSolutionStepValue(SOLID_NODAL_VOLUMETRIC_DEF_RATE) = DefVol;

				double nodalSigmaTot_xx = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0];
				double nodalSigmaTot_yy = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1];
				double nodalSigmaTot_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2];

				double nodalSigmaDev_xx = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] - DefVol / 3.0);
				double nodalSigmaDev_yy = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] - DefVol / 3.0);
				double nodalSigmaDev_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2];

				if (itNode->Is(SOLID))
				{
					nodalSigmaTot_xx += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[0];
					nodalSigmaTot_yy += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[1];
					nodalSigmaTot_xy += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[2];

					nodalSigmaDev_xx += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[0];
					nodalSigmaDev_yy += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[1];
					nodalSigmaDev_xy += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[2];
				}

				auto &r_stress_tensor2D = itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 0);
				r_stress_tensor2D[0] = nodalSigmaTot_xx;
				r_stress_tensor2D[1] = nodalSigmaTot_yy;
				r_stress_tensor2D[2] = nodalSigmaTot_xy;

				auto &r_dev_stress_tensor2D = itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 0);
				r_dev_stress_tensor2D[0] = nodalSigmaDev_xx;
				r_dev_stress_tensor2D[1] = nodalSigmaDev_yy;
				r_dev_stress_tensor2D[2] = nodalSigmaDev_xy;
			}
			else if (dimension == 3)
			{

				auto &r_stain_tensor3D = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE);
				r_stain_tensor3D[0] = SpatialVelocityGrad(0, 0);
				r_stain_tensor3D[1] = SpatialVelocityGrad(1, 1);
				r_stain_tensor3D[2] = SpatialVelocityGrad(2, 2);
				r_stain_tensor3D[3] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));
				r_stain_tensor3D[4] = 0.5 * (SpatialVelocityGrad(2, 0) + SpatialVelocityGrad(0, 2));
				r_stain_tensor3D[5] = 0.5 * (SpatialVelocityGrad(2, 1) + SpatialVelocityGrad(1, 2));

				const double DefVol = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] + itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] + itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2];

				itNode->GetSolutionStepValue(SOLID_NODAL_VOLUMETRIC_DEF_RATE) = DefVol;

				double nodalSigmaTot_xx = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0];
				double nodalSigmaTot_yy = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1];
				double nodalSigmaTot_zz = currFirstLame * DefVol + 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2];
				double nodalSigmaTot_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[3];
				double nodalSigmaTot_xz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[4];
				double nodalSigmaTot_yz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[5];

				double nodalSigmaDev_xx = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] - DefVol / 3.0);
				double nodalSigmaDev_yy = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] - DefVol / 3.0);
				double nodalSigmaDev_zz = 2.0 * deviatoricCoeff * (itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2] - DefVol / 3.0);
				double nodalSigmaDev_xy = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[3];
				double nodalSigmaDev_xz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[4];
				double nodalSigmaDev_yz = 2.0 * deviatoricCoeff * itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[5];

				if (itNode->Is(SOLID))
				{
					nodalSigmaTot_xx += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[0];
					nodalSigmaTot_yy += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[1];
					nodalSigmaTot_zz += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[2];
					nodalSigmaTot_xy += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[3];
					nodalSigmaTot_xz += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[4];
					nodalSigmaTot_yz += itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 1)[5];

					nodalSigmaDev_xx += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[0];
					nodalSigmaDev_yy += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[1];
					nodalSigmaDev_zz += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[2];
					nodalSigmaDev_xy += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[3];
					nodalSigmaDev_xz += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[4];
					nodalSigmaDev_yz += itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 1)[5];
				}

				auto &r_tensor3D = itNode->GetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS, 0);
				r_tensor3D[0] = nodalSigmaTot_xx;
				r_tensor3D[1] = nodalSigmaTot_yy;
				r_tensor3D[2] = nodalSigmaTot_zz;
				r_tensor3D[3] = nodalSigmaTot_xy;
				r_tensor3D[4] = nodalSigmaTot_xz;
				r_tensor3D[5] = nodalSigmaTot_yz;

				auto &r_dev_tensor3D = itNode->GetSolutionStepValue(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS, 0);
				r_dev_tensor3D[0] = nodalSigmaDev_xx;
				r_dev_tensor3D[1] = nodalSigmaDev_yy;
				r_dev_tensor3D[2] = nodalSigmaDev_zz;
				r_dev_tensor3D[3] = nodalSigmaDev_xy;
				r_dev_tensor3D[4] = nodalSigmaDev_xz;
				r_dev_tensor3D[5] = nodalSigmaDev_yz;
			}
		}

		void CalcNodalStrainsForSolidNode(ModelPart::NodeIterator itNode)
		{

			/* std::cout << "Calc Nodal Strains  " << std::endl; */
			ModelPart &rModelPart = BaseType::GetModelPart();

			const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			double detFgrad = 1.0;
			Matrix nodalFgrad = ZeroMatrix(dimension, dimension);
			Matrix FgradVel = ZeroMatrix(dimension, dimension);
			Matrix InvFgrad = ZeroMatrix(dimension, dimension);
			Matrix SpatialVelocityGrad = ZeroMatrix(dimension, dimension);

			nodalFgrad = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD);
			FgradVel = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL);

			// Inverse

			if (dimension == 2)
			{
				MathUtils<double>::InvertMatrix2(nodalFgrad, InvFgrad, detFgrad);
			}
			else if (dimension == 3)
			{
				MathUtils<double>::InvertMatrix3(nodalFgrad, InvFgrad, detFgrad);
			}

			// it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
			SpatialVelocityGrad = prod(FgradVel, InvFgrad);

			if (dimension == 2)
			{

				itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] = SpatialVelocityGrad(0, 0);
				itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] = SpatialVelocityGrad(1, 1);
				itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));

				itNode->FastGetSolutionStepValue(SOLID_NODAL_EQUIVALENT_STRAIN_RATE) = sqrt((2.0 * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] +
																							 2.0 * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] +
																							 4.0 * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2] * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2]));

				const double DefX = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0];
				const double DefY = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1];

				const double DefVol = DefX + DefY;

				itNode->GetSolutionStepValue(SOLID_NODAL_VOLUMETRIC_DEF_RATE) = DefVol;
			}
			else if (dimension == 3)
			{

				itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] = SpatialVelocityGrad(0, 0);
				itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] = SpatialVelocityGrad(1, 1);
				itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2] = SpatialVelocityGrad(2, 2);
				itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[3] = 0.5 * (SpatialVelocityGrad(1, 0) + SpatialVelocityGrad(0, 1));
				itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[4] = 0.5 * (SpatialVelocityGrad(2, 0) + SpatialVelocityGrad(0, 2));
				itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[5] = 0.5 * (SpatialVelocityGrad(2, 1) + SpatialVelocityGrad(1, 2));

				itNode->FastGetSolutionStepValue(SOLID_NODAL_EQUIVALENT_STRAIN_RATE) = sqrt(2.0 * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0] +
																							2.0 * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1] +
																							2.0 * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2] * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2] +
																							4.0 * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[3] * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[3] +
																							4.0 * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[4] * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[4] +
																							4.0 * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[5] * itNode->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[5]);

				const double DefX = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[0];
				const double DefY = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[1];
				const double DefZ = itNode->GetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE)[2];

				const double DefVol = DefX + DefY + DefZ;

				itNode->GetSolutionStepValue(SOLID_NODAL_VOLUMETRIC_DEF_RATE) = DefVol;
			}
		}

		void CalcNodalStrains()
		{

			/* std::cout << "Calc Nodal Strains  " << std::endl; */
			ModelPart &rModelPart = BaseType::GetModelPart();

			const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			// #pragma omp parallel
			//     {
			ModelPart::NodeIterator NodesBegin;
			ModelPart::NodeIterator NodesEnd;
			OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

			for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
			{

				const double nodalVolume = itNode->FastGetSolutionStepValue(NODAL_VOLUME);
				const double solidNodalVolume = itNode->FastGetSolutionStepValue(SOLID_NODAL_VOLUME);

				double theta = 1.0;

				if (itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true)
				{

					if (nodalVolume > 0)
					{
						// I have to compute the strains two times because one time is for the solid and the other for the fluid
						Vector nodalSFDneighboursId = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
						Vector rNodalSFDneigh = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);

						Matrix &interfaceFgrad = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
						Matrix &interfaceFgradVel = itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);

						if (interfaceFgrad.size1() != dimension)
							interfaceFgrad.resize(dimension, dimension, false);

						if (interfaceFgradVel.size1() != dimension)
							interfaceFgradVel.resize(dimension, dimension, false);

						noalias(interfaceFgrad) = ZeroMatrix(dimension, dimension);
						noalias(interfaceFgradVel) = ZeroMatrix(dimension, dimension);

						// Matrix interfaceFgrad    = ZeroMatrix(dimension,dimension);
						// Matrix interfaceFgradVel = ZeroMatrix(dimension,dimension);
						// the following function is more expensive than the general one because there is one loop more over neighbour nodes. This is why I do it here also for fluid interface nodes.
						ComputeAndStoreNodalDeformationGradientForInterfaceNode(itNode, nodalSFDneighboursId, rNodalSFDneigh, theta, interfaceFgrad, interfaceFgradVel);
						// itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)=interfaceFgrad;
						// itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL)=interfaceFgradVel;
						this->CalcNodalStrainsForNode(itNode);
					}

					if (solidNodalVolume > 0)
					{
						Vector solidNodalSFDneighboursId = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
						Vector rSolidNodalSFDneigh = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS);

						Matrix &solidInterfaceFgrad = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD);
						Matrix &solidInterfaceFgradVel = itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL);

						if (solidInterfaceFgrad.size1() != dimension)
							solidInterfaceFgrad.resize(dimension, dimension, false);

						if (solidInterfaceFgradVel.size1() != dimension)
							solidInterfaceFgradVel.resize(dimension, dimension, false);

						noalias(solidInterfaceFgrad) = ZeroMatrix(dimension, dimension);
						noalias(solidInterfaceFgradVel) = ZeroMatrix(dimension, dimension);

						// Matrix solidInterfaceFgrad       = ZeroMatrix(dimension,dimension);
						// Matrix solidInterfaceFgradVel    = ZeroMatrix(dimension,dimension);
						ComputeAndStoreNodalDeformationGradientForInterfaceNode(itNode, solidNodalSFDneighboursId, rSolidNodalSFDneigh, theta, solidInterfaceFgrad, solidInterfaceFgradVel);
						// itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD)=solidInterfaceFgrad;
						// itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL)=solidInterfaceFgradVel;
						// CalcNodalStrainsForInterfaceSolidNode(itNode);
						CalcNodalStrainsForSolidNode(itNode);
					}
				}
				else
				{
					if (itNode->Is(SOLID) && solidNodalVolume > 0)
					{
						ComputeAndStoreNodalDeformationGradientForSolidNode(itNode, theta);
						CalcNodalStrainsForSolidNode(itNode);
					}
					else if (nodalVolume > 0)
					{
						this->ComputeAndStoreNodalDeformationGradient(itNode, theta);
						this->CalcNodalStrainsForNode(itNode);
					}
				}
				if (nodalVolume == 0 && solidNodalVolume == 0)
				{ // if nodalVolume==0
					this->InitializeNodalVariablesForRemeshedDomain(itNode);
					InitializeNodalVariablesForSolidRemeshedDomain(itNode);
				}
			}
			// }
		}

		void ComputeAndStoreNodalDeformationGradientForSolidNode(ModelPart::NodeIterator itNode, double theta)
		{

			KRATOS_TRY;

			ModelPart &rModelPart = BaseType::GetModelPart();
			const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
			Vector nodalSFDneighboursId = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
			Vector rNodalSFDneigh = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS);
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
					for (unsigned int i = 0; i < neighSize - 1; i++) // neigh_nodes has one cell less than nodalSFDneighboursId becuase this has also the considered node ID at the beginning
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

			itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD) = Fgrad;
			itNode->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL) = FgradVel;
			KRATOS_CATCH("");
		}

		void ComputeAndStoreNodalDeformationGradientForInterfaceNode(ModelPart::NodeIterator itNode, Vector nodalSFDneighboursId, Vector rNodalSFDneigh, double theta, Matrix &Fgrad, Matrix &FgradVel)
		{

			KRATOS_TRY;

			ModelPart &rModelPart = BaseType::GetModelPart();
			const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
			const unsigned int neighSize = nodalSFDneighboursId.size();
			noalias(Fgrad) = ZeroMatrix(dimension, dimension);
			noalias(FgradVel) = ZeroMatrix(dimension, dimension);
			NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
			const unsigned int neighNodesSize = neighb_nodes.size();

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
					for (unsigned int i = 0; i < neighSize - 1; i++) // neigh_nodes has one cell less than nodalSFDneighboursId becuase this has also the considered node ID at the beginning
					{
						unsigned int other_neigh_nodes_id = nodalSFDneighboursId[i + 1];
						for (unsigned int k = 0; k < neighNodesSize; k++)
						{
							unsigned int neigh_nodes_id = neighb_nodes[k].Id();
							if (neigh_nodes_id == other_neigh_nodes_id)
							{
								dNdXi = rNodalSFDneigh[firstRow];
								dNdYi = rNodalSFDneigh[firstRow + 1];
								Fgrad(0, 0) += dNdXi * neighb_nodes[k].X();
								Fgrad(0, 1) += dNdYi * neighb_nodes[k].X();
								Fgrad(1, 0) += dNdXi * neighb_nodes[k].Y();
								Fgrad(1, 1) += dNdYi * neighb_nodes[k].Y();

								VelocityX = neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_X, 0) * theta + neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_X, 1) * (1 - theta);
								VelocityY = neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_Y, 0) * theta + neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_Y, 1) * (1 - theta);

								FgradVel(0, 0) += dNdXi * VelocityX;
								FgradVel(0, 1) += dNdYi * VelocityX;
								FgradVel(1, 0) += dNdXi * VelocityY;
								FgradVel(1, 1) += dNdYi * VelocityY;

								firstRow += 2;
								break;
							}
						}
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
						unsigned int other_neigh_nodes_id = nodalSFDneighboursId[i + 1];
						for (unsigned int k = 0; k < neighNodesSize; k++)
						{
							unsigned int neigh_nodes_id = neighb_nodes[k].Id();
							if (neigh_nodes_id == other_neigh_nodes_id)
							{

								dNdXi = rNodalSFDneigh[firstRow];
								dNdYi = rNodalSFDneigh[firstRow + 1];
								dNdZi = rNodalSFDneigh[firstRow + 2];

								VelocityX = neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_X, 0) * theta + neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_X, 1) * (1 - theta);
								VelocityY = neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_Y, 0) * theta + neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_Y, 1) * (1 - theta);
								VelocityZ = neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_Z, 0) * theta + neighb_nodes[k].FastGetSolutionStepValue(VELOCITY_Z, 1) * (1 - theta);

								Fgrad(0, 0) += dNdXi * neighb_nodes[k].X();
								Fgrad(0, 1) += dNdYi * neighb_nodes[k].X();
								Fgrad(0, 2) += dNdZi * neighb_nodes[k].X();

								Fgrad(1, 0) += dNdXi * neighb_nodes[k].Y();
								Fgrad(1, 1) += dNdYi * neighb_nodes[k].Y();
								Fgrad(1, 2) += dNdZi * neighb_nodes[k].Y();

								Fgrad(2, 0) += dNdXi * neighb_nodes[k].Z();
								Fgrad(2, 1) += dNdYi * neighb_nodes[k].Z();
								Fgrad(2, 2) += dNdZi * neighb_nodes[k].Z();

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
								break;
							}
						}
					}
				}
			}

			KRATOS_CATCH("");
		}

		void UpdateTopology(ModelPart &rModelPart, unsigned int echoLevel)
		{
			KRATOS_TRY;

			CalculateDisplacementsAndResetNodalVariables();
			BaseType::MoveMesh();
			BoundaryNormalsCalculationUtilities BoundaryComputation;
			BoundaryComputation.CalculateUnitBoundaryNormals(rModelPart, echoLevel);

			KRATOS_CATCH("");
		}

		void CalculateDisplacementsAndResetNodalVariables()
		{
			ModelPart &rModelPart = BaseType::GetModelPart();
			const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
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

				Vector &rSolidNodalSFDneighbours = i->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS);
				unsigned int solidSizeSDFNeigh = rSolidNodalSFDneighbours.size();

				i->FastGetSolutionStepValue(SOLID_NODAL_VOLUME) = 0;
				i->FastGetSolutionStepValue(SOLID_NODAL_MEAN_MESH_SIZE) = 0;
				i->FastGetSolutionStepValue(SOLID_NODAL_FREESURFACE_AREA) = 0;
				i->FastGetSolutionStepValue(SOLID_NODAL_VOLUMETRIC_DEF_RATE) = 0;
				i->FastGetSolutionStepValue(SOLID_NODAL_EQUIVALENT_STRAIN_RATE) = 0;

				noalias(rSolidNodalSFDneighbours) = ZeroVector(solidSizeSDFNeigh);

				Vector &rSolidSpatialDefRate = i->FastGetSolutionStepValue(SOLID_NODAL_SPATIAL_DEF_RATE);
				noalias(rSolidSpatialDefRate) = ZeroVector(sizeStrains);

				Matrix &rSolidFgrad = i->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD);
				noalias(rSolidFgrad) = ZeroMatrix(dimension, dimension);

				Matrix &rSolidFgradVel = i->FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD_VEL);
				noalias(rSolidFgradVel) = ZeroMatrix(dimension, dimension);
			}
			//  }
		}

		/// Turn back information as a string.
		std::string Info() const override
		{
			std::stringstream buffer;
			buffer << "NodalTwoStepVPStrategyForFSI";
			return buffer.str();
		}

		/// Print information about this object.
		void PrintInfo(std::ostream &rOStream) const override
		{
			rOStream << "NodalTwoStepVPStrategyForFSI";
		}

		// /// Print object's data.
		// void PrintData(std::ostream& rOStream) const override
		// {

		// }

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

		///@}
		///@name Private  Access
		///@{

		///@}
		///@name Private Inquiry
		///@{

		///@}
		///@name Un accessible methods
		///@{

	private:
	
		void FillNodalSFDVector()
		{

			ModelPart &rModelPart = BaseType::GetModelPart();

			for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
			{
				this->InitializeNodalVariablesForRemeshedDomain(itNode);

				InitializeNodalVariablesForSolidRemeshedDomain(itNode);

				if (itNode->FastGetSolutionStepValue(INTERFACE_NODE) == false)
				{
					this->SetNeighboursOrderToNode(itNode); // it assigns neighbours to inner nodes, filling NODAL_SFD_NEIGHBOURS_ORDER
					if (itNode->Is(SOLID))
					{
						SetNeighboursOrderToSolidNode(itNode); // it assigns neighbours to solid inner nodes, filling SOLID_NODAL_SFD_NEIGHBOURS_ORDER
					}
				}
				else
				{
					SetNeighboursOrderToInterfaceNode(itNode); // it assigns neighbours to interface nodes, filling SOLID_NODAL_SFD_NEIGHBOURS_ORDER for solids and NODAL_SFD_NEIGHBOURS_ORDER for fluids
				}
			}
		}

		/// Assignment operator.
		NodalTwoStepVPStrategyForFSI &operator=(NodalTwoStepVPStrategyForFSI const &rOther) {}

		/// Copy constructor.
		NodalTwoStepVPStrategyForFSI(NodalTwoStepVPStrategyForFSI const &rOther) {}

		///@}

	}; /// Class NodalTwoStepVPStrategyForFSI

	///@}
	///@name Type Definitions
	///@{

	///@}

	///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_NODAL_TWO_STEP_V_P_STRATEGY_H
