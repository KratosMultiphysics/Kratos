//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "custom_constitutive/custom_flow_rules/borja_casm_explicit_plastic_flow_rule.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "pfem_solid_mechanics_application_variables.h"

#include "custom_utilities/stress_invariants_utilities.hpp"

namespace Kratos
{
	// Version of Casm-like hiperelastic model. This is the Hyperelsatic model of Borja, that includes the Houlsby and the constant shear modulus as special cases.

	//************ CONSTRUCTOR ***********
	BorjaCasmExplicitFlowRule::BorjaCasmExplicitFlowRule()
		:NonAssociativeExplicitPlasticFlowRule()
	{
	}

	//*****************************INITIALIZATION CONSTRUCTOR*****************************
	//************************************************************************************

	BorjaCasmExplicitFlowRule::BorjaCasmExplicitFlowRule(YieldCriterionPointer pYieldCriterion)
		:NonAssociativeExplicitPlasticFlowRule(pYieldCriterion)
	{
		std::cout<<"   CASM FLOW RULE constructed "<<std::endl;
	}

  //********* ASSIGMENT OPERATOR
  BorjaCasmExplicitFlowRule& BorjaCasmExplicitFlowRule::operator=(BorjaCasmExplicitFlowRule const& rOther)
  {
		NonAssociativeExplicitPlasticFlowRule::operator=(rOther);
		mPlasticVariables = rOther.mPlasticVariables;
    return *this;
  }

  //********** COPY CONSTRUCTOR *********
  BorjaCasmExplicitFlowRule::BorjaCasmExplicitFlowRule(BorjaCasmExplicitFlowRule const& rOther)
		:NonAssociativeExplicitPlasticFlowRule(rOther),
		mPlasticVariables(rOther.mPlasticVariables)
  {
  }

  //*******   CLONE ********
  FlowRule::Pointer BorjaCasmExplicitFlowRule::Clone() const
  {
		FlowRule::Pointer p_clone(new BorjaCasmExplicitFlowRule(*this));
			return p_clone;
  }

  // ********** DESTRUCTOR **************
  BorjaCasmExplicitFlowRule::~BorjaCasmExplicitFlowRule()
  {
  }

  // ***** Initialize Material (due to the new variable) *****
  void BorjaCasmExplicitFlowRule::InitializeMaterial( YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties) //INHERIT
  {
		mPlasticVariables.clear();
		NonAssociativeExplicitPlasticFlowRule::InitializeMaterial( pYieldCriterion, pHardeningLaw, rMaterialProperties);
	}
	
	void BorjaCasmExplicitFlowRule::InitializeMaterial( const Properties& rMaterialProperties) //INHERIT
	{
		mPlasticVariables.clear();
		NonAssociativeExplicitPlasticFlowRule::InitializeMaterial( rMaterialProperties);
	}


	// *****************************************************************
	// ***** EXPLICIT RETURN MAPPING ALGORITHM (Sloan et al. 2002) *****
	// *****************************************************************
  bool BorjaCasmExplicitFlowRule::CalculateReturnMappingExpl(RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
  {
		//std::cout<< "  " << std::endl;
		//std::cout<< "  ............................ " << std::endl;
		//std::cout<< "  Explicit return mapping CASM " << std::endl;
		//std::cout<< "  ............................ " << std::endl;
		//std::cout<< "    F_incr:        " << rIncrementalDeformationGradient << std::endl;
		//std::cout<< "    b_n+1:         " << rNewElasticLeftCauchyGreen << std::endl;
		//std::cout<< "    StressMatrix:  " << rStressMatrix << std::endl;
		//std::cout<< "    ReturnMappingVariables: " << rReturnMappingVariables << std::endl;
		
    // I. Initialize some variables
//std::cout<< "  ... I. Initialization ... " << std::endl;
    bool PlasticityActive = false;
    rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);
		//get internal plastic variables
		PlasticVariablesType IntPlasticVariables = mPlasticVariables;
//mPlasticVariables.print();

		//reset IncrementalPlasticStrain and Total VolumetricStrain
		rReturnMappingVariables.PlasticVariables.EquivalentPlasticStrain		= IntPlasticVariables.EquivalentPlasticStrain;
		rReturnMappingVariables.PlasticVariables.DeltaEqPlasticStrain				= 0;
		rReturnMappingVariables.PlasticVariables.PlasticShearStrain					= IntPlasticVariables.PlasticShearStrain;
		rReturnMappingVariables.PlasticVariables.DeltaPlasticShearStrain		= 0;
		rReturnMappingVariables.PlasticVariables.PreconsolidationPressure		= IntPlasticVariables.PreconsolidationPressure;
		rReturnMappingVariables.PlasticVariables.Bonding										= 0; //CHANGED

		Matrix PreviousElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
		Vector PreviousStressVector;
		Vector NewStressVector;
		double Tolerance = 1.0E-4; // LMV
		
//std::cout<< "  ... I. Initialization completed ... " << std::endl<<std::endl;

		// II. Check for Yield Condition (at the beginning) and if necessary do yield surface correction
std::cout<< "  ... II. Initial yield condition ... " << std::endl;
		this->CalculateKirchhoffStressVector( PreviousElasticLeftCauchyGreen, PreviousStressVector);
std::cout<< "    b_n:   " << PreviousElasticLeftCauchyGreen << std::endl;
std::cout<< "    tau_n: " << PreviousStressVector << std::endl;
		rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PreviousStressVector, rReturnMappingVariables.PlasticVariables);

		if (rReturnMappingVariables.TrialStateFunction > Tolerance) {
			this->ReturnStressToYieldSurface( rReturnMappingVariables, PreviousElasticLeftCauchyGreen, PreviousStressVector, rReturnMappingVariables.TrialStateFunction, Tolerance);
			//this->ReturnStressToYieldSurfaceNormal( rReturnMappingVariables, PreviousElasticLeftCauchyGreen, PreviousStressVector, rReturnMappingVariables.TrialStateFunction, Tolerance);
			//THE INITIAL TENSIONAL STATE MAY BE OUTSIDE THE YIELD SURFACE?? (yes and no: no if you don't do nothing strange, however, if you apply LSInterpolation or change constant parameters during the simulation YES)
		}
		rReturnMappingVariables.DeltaBeta = 0;
		
std::cout<< "  ... II. Initial yield condition checked ... " << std::endl<<std::endl;


		// III. Compute an ELASTIC TRIAL STEP
std::cout<< "  ... III. Elastic trial step ... " << std::endl;
		ExplicitStressUpdateInformation  StressUpdateInformation;
		this->CalculateOneExplicitStep(rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, false, StressUpdateInformation);
		
		double ElasticTrialStateFunction;
//IntPlasticVariables.print();
		ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, NewStressVector, IntPlasticVariables);
std::cout<< "  ... III. Elastic trial step computed: tau_trial = "<< NewStressVector << " f = " << ElasticTrialStateFunction <<" ... "<< std::endl << std::endl;
		

		// IVa. PURELY ELASTIC STEP
		if ( ElasticTrialStateFunction <= Tolerance)   {
			PlasticityActive = false;
			
std::cout<< "  ... IVa. Purely elastic step computed " << std::endl;
		}
		// IVb. SOME PART OF THE STEP IS PLASTIC
		else {
			PlasticityActive = true; 

			if ( ( rReturnMappingVariables.TrialStateFunction < -Tolerance) && (ElasticTrialStateFunction > Tolerance) ) {
				//2b.1 ELASTIC & ELASTO-PLASTIC STEP 
//std::cout<< "  ... IVb.1 Elastic & elasto-plastic step " << std::endl;
				this->CalculateExplicitSolutionWithChange( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, Tolerance);
//rReturnMappingVariables.PlasticVariables.print();
std::cout<< "  ... IVb.1 Elastic & elasto-plastic step computed " << std::endl<<std::endl;
			}
			else {         
				bool UnloadingCondition;
				UnloadingCondition =  this->EvaluateElastoPlasticUnloadingCondition( UnloadingCondition, PreviousElasticLeftCauchyGreen, rIncrementalDeformationGradient, IntPlasticVariables, Tolerance); 

				if (UnloadingCondition ) {
					//2c.1 ELASTO-PLASTIC UNLOADING
//std::cout<< "  ... IVb.2 Elasto-plastic unloading step " << std::endl;
					this->CalculateExplicitSolutionWithChange( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, Tolerance);
//rReturnMappingVariables.PlasticVariables.print();
std::cout<< "  ... IVb.2 Elasto-plastic unloading step computed " << std::endl<<std::endl;
				}
				else {
					//2c. 2 ELASTO-PLASTIC STEP
//std::cout<< "  ... IVb.3 Elasto-plastic loading step " << std::endl;
					this->CalculateExplicitSolution( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, true, Tolerance);
//rReturnMappingVariables.PlasticVariables.print();
std::cout<< "  ... IVb.3 Elasto-plastic loading step computed " << std::endl<<std::endl;
				}
			}
		}


		// V. YIELD SURFACE DRIFT CORRECTION AFTER STRESS INTEGRATION
		if (PlasticityActive) {

			double DriftViolation;

			DriftViolation = mpYieldCriterion->CalculateYieldCondition(DriftViolation, NewStressVector, rReturnMappingVariables.PlasticVariables);
			//rReturnMappingVariables.EigenVectors = NewElasticLeftCauchyGreen; 
			double Beta = rReturnMappingVariables.DeltaBeta;
			if ( fabs(DriftViolation) > Tolerance ) {
//std::cout<< "  ... IVb. Yield surface drift violation: "<< fabs(DriftViolation) << std::endl;
				this->ReturnStressToYieldSurface( rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, DriftViolation, Tolerance);
//rReturnMappingVariables.PlasticVariables.print();
			}
			// LMV
			rReturnMappingVariables.DeltaBeta = Beta;
			
std::cout<< "  ... IVb. Yield surface drift corrected " << std::endl;
		}

		//set flags and compute rStressMatrix
		rStressMatrix = MathUtils<double>::StressVectorToTensor(NewStressVector);

		rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
		rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

		return PlasticityActive;
  }


  // ************************************************************************************
  // ***** EXPLICIT PART: Elasto-plastic time integration with adaptive substepping ***** INHERIT
  // ************************************************************************************
  void BorjaCasmExplicitFlowRule::CalculateExplicitSolution( const Matrix& rIncrementalDeformationGradient, const Matrix& rPreviousElasticCauchyGreen, RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, const double& rTolerance)
  {
		KRATOS_TRY

		double TimeStep = 0.5;
		double MinTimeStep = 1.0e-4;
		double DoneTimeStep = 0.0;
		double MaxTimeStep = 0.5;
		TimeStep = MaxTimeStep;

		Matrix ActualElasticLeftCauchyGreen = rPreviousElasticCauchyGreen;
		Matrix SubstepDeformationGradient; 

		//StressUpdateInformation to store data during time integration
		ExplicitStressUpdateInformation  StressUpdateInformation;
		StressUpdateInformation.clear();
		bool MayBeLast = false;

		//TIME INTEGRATION
		//with adaptive substepping: increment from n+alpha+doneStep to n+alpha+doneStep+TimeStep until n+1
		while (DoneTimeStep < 1.0)  {
//std::cout << " DONE TIME: " << DoneTimeStep << ", DELTA TIME: " << TimeStep  << std::endl;
			if (DoneTimeStep + TimeStep >= 1.0) {
				TimeStep = 1.0 - DoneTimeStep;
				MayBeLast = true;
			}

			this->ComputeSubstepIncrementalDeformationGradient( rIncrementalDeformationGradient, DoneTimeStep, DoneTimeStep + TimeStep, SubstepDeformationGradient);

			this->CalculateOneExplicitStep( SubstepDeformationGradient, ActualElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, rNewStressVector,  rElastoPlasticBool, StressUpdateInformation);

			//substep converged: save data & go to next pseudo-time step
			if ( StressUpdateInformation.StressErrorMeasure < rTolerance ) {
//std::cout << "   SUbstepping " << DoneTimeStep << " dt " << TimeStep << "Error " << StressUpdateInformation.StressErrorMeasure << std::endl;
					
				ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
				this->UpdateRadialReturnVariables( rReturnMappingVariables, StressUpdateInformation);

				DoneTimeStep += TimeStep;
				//return is final increment from n+1-TimeStep to n+1
				if (MayBeLast == true)
					return;
			}
			//if minimum adaptive time step is reached & not converged: throw error
			else {
				if (TimeStep <= MinTimeStep) {
					if ( StressUpdateInformation.StressErrorMeasure > rTolerance) {
						std::cout << "ExplicitStressIntegrationDidNotConverge: StressError: " << StressUpdateInformation.StressErrorMeasure << " dt " << TimeStep <<  " meanStress " << (rNewStressVector(0) + rNewStressVector(1) + rNewStressVector(2) ) / 3.0 << std::endl ;
						std::cout << "  rIncrementalDeformationGradient: " << rIncrementalDeformationGradient << std::endl;
						std::cout << "  rPreviousElasticCauchyGreen: " << rPreviousElasticCauchyGreen << std::endl;
						std::cout << "  rNewElasticLeftCauchyGreen: " << rNewElasticLeftCauchyGreen << std::endl;
						std::cout << "  rNewStressVector: " << rNewStressVector << std::endl;
						double Hola;
						Hola = mpYieldCriterion->CalculateYieldCondition(Hola, rNewStressVector, StressUpdateInformation.PlasticVariablesUpdate);
						std::cout << "  rYieldValue: " << Hola << std::endl;
					}

					DoneTimeStep += TimeStep;

					ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
					this->UpdateRadialReturnVariables( rReturnMappingVariables, StressUpdateInformation);
					if (MayBeLast == true )
						return;
				}
				else if ( std::isnan(TimeStep) ) {

						std::cout << "ExplicitStressIntegrationDidNotConverge: StressError: " << StressUpdateInformation.StressErrorMeasure << " dt " << TimeStep <<  " meanStress " << (rNewStressVector(0) + rNewStressVector(1) + rNewStressVector(2) ) / 3.0 << std::endl ;
						std::cout << "  rIncrementalDeformationGradient: " << rIncrementalDeformationGradient << std::endl;
						std::cout << "  rPreviousElasticCauchyGreen: " << rPreviousElasticCauchyGreen << std::endl;
						std::cout << "  rNewElasticLeftCauchyGreen: " << rNewElasticLeftCauchyGreen << std::endl;
						std::cout << "  StressUpdateInformation: " << StressUpdateInformation.DeltaDeltaPlastic << std::endl;
						std::cout << "  rNewStressVector: " << rNewStressVector << std::endl;
						double Hola;
						Hola = mpYieldCriterion->CalculateYieldCondition(Hola, rNewStressVector, StressUpdateInformation.PlasticVariablesUpdate);
						std::cout << "  rYieldValue: " << Hola << std::endl;

						return;
				}
			}

			//compute new adaptive step increment depending on StressErrorMeasure
			TimeStep *= 0.9* (pow( rTolerance / (StressUpdateInformation.StressErrorMeasure+ 1e-8), 0.5 ));
			TimeStep = std::max(TimeStep, MinTimeStep);
			TimeStep = std::min(TimeStep, MaxTimeStep);
				         
//std::cout<<"rTolerance: "<<rTolerance<<std::endl;
//std::cout<<"StressErrorMeasure: "<<StressUpdateInformation.StressErrorMeasure<<std::endl;
//std::cout<<"actual time step: "<<TimeStep<<std::endl;
				
		}
		
    KRATOS_CATCH("")
  }


  // ******************************************************************************************
	// ***** EXPLICIT PART: Find yield surface intersection and compute elasto-plastic step ***** INHERIT
	// ******************************************************************************************
  void BorjaCasmExplicitFlowRule::CalculateExplicitSolutionWithChange(const Matrix& rDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector,  const double& rTolerance)
  {
//std::cout<<"    ... Calculating explicit solution with change ... "<<std::endl;
		//bisection method to find yield surface intersection
		double InitialPosition = 0.0;
		double EndPosition = 1.0;
		double HalfPosition;

		double InitialStateFunction = -1.0;
		double EndStateFunction = 1.0;
		double HalfPositionStateFunction;

		Matrix HalfPositionDeformationGradient;
		Matrix HalfPositionElasticCauchyGreen;

		ExplicitStressUpdateInformation StressUpdateInformation;

		double ErrorMeasure1;
		double ErrorMeasure2;

		for (unsigned int i = 0 ; i< 150; ++i) {

			HalfPosition = 0.5*(InitialPosition + EndPosition);

			this->ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.0, HalfPosition, HalfPositionDeformationGradient);

			this->CalculateOneExplicitStep( HalfPositionDeformationGradient, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables, HalfPositionElasticCauchyGreen, rNewStressVector, false, StressUpdateInformation);  
			HalfPositionStateFunction = mpYieldCriterion->CalculateYieldCondition( HalfPositionStateFunction, rNewStressVector, rReturnMappingVariables.PlasticVariables);

			if ( HalfPositionStateFunction < 0.0)  {
				InitialStateFunction = HalfPositionStateFunction;
				InitialPosition = HalfPosition;
			}
			else {
				EndPosition = HalfPosition;
				EndStateFunction = HalfPositionStateFunction;
			}
			ErrorMeasure1 = fabs(InitialStateFunction-EndStateFunction);
			ErrorMeasure2 = fabs(InitialPosition - EndPosition);

			if ( (ErrorMeasure1 < rTolerance ) && (ErrorMeasure2 < rTolerance)) {
				break;
			}
		}
		
//std::cout<<"    ... Yield surface intersection found ... "<<std::endl;

		// COMPUTE ELASTIC STEP;
		this->ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.0, HalfPosition, HalfPositionDeformationGradient);
		
		this->CalculateOneExplicitStep( HalfPositionDeformationGradient, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables, HalfPositionElasticCauchyGreen, rNewStressVector, false, StressUpdateInformation);

//std::cout<<"    ... Elastic step computed: F_el = "<<HalfPositionDeformationGradient<<", F = "<<rDeformationGradient<<" ... "<<std::endl;

		// COMPUTE ELASTO-PLASTIC STEP
		this->ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, HalfPosition, 1, HalfPositionDeformationGradient); 
		
//std::cout<<"      ... Plastic sub gradient F_pl = "<<HalfPositionDeformationGradient<<" ... "<<std::endl;

		this->CalculateExplicitSolution( HalfPositionDeformationGradient, HalfPositionElasticCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, rNewStressVector,  true, rTolerance);
  
//std::cout<<"    ... Elasto-plastic step computed ... "<<std::endl;
  }


  // *********************************************************************
  // ***** EXPLICIT PART: Update one elasto-plastic strain increment ***** INHERIT
  // *********************************************************************
  void BorjaCasmExplicitFlowRule::CalculateOneExplicitPlasticStep(const Matrix& rDeltaDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen,
   const PlasticVariablesType& rPreviousPlasticVariables, Matrix& rNewElasticLeftCauchyGreen, PlasticVariablesType& rNewPlasticVariables, double& rDeltaPlastic)
  {
    Matrix ElasticMatrix;
    Vector ElasticStrainVector;
	
		//calculate elastic matrix De
    ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rPreviousElasticLeftCauchyGreen);
    this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);

		//calculate yield function and plastic potential derivatives
    AuxiliarDerivativesStructure AuxiliarDerivatives;
    this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rPreviousPlasticVariables);

		//calculate hardening modulus H
    double H;
    this->ComputePlasticHardeningParameter(ElasticStrainVector, rPreviousPlasticVariables, H);

		//calculate incremtental Hencky strain vector eps
    ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient)) );

		//calculate plastic multiplier DeltaGamma
    Vector auxVector;
    auxVector = prod(ElasticMatrix, ElasticStrainVector);
    double DeltaGamma = 0.0;
    for (unsigned int i = 0; i<6; ++i)
      DeltaGamma += auxVector(i)*AuxiliarDerivatives.YieldFunctionD(i);

    double auxDenominador = H + MathUtils<double>::Dot( AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD));
//std::cout<<"---------------------- df: "<<AuxiliarDerivatives.YieldFunctionD<<std::endl;
//std::cout<<"---------------------- dg: "<<AuxiliarDerivatives.PlasticPotentialD<<std::endl;
//std::cout<<"---------------------- H: "<<H<<", K_p: "<<auxDenominador<<"    nom: "<<DeltaGamma<<", gamma: "<<DeltaGamma/auxDenominador<<std::endl<<std::endl;
    //DeltaGamma /= auxDenominador;

    if ( std::isnan(DeltaGamma/auxDenominador) )
    {
    	std::cout<<std::endl<<"NaN in plastic multiplier: "<<std::endl;
    	std::cout<<"  df: "<< AuxiliarDerivatives.YieldFunctionD <<std::endl;
    	std::cout<<"  dg: "<< AuxiliarDerivatives.PlasticPotentialD <<std::endl;
    	std::cout<<"  H: "<<H<<", K_p: "<<auxDenominador<<", nom: "<<DeltaGamma<<", gamma: "<<DeltaGamma/auxDenominador<<std::endl<<std::endl;
    }

    DeltaGamma /= auxDenominador;

    if (DeltaGamma < 0)// || std::isnan(DeltaGamma) )
      DeltaGamma = 0;

    rDeltaPlastic = DeltaGamma;

	//calculate rNewElasticLeftCauchyGreen based on exponential variation of the plastic deformation gradient (Simo 1998)
    Vector MuyAuxiliar;
    MuyAuxiliar = -DeltaGamma * AuxiliarDerivatives.PlasticPotentialD / 2.0; 

    Matrix UpdateMatrix;
    UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor (MuyAuxiliar);
    UpdateMatrix = prod( rDeltaDeformationGradient, UpdateMatrix);

    rNewElasticLeftCauchyGreen = prod(UpdateMatrix, rPreviousElasticLeftCauchyGreen);
    rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans(UpdateMatrix));

		//calculate incremental equivalent plastic strain
    double DeltaAlpha = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
      DeltaAlpha += DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i);

		//calculate incremental plastic shear strain
    double DeltaBeta = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
      DeltaBeta += pow( DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i) - DeltaAlpha/3.0, 2);

    for (unsigned int i = 3; i < 6; ++i)
      DeltaBeta += 2.0*pow( DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i)/2.0, 2);
    DeltaBeta = sqrt( 2.0/3.0*DeltaBeta );
    
    //update strain variables
    rNewPlasticVariables = rPreviousPlasticVariables;
		rNewPlasticVariables.EquivalentPlasticStrain 	  += DeltaAlpha;
		rNewPlasticVariables.DeltaEqPlasticStrain 			+= DeltaAlpha;
		rNewPlasticVariables.PlasticShearStrain 				+= DeltaBeta;
		rNewPlasticVariables.DeltaPlasticShearStrain 		+= DeltaBeta;
		
		//updated hardening parameters
		mpYieldCriterion->GetHardeningLaw().CalculateHardening(rNewPlasticVariables, DeltaAlpha, DeltaBeta);

		if (rNewPlasticVariables.PreconsolidationPressure > 0 ){
			std::cout<<std::endl<<"P0 > 0 in CalculateOneExplicitPlasticStep"<<std::endl;
		}
  }


  // **************************************************************************
  // ***** EXPLICIT PART: Calculate single elasto-plastic or elastic step ***** INHERIT
  // **************************************************************************
  void BorjaCasmExplicitFlowRule::CalculateOneExplicitStep(const Matrix& rDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, const RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, ExplicitStressUpdateInformation& rStressUpdateInformation)
  {
		//calculate elasto-plastic step with two different discretizations and computes the error
		if ( rElastoPlasticBool)  {
		
			PlasticVariablesType NewPlasticVariables, IntermediatePlasticVariables;
			double DeltaPlastic; 
			Vector FirstApproxStressVector;
 
			//ONE STEP UPDATE -> FirstApproxStressVector	
			this->CalculateOneExplicitPlasticStep( rDeformationGradient, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables.PlasticVariables, rNewElasticLeftCauchyGreen, NewPlasticVariables, DeltaPlastic);
			this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, FirstApproxStressVector); 


			//TWO STEP UPDATE
			Matrix HalfStepIncrementalDeformation;
 
			//first half step
			ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.0, 0.5, HalfStepIncrementalDeformation);
			Matrix SecondApproxLeftCauchyGreen; 
			this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables.PlasticVariables,
			 SecondApproxLeftCauchyGreen, IntermediatePlasticVariables, DeltaPlastic);
			
			rStressUpdateInformation.DeltaDeltaPlastic 					 	= DeltaPlastic;
		
			//second half step -> rNewStressVector
			ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.5, 1.0, HalfStepIncrementalDeformation);
			this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, SecondApproxLeftCauchyGreen, IntermediatePlasticVariables,
			 rNewElasticLeftCauchyGreen, NewPlasticVariables, DeltaPlastic);
			this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector);
			
			rStressUpdateInformation.DeltaDeltaPlastic 					 += DeltaPlastic;
			rStressUpdateInformation.PlasticVariablesUpdate				= NewPlasticVariables;

			//error measure calculated from the two obtained stress vectors (FirstApproxStressVector & rNewStressVector)
			rStressUpdateInformation.StressErrorMeasure = 0.0;
			double Denominador = 0.0;
			for (unsigned int i = 0; i<6; ++i) {
				rStressUpdateInformation.StressErrorMeasure += pow( rNewStressVector(i) - FirstApproxStressVector(i), 2);
				Denominador  += pow(rNewStressVector(i), 2);
			}
			rStressUpdateInformation.StressErrorMeasure = 1.0* pow( rStressUpdateInformation.StressErrorMeasure/Denominador, 0.5);
		}
		//or calculate purely elastic step
		else  {
			rNewElasticLeftCauchyGreen = prod( rPreviousElasticLeftCauchyGreen, trans( rDeformationGradient) );
			rNewElasticLeftCauchyGreen = prod( rDeformationGradient, rNewElasticLeftCauchyGreen );

			rStressUpdateInformation.StressErrorMeasure = 0.0;

			this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector);
		}
  }


  // ******************************************************
	// ***** DRIFT CORRECTION : CUTTING PLANE ALGORITHM ***** INHERIT
  // ******************************************************
  void BorjaCasmExplicitFlowRule::ReturnStressToYieldSurface( RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rStressVector, double& rDrift, const double& rTolerance)
  {
		AuxiliarDerivativesStructure  AuxiliarDerivatives;
		Matrix ElasticMatrix;
		double H;
		double DeltaGamma;
		Vector ElasticCorrection;
		Matrix CorrectedLeftCauchyGreen;
		Matrix UpdateMatrix;
		Vector ActualElasticHenckyStrain;
		Vector StressVector;
		rReturnMappingVariables.DeltaBeta = 0.0;

		for (unsigned int i = 0; i < 150; ++i) {

			ActualElasticHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain( rNewElasticLeftCauchyGreen);

			this->UpdateDerivatives(ActualElasticHenckyStrain, AuxiliarDerivatives, rReturnMappingVariables.PlasticVariables);
			this->ComputeElasticMatrix(ActualElasticHenckyStrain, ElasticMatrix);
			this->ComputePlasticHardeningParameter(ActualElasticHenckyStrain, rReturnMappingVariables.PlasticVariables, H);
			
//std::cout<<std::endl<< AuxiliarDerivatives.YieldFunctionD<<std::endl;
//std::cout<< AuxiliarDerivatives.PlasticPotentialD<<std::endl;
//rReturnMappingVariables.PlasticVariables.print();

			//correction multiplier so that total strain increment stays unchanged
			DeltaGamma = rDrift;
			DeltaGamma /= ( H + MathUtils<double>::Dot(AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD)));

			//calculate new elastic Left Cauchy Green tensor
			Vector MuyAuxiliar = -DeltaGamma*AuxiliarDerivatives.PlasticPotentialD/ 2.0;
			UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor( MuyAuxiliar);
			CorrectedLeftCauchyGreen =  prod((UpdateMatrix), rNewElasticLeftCauchyGreen);
			CorrectedLeftCauchyGreen =  prod( CorrectedLeftCauchyGreen, trans(UpdateMatrix));

			//calculate new incremental plastic strains
			double AlphaUpdate = 0.0;
			for (unsigned int j = 0; j < 3; ++j)
				AlphaUpdate += DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j);

			double BetaUpdate = 0.0;
			for (unsigned int j = 0; j < 3; ++j)
				BetaUpdate += pow(DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j) - AlphaUpdate/3.0, 2);

			for (unsigned int j = 3; j < 6; ++j)
				BetaUpdate += 2.0*pow(DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j) / 2.0, 2);
			BetaUpdate = sqrt(2.0/3.0*BetaUpdate);
//std::cout<<" AlphaUpdate: "<<AlphaUpdate<<",  BetaUpdate: "<<BetaUpdate<<",  DeltaGamma: "<<DeltaGamma<< ",  H: "<<H<<std::endl;
			//update strain variables
			rReturnMappingVariables.PlasticVariables.EquivalentPlasticStrain += AlphaUpdate;
			rReturnMappingVariables.PlasticVariables.DeltaEqPlasticStrain += AlphaUpdate;
			rReturnMappingVariables.PlasticVariables.PlasticShearStrain += BetaUpdate;
			rReturnMappingVariables.PlasticVariables.DeltaPlasticShearStrain += BetaUpdate;
//rReturnMappingVariables.PlasticVariables.print();
			
			//calculate updated hardening parameters
			//this->CalculateHardening(rReturnMappingVariables.PlasticVariables, AlphaUpdate, BetaUpdate);
			mpYieldCriterion->GetHardeningLaw().CalculateHardening(rReturnMappingVariables.PlasticVariables, AlphaUpdate, BetaUpdate);

			if (rReturnMappingVariables.PlasticVariables.PreconsolidationPressure > 0 ){
			std::cout<<std::endl<<"P0 > 0 in ReturnStressToYieldSurface"<<std::endl;
		}
			
			rNewElasticLeftCauchyGreen = CorrectedLeftCauchyGreen;

			//check yield surface violation
			this->CalculateKirchhoffStressVector( CorrectedLeftCauchyGreen, rStressVector);
			rDrift = mpYieldCriterion->CalculateYieldCondition( rDrift, rStressVector, rReturnMappingVariables.PlasticVariables);

//std::cout<<rDrift<<std::endl;

			if( fabs(rDrift) < 1.0*rTolerance) {
				return;
			}
		}
		
		//leaving drift correction after 150 runs
		if ( fabs(rDrift) > 150.0*rTolerance) {
			std::cout << " " << std::endl;
			std::cout << "Leaving consistent drift correction WITHOUT converging"<< std::endl;
			std::cout << "  Drift: " 				<< rDrift << std::endl;
			std::cout << "  StressVector: " << rStressVector << std::endl;
			std::cout << "  epsVol: " 			<< rReturnMappingVariables.PlasticVariables.EquivalentPlasticStrain << std::endl;
			std::cout << "  epsShear: " 		<< rReturnMappingVariables.PlasticVariables.PlasticShearStrain << std::endl;
			
			//compute Lode angle to see what happens
			Matrix SM = ZeroMatrix(3,3);
			SM = MathUtils<double>::StressVectorToTensor( rStressVector);
			double p = 0.0;
			for (unsigned int i = 0; i < 3; ++i)
				p += SM(i,i) / 3.0;
			double J2 = 0.0;
			for (unsigned int i = 0; i < 3; ++i) {
				SM(i,i) -= p;     
				J2 += pow( SM(i,i)- p, 2);
			}
			for (unsigned int i = 0; i < 6; ++i)
				J2 += 2.0*pow( rStressVector(i), 2);
			J2 = sqrt(J2 / 2.0);

			double Lode = MathUtils<double>::Det( SM);
			Lode *= 3.0*sqrt(3.0) / 2.0;
			Lode /= pow( J2, 3);

			std::cout << "  LODE: " << std::asin( Lode) / 3.0 / 3.14159265359*180.0 << std::endl;
		}
  }


  // *************************************************
	// ***** COMPUTE ELASTO-PLASTIC TANGENT MATRIX ***** INHERIT
	// *************************************************
  void BorjaCasmExplicitFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rLeftCauchyGreenMatrix, const double& rAlpha, Matrix& rElasticMatrix)
  {

		int StressIntTechnique = 0; //mpYieldCriterion->GetHardeningLaw().GetProperties()[STRESS_INTEGRATION_STRATEGY];

		//compute Hencky strain from Left Cauchy Green
		Vector ElasticStrainVector = ZeroVector(6);
		ElasticStrainVector = this->ConvertCauchyGreenTensorToHenckyStrain(rLeftCauchyGreenMatrix);

		//compute elastic matrix
		this->ComputeElasticMatrix(ElasticStrainVector, rElasticMatrix);

		//IMPLEX
		//if ( StressIntTechnique == 3) // mIMPLEX // this case allways goes with the elastic matrix.
		if ( rReturnMappingVariables.Options.Is(IMPLEX_ACTIVE) ) {
			 return; 
		}
		else if ( StressIntTechnique == 4) // IMPLEX // this case may go with elastic or the "second derivative plus"
		{
			if ( fabs( rReturnMappingVariables.NormIsochoricStress) < 1e-8 ) // ONLY ELASTIC. No Plastic part
			{
				return;
			}

			AuxiliarDerivativesStructure AuxiliarDerivatives;
			this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rReturnMappingVariables.PlasticVariables);

			Matrix SystemMatrix = rElasticMatrix;

			bool singular = SolidMechanicsMathUtilities<double>::InvertMatrix( rElasticMatrix, SystemMatrix);
			if (  ! singular)
			{
				SystemMatrix += rReturnMappingVariables.NormIsochoricStress *  AuxiliarDerivatives.PlasticPotentialDD;

				singular = SolidMechanicsMathUtilities<double>::InvertMatrix( SystemMatrix, SystemMatrix);
			}
			if ( !singular)
				rElasticMatrix = SystemMatrix;

			return;
		}

		//IMPLICIT OR EXPLICIT
		if ( rReturnMappingVariables.Options.Is( FlowRule::PLASTIC_REGION ) )
		{
			AuxiliarDerivativesStructure AuxiliarDerivatives;
			this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rReturnMappingVariables.PlasticVariables);

			//IMPLICIT CASE: add second derivative to the Elastic Matrix
			if ( StressIntTechnique == 2) // IMPLICIT CASE
			{
				Matrix SystemMatrix = rElasticMatrix;

				bool singular = SolidMechanicsMathUtilities<double>::InvertMatrix( rElasticMatrix, SystemMatrix);
				if ( ! singular) {
					SystemMatrix += rReturnMappingVariables.DeltaBeta *  AuxiliarDerivatives.PlasticPotentialDD;

					singular = SolidMechanicsMathUtilities<double>::InvertMatrix( SystemMatrix, SystemMatrix);
				}
				if ( !singular)
				{
					rElasticMatrix = SystemMatrix;
				}
			}

			//compute elasto-plastic matrix
			double H;
			this->ComputePlasticHardeningParameter(ElasticStrainVector, rReturnMappingVariables.PlasticVariables, H);

			Vector AuxVectorF;
			Vector AuxVectorG;

			AuxVectorF = prod(trans(AuxiliarDerivatives.YieldFunctionD), rElasticMatrix);
			AuxVectorG = prod( rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD);

			Matrix PlasticUpdate;
			PlasticUpdate = MyCrossProduct(rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD, AuxiliarDerivatives.YieldFunctionD);

			rElasticMatrix -= 1.0*PlasticUpdate / ( H + MathUtils<double>::Dot(AuxVectorF, AuxiliarDerivatives.PlasticPotentialD));
		}
  }


  // ***** Check unloading condition ***** INHERIT
	bool& BorjaCasmExplicitFlowRule::EvaluateElastoPlasticUnloadingCondition( bool& rUnloadingCondition, const Matrix& rElasticLeftCauchyGreen, const Matrix& rDeltaDeformationGradient, const PlasticVariablesType& rPlasticVariables, const double& rTolerance)
  {

		rUnloadingCondition = false;

		Vector StressVector;
		Vector YieldFunctionD;
		this->CalculateKirchhoffStressVector( rElasticLeftCauchyGreen, StressVector);
		mpYieldCriterion->CalculateYieldFunctionDerivative(StressVector, YieldFunctionD, rPlasticVariables);

		Matrix ElasticMatrix;
		Vector ElasticStrainVector ;
		Vector IncrementalHencky;

		ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rElasticLeftCauchyGreen);
		IncrementalHencky = ConvertCauchyGreenTensorToHenckyStrain( prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient) ));
		this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);

		Vector IncrementalElasticStress;
		IncrementalElasticStress = prod(ElasticMatrix, IncrementalHencky);

		double Numerador = 0.0;
		double Denom1 = 0.0;
		double Denom2 = 0.0;

		for (unsigned int i = 0; i < 6; ++i) {
			 Numerador += YieldFunctionD(i)*IncrementalElasticStress(i);
			 Denom1 += pow( YieldFunctionD(i), 2);
			 Denom2 += pow( IncrementalElasticStress(i), 2);
		}

		Denom1 = pow( Denom1*Denom2, 0.5);

		Numerador /= Denom1;

		if (Numerador < rTolerance) 
			 rUnloadingCondition = true;

		return rUnloadingCondition;

  }


  // ***********************************************
  // ***** Compute plastic hardening modulus H ***** MODIFICATION
  // ***********************************************
/*	void BorjaCasmExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const PlasticVariablesType& rPlasticVariables, double& rH)
	{
		//std::cout<<"BorjaCasmExplicitFlowRule::ComputePlasticHardeningParameter"<<std::endl;

		//get material constants
		const double ShearM 			= mpYieldCriterion->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		const double SpacingR 			= mpYieldCriterion->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN 			= mpYieldCriterion->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];
		const double SwellingSlope		= mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
		const double OtherSlope      	= mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];

		//calculate hardening parameters
		double P0 = rPlasticVariables.PreconsolidationPressure;

		//calcualte Kirchhoff stress vector and invariants
		Vector StressVector = ZeroVector(6);
		this->CalculateKirchhoffStressVector(rHenckyStrainVector, StressVector);
		double MeanStress, J2, LodeAngle;
		StressInvariantsUtilities::CalculateStressInvariants( StressVector, MeanStress, J2, LodeAngle);
		
		//calculate third invariant effect 
		double ThirdInvariantEffect = mpYieldCriterion->EvaluateThirdInvariantEffectMC( LodeAngle);
		
		//calculate derivatives with respect to P & J2
		double FirstDerivativeP = 0;
		double FirstDerivativeJ2 = 0;
		CalculatePlasticPotentialDerivativesPJ2(StressVector, FirstDerivativeP, FirstDerivativeJ2, rPlasticVariables);	
		
		//calculate d_P0/d_gamma
		double dP0dGamma = ( FirstDerivativeP ) * P0/(OtherSlope-SwellingSlope);
		
		//calculate hardening modulus H
		rH = 1/(P0*log(SpacingR)) * ( dP0dGamma );
		rH -= ( ShapeN*pow(J2*pow(3.0,0.5),ShapeN) )/( pow(ShearM/ThirdInvariantEffect,ShapeN)*pow(-P0,ShapeN+1) );

		//std::cout<<"rh: "<<rH<<std::endl;

	}
*/


	// COMPUTE THE PLASTIC HARDENING PARAMETER - COMPARE !!!!!!!!!!!
	void BorjaCasmExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const PlasticVariablesType& rPlasticVariables, double& rH)
	{
//std::cout<<"BorjaCasmExplicitFlowRule::ComputePlasticHardeningParameter"<<std::endl;

		//get material constants
		const double SpacingR 			= mpYieldCriterion->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double SwellingSlope		= mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
  		const double OtherSlope      	= mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];

  		//get hardening variable
  		double P0 = rPlasticVariables.PreconsolidationPressure;

		//calcualte Kirchhoff stress vector
		Vector StressVector = ZeroVector(6);
		this->CalculateKirchhoffStressVector(rHenckyStrainVector, StressVector);

		//calculate d_h/d_eps^p
		Vector PreconDerivativeEpsVol = ZeroVector(6);
		for (unsigned int i = 0; i<3; ++i)
			PreconDerivativeEpsVol(i) = 1.0;
		PreconDerivativeEpsVol *= -P0/(OtherSlope - SwellingSlope);
		
		//calculate d_g/d_sig
		Vector PlasticPotentialD = ZeroVector(6);
		Matrix PlasticPotentialDD = ZeroMatrix(1,1);
		this->CalculatePlasticPotentialDerivatives(StressVector, PlasticPotentialD, PlasticPotentialDD, rPlasticVariables);
		
		//compute H = - d_f/d_h * < d_h/d_eps^p, d_g/d_sig >
		rH = 1.0/(P0*log(SpacingR));
		rH *= MathUtils<double>::Dot( PreconDerivativeEpsVol, PlasticPotentialD );

std::cout<<"casm rh: "<<rH<<std::endl;

	}


	// *******************************************************************************
	// ***** Calculate plastic potential derivative with respect to sigma tensor ***** INHERIT
	// *******************************************************************************
	void BorjaCasmExplicitFlowRule::CalculatePlasticPotentialDerivatives( const Vector& rStressVector, Vector& rFirstDerivative, Matrix & rSecondDerivative, const PlasticVariablesType& rPlasticVariables )
	{
		// set second derivatives to zero as not needed for explicit method
		rSecondDerivative = ZeroMatrix(1,1);
		
		//calculate derivatives with respect to P & J2
		double FirstDerivativeP = 0;
		double FirstDerivativeJ2 = 0;
		this->CalculatePlasticPotentialDerivativesPJ2(rStressVector, FirstDerivativeP, FirstDerivativeJ2, rPlasticVariables);
		
		//stress invariants derivatives
		Vector V1, V2;
		StressInvariantsUtilities::CalculateDerivativeVectors( rStressVector, V1, V2);
		
		// calculate d_g/d_sig = d_g/d_p * d_p/d_sig + d_g/d_J * d_J/d_sig + 0
		rFirstDerivative = FirstDerivativeP * V1 + FirstDerivativeJ2 * V2;
		
//std::cout<<"++++++++++ dP: "<<FirstDerivativeP<<", dJ: "<<FirstDerivativeJ2<<std::endl;
//std::cout<<"++++++++++ V1: "<<V1<<", V2: "<<V2<<std::endl<<std::endl;
	}


	// *************************************************************************
	// ***** Calculate plastic potential derivative with respect to P & J2 ***** MODIFICATION
	// *************************************************************************
	void BorjaCasmExplicitFlowRule::CalculatePlasticPotentialDerivativesPJ2( const Vector& rStressVector, double& rFirstDerivativeP, double& rFirstDerivativeJ2, const PlasticVariablesType& rPlasticVariables)
	{
		//get inclination of CS-line
		const double ShearM = mpYieldCriterion->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		const double Lambda = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];
		const double Kappa = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
		const double SpacingR = mpYieldCriterion->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN = mpYieldCriterion->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];
		
		//calcualte m-factor
		const double FactorM = 2.0/3.0 * (pow(ShearM*(6-ShearM),ShapeN) - pow(3.0*ShearM,ShapeN))/( (1-Kappa/Lambda)*(6.0-ShearM)*pow(3.0*ShearM,ShapeN-1) );
		
		//stress invariants & invariants derivatives
		double MeanStress, J2, LodeAngle;
		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2, LodeAngle);
		
		//calculate third invariant effect 
		double ThirdInvariantEffect = 1.0;
		ThirdInvariantEffect = mpYieldCriterion->EvaluateThirdInvariantEffectMC( LodeAngle);
		
		const double flowID = 5;
		if (flowID == 1)
		{
			//Rowe
			rFirstDerivativeP = 27.0*( ShearM/ThirdInvariantEffect*MeanStress + pow(3.0,0.5)*J2 ) / ( (3.0*MeanStress-2.0*pow(3.0,0.5)*J2) * (3.0*MeanStress+pow(3.0,0.5)*J2) );
			rFirstDerivativeJ2 = 3.0*pow(3.0,0.5)*(-9.0*MeanStress - ShearM/ThirdInvariantEffect*(3.0*MeanStress + 2.0*pow(3.0,0.5)*J2)) / ( (3.0*MeanStress-2.0*pow(3.0,0.5)*J2) * (3.0*MeanStress+pow(3.0,0.5)*J2) );
std::cout<<"casm ROWE: dp = "<<rFirstDerivativeP<<" dJ = "<<rFirstDerivativeJ2<<std::endl;
		}
		else if (flowID == 2)
		{
			//Yu
			rFirstDerivativeP = ( ShapeN*(FactorM-1.0)/MeanStress + ( ShapeN*FactorM*(FactorM-1.0)*pow(pow(3.0,0.5)*J2/(ShearM/ThirdInvariantEffect*(-MeanStress)),ShapeN) )/( -MeanStress*(1.0+(FactorM-1.0)* pow(pow(3.0,0.5)*J2/(ShearM/ThirdInvariantEffect*(-MeanStress)),ShapeN)) ) );
			rFirstDerivativeJ2 = ( ShapeN*FactorM*(FactorM-1.0)*pow(pow(3.0,0.5)*J2/(ShearM/ThirdInvariantEffect*(-MeanStress)),ShapeN) )/( J2*(1.0+(FactorM-1.0)* pow(pow(3.0,0.5)*J2/(ShearM/ThirdInvariantEffect*(-MeanStress)),ShapeN)) );
		}
		else if (flowID == 3)
		{
			//alternative
			const double Beta = 1.00;
			rFirstDerivativeP = 1/( MeanStress * log(SpacingR) ) + ( ShapeN * pow( Beta*pow(3,1/2)*J2 , ShapeN) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-MeanStress,ShapeN+1) );
			rFirstDerivativeJ2 = ( ShapeN * pow(3*Beta,ShapeN/2) * pow(J2,ShapeN-1) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-MeanStress,ShapeN) );
		}
		else if (flowID == 4)
		{
			//alternative
			rFirstDerivativeP = 1/( MeanStress * log(SpacingR) ) + ( ShapeN * pow( pow(3,1/2)*J2 , ShapeN) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-MeanStress,ShapeN+1) );
			rFirstDerivativeJ2 = ( ShapeN * pow(3,ShapeN/2) * pow(J2,ShapeN-1) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-MeanStress,ShapeN) );
		}
		else if (flowID == 5)
		{
			//MCC
			double chiMCC = MeanStress + 3.0/MeanStress*pow( J2*ThirdInvariantEffect/ShearM, 2 );
			rFirstDerivativeP = 2.0*MeanStress - chiMCC;
			rFirstDerivativeJ2 = 2.0 * 3.0 * pow( ThirdInvariantEffect / ShearM, 2) * J2;
std::cout<<"casm MCC: dp = "<<rFirstDerivativeP<<" dJ = "<<rFirstDerivativeJ2<<std::endl;
		}
	}


	// **********************************************************************
	// ***** Calculate yield function and plastic potential derivatives ***** INHERIT
	// **********************************************************************
	void BorjaCasmExplicitFlowRule::UpdateDerivatives(const Vector& rHenckyStrain, AuxiliarDerivativesStructure& rAuxiliarDerivatives, const PlasticVariablesType& rPlasticVariables)
  {
		Vector StressVector = ZeroVector(6);

		this->CalculateKirchhoffStressVector( rHenckyStrain, StressVector);

		this->CalculatePlasticPotentialDerivatives(StressVector, rAuxiliarDerivatives.PlasticPotentialD, rAuxiliarDerivatives.PlasticPotentialDD, rPlasticVariables);
//rAuxiliarDerivatives.PlasticPotentialD = ZeroVector(1);

		mpYieldCriterion->CalculateYieldFunctionDerivative(StressVector, rAuxiliarDerivatives.YieldFunctionD, rPlasticVariables);

		// set plastic potential derivatives equal to yield surface derivatives if zero matrices are obtained
		if ( rAuxiliarDerivatives.PlasticPotentialD.size() == 1 ) {
			rAuxiliarDerivatives.PlasticPotentialD = rAuxiliarDerivatives.YieldFunctionD;
			if ( rAuxiliarDerivatives.PlasticPotentialDD.size1() == 1) {
					//mpYieldCriterion->CalculateYieldFunctionSecondDerivative( StressVector, rAuxiliarDerivatives.PlasticPotentialDD, rAlpha);

				if ( rAuxiliarDerivatives.PlasticPotentialDD.size1() == 1) {
					rAuxiliarDerivatives.PlasticPotentialDD = ZeroMatrix(6,6);
				}
			}
		}
		else
		{
		} 
  }

  // ***********************************************************************
  // ***** Calculate Kirchhoff stress vector from HENCKY STRAIN VECTOR ***** INHERIT
  // ***********************************************************************
  void BorjaCasmExplicitFlowRule::CalculateKirchhoffStressVector(const Vector& rHenckyStrainVector, Vector& rKirchhoffStressVector)
  {
		Vector IdentityVector = ZeroVector(6);

		for (unsigned int i = 0; i < 3; ++i)
			IdentityVector(i) = 1.0;

		//calcualte volumetric Hencky strain
		double MeanStress; 
		double VolumetricStrain = MathUtils<double>::Dot( trans(rHenckyStrainVector), IdentityVector);
		//calculate deviatoric Hencky strain vector
		Vector DeviatoricStrainVector; 
		DeviatoricStrainVector = rHenckyStrainVector -  (VolumetricStrain/3.0)*IdentityVector;
		
		//evaluate mean stress
		EvaluateMeanStress(VolumetricStrain, DeviatoricStrainVector, MeanStress);
		//evaluate deviatoric stress
		EvaluateDeviatoricStress(VolumetricStrain, DeviatoricStrainVector, rKirchhoffStressVector);
		//calculate total Kirchhoff stress vector
		rKirchhoffStressVector += MeanStress*IdentityVector;
  }


  // **********************************************************************************
  // ***** Calculate Kirchhoff stress vector from LEFT CAUCHY GREEN STRAIN TENSOR ***** INHERIT
  // **********************************************************************************
  void BorjaCasmExplicitFlowRule::CalculateKirchhoffStressVector(const Matrix& rElasticLeftCauchyGreen, Vector& rNewStressVector)
  {
    Vector ElasticHenckyStrainVector;
    ElasticHenckyStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rElasticLeftCauchyGreen); 
    this->CalculateKirchhoffStressVector(ElasticHenckyStrainVector, rNewStressVector);
  }


  // *********************************************************************************************
  // ***** Calculate the mean stress based on volumetric strain and deviatoric strain vector ***** INHERIT
  // *********************************************************************************************
  void BorjaCasmExplicitFlowRule::EvaluateMeanStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, double& rMeanStress)
  {
	double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;    

    double DeviatoricStrain2Norm = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
			DeviatoricStrain2Norm += pow(rDeviatoricStrainVector(i), 2);

    for (unsigned int i = 3; i < 6; ++i)
      DeviatoricStrain2Norm += 2.0*pow(rDeviatoricStrainVector(i)/2.0, 2);

    rMeanStress = -ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope) * (1.0 + 1.0*AlphaShear*DeviatoricStrain2Norm/SwellingSlope);
  }

  // ***************************************************************
  // ***** Calculate the mean stress from Hencky strain vector ***** INHERIT
  // ***************************************************************
  void BorjaCasmExplicitFlowRule::EvaluateMeanStress(const Vector& rHenckyStrainVector, double& rMeanStress)
  {
		//calculate volumetric Hencky strain
    double VolumetricStrain = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
			VolumetricStrain += rHenckyStrainVector(i);
		//calculate deviatoric Hencky strain vector
    Vector DeviatoricStrain = rHenckyStrainVector;
    for (unsigned int i = 0; i < 3; ++i)
			DeviatoricStrain(i) -= VolumetricStrain/3.0;
		//calculate the mean stress based on eps_vol and eps_dev-vector
    this->EvaluateMeanStress(VolumetricStrain, DeviatoricStrain, rMeanStress);
  }


  // ******************************************************************
	// ***** Evaluate the deviatoric part of the hyperelastic model ***** INHERIT
	// ******************************************************************
  void BorjaCasmExplicitFlowRule::EvaluateDeviatoricStress(const double& rVolumetricStrain, const Vector & rDeviatoricStrainVector, Vector& rDeviatoricStress)
  {
		//calculate reference pressure p_ref = p_c0/OCR
    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;
    //get kappa, alpha shear, G0
    double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];
    double ConstantShearModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];

		//
    rDeviatoricStress = rDeviatoricStrainVector;
    double ShearModulus = AlphaShear*ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope);
    rDeviatoricStress *= 2.0*( ShearModulus + ConstantShearModulus);
		// BECAUSE OF VOIGT NOTATION: e = [e11 e22 e33 2*e23 2*e13 2*e12]'
    for (unsigned int i = 3; i<6; ++i){
      rDeviatoricStress(i) /= 2.0;  
    }
  }


  // ***********************************************
	// ***** EVALUATE THE TANGENT ELASTIC MATRIX ***** INHERIT
	// ***********************************************
	void BorjaCasmExplicitFlowRule::ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix )
	{
//std::cout<<"BorjaCasmExplicitFlowRule::ComputeElasticMatrix"<<std::endl;
		Matrix FourthOrderIdentity = ZeroMatrix(6,6);
		for (unsigned int i = 0; i<3; ++i)
			 FourthOrderIdentity(i,i) = 1.0;

		for (unsigned int i = 3; i<6; ++i)
			 FourthOrderIdentity(i,i) = 0.50;
		// VOIGT NOTATION AND NOT KELVIN

		Matrix IdentityCross = ZeroMatrix(6,6);
		for (unsigned int i = 0; i<3; ++i) {
			 for (unsigned int j = 0; j<3; ++j) {
					IdentityCross(i,j) = 1.0;
			 }
		}


		Vector StressVector = ZeroVector(6);
		this->CalculateKirchhoffStressVector(rElasticStrainVector, StressVector);


		double MeanStress = 0.0;
		double VolumetricStrain = 0.0;
		for (unsigned int i = 0; i<3; i++) {
			 MeanStress += StressVector(i);
			 VolumetricStrain += rElasticStrainVector(i);
		}
		MeanStress /= 3.0;


		double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
		double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
		ReferencePressure /= OCR;    

		double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
		double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

		double ConstantShearModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];


		rElasticMatrix  = (-1.0/SwellingSlope)*MeanStress*IdentityCross;
		rElasticMatrix += 2.0*AlphaShear*ReferencePressure*std::exp(-VolumetricStrain/SwellingSlope)*(FourthOrderIdentity - (1.0/3.0)*IdentityCross);



		Vector StrainVector = rElasticStrainVector; 
		for (unsigned int i = 0; i < 3; ++i)
			 StrainVector(i) -= VolumetricStrain / 3.0;

		double Modulus = 2.0 * ReferencePressure * exp( - VolumetricStrain/SwellingSlope) * ( AlphaShear / SwellingSlope);

		for (unsigned int i = 3; i < 6 ; ++i)
			 StrainVector(i) /= 2.0;


		// PARTE ASQUEROSA
		for (unsigned int i = 0; i<3; ++i) {
			 for (unsigned int j = 0; j<3; ++j) {
					rElasticMatrix(i,j) -= Modulus * (StrainVector(i) ); //-MeanStress);
					rElasticMatrix(i,j) -= Modulus * (StrainVector(j) ); //-MeanStress);
			 }
		}

		for (unsigned int i = 0; i<3; ++i) {
			 for (unsigned int j = 3; j < 6; ++j) {
					rElasticMatrix(i,j) -= Modulus*(StrainVector(j));///2.0;
			 }
		}

		for (unsigned int i = 3; i<6; ++i) {
			 for (unsigned int j = 0; j<3; ++j) {
					rElasticMatrix(i,j) -= Modulus*(StrainVector(i));///2.0;
			 }
		}

		// AND THE PART DUE TO THE INITIAL SHEAR MODULUS

		rElasticMatrix +=  2.0*ConstantShearModulus * ( FourthOrderIdentity - (1.0/3.0)*IdentityCross );

	}

	// *************************************************************
	// ***** Set initial bonding and preconsolidation pressure *****
	//
	void BorjaCasmExplicitFlowRule::SetPlasticVariables( const double& rInitialPreconPressure, const double& rInitialBonding)
	{
		//mPlasticVariables.Bonding = rInitialBonding;
		mPlasticVariables.PreconsolidationPressure = rInitialPreconPressure;
	}


	// *************************************************
	// ***** Set initial preconsolidation pressure ***** INHERIT
	//
	void BorjaCasmExplicitFlowRule::SetPreconsolidation( const double& rInitialPreconPressure)
	{
		mPlasticVariables.PreconsolidationPressure = rInitialPreconPressure;
	}


	// *********************************************
	// ***** Upadte internal plastic variables ***** INHERIT
	//
	bool BorjaCasmExplicitFlowRule::UpdateInternalVariables( RadialReturnVariables & rReturnMappingVariables)
	{
		mPlasticVariables												= rReturnMappingVariables.PlasticVariables;

		mPlasticMultiplierVelocity 							= rReturnMappingVariables.DeltaBeta / rReturnMappingVariables.DeltaTime;
		
		return 0;
 	}


	// *****************************************************************
	// ***** Update ReturnMappingVariables during time integration ***** INHERIT
 	//
  void BorjaCasmExplicitFlowRule::UpdateRadialReturnVariables( RadialReturnVariables& rReturnMappingVariables, const ExplicitStressUpdateInformation& rStressUpdateInformation)
  {
		rReturnMappingVariables.PlasticVariables 								= rStressUpdateInformation.PlasticVariablesUpdate;
		
		rReturnMappingVariables.DeltaBeta 											 += rStressUpdateInformation.DeltaDeltaPlastic;	
  }


  void BorjaCasmExplicitFlowRule::save( Serializer& rSerializer) const 
  {
		KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonAssociativeExplicitPlasticFlowRule )
  }

  void BorjaCasmExplicitFlowRule::load( Serializer& rSerializer)
  {
		KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonAssociativeExplicitPlasticFlowRule )
  }

} //end namespace kratos
