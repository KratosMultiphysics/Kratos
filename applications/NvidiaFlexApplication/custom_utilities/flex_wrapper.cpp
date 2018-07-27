// Author: Salva Latorre

#include "flex_wrapper.h"

namespace Kratos {

    FlexWrapper::FlexWrapper() {}

    FlexWrapper::~FlexWrapper() {}
    
    void FlexWrapper::Initialize() {
           
        mFlexLibrary = NvFlexInit();
        
        InitializeNvFlexSolverDescParams(mSolverDescriptor);
        
        mFlexSolver = NvFlexCreateSolver(mFlexLibrary, &mSolverDescriptor);
    }
    
    void FlexWrapper::UpdateFlex() {
    
        TransferDataFromKratosToFlex();
        
        UpdateFlexParameters();
    }
    
    void FlexWrapper::TransferDataFromKratosToFlex() {
    
        NvFlexVector<Vec4> mFlexPositions(mFlexLibrary);
        //positions.map();        

	for (unsigned int i = 0; i < mNumberOfParticles; ++i) {
            Vec4 positions_data = { 0.5, 0.5, 0.5, 1};
            mFlexPositions.push_back(positions_data);
	}
        
        NvFlexVector<Vec3> mFlexVelocities(mFlexLibrary);
        //velocities.map();        

	for (unsigned int i = 0; i < mNumberOfParticles; ++i) {
            Vec3 velocities_data = { 0.5, 0.5, 0.5};
            mFlexVelocities.push_back(velocities_data);
	}
        
        NvFlexVector<int> mFlexPhases(mFlexLibrary);
        //phases.map();        

	for (unsigned int i = 0; i < mNumberOfParticles; ++i) {
            mFlexPhases.push_back(mPhaseType);
	}
    }
    
    void FlexWrapper::UpdateFlexParameters() {
    
        InitializeNvFlexParams(mFlexParameters);
        
        InitializeNvFlexCopyDescParams(mFlexCopyDescriptor);
    }
    
    void FlexWrapper::Solve() {
            
        // map buffers for reading / writing

        // unmap buffers
        NvFlexUnmap(mFlexPositions->buffer);
        NvFlexUnmap(mFlexVelocities->buffer);
        NvFlexUnmap(mFlexPhases->buffer);

        // write to device (async)
        NvFlexSetParams(mFlexSolver, &mFlexParameters);

        NvFlexSetParticles(mFlexSolver, mFlexPositions->buffer, &mFlexCopyDescriptor);
        NvFlexSetVelocities(mFlexSolver, mFlexVelocities->buffer, &mFlexCopyDescriptor);
        //NvFlexSetNormals(solver, mNormals, &copyDesc);
        NvFlexSetPhases(mFlexSolver, mFlexPhases->buffer, &mFlexCopyDescriptor);
        //NvFlexSetRestParticles(solver, mRestPositions, &copyDesc);
        //NvFlexSetActive(solver, mActiveIndices, &copyDesc);
        NvFlexSetActiveCount(mFlexSolver, mNumberOfParticles);
        //NvFlexSetRigids(g_solver, g_buffers->rigidOffsets->buffer, g_buffers->rigidIndices->buffer, g_buffers->rigidLocalPositions->buffer, g_buffers->rigidLocalNormals->buffer, g_buffers->rigidCoefficients->buffer,
        //g_buffers->rigidPlasticThresholds->buffer, g_buffers->rigidPlasticCreeps->buffer, g_buffers->rigidRotations->buffer, g_buffers->rigidTranslations->buffer, g_buffers->rigidOffsets.size() - 1, g_buffers->rigidIndices.size());

        NvFlexUpdateSolver(mFlexSolver, mDeltaTime, 1, false);
        
        TransferDataFromFlexToKratos();
    }
    
    void FlexWrapper::RunSimulation() {
                               
        Initialize();
        
        UpdateFlex();
        
        Solve();
        
        Finalize();        
    }
    
    void FlexWrapper::Finalize() {
    
        NvFlexFreeBuffer(mFlexPositions->buffer);
        NvFlexFreeBuffer(mFlexVelocities->buffer);
        NvFlexFreeBuffer(mFlexPhases->buffer);
        NvFlexDestroySolver(mFlexSolver);
        NvFlexShutdown(mFlexLibrary);
    }
    
    void FlexWrapper::TransferDataFromFlexToKratos() {
    
        // read back (async)
        NvFlexGetParticles(mFlexSolver, mFlexPositions->buffer, &mFlexCopyDescriptor);
        NvFlexGetVelocities(mFlexSolver, mFlexVelocities->buffer, &mFlexCopyDescriptor);
        NvFlexGetPhases(mFlexSolver, mFlexPhases->buffer, &mFlexCopyDescriptor);
    }
    
    void FlexWrapper::InitializeNvFlexSolverDescParams(NvFlexSolverDesc& mSolverDescriptor) {
                
        mSolverDescriptor.featureMode = eNvFlexFeatureModeSimpleSolids;
        mSolverDescriptor.maxParticles = 10;
	mSolverDescriptor.maxDiffuseParticles = 10;
	mSolverDescriptor.maxNeighborsPerParticle = 32;
	mSolverDescriptor.maxContactsPerParticle = 6;
    }
    
    void FlexWrapper::InitializeNvFlexParams(NvFlexParams& mFlexParameters) {
                        
        mFlexParameters.gravity[0] = 0.0f;
	mFlexParameters.gravity[1] = -9.8f;
	mFlexParameters.gravity[2] = 0.0f;

	mFlexParameters.wind[0] = 0.0f;
	mFlexParameters.wind[1] = 0.0f;
	mFlexParameters.wind[2] = 0.0f;

	mFlexParameters.radius = 0.15f;
	mFlexParameters.viscosity = 0.0f;
	mFlexParameters.dynamicFriction = 0.0f;
	mFlexParameters.staticFriction = 0.0f;
	mFlexParameters.particleFriction = 0.0f;
	mFlexParameters.freeSurfaceDrag = 0.0f;
	mFlexParameters.drag = 0.0f;
	mFlexParameters.lift = 0.0f;
	mFlexParameters.numIterations = 3;
	mFlexParameters.fluidRestDistance = 0.0f;
	mFlexParameters.solidRestDistance = 0.0f;

	mFlexParameters.anisotropyScale = 1.0f;
	mFlexParameters.anisotropyMin = 0.1f;
	mFlexParameters.anisotropyMax = 2.0f;
	mFlexParameters.smoothing = 1.0f;

	mFlexParameters.dissipation = 0.0f;
	mFlexParameters.damping = 0.0f;
	mFlexParameters.particleCollisionMargin = 0.0f;
	mFlexParameters.shapeCollisionMargin = 0.0f;
	mFlexParameters.collisionDistance = 0.0f;
	mFlexParameters.sleepThreshold = 0.0f;
	mFlexParameters.shockPropagation = 0.0f;
	mFlexParameters.restitution = 0.0f;

	mFlexParameters.maxSpeed = 1000.0f; //FLT_MAX;
	mFlexParameters.maxAcceleration = 100.0f;	// approximately 10x gravity

	mFlexParameters.relaxationMode = eNvFlexRelaxationLocal;
	mFlexParameters.relaxationFactor = 1.0f;
	mFlexParameters.solidPressure = 1.0f;
	mFlexParameters.adhesion = 0.0f;
	mFlexParameters.cohesion = 0.025f;
	mFlexParameters.surfaceTension = 0.0f;
	mFlexParameters.vorticityConfinement = 0.0f;
	mFlexParameters.buoyancy = 1.0f;
	mFlexParameters.diffuseThreshold = 100.0f;
	mFlexParameters.diffuseBuoyancy = 1.0f;
	mFlexParameters.diffuseDrag = 0.8f;
	mFlexParameters.diffuseBallistic = 16;
	mFlexParameters.diffuseLifetime = 2.0f;

	// planes created after particles
	mFlexParameters.numPlanes = 1;
	
	if (mFlexParameters.solidRestDistance == 0.0f) mFlexParameters.solidRestDistance = mFlexParameters.radius;

	// if fluid present then we assume solid particles have the same radius
	if (mFlexParameters.fluidRestDistance > 0.0f) mFlexParameters.solidRestDistance = mFlexParameters.fluidRestDistance;

	// set collision distance automatically based on rest distance if not already set
	if (mFlexParameters.collisionDistance == 0.0f) mFlexParameters.collisionDistance = std::max(mFlexParameters.solidRestDistance, mFlexParameters.fluidRestDistance)*0.5f; //Max(g_params.solidRestDistance, g_params.fluidRestDistance)*0.5f;

	// default particle friction to 10% of shape friction
	if (mFlexParameters.particleFriction == 0.0f) mFlexParameters.particleFriction = mFlexParameters.dynamicFriction*0.1f;

	// add a margin for detecting contacts between particles and shapes
	if (mFlexParameters.shapeCollisionMargin == 0.0f) mFlexParameters.shapeCollisionMargin = mFlexParameters.collisionDistance*0.5f;        
    }
    
    void FlexWrapper::InitializeNvFlexCopyDescParams(NvFlexCopyDesc& mFlexCopyDescriptor) {
        
        mFlexCopyDescriptor.dstOffset = 0;
	mFlexCopyDescriptor.srcOffset = 0;
	mFlexCopyDescriptor.elementCount = mNumberOfParticles;
    }

    std::string FlexWrapper::Info() const {
        
        std::stringstream buffer;
        buffer << "FlexWrapper" ;
        return buffer.str();
    }

    void FlexWrapper::PrintInfo(std::ostream& rOStream) const {rOStream << "FlexWrapper";}

    void FlexWrapper::PrintData(std::ostream& rOStream) const {}
} // namespace Kratos
