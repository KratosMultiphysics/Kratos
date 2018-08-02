// Author: Salva Latorre
#include "includes/define.h"
#include "flex_wrapper.h"
#include "input_output/logger.h"
#include "includes/variables.h"
#include "includes/dem_variables.h"

namespace Kratos {

    FlexWrapper::FlexWrapper(ModelPart& rSpheresModelPart):mrSpheresModelPart(rSpheresModelPart) {
        mFlexLibrary = NvFlexInit();
        if (mFlexLibrary == NULL) {
		    KRATOS_ERROR << "Could not initialize Nvidia Flex, exiting." << std::endl;
	    }
        else {
            KRATOS_INFO("Flex: ") << "Nvidia Flex Initialized correctly" << std::endl;
            char device_name[256];
            strcpy(device_name, NvFlexGetDeviceName(mFlexLibrary));
	        KRATOS_INFO("Flex: ") << "Computing Device: "<< device_name << std::endl;
        }
        mFlexPositions = new NvFlexVector<Vec4>(mFlexLibrary);
        mFlexVelocities = new NvFlexVector<Vec3>(mFlexLibrary);
        mFlexPhases = new NvFlexVector<int>(mFlexLibrary);
        mActiveIndices = new NvFlexVector<int>(mFlexLibrary);

        NvFlexSetSolverDescDefaults(&mSolverDescriptor);
        InitializeNvFlexSolverDescParams(mSolverDescriptor);
        mFlexSolver = NvFlexCreateSolver(mFlexLibrary, &mSolverDescriptor);
    }

    FlexWrapper::~FlexWrapper() {
        NvFlexDestroySolver(mFlexSolver);
		mFlexSolver = NULL;
        delete mFlexPositions;
        delete mFlexVelocities;
        delete mFlexPhases;
        delete mActiveIndices;
    }

    void FlexWrapper::UpdateFlex() {

        TransferDataFromKratosToFlex();

        UpdateFlexParameters();
    }

    void FlexWrapper::TransferDataFromKratosToFlex() {


        const size_t number_of_nodes = mrSpheresModelPart.Nodes().size();
        mSolverDescriptor.maxParticles = number_of_nodes;
        mSolverDescriptor.maxDiffuseParticles = number_of_nodes;
        NvFlexDestroySolver(mFlexSolver);
        mFlexSolver = NvFlexCreateSolver(mFlexLibrary, &mSolverDescriptor);
        NvFlexSetParams(mFlexSolver, &mFlexParameters);
        mFlexPositions->resize(number_of_nodes);
        mFlexVelocities->resize(number_of_nodes);
        mFlexPhases->resize(number_of_nodes);
        mActiveIndices->resize(number_of_nodes);

        mFlexPositions->map();
        mFlexVelocities->map();
        mFlexPhases->map();
        mActiveIndices->map();

        Vec4 aux_vec4;
        Vec3 aux_vec3;
        int phase = NvFlexMakePhase(0, eNvFlexPhaseSelfCollide | eNvFlexPhaseFluid);

        for (size_t i=0; i< number_of_nodes; i++) {
            const auto node_it = mrSpheresModelPart.Nodes().begin() + i;

            //positions
            const auto& coords = node_it->Coordinates();
            aux_vec4[0] = coords[0];
            aux_vec4[1] = coords[1];
            aux_vec4[2] = coords[2];
            aux_vec4[3] = 1.0 / node_it->FastGetSolutionStepValue(NODAL_MASS);
            NvFlexVector<Vec4>& array_of_positions = *mFlexPositions;
            array_of_positions[i] = aux_vec4;

            //velocities
            const auto& vel = node_it->FastGetSolutionStepValue(VELOCITY);
            aux_vec3[0] = vel[0];
            aux_vec3[1] = vel[1];
            aux_vec3[2] = vel[2];
            NvFlexVector<Vec3>& array_of_velocities = *mFlexVelocities;
            array_of_velocities[i] = aux_vec3;

            //phases
            NvFlexVector<int>& array_of_phases = *mFlexPhases;
            array_of_phases[i] = phase;

            //indices
            NvFlexVector<int>& array_of_indices = *mActiveIndices;
            array_of_indices[i] = node_it->Id();
        }

        mFlexPositions->unmap();
        mFlexVelocities->unmap();
        mFlexPhases->unmap();
        mActiveIndices->unmap();

        // send any particle updates to the solver
        NvFlexSetParticles(mFlexSolver, mFlexPositions->buffer, NULL);
        NvFlexSetVelocities(mFlexSolver, mFlexVelocities->buffer, NULL);
        NvFlexSetPhases(mFlexSolver, mFlexPhases->buffer, NULL);
        NvFlexSetActive(mFlexSolver, mActiveIndices->buffer, NULL);
        NvFlexSetActiveCount(mFlexSolver, mActiveIndices->size());
    }

    void FlexWrapper::SolveTimeSteps(double dt, int number_of_substeps) {
        NvFlexUpdateSolver(mFlexSolver, (float)dt, number_of_substeps, false);
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

        NvFlexGetParticles(mFlexSolver, mFlexPositions->buffer, &mFlexCopyDescriptor);
        NvFlexGetVelocities(mFlexSolver, mFlexVelocities->buffer, &mFlexCopyDescriptor);

        mFlexPositions->map();
        mFlexVelocities->map();
        mFlexPhases->map();
        mActiveIndices->map();

        const size_t number_of_nodes = mrSpheresModelPart.Nodes().size();

        for (size_t i=0; i< number_of_nodes; i++) {
            const auto node_it = mrSpheresModelPart.Nodes().begin() + i;

            //positions
            auto& coords = node_it->Coordinates();
            NvFlexVector<Vec4>& array_of_positions = *mFlexPositions;
            coords[0] = array_of_positions[i][0];
            coords[1] = array_of_positions[i][1];
            coords[2] = array_of_positions[i][2];

            //velocities
            auto& vel = node_it->FastGetSolutionStepValue(VELOCITY);
            NvFlexVector<Vec3>& array_of_velocities = *mFlexVelocities;
            vel[0] = array_of_velocities[i][0];
            vel[1] = array_of_velocities[i][1];
            vel[2] = array_of_velocities[i][2];
        }



    }

    void FlexWrapper::InitializeNvFlexSolverDescParams(NvFlexSolverDesc& mSolverDescriptor) {

        mSolverDescriptor.featureMode = eNvFlexFeatureModeSimpleSolids;
        mSolverDescriptor.maxParticles = mrSpheresModelPart.NumberOfElements();
	    mSolverDescriptor.maxDiffuseParticles = mrSpheresModelPart.NumberOfElements();;
	    mSolverDescriptor.maxNeighborsPerParticle = 32;
	    mSolverDescriptor.maxContactsPerParticle = 10;
    }

    void FlexWrapper::InitializeNvFlexParams(NvFlexParams& mFlexParameters) {

        const array_1d<double,3>& gravity = mrSpheresModelPart.GetProcessInfo()[GRAVITY];
        mFlexParameters.gravity[0] = (float)gravity[0];
        mFlexParameters.gravity[1] = (float)gravity[1];
        mFlexParameters.gravity[2] = (float)gravity[2];

        mFlexParameters.wind[0] = 0.0f;
        mFlexParameters.wind[1] = 0.0f;
        mFlexParameters.wind[2] = 0.0f;

        mFlexParameters.radius = 0.15f;
        const size_t number_of_nodes = mrSpheresModelPart.Nodes().size();
        for (size_t i=0; i< number_of_nodes; i++) {
            const auto node_it = mrSpheresModelPart.Nodes().begin() + i;
            mFlexParameters.radius = node_it->FastGetSolutionStepValue(RADIUS);
            break;
        }
        //TODO: MA: check that all radii are the same!!

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

        const size_t number_of_nodes = mrSpheresModelPart.Nodes().size();
	    mFlexCopyDescriptor.elementCount = number_of_nodes;
    }

    std::string FlexWrapper::Info() const {

        std::stringstream buffer;
        buffer << "FlexWrapper" ;
        return buffer.str();
    }

    void FlexWrapper::PrintInfo(std::ostream& rOStream) const {rOStream << "FlexWrapper";}

    void FlexWrapper::PrintData(std::ostream& rOStream) const {}
} // namespace Kratos
