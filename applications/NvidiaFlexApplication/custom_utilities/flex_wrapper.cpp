// Author: Salva Latorre
#include "includes/define.h"
#include "flex_wrapper.h"
#include "input_output/logger.h"
#include "includes/variables.h"
#include "includes/dem_variables.h"

#define NV_FLEX_VERSION 120
bool NvidiaFlexError = false;

namespace Kratos {

    static void NvidiaFlexErrorCallback(NvFlexErrorSeverity severity, const char* msg, const char* file, int line)
    {
        printf("Flex: %s - %s:%d\n", msg, file, line);
        NvidiaFlexError = (severity == eNvFlexLogError); //eNvFlexLogWarning eNvFlexLogInfo eNvFlexLogDebug eNvFlexLogAll
    }

    FlexWrapper::FlexWrapper(ModelPart& rSpheresModelPart, ParticleCreatorDestructor& rParticleCreatorDestructor):mFlexSolver(NULL), mrSpheresModelPart(rSpheresModelPart), mrParticleCreatorDestructor(rParticleCreatorDestructor) {

        mMaxparticles = 10000;
        // a setting of -1 means Flex will use the device specified in the NVIDIA control panel
	    mInitDesc.deviceIndex = -1;
	    mInitDesc.enableExtensions = true;
	    mInitDesc.renderDevice = 0;
	    mInitDesc.renderContext = 0;
	    mInitDesc.computeContext = 0;
	    mInitDesc.computeType = eNvFlexCUDA;
        mInitDesc.runOnRenderContext = false;

        mFlexLibrary = NvFlexInit(NV_FLEX_VERSION, NvidiaFlexErrorCallback, &mInitDesc);
        if (NvidiaFlexError || mFlexLibrary == NULL) {
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
        //mFlexRestPositions = new NvFlexVector<Vec4>(mFlexLibrary);

        NvFlexSetSolverDescDefaults(&mSolverDescriptor);

        //mSolverDescriptor.featureMode = eNvFlexFeatureModeSimpleSolids;
        mSolverDescriptor.maxNeighborsPerParticle = 32;
	    mSolverDescriptor.maxContactsPerParticle = 10;
        mSolverDescriptor.maxParticles = mMaxparticles;
        mSolverDescriptor.maxDiffuseParticles = 0;
        mFlexSolver = NvFlexCreateSolver(mFlexLibrary, &mSolverDescriptor);
    }

    FlexWrapper::~FlexWrapper() {
        NvFlexDestroySolver(mFlexSolver);
		mFlexSolver = NULL;
        delete mFlexPositions;
        delete mFlexVelocities;
        delete mFlexPhases;
        delete mActiveIndices;
        //delete mFlexRestPositions;
    }

    void FlexWrapper::SetNvFlexCopyDescParams(NvFlexCopyDesc& mFlexCopyDescriptor) {
        mFlexCopyDescriptor.dstOffset = 0;
	    mFlexCopyDescriptor.srcOffset = 0;

        mFlexCopyDescriptor.elementCount = mrSpheresModelPart.Nodes().size();
	    // mFlexCopyDescriptor.elementCount = mMaxparticles;
    }

    void FlexWrapper::UpdateFlex() {
        SetNvFlexCopyDescParams(mFlexCopyDescriptor);
        SetNvFlexParams(mFlexParameters);
        TransferDataFromKratosToFlex();
    }

    void FlexWrapper::TransferDataFromKratosToFlex() {

        NvFlexSetParams(mFlexSolver, &mFlexParameters);

        const size_t number_of_nodes = mrSpheresModelPart.Nodes().size();
        mFlexPositions->resize(number_of_nodes);
        mFlexVelocities->resize(number_of_nodes);
        mFlexPhases->resize(number_of_nodes);
        mActiveIndices->resize(number_of_nodes);
        //mFlexRestPositions->resize(0);

        //Vec4 aux_vec4;
        // Vec3 aux_vec3;
        int phase = NvFlexMakePhase(0, eNvFlexPhaseSelfCollide);

        for (size_t i=0; i< number_of_nodes; i++) {
            const auto node_it = mrSpheresModelPart.Nodes().begin() + i;

            //positions
            const auto& coords = node_it->Coordinates();
            (*mFlexPositions)[i] = Vec4((float)coords[0], (float)coords[1], (float)coords[2], (float)(1.0 / node_it->FastGetSolutionStepValue(NODAL_MASS)));

            //velocities
            const auto& vel = node_it->FastGetSolutionStepValue(VELOCITY);
            (*mFlexVelocities)[i] = Vec3((float)vel[0], (float)vel[1], (float)vel[2]);

            //phases
            (*mFlexPhases)[i] = phase;

            //indices
            (*mActiveIndices)[i] = i; //node_it->Id();

            //(*mFlexRestPositions)[i] = (*mFlexPositions)[i];
        }

        mFlexPositions->unmap();
        mFlexVelocities->unmap();
        mFlexPhases->unmap();
        mActiveIndices->unmap();
        //mFlexRestPositions->unmap();

        // send any particle updates to the solver
        NvFlexSetParticles(mFlexSolver, mFlexPositions->buffer, &mFlexCopyDescriptor);
        NvFlexSetVelocities(mFlexSolver, mFlexVelocities->buffer, &mFlexCopyDescriptor);
        NvFlexSetPhases(mFlexSolver, mFlexPhases->buffer, &mFlexCopyDescriptor);
        NvFlexSetActive(mFlexSolver, mActiveIndices->buffer, &mFlexCopyDescriptor);
        NvFlexSetActiveCount(mFlexSolver, mActiveIndices->size());
        //NvFlexSetRestParticles(mFlexSolver, mFlexRestPositions->buffer, &mFlexCopyDescriptor);
    }

    void FlexWrapper::SetNvFlexParams(NvFlexParams& mFlexParameters) {

        const array_1d<double,3>& gravity = mrSpheresModelPart.GetProcessInfo()[GRAVITY];
        mFlexParameters.gravity[0] = (float)gravity[0];
        mFlexParameters.gravity[1] = (float)gravity[1];
        mFlexParameters.gravity[2] = (float)gravity[2];

        mFlexParameters.wind[0] = 0.0f;
        mFlexParameters.wind[1] = 0.0f;
        mFlexParameters.wind[2] = 0.0f;

        mFlexParameters.radius = 0.00015f;
        /*const size_t number_of_nodes = mrSpheresModelPart.Nodes().size();
        for (size_t i=0; i< number_of_nodes; i++) {
            const auto node_it = mrSpheresModelPart.Nodes().begin() + i;
            mFlexParameters.radius = node_it->FastGetSolutionStepValue(RADIUS);
            break;
        }*/
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
        //mFlexParameters.solidRestDistance = 0.0f;
        mFlexParameters.solidRestDistance = mFlexParameters.radius;

        mFlexParameters.anisotropyScale = 1.0f;
        mFlexParameters.anisotropyMin = 0.1f;
        mFlexParameters.anisotropyMax = 2.0f;
        mFlexParameters.smoothing = 1.0f;

        mFlexParameters.dissipation = 0.0f;
        mFlexParameters.damping = 0.0f;
        mFlexParameters.particleCollisionMargin = 0.0f;
        mFlexParameters.shapeCollisionMargin = 0.0f;
        //mFlexParameters.collisionDistance = 0.0f;
        mFlexParameters.collisionDistance = mFlexParameters.radius;
        mFlexParameters.sleepThreshold = 0.0f;
        mFlexParameters.shockPropagation = 0.0f;
        mFlexParameters.restitution = 0.5f;

        mFlexParameters.maxSpeed = 100.0f; //FLT_MAX;
        mFlexParameters.maxAcceleration = 100.0f;	// approximately 10x gravity

        mFlexParameters.relaxationMode = eNvFlexRelaxationLocal;
        mFlexParameters.relaxationFactor = 1.0f;
        mFlexParameters.solidPressure = 1.0f;
        mFlexParameters.adhesion = 0.0f;
        mFlexParameters.cohesion = 0.025f;
        mFlexParameters.surfaceTension = 0.0f;
        mFlexParameters.vorticityConfinement = 0.0f;
        mFlexParameters.buoyancy = 1.0f;
        mFlexParameters.diffuseThreshold = 1000.0f;
        mFlexParameters.diffuseBuoyancy = 0.0f;
        mFlexParameters.diffuseDrag = 0.0f;
        mFlexParameters.diffuseBallistic = 16;
        mFlexParameters.diffuseLifetime = 0.0f;

        // planes created after particles
        //MAC: ALL of this can be done a lot simpler if high node and low node are in some Kratos processinfo.
        const array_1d<double, 3 >& low_point = mrParticleCreatorDestructor.GetLowNode();
        const array_1d<double, 3 >& high_point = mrParticleCreatorDestructor.GetHighNode();

        mFlexParameters.numPlanes = 6;
        (Vec4&)mFlexParameters.planes[0] = Vec4(0.0f, 1.0f, 0.0f, (float)-low_point[1]);
        (Vec4&)mFlexParameters.planes[1] = Vec4(0.0f, 0.0f, 1.0f, (float)-low_point[2]);
        (Vec4&)mFlexParameters.planes[2] = Vec4(1.0f, 0.0f, 0.0f, (float)-low_point[0]);
        (Vec4&)mFlexParameters.planes[3] = Vec4(-1.0f, 0.0f, 0.0f, (float)high_point[0]);
        (Vec4&)mFlexParameters.planes[4] = Vec4(0.0f, 0.0f, -1.0f, (float)high_point[2]);
        (Vec4&)mFlexParameters.planes[5] = Vec4(0.0f, -1.0f, 0.0f, (float)high_point[1]);
        (Vec4&)mFlexParameters.planes[6] = Vec4(0.0f, -1.0f, 0.0f, (float)high_point[1]);
        (Vec4&)mFlexParameters.planes[7] = Vec4(0.0f, -1.0f, 0.0f, (float)high_point[1]);

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

    void FlexWrapper::SolveTimeSteps(double dt, int number_of_substeps) {
        NvFlexUpdateSolver(mFlexSolver, (float)dt, number_of_substeps, false);
    }

    void FlexWrapper::TransferDataFromFlexToKratos() {

        NvFlexGetParticles(mFlexSolver, mFlexPositions->buffer, &mFlexCopyDescriptor);
        NvFlexGetVelocities(mFlexSolver, mFlexVelocities->buffer, &mFlexCopyDescriptor);

        mFlexPositions->map();
        mFlexVelocities->map();
        //mFlexPhases->map();
        //mActiveIndices->map();
        //mFlexRestPositions->map();

        const size_t number_of_nodes = mrSpheresModelPart.Nodes().size();
        for (size_t i=0; i< number_of_nodes; i++) {
            const auto node_it = mrSpheresModelPart.Nodes().begin() + i;

            //positions
            auto& coords = node_it->Coordinates();
            coords[0] = (*mFlexPositions)[i][0];
            coords[1] = (*mFlexPositions)[i][1];
            coords[2] = (*mFlexPositions)[i][2];
            auto& disp = node_it->FastGetSolutionStepValue(DISPLACEMENT);
            disp[0] = coords[0] - node_it->X0();
            disp[1] = coords[1] - node_it->Y0();
            disp[2] = coords[2] - node_it->Z0();

            //velocities
            auto& vel = node_it->FastGetSolutionStepValue(VELOCITY);
            vel[0] = (*mFlexVelocities)[i][0];
            vel[1] = (*mFlexVelocities)[i][1];
            vel[2] = (*mFlexVelocities)[i][2];
        }

        mFlexPositions->unmap();
        mFlexVelocities->unmap();
        //mFlexPhases->unmap();
        //mActiveIndices->unmap();
    }

    void FlexWrapper::Finalize() {

        NvFlexFreeBuffer(mFlexPositions->buffer);
        NvFlexFreeBuffer(mFlexVelocities->buffer);
        NvFlexFreeBuffer(mFlexPhases->buffer);
        NvFlexFreeBuffer(mActiveIndices->buffer);
        NvFlexDestroySolver(mFlexSolver);
        NvFlexShutdown(mFlexLibrary);
        //NvFlexFreeBuffer(mFlexRestPositions->buffer);
    }

    std::string FlexWrapper::Info() const {

        std::stringstream buffer;
        buffer << "FlexWrapper" ;
        return buffer.str();
    }

    void FlexWrapper::PrintInfo(std::ostream& rOStream) const {rOStream << "FlexWrapper";}

    void FlexWrapper::PrintData(std::ostream& rOStream) const {}
} // namespace Kratos
