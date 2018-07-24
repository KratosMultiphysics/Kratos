// This code contains NVIDIA Confidential Information and is disclosed to you
// under a form of NVIDIA software license agreement provided separately to you.
//
// Notice
// NVIDIA Corporation and its licensors retain all intellectual property and
// proprietary rights in and to this software and related documentation and
// any modifications thereto. Any use, reproduction, disclosure, or
// distribution of this software and related documentation without an express
// license agreement from NVIDIA Corporation is strictly prohibited.
//
// ALL NVIDIA DESIGN SPECIFICATIONS, CODE ARE PROVIDED "AS IS.". NVIDIA MAKES
// NO WARRANTIES, EXPRESSED, IMPLIED, STATUTORY, OR OTHERWISE WITH RESPECT TO
// THE MATERIALS, AND EXPRESSLY DISCLAIMS ALL IMPLIED WARRANTIES OF NONINFRINGEMENT,
// MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Information and code furnished is believed to be accurate and reliable.
// However, NVIDIA Corporation assumes no responsibility for the consequences of use of such
// information or for any infringement of patents or other rights of third parties that may
// result from its use. No license is granted by implication or otherwise under any patent
// or patent rights of NVIDIA Corporation. Details are subject to change without notice.
// This code supersedes and replaces all information previously supplied.
// NVIDIA Corporation products are not authorized for use as critical
// components in life support devices or systems without express written approval of
// NVIDIA Corporation.
//
// Copyright (c) 2013-2017 NVIDIA Corporation. All rights reserved.

#ifndef NV_FLEX_EXT_H
#define NV_FLEX_EXT_H

/** \file NvFlexExt.h 
 * The main include file for the Flex extensions API, this is a collection of helper functions
 * for asset creation, scene management, and sample code that builds on the Flex core API.
 */

#include "NvFlex.h"

#include <cassert>
#include <cstddef>

// A vector type that wraps a NvFlexBuffer, behaves like a standard vector for POD types (no construction)
// The vector must be mapped using map() before any read/write access to elements or resize operation

template <typename T>
struct NvFlexVector
{
	NvFlexVector(NvFlexLibrary* l, int size = 0, NvFlexBufferType type = eNvFlexBufferHost) : lib(l), buffer(NULL), mappedPtr(NULL), count(0), capacity(0), type(type)
	{
		if (size)
		{
			resize(size);

			// resize implicitly maps, unmap initial allocation
			unmap();
		}		
	}
	
	NvFlexVector(NvFlexLibrary* l, const T* ptr, int size, NvFlexBufferType type = eNvFlexBufferHost) : lib(l), buffer(NULL), mappedPtr(NULL), count(0), capacity(0), type(type)
	{
		assign(ptr, size);
		unmap();
	}
	

	~NvFlexVector()	
	{
		destroy();
	}

	NvFlexLibrary* lib;
	NvFlexBuffer* buffer;

	T* mappedPtr;
	int count;
	int capacity;
	NvFlexBufferType type;

	// reinitialize the vector leaving it unmapped
	void init(int size)
	{
		destroy();
		resize(size);
		unmap();
	}

	void destroy()
	{	
		if (mappedPtr)
			NvFlexUnmap(buffer);

		if (buffer)
			NvFlexFreeBuffer(buffer);

		mappedPtr = NULL;
		buffer = NULL;
		capacity = 0;		
		count = 0;
	}

	void map(int flags=eNvFlexMapWait)
	{
		if (!buffer)
			return;

		assert(!mappedPtr);
		mappedPtr = (T*)NvFlexMap(buffer, flags);
	}
	
	void unmap()
	{
		if (!buffer)
			return;

		assert(mappedPtr);

		NvFlexUnmap(buffer);
		mappedPtr = 0;
	}

	const T& operator[](int index) const
	{	
		assert(mappedPtr);
		assert(index < count);

		return mappedPtr[index];
	}

	T& operator[](int index)
	{
		assert(mappedPtr);
		assert(index < count);

		return mappedPtr[index];
	}

	void push_back(const T& t)
	{
		assert(mappedPtr || !buffer);

		reserve(count+1);

		// copy element
		mappedPtr[count++] = t;		
	}

	void assign(const T* srcPtr, int newCount)
	{
		assert(mappedPtr || !buffer);

		resize(newCount);

		memcpy(mappedPtr, srcPtr, newCount*sizeof(T));
	}

	void copyto(T* dest, int count)	
	{
		assert(mappedPtr);

		memcpy(dest, mappedPtr, sizeof(T)*count);
	}

	int size() const { return count; }

	bool empty() const { return size() == 0; }

	const T& back() const
	{
		assert(mappedPtr);
		assert(!empty());

		return mappedPtr[count-1];
	}

	void reserve(int minCapacity)
	{
		if (minCapacity > capacity)
		{
			// growth factor of 1.5
			const int newCapacity = minCapacity*3/2;

			NvFlexBuffer* newBuf = NvFlexAllocBuffer(lib, newCapacity, sizeof(T), type);

			// copy contents to new buffer			
			void* newPtr = NvFlexMap(newBuf, eNvFlexMapWait);
			memcpy(newPtr, mappedPtr, count*sizeof(T));

			// unmap old buffer, but leave new buffer mapped
			unmap();
			
			if (buffer)
				NvFlexFreeBuffer(buffer);

			// swap
			buffer = newBuf;
			mappedPtr = (T*)newPtr;
			capacity = newCapacity;			
		}
	}

	// resizes mapped buffer and leaves new buffer mapped 
	void resize(int newCount)
	{
		assert(mappedPtr || !buffer);

		reserve(newCount);	

		// resize but do not initialize new entries
		count = newCount;
	}

	void resize(int newCount, const T& val)
	{
		assert(mappedPtr || !buffer);

		const int startInit = count;
		const int endInit = newCount;

		resize(newCount);

		// init any new entries
		for (int i=startInit; i < endInit; ++i)
			mappedPtr[i] = val;
	}
};

extern "C" {

/**
 * Helper struct for storing the state of a moving frame, see NvFlexExtMovingFrameInit()
 */
struct NvFlexExtMovingFrame
{
	float position[3];
	float rotation[4];

	float velocity[3];
	float omega[3];

	float acceleration[3];
	float tau[3];

	float delta[4][4];
};

/**
 * Creates a new moving frame struct. This helper method is used to calculate inertial forces for particles
 * inside an attached parent frame. For example, when simulating cloth attached to the character, we would like to perform
 * a local space simulation of the cloth to avoid excessive stretching and collision issues during fast animations.
 * However, we would still like the cloth to respond to character movements in at least a limited, or controlled fashion.
 * The NvFlexExtMovingFrame provides a way to include or remove these inertial forces. The basic usage is as follows:
 *  
 \code{.c}

	NvFlexExtMovingFrame frame;
	NvFlexExtMovingFrameInit(&frame, initialTranslation, initialRotation);

	const linearInertiaScale = 0.25f;
	const angularInertiaScale 0.5;

	while(simulating)
	{
		float3 newPosition;
		float4 newRotation;
		
		// move parent frame (character / emitter) according to application's animation system
		Animate(newPosition, newRotation);
		
		// update the frame
		NvFlexExtMovingFrameUpdate(frame, newPosition, newRotation, dt);

		// apply inertial forces and update particles
		NvFlexExtMovingFrameApply(frame, particlePositions, particleVelocities, numParticles, linearInertiaScale, angularInertiaScale, dt);
	}

 \endcode
 
 * @param[in] frame A pointer to a user-allocated NvFlexExtMovingFrame struct
 * @param[in] worldTranslation A pointer to a vec3 storing the frame's initial translation in world space
 * @param[in] worldRotation A pointer to a quaternion storing the frame's initial rotation in world space
 */
NV_FLEX_API void NvFlexExtMovingFrameInit(NvFlexExtMovingFrame* frame, const float* worldTranslation, const float* worldRotation);

/* Update a frame to a new position, this will automatically update the velocity and acceleration of
 * the frame, which can then be used to calculate inertial forces. This should be called once per-frame
 * with the new position and time-step used when moving the frame.
 *
 * @param[in] frame A pointer to a user-allocated NvFlexExtMovingFrame struct
 * @param[in] worldTranslation A pointer to a vec3 storing the frame's initial translation in world space
 * @param[in] worldRotation A pointer to a quaternion storing the frame's initial rotation in world space
 * @param[in] dt The time that elapsed since the last call to the frame update
 */
NV_FLEX_API void NvFlexExtMovingFrameUpdate(NvFlexExtMovingFrame* frame, const float* worldTranslation, const float* worldRotation, float dt);

/* Teleport particles to the frame's new position and apply the inertial forces
 *
 * @param[in] frame A pointer to a user-allocated NvFlexExtMovingFrame struct
 * @param[in] positions A pointer to an array of particle positions in (x, y, z, 1/m) format
 * @param[in] velocities A pointer to an array of particle velocities in (vx, vy, vz) format
 * @param[in] numParticles The number of particles to update
 * @param[in] linearScale How strongly the translational inertial forces should be applied, 0.0 corresponds to a purely local space simulation removing all inertial forces, 1.0 corresponds to no inertial damping and has no benefit over regular world space simulation
 * @param[in] angularScale How strongly the angular inertial forces should be applied, 0.0 corresponds to a purely local space simulation, 1.0 corresponds to no inertial damping
 * @param[in] dt The time that elapsed since the last call to the frame update, should match the value passed to NvFlexExtMovingFrameUpdate()
 */
NV_FLEX_API void NvFlexExtMovingFrameApply(NvFlexExtMovingFrame* frame, float* positions, float* velocities, int numParticles, float linearScale, float angularScale, float dt);


/** 
 * Represents a group of particles and constraints, each asset 
 * can be instanced into a container using NvFlexExtCreateInstance()
 */
struct NvFlexExtAsset
{	
	// particles
	float* particles;				//!< Local space particle positions, x,y,z,1/mass
	int numParticles;				//!< Number of particles
	int maxParticles;				//!< Maximum number of particles, allows extra space for tearable assets which duplicate particles

	// springs
	int* springIndices;				//!< Spring indices
	float* springCoefficients;		//!< Spring coefficients
	float* springRestLengths;		//!< Spring rest-lengths
	int numSprings;					//!< Number of springs

	// shapes
	int* shapeIndices;				//!< The indices of the shape matching constraints
	int numShapeIndices;			//!< Total number of indices for shape constraints	
	int* shapeOffsets;				//!< Each entry stores the end of the shape's indices in the indices array (exclusive prefix sum of shape lengths)
	float* shapeCoefficients;		//!< The stiffness coefficient for each shape
	float* shapeCenters;			//!< The position of the center of mass of each shape, an array of vec3s mNumShapes in length
	int numShapes;					//!< The number of shape matching constraints

	// plastic deformation
	float* shapePlasticThresholds;	//!< The plastic threshold coefficient for each shape
	float* shapePlasticCreeps;		//!< The plastic creep coefficient for each shape

	// faces for cloth
	int* triangleIndices;			//!< Indexed triangle mesh indices for clothing
	int numTriangles;				//!< Number of triangles

	// inflatable params
	bool inflatable;				//!< Whether an inflatable constraint should be added
	float inflatableVolume;			//!< The rest volume for the inflatable constraint
	float inflatablePressure;		//!< How much over the rest volume the inflatable should attempt to maintain
	float inflatableStiffness;		//!< How stiff the inflatable is
};

/** 
 * Represents an instance of a FlexAsset in a container
 */
struct NvFlexExtInstance
{
	int* particleIndices;			//!< Simulation particle indices
	int numParticles;				//!< Number of simulation particles
	
	int triangleIndex;				//!< Index in the container's triangle array
	int shapeIndex;					//!< Index in the container's shape body constraints array	
	int inflatableIndex;			//!< Index in the container's inflatables array

	float* shapeTranslations;		//!< Shape matching group translations (vec3s)
	float* shapeRotations;			//!< Shape matching group rotations (quaternions)

	const NvFlexExtAsset* asset;	//!< Source asset used to create this instance
	
	void* userData;					//!< User data pointer
};

/**
* Represents a soft joint with a radius overlapping different flex objects
* Each soft joint can be spawned into a container using NvFlexExtCreateSoftJoint()
*/
struct NvFlexExtSoftJoint
{
	int* particleIndices;			//!< Global indices
	float* particleLocalPositions;  //!< Relative offsets from the particles of the joint to the center
	int shapeIndex;					//!< Index in the container's shape body constraints array	
	int numParticles;				//!< Number of particles in the joint

	float shapeTranslations[3];		//!< Joint shape matching group translations (vec3s)
	float shapeRotations[4];		//!< Joint shape matching group rotations (quaternions)

	float stiffness;				//!< Joint stiffness

	bool initialized;				//!< Joint status flag
};

/** 
 * Opaque type representing a simulation
 */
typedef struct NvFlexExtContainer NvFlexExtContainer;

/**
 * Create an index buffer of unique vertices in the mesh (collapses vertices in the same position even if they have different normals / texcoords).
 * This can be used to create simulation meshes from render meshes, and is typically done as a pre-pass before calling NvFlexExtCreateClothFromMesh().
 *
 * @param[in] vertices A pointer to an array of float3 positions
 * @param[in] numVertices The number of vertices in the mesh
 * @param[out] uniqueVerts A list of unique mesh vertex indices, should be numVertices in length (worst case all verts are unique)
 * @param[out] originalToUniqueMap Mapping from the original vertex index to the unique vertex index, should be numVertices in length
 * @param[in] threshold The distance below which two vertices are considered duplicates
 * @return The number of unique vertices in the mesh
 */
NV_FLEX_API int NvFlexExtCreateWeldedMeshIndices(const float* vertices, int numVertices, int* uniqueVerts, int* originalToUniqueMap, float threshold);

/**
 * Create a cloth asset consisting of stretch and bend distance constraints given an indexed triangle mesh. Stretch constraints will be placed along
 * triangle edges, while bending constraints are placed over two edges.
 *
 * @param[in] particles Positions and masses of the particles in the format [x, y, z, 1/m]
 * @param[in] numParticles The number of particles
 * @param[in] indices The triangle indices, these should be 'welded' using NvFlexExtCreateWeldedMeshIndices() first
 * @param[in] numTriangles The number of triangles
 * @param[in] stretchStiffness The stiffness coefficient for stretch constraints
 * @param[in] bendStiffness The stiffness coefficient used for bending constraints
 * @param[in] tetherStiffness If > 0.0f then the function will create tethers attached to particles with zero inverse mass. These are unilateral, long-range attachments, which can greatly reduce stretching even at low iteration counts.
 * @param[in] tetherGive Because tether constraints are so effective at reducing stiffness, it can be useful to allow a small amount of extension before the constraint activates.
 * @param[in] pressure If > 0.0f then a volume (pressure) constraint will also be added to the asset, the rest volume and stiffness will be automatically computed by this function
 * @return A pointer to an asset structure holding the particles and constraints
 */
NV_FLEX_API NvFlexExtAsset* NvFlexExtCreateClothFromMesh(const float* particles, int numParticles, const int* indices, int numTriangles, float stretchStiffness, float bendStiffness, float tetherStiffness, float tetherGive, float pressure);

/**
 * Create a cloth asset consisting of stretch and bend distance constraints given an indexed triangle mesh. This creates an asset with the same
 * structure as NvFlexExtCreateClothFromMesh(), however tether constraints are not supported, and additional information regarding mesh topology
 * will be stored with the asset to allow tearing.
 *
 * @note: Typically each instance of a tearable cloth mesh will have it's own asset. This is because the asset holds the topology of the mesh which is
 * unique for each instance.
 *
 * @param[in] particles Positions and masses of the particles in the format [x, y, z, 1/m]
 * @param[in] numParticles The number of particles
 * @param[in] maxParticles The maximum number of particles for this asset, this will limit the amount of tearing that can be performed.
 * @param[in] indices The triangle indices, these should be 'welded' using NvFlexExtCreateWeldedMeshIndices() first
 * @param[in] numTriangles The number of triangles
 * @param[in] stretchStiffness The stiffness coefficient for stretch constraints
 * @param[in] bendStiffness The stiffness coefficient used for bending constraints
 * @param[in] pressure If > 0.0f then a volume (pressure) constraint will also be added to the asset, the rest volume and stiffness will be automatically computed by this function
 * @return A pointer to an asset structure holding the particles and constraints
 */
NV_FLEX_API NvFlexExtAsset* NvFlexExtCreateTearingClothFromMesh(const float* particles, int numParticles, int maxParticles, const int* indices, int numTriangles, float stretchStiffness, float bendStiffness, float pressure);

/**
 * Destroy an asset created with NvFlexExtCreateTearingClothFromMesh()
 * @param[in] asset The asset to be destroyed.
 */
NV_FLEX_API void NvFlexExtDestroyTearingCloth(NvFlexExtAsset* asset);

/**
 * Particles and vertices may need to be copied during tearing. Because the user may maintain particle data 
 * outside of Flex, this structure describes how to update the particle data. The application should copy each 
 * existing particle given by srcIndex (in the asset's particle array) to the destIndex (also in the asset's array).
 */
struct NvFlexExtTearingParticleClone
{
	int srcIndex;	
	int destIndex;
};

/**
 * The mesh topology may need to be updated to reference new particles during tearing. Because the user may maintain
 * mesh topology outside of Flex, this structure describes the necessary updates that should be performed on the mesh.
 * The triIndex member is the index of the index to be updated, e.g.:
 * a triIndex value of 4 refers to the index 1 vertex (4%3) of the index 1 triangle (4/3). This entry in the indices
 * array should be updated to point to the newParticleIndex.
 */
struct NvFlexExtTearingMeshEdit
{
	int triIndex;			// index into the triangle indices array to update
	int newParticleIndex;	// new value for the index
};

/**
 * Perform cloth mesh tearing, this function will calculate the strain on each distance constraint and perform splits if it is
 * above a certain strain threshold (i.e.: length/restLength > maxStrain).
 *
 * @param[in] asset The asset describing the cloth constraint network, this must be created with NvFlexExtCreateTearingClothFromMesh()
 * @param[in] maxStrain The maximum allowable strain on each edge
 * @param[in] maxSplits The maximum number of constraint breaks that will be performed, this controls the 'rate' of mesh tearing
 * @param[in] particleCopies Pointer to an array of NvFlexExtTearingParticleClone structures that describe the particle copies that need to be performed
 * @param[in] numParticleCopies Pointer to an integer that will have the number of copies performed written to it
 * @param[in] maxCopies The maximum number of particle copies that will be performed, multiple particles copies may be performed in response to one split
 * @param[in] triangleEdits Pointer to an array of NvFlexExtTearingMeshEdit structures that describe the topology updates that need to be performed
 * @param[in] numTriangleEdits Pointer to an integer that will have the number of topology updates written to it
 * @param[in] maxEdits The maximum number of index buffer edits that will be output
 */
NV_FLEX_API void NvFlexExtTearClothMesh(NvFlexExtAsset* asset, float maxStrain,  int maxSplits, NvFlexExtTearingParticleClone* particleCopies, int* numParticleCopies, int maxCopies, NvFlexExtTearingMeshEdit* triangleEdits, int* numTriangleEdits, int maxEdits);

/**
 * Create a shape body asset from a closed triangle mesh. The mesh is first voxelized at a spacing specified by the radius, and particles are placed at occupied voxels.
 *
 * @param[in] vertices Vertices of the triangle mesh
 * @param[in] numVertices The number of vertices
 * @param[in] indices The triangle indices
 * @param[in] numTriangleIndices The number of triangles indices (triangles*3)
 * @param[in] radius The spacing used for voxelization, note that the number of voxels grows proportional to the inverse cube of radius, currently this method limits construction to resolutions < 64^3
 * @param[in] expand Particles will be moved inwards (if negative) or outwards (if positive) from the surface of the mesh according to this factor
 * @return A pointer to an asset structure holding the particles and constraints
 */
NV_FLEX_API NvFlexExtAsset* NvFlexExtCreateRigidFromMesh(const float* vertices, int numVertices, const int* indices, int numTriangleIndices, float radius, float expand);

/**
* Create a shape body asset from a closed triangle mesh. The mesh is first voxelized at a spacing specified by the radius, and particles are placed at occupied voxels.
*
* @param[in] vertices Vertices of the triangle mesh
* @param[in] numVertices The number of vertices
* @param[in] indices The triangle indices
* @param[in] numTriangleIndices The number of triangles indices (triangles*3)
* @param[in] particleSpacing The spacing to use when creating particles
* @param[in] volumeSampling Control the resolution the mesh is voxelized at in order to generate interior sampling, if the mesh is not closed then this should be set to zero and surface sampling should be used instead
* @param[in] surfaceSampling Controls how many samples are taken of the mesh surface, this is useful to ensure fine features of the mesh are represented by particles, or if the mesh is not closed 
* @param[in] clusterSpacing The spacing for shape-matching clusters, should be at least the particle spacing
* @param[in] clusterRadius Controls the overall size of the clusters, this controls how much overlap  the clusters have which affects how smooth the final deformation is, if parts of the body are detaching then it means the clusters are not overlapping sufficiently to form a fully connected set of clusters
* @param[in] clusterStiffness Controls the stiffness of the resulting clusters
* @param[in] linkRadius Any particles below this distance will have additional distance constraints created between them
* @param[in] linkStiffness The stiffness of distance links
* @param[in] globalStiffness If this parameter is > 0.0f, adds an additional global cluster that consists of all particles in the shape. The stiffness of this cluster is the globalStiffness.
* @param[in] clusterPlasticThreshold Particles belonging to rigid shapes that move with a position delta magnitude > threshold will be permanently deformed in the rest pose, if clusterPlasticCreep > 0.0f
* @param[in] clusterPlasticCreep Controls the rate at which particles in the rest pose are deformed for particles passing the deformation threshold
* @return A pointer to an asset structure holding the particles and constraints
*/
NV_FLEX_API NvFlexExtAsset* NvFlexExtCreateSoftFromMesh(const float* vertices, int numVertices, const int* indices, int numTriangleIndices, float particleSpacing, float volumeSampling, float surfaceSampling, float clusterSpacing, float clusterRadius, float clusterStiffness, float linkRadius, float linkStiffness, float globalStiffness, float clusterPlasticThreshold, float clusterPlasticCreep);

/**
 * Frees all memory associated with an asset created by one of the creation methods
 * param[in] asset The asset to destroy.
 */
NV_FLEX_API void NvFlexExtDestroyAsset(NvFlexExtAsset* asset);

/**
* Creates information for linear blend skining a graphics mesh to a set of transforms (bones)
*
* @param[in] vertices Vertices of the triangle mesh
* @param[in] numVertices The number of vertices
* @param[in] bones Pointer to an array of vec3 positions representing the bone positions
* @param[in] numBones Then number of bones
* @param[in] falloff The speed at which the bone's influence on a vertex falls off with distance
* @param[in] maxDistance The maximum distance a bone can be from a vertex before it will not influence it any more
* @param[out] skinningWeights The normalized weights for each bone, there are up to 4 weights per-vertex so this should be numVertices*4 in length
* @param[out] skinningIndices The indices of each bone corresponding to the skinning weight, will be -1 if this weight is not used
*/
NV_FLEX_API void NvFlexExtCreateSoftMeshSkinning(const float* vertices, int numVertices, const float* bones, int numBones, float falloff, float maxDistance, float* skinningWeights, int* skinningIndices);

/**
 * Creates a wrapper object around a Flex solver that can hold assets / instances, the container manages sending and retrieving partical data from the solver
 *
 * @param[in] lib The library instance to use
 * @param[in] solver The solver to wrap
 * @param[in] maxParticles The maximum number of particles to manage
 * @return A pointer to the new container
 */
NV_FLEX_API NvFlexExtContainer* NvFlexExtCreateContainer(NvFlexLibrary* lib, NvFlexSolver* solver, int maxParticles);

/**
 * Frees all memory associated with a container
 *
 * @param[in] container The container to destroy
 */
NV_FLEX_API void NvFlexExtDestroyContainer(NvFlexExtContainer* container);

/**
 * Allocates particles in the container.
 *
 * @param[in] container The container to allocate out of
 * @param[in] n The number of particles to allocate
 * @param[out] indices An n-length array of ints that will store the indices to the allocated particles
 */
NV_FLEX_API int  NvFlexExtAllocParticles(NvFlexExtContainer* container, int n, int* indices);

/**
 * Free allocated particles
 *
 * @param[in] container The container to free from
 * @param[in] n The number of particles to free
 * @param[in] indices The indices of the particles to free
 */
NV_FLEX_API void NvFlexExtFreeParticles(NvFlexExtContainer* container, int n, const int* indices);


/**
 * Retrives the indices of all active particles
 *
 * @param[in] container The container to free from
 * @param[out] indices Returns the number of active particles
 * @return The number of active particles
 */
NV_FLEX_API int NvFlexExtGetActiveList(NvFlexExtContainer* container, int* indices);


struct NvFlexExtParticleData
{
	float* particles;		//!< Receives a pointer to the particle position / mass data								
	float* restParticles;	//!< Receives a pointer to the particle's rest position (used for self collision culling)
	float* velocities;		//!< Receives a pointer to the particle velocity data
	int* phases;			//!< Receives a pointer to the particle phase data
	float* normals;			//!< Receives a pointer to the particle normal data with 16 byte stride in format [nx, ny, nz, nw]

	const float* lower;		//!< Receive a pointer to the particle lower bounds [x, y, z]
	const float* upper;		//!< Receive a pointer to the particle upper bounds [x, y, z]
};

/** 
 * Returns pointers to the internal data stored by the container. These are host-memory pointers, and will 
 * remain valid NvFlexExtUnmapParticleData() is called.
 *
  @param container The container whose data should be accessed
 */
NV_FLEX_API NvFlexExtParticleData NvFlexExtMapParticleData(NvFlexExtContainer* container);
NV_FLEX_API void NvFlexExtUnmapParticleData(NvFlexExtContainer* container);

struct NvFlexExtTriangleData
{
	int* indices;		//!< Receives a pointer to the array of triangle index data
	float* normals;		//!< Receives a pointer to an array of triangle normal data stored with 16 byte stride, i.e.: [nx, ny, nz]
};

/** 
 * Access triangle constraint data, see NvFlexExtGetParticleData() for notes on ownership.
 *
 * @param container The container to retrive from
 */
NV_FLEX_API NvFlexExtTriangleData NvFlexExtMapTriangleData(NvFlexExtContainer* container);

/** 
 * Unmap triangle data, see NvFlexExtMapTriangleData()
 */
NV_FLEX_API void NvFlexExtUnmapTriangleData(NvFlexExtContainer* container);

struct NvFlexExtShapeData
{
	float* rotations;	//!< Receives a pointer to the array quaternion rotation data in [x, y z, w] format
	float* positions;	//!< Receives a pointer to an array of shape body translations in [x, y, z] format
	int n;				//!< Number of valid tranforms
};

/** 
 * Access shape body constraint data, see NvFlexExtGetParticleData() for notes on ownership.
 *
 * @param container The container to retrive from
 */
NV_FLEX_API NvFlexExtShapeData NvFlexExtMapShapeData(NvFlexExtContainer* container);

/** 
 * Unmap shape transform data, see NvFlexExtMapShapeData()
 */
NV_FLEX_API void NvFlexExtUnmapShapeData(NvFlexExtContainer* container);

/**
 * Creates an instance of an asset, the container will internally store a reference to the asset so it should remain valid for the instance lifetime. This
 * method will allocate particles for the asset, assign their initial positions, velocity and phase.
 *
 * @param[in] container The container to spawn into
 * @param[in] particleData Pointer to a mapped particle data struct, returned from NvFlexExtMapParticleData()
 * @param[in] asset The asset to be spawned
 * @param[in] transform A pointer to a 4x4 column major, column vector transform that specifies the initial world space configuration of the particles
 * @param[in] vx The velocity of the particles along the x axis
 * @param[in] vy The velocity of the particles along the y axis
 * @param[in] vz The velocity of the particles along the z axis
 * @param[in] phase The phase used for the particles
 * @param[in] invMassScale A factor applied to the per particle inverse mass
 * @return A pointer to the instance of the asset
 */
NV_FLEX_API NvFlexExtInstance* NvFlexExtCreateInstance(NvFlexExtContainer* container,  NvFlexExtParticleData* particleData, const NvFlexExtAsset* asset, const float* transform, float vx, float vy, float vz, int phase, float invMassScale);

/** Destoy an instance of an asset
 *
 * @param[in] container The container the instance belongs to
 * @param[in] instance The instance to destroy
 */
NV_FLEX_API void NvFlexExtDestroyInstance(NvFlexExtContainer* container, const NvFlexExtInstance* instance);

/** Notifies the container that asset data has changed and needs to be sent to the GPU
 *  this should be called if the constraints for an existing asset are modified by the user
 *
 * @param[in] container The container the instance referencing the asset belongs to
 * @param[in] asset The asset which was modified (can be NULL)
 */
NV_FLEX_API void NvFlexExtNotifyAssetChanged(NvFlexExtContainer* container, const NvFlexExtAsset* asset);

/**
 * Updates the container, applies force fields, steps the solver forward in time, updates the host with the results synchronously.
 * This is a helper function which performs a synchronous update using the following flow.
 *
    \code{.c}
		// async update GPU data
		NvFlexExtPushToDevice(container);

		// update solver
		NvFlexUpdateSolver(container, dt, iterations);

		// async read data back to CPU
		NvFlexExtPullFromDevice(container);

		// read / write particle data on CPU
		NvFlexExtParticleData data = NvFlexExtMapParticleData();

		// CPU particle processing
		ProcessParticles(data);

		// unmap data
		NvFlexExtUnmapParticleData();

  \endcode
  @param[in] container The container to update
  @param[in] dt The time-step in seconds
  @param[in] numSubsteps The number of substeps to perform
  @param[in] enableTimers Whether to record detailed timers, see NvFlexUpdateSolver()
 */
NV_FLEX_API void NvFlexExtTickContainer(NvFlexExtContainer* container, float dt, int numSubsteps, bool enableTimers=false);

/**
 * Updates the device asynchronously, transfers any particle and constraint changes to the flex solver, 
 * expected to be called in the following sequence: NvFlexExtPushToDevice, NvFlexUpdateSolver, NvFlexExtPullFromDevice, flexSynchronize
 * @param[in] container The container to update
 */
NV_FLEX_API void NvFlexExtPushToDevice(NvFlexExtContainer* container);

/**
 * Updates the host asynchronously, transfers particle and constraint data back to he host, 
 * expected to be called in the following sequence: NvFlexExtPushToDevice, NvFlexUpdateSolver, NvFlexExtPullFromDevice
 * @param[in] container The container to update
 */ 
NV_FLEX_API void NvFlexExtPullFromDevice(NvFlexExtContainer* container);

/**
 * Synchronizes the per-instance data with the container's data, should be called after the synchronization with the solver read backs are complete
 *
 * @param[in] container The instances belonging to this container will be updated
 */ 
NV_FLEX_API void NvFlexExtUpdateInstances(NvFlexExtContainer* container);


/** 
 * Controls the way that force fields affect particles
 */
enum NvFlexExtForceMode
{
	//! Apply field value as a force. 
    eNvFlexExtModeForce				=      0,

	//! Apply field value as an impulse. 
    eNvFlexExtModeImpulse			=      1,

	//! Apply field value as a velocity change. 
    eNvFlexExtModeVelocityChange	=      2,
};

/** 
 * Force field data, currently just supports radial fields
 */
struct NvFlexExtForceField
{
	float mPosition[3];		//!< Center of force field
	float mRadius;			//!< Radius of the force field
	float mStrength;		//!< Strength of the force field
	NvFlexExtForceMode mMode;	//!< Mode of field application
	bool mLinearFalloff;	//!< Linear or no falloff 
};

/** 
 * Opaque type representing a force field callback structure that ecapsulates 
 * the force field kernels and associated data applied as a callback during the Flex update
 */
typedef struct NvFlexExtForceFieldCallback NvFlexExtForceFieldCallback;

/**
 * Create a NvFlexExtForceFieldCallback structure, each callback is associated with the
 * passed in solver once the NvFlexExtSetForceFields() is called.
 *
 * @param[in] solver A valid solver created with NvFlexCreateSolver()
 * @return A pointer to a callback structure
 */
NV_FLEX_API NvFlexExtForceFieldCallback* NvFlexExtCreateForceFieldCallback(NvFlexSolver* solver);

/**
 * Destroy the force field callback
 *
 * @param[in] callback A valid solver created with NvFlexExtCreateForceFieldCallback()
 */
NV_FLEX_API void NvFlexExtDestroyForceFieldCallback(NvFlexExtForceFieldCallback* callback);

/**
 * Set force fields on the container, these will be applied during the Flex update
 *
 * @param[in] callback The callback to update
 * @param[in] forceFields A pointer to an array of force field data, may be host or GPU memory
 * @param[in] numForceFields The number of force fields to send to the device
 */
NV_FLEX_API void NvFlexExtSetForceFields(NvFlexExtForceFieldCallback* callback, const NvFlexExtForceField* forceFields, int numForceFields);

/**
* Create a soft joint, the container will internally store a reference to the joint array
*
* @param[in] container The container to spawn into
* @param[in] particleIndices A pointer to an array of particle indices
* @param[in] particleLocalPositions A pointer to an array of particle local positions
* @param[in] numJointParticles The number of particles in the joint
* @param[in] stiffness The stiffness of the joint
* @return A pointer to the soft joint
*/
NV_FLEX_API NvFlexExtSoftJoint* NvFlexExtCreateSoftJoint(NvFlexExtContainer* container, const int* particleIndices, const float* particleLocalPositions, const int numJointParticles, const float stiffness);

/** Destroy a soft joint
*
* @param[in] container The container the joint belongs to
* @param[in] joint The soft joint to destroy
*/
NV_FLEX_API void NvFlexExtDestroySoftJoint(NvFlexExtContainer* container, NvFlexExtSoftJoint* joint);

/** Transform all the local particles of the soft joint
*
* @param[in] container The container to spawn into
* @param[in] joint The soft joint to destroy
* @param[in] position A pointer to a vec3 storing the soft joint new position
* @param[in] rotation A pointer to a quaternion storing the soft joint new rotation
*/
NV_FLEX_API void NvFlexExtSoftJointSetTransform(NvFlexExtContainer* container, NvFlexExtSoftJoint* joint, const float* position, const float* rotation);

} // extern "C"

#endif // NV_FLEX_EXT_H

