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

#ifndef NV_FLEX_DEVICE_H
#define NV_FLEX_DEVICE_H

//! \cond HIDDEN_SYMBOLS
#ifndef NV_FLEX_API
#if _WIN32
#define NV_FLEX_API __declspec(dllexport)
#else
#define NV_FLEX_API
#endif
#endif
//! \endcond

/** \file 
 * NvFlexDevice is an optional helper library that performs some
 * initialization tasks related to GPU device management.
 * The library can be used to query the NVIDIA PhysX control panel for
 * the selected "PhysX" GPU, and to create an optimized CUDA context.
 * Currently the library is a closed source component but is purely optional.
 * See the FlexDemo for an example of how to use the device API.
 */

/**
 * Returns the CUDA ordinal of the GPU selected as "PhysX" in the NVIDIA control panel.
 * Returns -1 if there is no NVIDIA CUDA device available.
 *
 * @note The returned ordinal is a CUDA ordinal and does not correspond to the DXGI
 * ordinal. D3D users should use their own device selection method and pass the appropriate
 * DXGI device index or custom D3D devices to NvFlexInit().
 */
NV_FLEX_API int NvFlexDeviceGetSuggestedOrdinal();

/**
 * Creates a CUDA context optimized for Flex, returns true on success and sets the context 
 * as current on the calling thread. If using this method to initialize CUDA then you should 
 * ensure that no prior CUDA calls are made prior to avoid creating multiple contexts.
 *
 * @param[in] ordinal The CUDA ordinal of the GPU to create the context on, this can be the suggested ordinal (see flexDeviceGetSuggestedOrdinal()), or a manually selected ordinal.
 */
NV_FLEX_API bool NvFlexDeviceCreateCudaContext(int ordinal);

/**
 * Destroy the context associated with the current thread, can be used to destroy the CUDA context created by flexDeviceCreateCudaContext().
 */
NV_FLEX_API void NvFlexDeviceDestroyCudaContext();


#endif // NV_FLEX_DEVICE_H