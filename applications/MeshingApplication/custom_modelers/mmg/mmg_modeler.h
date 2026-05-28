// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "modeler/modeler.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/mmg/mmg_utilities.h"

// NOTE: The following contains the license of the MMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

namespace Kratos
{

/**
 * @class MmgModeler
 * @ingroup MeshingApplication
 * @brief Modeler that wraps MmgProcess for use in the Kratos modeler pipeline.
 * @details This modeler delegates remeshing to MmgProcess and is invoked during
 * the SetupModelPart stage of the modeler pipeline. It is templated on the MMG
 * library variant (MMG2D, MMG3D, or MMGS).
 * @author Vicente Mataix Ferrandiz
 */
template<MMGLibrary TMMGLibrary>
class KRATOS_API(MESHING_APPLICATION) MmgModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(MmgModeler);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor (used for registry prototype only).
     */
    MmgModeler() = default;

    /**
     * @brief Constructor with model and parameters.
     * @param rModel The model owning the model part to remesh.
     * @param ModelerParameters Configuration parameters (must include "model_part_name").
     */
    MmgModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters()
        );

    ~MmgModeler() override = default;

    /**
     * @brief Factory method for the Modeler registry.
     */
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters
        ) const override;

    ///@}
    ///@name Modeler Stages
    ///@{

    /**
     * @brief Runs MMG remeshing on the configured ModelPart.
     */
    void SetupModelPart() override;

    /**
     * @brief Returns the default parameters.
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override { return "MmgModeler"; }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;

    ///@}

}; // class MmgModeler

} // namespace Kratos
