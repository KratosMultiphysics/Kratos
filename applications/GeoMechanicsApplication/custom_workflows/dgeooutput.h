// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

// System includes
#include <includes/gid_io.h>

class NodeOperation
{
public:
    virtual ~NodeOperation() = default;
    virtual void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part);
};
class NodeDISPLACEMENT : public NodeOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class NodeTOTAL_DISPLACEMENT : public NodeOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class NodeWATER_PRESSURE : public NodeOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class NodeNORMAL_FLUID_FLUX : public NodeOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class NodeVOLUME_ACCELERATION : public NodeOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class NodeHYDRAULIC_DISCHARGE : public NodeOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class NodeHYDRAULIC_HEAD : public NodeOperation
{
public:
    void write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) override;
};

class GaussOperation
{
public:
    virtual ~GaussOperation() = default;
    virtual void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part);
};

class GaussFLUID_FLUX_VECTOR : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class GaussHYDRAULIC_HEAD : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class GaussLOCAL_FLUID_FLUX_VECTOR : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class GaussLOCAL_PERMEABILITY_MATRIX : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class GaussPERMEABILITY_MATRIX : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class GaussDEGREE_OF_SATURATION : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class GaussDERIVATIVE_OF_SATURATION : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class GaussRELATIVE_PERMEABILITY : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class GaussPIPE_ACTIVE : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};
class GaussPIPE_HEIGHT : public GaussOperation
{
public:
    void write(Kratos::GidIO<> &gid_io, Kratos::ModelPart &model_part) override;
};

namespace Kratos
{
    class KratosGeoOutput
    {
    public:
        KratosGeoOutput() = delete;
        ~KratosGeoOutput() {};
        static void outputGiD(Model& model, ModelPart& model_part, Parameters parameters, std::string workingDirectory);
    };
}