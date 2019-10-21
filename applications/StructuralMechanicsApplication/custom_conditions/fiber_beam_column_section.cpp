// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

// System includes
#include<iostream>

// External includes

// Project includes
#include "custom_conditions/fiber_beam_column_section.hpp"

namespace Kratos {

FiberBeamColumnSection::FiberBeamColumnSection(IndexType NewId)
    : mId(NewId) {}

FiberBeamColumnSection::FiberBeamColumnSection(
    IndexType NewId,
    double position,
    double weight,
    int number_of_fibers_y,
    int number_of_fibers_z,
    double width,
    double height)
: mId(NewId), mPosition(position), mWeight(weight)
{
    const double area = (width / number_of_fibers_y) * (height / number_of_fibers_z);
    // for (int i = 0; i < number_of_fibers_y; ++i){
    //     double y = width / number_of_fibers_y * (i + 0.5);
    //     for (int j = 0; j < number_of_fibers_z; j++)
    //     {
    //         double z = height / number_of_fibers_z * (j + 0.5);
    //         mFibers[i*number_of_fibers_y + j] = FiberBeamColumnUniaxialFiber(
    //             y, z, area //Material law
    //         );
    //     }

    // }
}

// FiberBeamColumnSection::FiberBeamColumnSection(FiberBeamColumnSection const& rOther)
//     : mId(rOther.mId), mPosition(rOther.mPosition), mWeight(rOther.mWeight) {}

// FiberBeamColumnSection::~FiberBeamColumnSection() {}

// FiberBeamColumnSection& FiberBeamColumnSection::operator=(FiberBeamColumnSection const& rOther)
// {
//     return *this;
// }

void FiberBeamColumnSection::Initialize()
{
}


std::string FiberBeamColumnSection::Info() const {
    std::stringstream buffer;
    buffer << "FiberBeamColumnSection #" << Id();
    return buffer.str();
}

void FiberBeamColumnSection::PrintInfo(std::ostream& rOStream) const {
    rOStream << "FiberBeamColumnSection #" << Id();
}

void FiberBeamColumnSection::PrintData(std::ostream& rOStream) const {
    rOStream << "    Position : " << mPosition << std::endl;
    rOStream << "    Weight   : " << mWeight;
}

void FiberBeamColumnSection::save(Serializer& rSerializer) const
{
    rSerializer.save("mId", mId);
    rSerializer.save("mPosition", mPosition);
    rSerializer.save("mWeight", mWeight);
}
void FiberBeamColumnSection::load(Serializer& rSerializer)
{
    rSerializer.load("mId", mId);
    rSerializer.load("mPosition", mPosition);
    rSerializer.load("mWeight", mWeight);
}

} // namespace Kratos.
