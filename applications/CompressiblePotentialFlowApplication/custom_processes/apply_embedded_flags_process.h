//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez
//

#ifndef KRATOS_APPLY_EMBEDDED_FLAGS_PROCESS_H
#define KRATOS_APPLY_EMBEDDED_FLAGS_PROCESS_H

#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{

class ApplyEmbeddedFlagsProcess: public Process
{
public:

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(ApplyEmbeddedFlagsProcess);

    /// Constructor for ApplyEmbeddedFlagsProcess Process
    ApplyEmbeddedFlagsProcess(ModelPart& rModelPart
                    ):
        Process(),
        mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~ApplyEmbeddedFlagsProcess() override {}

    void operator()()
    {
        Execute();
    }


    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY;

        for(auto it=mrModelPart.ElementsBegin(); it!=mrModelPart.ElementsEnd(); ++it)
        {
            it->Set(BOUNDARY,false);
            it->Set(ACTIVE,true);

            auto geom = it->GetGeometry();
            bool has_positive_values = false;
            bool has_negative_values = false;
            for(unsigned int i=0; i<geom.size(); ++i)
            {
                double distance = geom[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE);
                if (distance > 0.0)
                    has_positive_values = true;
                else
                    has_negative_values = true;
            }

            if (has_positive_values && has_negative_values)
                it->Set(BOUNDARY,true);
            else if (has_negative_values)
                it->Set(ACTIVE,false);
        }

        KRATOS_CATCH("");
    }


    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyEmbeddedFlagsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyEmbeddedFlagsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


private:

    ModelPart& mrModelPart;

    ApplyEmbeddedFlagsProcess& operator=(ApplyEmbeddedFlagsProcess const& rOther);

    /// Copy constructor.
    ApplyEmbeddedFlagsProcess(ApplyEmbeddedFlagsProcess const& rOther);

}; // Class Process

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyEmbeddedFlagsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyEmbeddedFlagsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos


#endif // KRATOS_APPLY_FLAGS_PROCESS_H
