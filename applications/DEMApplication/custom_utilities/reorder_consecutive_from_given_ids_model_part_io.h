//
//  Author:    Miguel Angel Celigueta
//

#if !defined(KRATOS_REORDER_CONSECUTIVE_FROM_GIVEN_IDS_MODEL_PART_IO_H_INCLUDED )
#define  KRATOS_REORDER_CONSECUTIVE_FROM_GIVEN_IDS_MODEL_PART_IO_H_INCLUDED


#include "includes/reorder_consecutive_model_part_io.h"

namespace Kratos
{


class ReorderConsecutiveFromGivenIdsModelPartIO : public ReorderConsecutiveModelPartIO
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ReorderConsecutiveFromGivenIdsModelPartIO);
    ReorderConsecutiveFromGivenIdsModelPartIO(std::string const& Filename, const int node_id=0, const int element_id=0, const int condition_id=0, const Flags Options = IO::READ|IO::NOT_IGNORE_VARIABLES_ERROR):
                ReorderConsecutiveModelPartIO(Filename, Options)
    {
        mNumberOfNodes = node_id;
        mNumberOfElements = element_id;
        mNumberOfConditions = condition_id;
    }

    virtual ~ReorderConsecutiveFromGivenIdsModelPartIO(){}

protected:

private:

}; // Class ReorderConsecutiveFromGivenIdsModelPartIO

}  // namespace Kratos.

#endif // KRATOS_REORDER_CONSECUTIVE_FROM_GIVEN_IDS_MODEL_PART_IO_H_INCLUDED  defined

