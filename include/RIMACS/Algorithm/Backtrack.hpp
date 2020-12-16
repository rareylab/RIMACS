
#pragma once

#include <RIMACS/Forward.hpp>

#include <RIMACS/Algorithm/SubgraphMapping.hpp>
#include <RIMACS/Algorithm/VertexMapping.hpp>


namespace RIMACS {

bool backtrack(
    SubgraphMapping<false>& mapping,
    VertexMapping<false>& vm);

bool backtrack(
    SubgraphMapping<true>& mapping,
    VertexMapping<true>& vm);

} // namespace RIMACS
