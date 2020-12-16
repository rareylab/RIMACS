
#pragma once

#include <limits>
#include <vector>
#include <memory>

namespace RIMACS {

using Index = size_t;
using MappingIndex = unsigned;
using ComponentIndex = int;

class Config;

enum class NodeMappingStatus : MappingIndex {
  Mappable=std::numeric_limits<MappingIndex>::max(),
  Unmappable = Mappable - 1
};

using EquivalenceClasses = std::vector<MappingIndex>;

using CompatibilityHint = std::vector<std::vector<MappingIndex>>;

class AdjacencyFunctor;

class NodeCompatibilityFunctor;

class EdgeCompatibilityFunctor;

class MCSResult;
using MCSResults = std::vector<MCSResult>;

using CompatibilityFunctors = std::pair<std::reference_wrapper<const NodeCompatibilityFunctor>,
                                        std::reference_wrapper<const EdgeCompatibilityFunctor>>;
using CompatibilityFunctorsPtr = std::pair<std::shared_ptr<const NodeCompatibilityFunctor>,
                                           std::shared_ptr<const EdgeCompatibilityFunctor>>;

struct MappingPair;
struct MappingTuple;
struct AdjacentNode;
template<bool use_weights>
using MappingTupleType = std::conditional_t<use_weights, MappingTuple, MappingPair>;

} // namespace RIMACS

// Tests

class VertexMappingTester;
class SubgraphMappingTester;
