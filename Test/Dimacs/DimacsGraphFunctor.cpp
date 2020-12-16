
#include "DimacsGraphFunctor.hpp"

namespace Dimacs {

BasicDimacsAdjacencyFunctor::~BasicDimacsAdjacencyFunctor() = default;
DimacsAdjacencyFunctor::~DimacsAdjacencyFunctor() = default;
DimacsNodeCompatibilityFunctor::~DimacsNodeCompatibilityFunctor() = default;
DimacsEdgeCompatibilityFunctor::~DimacsEdgeCompatibilityFunctor() = default;

double DimacsNodeCompatibilityFunctor::estimateNodeSimilarity(
    RIMACS::MappingIndex queryNode,
    RIMACS::MappingIndex targetNode) const
{
  // Default implementation, overload just as demo
  // If you have any information about your data which nodes are favorable to be mapped,
  // then overload estimateNodeSimilarity
  return nodesAreCompatible(queryNode, targetNode) ? 1.0 : 0.0;
}

} // namespace Dimacs