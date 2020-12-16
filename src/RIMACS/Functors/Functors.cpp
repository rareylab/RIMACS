
#include <queue>

#include <RIMACS/Functors/Functors.hpp>

namespace RIMACS {

// "Implement" the destructors in the cpp, such that at least one method is declared out of
// line for the functors.
// This way their virtual table is not emitted into every translation unit.
CachingRequiredMinimalAdjacencyFunctor::~CachingRequiredMinimalAdjacencyFunctor() = default;
AdjacencyFunctor::~AdjacencyFunctor() = default;
NodeCompatibilityFunctor::~NodeCompatibilityFunctor() = default;
EdgeCompatibilityFunctor::~EdgeCompatibilityFunctor() = default;


std::vector<ComponentIndex> AdjacencyFunctor::getNodeOrderSuggestionVector(
    double partial) const
{
  MappingIndex nofNodes = getNofNodes();
  std::vector<ComponentIndex> atomOrderSuggestion;
  atomOrderSuggestion.reserve(nofNodes + 1);
  std::vector<std::pair<MappingIndex, MappingIndex>> maxDistances;
  maxDistances.reserve(nofNodes);
  for(Index i = 0; i < nofNodes; ++i) {
    std::vector<bool> visited(nofNodes, false);
    std::queue<std::pair<MappingIndex, MappingIndex>> queue({{i, 0}});
    visited[i] = true;
    MappingIndex maxDist = 0;
    while(!queue.empty()) {
      std::pair<MappingIndex, MappingIndex> elem = queue.front();
      queue.pop();
      for(const MappingIndex& neighbour : nodeGetNeighbours(elem.first)) {
        if(!visited[neighbour]) {
          visited[neighbour] = true;
          maxDist = elem.second + 1;
          queue.emplace(neighbour, maxDist);
        }
      }
    }
    maxDistances.emplace_back(maxDist, i);
  }
  std::sort(maxDistances.begin(), maxDistances.end());
  Index lastDist = maxDistances.front().first;
  ComponentIndex nofDistances = 0;
  atomOrderSuggestion.resize(nofNodes);
  for(const std::pair<MappingIndex, MappingIndex>& distancePair: maxDistances) {
    atomOrderSuggestion[distancePair.second] = (nofDistances += static_cast<ComponentIndex>(lastDist != distancePair.first));
    lastDist = distancePair.first;
  }
  atomOrderSuggestion.push_back(nofDistances);
  ComponentIndex refDist = (nofDistances) * partial;
  for(ComponentIndex& dist : atomOrderSuggestion) {
    dist = std::abs(dist - refDist);
  }
  atomOrderSuggestion.back() = std::max(refDist, atomOrderSuggestion.back());
  return atomOrderSuggestion;
}

} // namespace RIMACS

