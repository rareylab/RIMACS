
#include <utility>
#include <functional>
#include <unordered_map>

#include <RIMACS/ASSERT.hpp>

#include <RIMACS/Algorithm/VertexMapping.hpp>
#include <RIMACS/Algorithm/SubgraphMapping.hpp>
#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Algorithm/Backtrack.hpp>
#include <RIMACS/Config.hpp>
#include <RIMACS/Functors/LineGraph.hpp>
#include <RIMACS/Functors/GraphCachingFunctor.hpp>
#include <RIMACS/Functors/CachingEdgeCompatibilityFunctor.hpp>

#include <RIMACS/MCSRunner.hpp>


namespace RIMACS {

namespace {

/**
 * @brief The SwappedNodeCompatibilityFunctor class
 *  Internal class swapping the input parameters
 *  for all calls to a base instance of NodeCompatibilityFunctor
 */
class SwappedNodeCompatibilityFunctor : public NodeCompatibilityFunctor {
public:
  explicit SwappedNodeCompatibilityFunctor(
      const NodeCompatibilityFunctor& base)
    : m_base(base)
  {}


  bool nodesAreCompatible(
      MappingIndex from,
      MappingIndex to) const override
  {
    return m_base.nodesAreCompatible(to, from);
  }


  double getWeight(
      MappingIndex from,
      MappingIndex to) const override
  {
    return m_base.getWeight(to, from);
  }

  double estimateNodeSimilarity(
      MappingIndex from,
      MappingIndex to) const override
  {
    return m_base.estimateNodeSimilarity(to, from);
  }


  bool mappingIsValid(
      typename std::vector<MappingPair>::const_iterator beginMapping,
      typename std::vector<MappingPair>::const_iterator endMapping) const override
  {
    return m_base.mappingIsValid(do_swap(beginMapping), do_swap(endMapping));
  }


  bool mappingIsValid(
      SwappedMappingIterator beginMapping,
      SwappedMappingIterator endMapping) const override
  {
    RIMACS_ASSERT(false && "There should never be a twice swapped node mapping");
    return true;
  }

  double estimateFinalWeight(
      typename std::vector<MappingPair>::const_iterator beginMapping,
      typename std::vector<MappingPair>::const_iterator endMapping,
      double currentWeight) const override
  {
    return m_base.estimateFinalWeight(do_swap(beginMapping), do_swap(endMapping),
                                      currentWeight);
  }
  double estimateFinalWeight(
      SwappedMappingIterator /*beginMapping*/,
      SwappedMappingIterator /*endMapping*/,
      double currentWeight) const override
  {
    RIMACS_ASSERT(false && "There should never be a twice swapped node mapping");
    return 0.0;
  }


private:
  const NodeCompatibilityFunctor& m_base;
};


/**
 * @brief The SwappedEdgeCompatibilityFunctor class
 *  Internal class swapping the input parameters
 *  for all calls to a base instance of EdgeCompatibilityFunctor
 */
class SwappedEdgeCompatibilityFunctor : public EdgeCompatibilityFunctor {
public:

  explicit SwappedEdgeCompatibilityFunctor(
      const EdgeCompatibilityFunctor& base)
    : m_base(base)
  {}


  bool edgesAreCompatible(
      MappingIndex from,
      MappingIndex to) const override
  {
    return m_base.edgesAreCompatible(to, from);
  }


private:
  const EdgeCompatibilityFunctor& m_base;
};



// Helper function extending an initial mapping.
bool test_extension_node_compatibility(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    std::unordered_map<MappingIndex, MappingIndex>& nodeMappingsG1,
    std::unordered_map<MappingIndex, MappingIndex>& nodeMappingsG2,
    const MappingPair& mp) {
  if(!nodeComp.nodesAreCompatible(mp.m_from, mp.m_to)) {
    return false;
  }
  if(nodeMappingsG1.empty()) {
    return true;
  }
  for(const MappingIndex& neighbour : query.nodeGetNeighbours(mp.m_from)) {
    auto it = nodeMappingsG1.find(neighbour);
    if(it != nodeMappingsG1.end()) {
      if(!edgeComp.edgesAreCompatible(query.getEdgeId(mp.m_from, neighbour),
                                      target.getEdgeId(mp.m_to, it->second))) {
        return false;
      }
    }
  }
  for(const MappingIndex& neighbour : target.nodeGetNeighbours(mp.m_to)) {
    auto it = nodeMappingsG2.find(neighbour);
    if(it != nodeMappingsG2.end()) {
      if(!edgeComp.edgesAreCompatible(query.getEdgeId(mp.m_from, it->second),
                                      target.getEdgeId(mp.m_to, neighbour))) {
        return false;
      }
    }
  }
  return true;
}


template<bool weighted>
bool handle_initial_mapping(
    SubgraphMapping<weighted>& mapping,
    std::unique_ptr<VertexMapping<weighted>>& vmPtr,
    const std::vector<MappingPair>& initialMapping,
    bool assumeInitialMappingAsCompatible)
{
  std::unordered_map<MappingIndex, MappingIndex> nodeMappingsG1;
  std::unordered_map<MappingIndex, MappingIndex> nodeMappingsG2;
  if(assumeInitialMappingAsCompatible) {
    RIMACS_ASSERT(initialMapping.size() == 1);
    const MappingPair& mp = initialMapping.front();
    mapping.updateConnectedComponents(mp.m_from, *vmPtr, true);
    vmPtr = std::make_unique<VertexMapping<weighted>>(mapping.handleInitialMapping(*vmPtr, mp));
  } else {
    for(const MappingPair& mp: initialMapping) {
      if(test_extension_node_compatibility(mapping.getQueryGraph(), mapping.getTargetGraph(),
                                           mapping.getNodeComp(), mapping.getEdgeComp(),
                                           nodeMappingsG1, nodeMappingsG2, mp)) {
        mapping.updateConnectedComponents(mp.m_from, *vmPtr, true);
        vmPtr = std::make_unique<VertexMapping<weighted>>(mapping.handleInitialMapping(*vmPtr, mp));
        nodeMappingsG1.emplace(mp.m_from, mp.m_to);
        nodeMappingsG2.emplace(mp.m_to, mp.m_from);
      } else {
        return false;
      }
    }
  }
  return true;
}

template<bool weighted>
std::pair<std::vector<MCSResult>, Index> run_mcs_impl(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    const std::vector<MappingPair>& initialMapping,
    const Config& config,
    const std::vector<MappingIndex>* queryEquivalenceClasses,
    const std::vector<MappingIndex>* targetEquivalenceClasses,
    bool ccFeature,
    bool assumeInitialMappingAsCompatible)
{
  RIMACS_ASSERT((!queryEquivalenceClasses || queryEquivalenceClasses->size() == query.getNofNodes())
                && "The provided equivalence classes must correspond to the graph, meaning "
                "that for each node of the graph, there must be an entry of an equivalence class");
  RIMACS_ASSERT((!targetEquivalenceClasses || targetEquivalenceClasses->size() == target.getNofNodes())
                && "The provided equivalence classes must correspond to the graph, meaning "
                "that for each node of the graph, there must be an entry of an equivalence class");
  SubgraphMapping<weighted> mapping(query, target, nodeComp, edgeComp, config,
                                    queryEquivalenceClasses, targetEquivalenceClasses);

  mapping.enableConnectedComponentTrackingFeature(ccFeature);
  std::unique_ptr<VertexMapping<weighted>> vmPtr;
  if(assumeInitialMappingAsCompatible && initialMapping.size() == 1) {
    vmPtr = std::make_unique<VertexMapping<weighted>>(query, target, nodeComp, initialMapping.front());
  } else {
    vmPtr = std::make_unique<VertexMapping<weighted>>(query, target, nodeComp);
  }
  mapping.markUnmappableNodes(*vmPtr);
  if(!initialMapping.empty()) {
    if(!handle_initial_mapping(mapping, vmPtr, initialMapping,
                               assumeInitialMappingAsCompatible)) {
      return std::make_pair(MCSResults(), Index{0});
    }
  }
  backtrack(mapping, *vmPtr);
  return std::make_pair(mapping.moveMCSResults(), mapping.getCount());
}


} // namespace


std::pair<std::vector<MCSResult>, Index> MCSRunner::runGenericMCSIncludingCount(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    const std::vector<MappingPair>& initialMapping,
    const Config& config,
    bool weighted,
    const std::vector<MappingIndex>* queryEquivalenceClasses,
    const std::vector<MappingIndex>* targetEquivalenceClasses,
    bool ccFeature,
    bool assumeInitialMappingAsCompatible)
{
  bool needSwap = target.getNofNodes() < query.getNofNodes();
  if(needSwap) {
    SwappedNodeCompatibilityFunctor swappedNodeComp(nodeComp);
    SwappedEdgeCompatibilityFunctor swappedEdgeComp(edgeComp);
    std::vector<MappingPair> ini = initialMapping;
    std::for_each(ini.begin(), ini.end(), [](MappingPair& mp) {std::swap(mp.m_from, mp.m_to);});
    std::pair<std::vector<MCSResult>, Index> tmpRes
        = runGenericMCSIncludingCount(target, query, swappedNodeComp, swappedEdgeComp, ini, config,
                                      weighted, targetEquivalenceClasses, queryEquivalenceClasses, ccFeature,
                                      assumeInitialMappingAsCompatible);
    for(MCSResult& res : tmpRes.first) {
      res.doSwap(initialMapping.size());
    }
    return tmpRes;
  }
  std::unique_ptr<AdjacencyFunctor> cachedQueryAdj;
  std::unique_ptr<AdjacencyFunctor> cachedTargetAdj;
  std::unique_ptr<EdgeCompatibilityFunctor> cachedEdgeComp;
  std::unique_ptr<std::vector<MappingIndex>> otherEquivalenceClasses;
  // Ensure that the graphs adjacency is cached
  if(!dynamic_cast<const CachingAdjacencyFunctor*>(&query)) {
    cachedQueryAdj = makeUniqueReferenceGraphCachingFunctor(query);
  }
  if(!dynamic_cast<const CachingAdjacencyFunctor*>(&target)) {
    cachedTargetAdj = makeUniqueReferenceGraphCachingFunctor(target);
  }
  if(!dynamic_cast<const CachingEdgeCompatibilityFunctor*>(&edgeComp) && !dynamic_cast<const LineGraph*>(&target)) {
    cachedEdgeComp = std::make_unique<CachingEdgeCompatibilityFunctor>(edgeComp,
                                                                       query.getNofEdges(), target.getNofEdges());
  }
  if(static_cast<bool>(queryEquivalenceClasses) != static_cast<bool>(targetEquivalenceClasses)) {
    if(!queryEquivalenceClasses) {
      otherEquivalenceClasses = std::make_unique<std::vector<MappingIndex>>(query.getNofNodes());
      queryEquivalenceClasses = otherEquivalenceClasses.get();
    } else {
      otherEquivalenceClasses = std::make_unique<std::vector<MappingIndex>>(target.getNofNodes());
      targetEquivalenceClasses = otherEquivalenceClasses.get();
    }
    std::iota(otherEquivalenceClasses->begin(), otherEquivalenceClasses->end(), 0);
  }
  if(weighted) {
    return run_mcs_impl<true>(cachedQueryAdj ? *cachedQueryAdj : query,
                              cachedTargetAdj ? *cachedTargetAdj : target,
                              nodeComp, cachedEdgeComp ? *cachedEdgeComp : edgeComp, initialMapping,
                              config, queryEquivalenceClasses, targetEquivalenceClasses, ccFeature,
                              assumeInitialMappingAsCompatible);
  }
  return run_mcs_impl<false>(cachedQueryAdj ? *cachedQueryAdj : query,
                             cachedTargetAdj ? *cachedTargetAdj : target,
                             nodeComp, cachedEdgeComp ? *cachedEdgeComp : edgeComp, initialMapping,
                             config, queryEquivalenceClasses, targetEquivalenceClasses, ccFeature,
                             assumeInitialMappingAsCompatible);
}


std::pair<std::vector<MCSResult>, Index> MCSRunner::runGenericMCSIncludingCount(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    const MappingPair& initialMapping,
    const Config& config,
    bool weighted,
    const std::vector<MappingIndex>* queryEquivalenceClasses,
    const std::vector<MappingIndex>* targetEquivalenceClasses,
    bool ccFeature,
    bool assumeInitialMappingAsCompatible)
{
  if(initialMapping.m_from != std::numeric_limits<MappingIndex>::max()
     && initialMapping.m_to != std::numeric_limits<MappingIndex>::max()) {
    return runGenericMCSIncludingCount(query, target, nodeComp, edgeComp,
                                       std::vector<MappingPair>({initialMapping}), config,
                                       weighted, queryEquivalenceClasses, targetEquivalenceClasses,
                                       ccFeature, assumeInitialMappingAsCompatible);
  }
  return runGenericMCSIncludingCount(query, target, nodeComp, edgeComp, std::vector<MappingPair>(), config,
                                     weighted, queryEquivalenceClasses, targetEquivalenceClasses,
                                     ccFeature, assumeInitialMappingAsCompatible);
}


std::pair<std::vector<MCSResult>, Index> MCSRunner::runGenericMCSIncludingCount(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    const Config& config,
    bool weighted,
    const std::vector<MappingIndex>* queryEquivalenceClasses,
    const std::vector<MappingIndex>* targetEquivalenceClasses,
    bool ccFeature)
{
  return runGenericMCSIncludingCount(query, target, nodeComp, edgeComp, std::vector<MappingPair>(),
                                     config, weighted, queryEquivalenceClasses, targetEquivalenceClasses,
                                     ccFeature);
}


std::pair<std::vector<MCSResult>, Index> MCSRunner::runGenericMCESIncludingCount(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    const Config& config,
    bool ccFeature)
{
  LineGraph queryLine(query);
  LineGraph targetLine(target);
  LineGraphNodeCompatibilityFunctor lineNodeComp(queryLine, targetLine, nodeComp, edgeComp);
  LineGraphEdgeCompatibilityFunctor lineEdgeComp(queryLine, targetLine, nodeComp);

  return runGenericMCSIncludingCount(queryLine, targetLine, lineNodeComp, lineEdgeComp, config, false,
                                     nullptr, nullptr, ccFeature);
}


std::tuple<std::vector<MCSResult>, Index, double> MCSRunner::runGenericNodeMappingMCES(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    const Config& config,
    bool ccFeature)
{
  LineGraph queryLine(query);
  LineGraph targetLine(target);
  PreferMaxNumberOfMappedNodesLineGraphNodeCompatibilityFunctor lineNodeComp(queryLine, targetLine,
                                                                             nodeComp, edgeComp, config);
  LineGraphEdgeCompatibilityFunctor lineEdgeComp(queryLine, targetLine, nodeComp);

  std::pair<std::vector<MCSResult>, Index> base_result
      = runGenericMCSIncludingCount(queryLine, targetLine, lineNodeComp, lineEdgeComp, config,
                                    false, nullptr, nullptr, ccFeature);
  std::tuple<std::vector<MCSResult>, Index, double> result;
  std::get<2>(result) = std::numeric_limits<double>::quiet_NaN();
  EdgeMappingConverter conv{queryLine, targetLine, nodeComp};
  if(!base_result.first.empty()) {
    double mappedEdges = base_result.first.front().size();
    std::get<2>(result) = mappedEdges / (queryLine.getNofNodes() + targetLine.getNofNodes() - mappedEdges);
  }
  if(config.resultType() == Config::ResultType::Maximum) {
    std::for_each(base_result.first.begin(), base_result.first.end(), [&conv](MCSResult& res) {
      std::vector<MappingPair> mapping;
      std::vector<MappingIndex> compSizes;
      std::tie(mapping, compSizes) = conv(res.getMappings(), res.getComponentSizes());
      res = MCSResult(std::move(mapping), std::move(compSizes), res.size());
    });
  } else {
    std::vector<MCSResult> convertedResults;
    convertedResults.reserve(base_result.first.size());
    std::for_each(base_result.first.begin(), base_result.first.end(),
                  [&conv, &convertedResults](const MCSResult& res) {
      auto multiConversionVector = conv(res.getMappings(), res.getComponentSizes(), true);
      for(auto& resultPair : multiConversionVector) {
        convertedResults.emplace_back(std::move(resultPair.first), std::move(resultPair.second),
                                      res.size());
      }
    });
    std::sort(convertedResults.begin(), convertedResults.end());
    auto newEndIt = std::unique(convertedResults.begin(), convertedResults.end());
    RIMACS_ASSERT(newEndIt == convertedResults.end());
    convertedResults.erase(newEndIt, convertedResults.end());
    base_result.first = std::move(convertedResults);
  }
  std::get<0>(result) = std::move(base_result.first);
  std::get<1>(result) = base_result.second;
  return result;
}


namespace {


template<bool weighted>
std::pair<MCSResults, Index> run_hinted_impl(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    const std::vector<std::vector<MappingIndex>>& compatibilityHint,
    const Config& config,
    const std::vector<MappingPair>& initialMapping,
    const std::vector<MappingIndex>* queryEquivalenceClasses,
    const std::vector<MappingIndex>* targetEquivalenceClasses)
{
  SubgraphMapping<weighted> mapping(query, target, nodeComp, edgeComp, config,
                                    queryEquivalenceClasses, targetEquivalenceClasses);
  auto vmPtr = std::make_unique<VertexMapping<weighted>>(query, target, nodeComp, compatibilityHint);
  if(!initialMapping.empty()) {
    if(!handle_initial_mapping(mapping, vmPtr, initialMapping, false)) {
      return std::make_pair(MCSResults(), Index{0});
    }
  }
  backtrack(mapping, *vmPtr);
  return std::make_pair(mapping.moveMCSResults(), mapping.getCount());
}


} // namespace


std::pair<MCSResults, Index> MCSRunner::runHintedMCS(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    const std::vector<std::vector<MappingIndex>>& compatibilityHint,
    const Config& config,
    const std::vector<MappingPair>& initialMapping,
    bool weighted,
    const std::vector<MappingIndex>* queryEquivalenceClasses,
    const std::vector<MappingIndex>* targetEquivalenceClasses)
{
  if(weighted) {
    return run_hinted_impl<true>(query, target, nodeComp, edgeComp, compatibilityHint, config,
                                 initialMapping, queryEquivalenceClasses, targetEquivalenceClasses);
  } else {
    return run_hinted_impl<false>(query, target, nodeComp, edgeComp, compatibilityHint, config,
                                  initialMapping, queryEquivalenceClasses, targetEquivalenceClasses);
  }
}

} // namespace RIMACS
