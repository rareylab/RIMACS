
#pragma once

#include <numeric>
#include <cmath>

#include <boost/range/adaptor/reversed.hpp>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Algorithm/VertexMapping.hpp>
#include <RIMACS/Config.hpp>
#include <RIMACS/MCSResult.hpp>
#include <RIMACS/Algorithm/CompareUtils.hpp>
#include <RIMACS/Utils/SortedVector.hpp>

#include <RIMACS/Algorithm/WeightTypeTraits.hpp>
#include <RIMACS/Algorithm/SubgraphMapping.hpp>

namespace RIMACS {

namespace {

template<bool weighted>
std::enable_if_t<!weighted, ComponentIndex> expand_to_connected_component_size(
    ComponentIndex nofNodes)
{
  return nofNodes;
}
template<bool weighted>
std::enable_if_t<weighted, std::pair<ComponentIndex, double>> expand_to_connected_component_size(
    ComponentIndex nofNodes)
{
  return std::pair<ComponentIndex, double>(nofNodes, nofNodes);
}

} // namespace

template<bool weighted>
SubgraphMapping<weighted>::SubgraphMapping(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const EdgeCompatibilityFunctor& edgeComp,
    const Config& config,
    const std::vector<MappingIndex>* queryEquivalenceClasses,
    const std::vector<MappingIndex>* targetEquivalenceClasses,
    Base::UserLog* info)
  : m_query(query)
  , m_target(target)
  , m_compatibility(nodeComp, edgeComp)
  , m_nodeStatus(query.getNofNodes(), static_cast<MappingIndex>(NodeMappingStatus::Mappable))
  , m_equivalentMcs(Compare::custom_compare_MCSResult_compare<>(queryEquivalenceClasses,
                                                                targetEquivalenceClasses))
  , m_config(config)
  , m_currentWeight(weighted, 0.0)
  , m_queryEquivalenceClasses(queryEquivalenceClasses)
  , m_targetEquivalenceClasses(targetEquivalenceClasses)
  , m_connectedComponentIds(std::vector<ComponentIndex>(query.getNofNodes(), -1))
  , m_connectedSizes(1, expand_to_connected_component_size<weighted>(
      -static_cast<ComponentIndex>(query.getNofNodes())))
  , m_maxConnectedComponentIndex(1)
  , m_use_cc_feature(getConnectedComponentFeatureSuggestion(config))
  , m_info(info)
{
  m_connectedComponentSizes.reserve(std::min(config.nofComponents(),
                                             m_query.get().getNofNodes() / 2 + 1));
  m_connectedComponentSizes.push_back(0);
  m_nodeMappings.reserve(query.getNofNodes());
}


template<bool weighted>
std::vector<MCSResult>&& SubgraphMapping<weighted>::moveMCSResults()
{
  if(!m_equivalentMcs.empty()) {
    m_mcs = m_equivalentMcs.moveOut();
    for(MCSResult& result : m_mcs) {
      if(result.m_nodeMapping.size() > m_initialMappingSize) {
        result.m_cumulatedComponentSizes.insert(result.m_cumulatedComponentSizes.begin(),
                                                m_initialMappingSize);
        result.sortResults<false>(m_initialMappingSize);
      }
    }
  }
  return std::move(m_mcs);
}


template<bool weighted>
void SubgraphMapping<weighted>::removeMapping(
    const std::vector<MappingIndex>& locallyUnmappableNodes)
{
  --m_connectedComponentSizes.front();
  for(const MappingIndex& node : locallyUnmappableNodes) {
    if(m_connectedComponentIds.at(node) < 0) {
      restoreConnectedComponents(node);
    }
  }
  m_nodeStatus.at(m_nodeMappings.back().m_from)
      = static_cast<MappingIndex>(NodeMappingStatus::Mappable);
  m_nodeMappings.pop_back();
  if(weighted) {
    m_currentWeight.pop_back();
  }
  m_currentExpansionIsValid = true;
}

// can't use namespace {} for the following inline functions.
// Not all cpp files including this template implementation header use all functions, such that
// the compiler assumes, they were defined but not used.

inline std::pair<ComponentIndex, double> add(
    const std::pair<ComponentIndex, double>& p1,
    const std::pair<ComponentIndex, double>& p2)
{
  return std::make_pair(p1.first + p2.first, p1.second + p2.second);
}


template<typename T>
inline T add(
    const T& lhs,
    const T& rhs)
{
  return lhs + rhs;
}


template<bool weighted>
typename SubgraphMapping<weighted>::ConnectionSize SubgraphMapping<weighted>::countExtendableNodes() const
{
  std::vector<ConnectionSize> selectableComponents = m_connectedSizes;
  // The first component contains all initially unmappable nodes
  getComponentSize(selectableComponents, 0) = 0;
  if(weighted) {
    selectableComponents.erase(std::remove_if(selectableComponents.begin() + 1,
                                              selectableComponents.end(),
                                              [](const ConnectionSize& cs) {
                                                  return getCompatibleNode(cs) == 0;
                                              }),
                               selectableComponents.end());
  }
  std::sort(selectableComponents.begin(), selectableComponents.end(), WeightedNodeCompare());
  auto componentSearcher = [](const ConnectionSize& cs) {
    return getCompatibleNode(cs) == ComponentIndex{0};
  };
  auto disconnectedComponentSearcher = [this](const ConnectionSize& cs) {
    return getComponentSize(cs) > -static_cast<ComponentIndex>(m_config.componentSize());
  };
  auto connectedNodeEnd = std::find_if(selectableComponents.begin(), selectableComponents.end(),
                                       componentSearcher);
  ConnectionSize connectedSize
      = std::accumulate(selectableComponents.begin(), connectedNodeEnd, ConnectionSize{},
                          ConnectionSizeAccumulator<true>());
  ConnectionSize disconnectedSize{};
  if(m_config.isDisconnected()) {
    ComponentIndex availableDisconnectedComponents
        = m_config.nofComponents() - (m_connectedComponentSizes.size());
    auto disconnectedNodeEnd = std::find_if(selectableComponents.rbegin(), selectableComponents.rend(),
                                            disconnectedComponentSearcher);
    if(std::distance(selectableComponents.rbegin(), disconnectedNodeEnd)
       > availableDisconnectedComponents) {
      disconnectedSize
          = std::accumulate(selectableComponents.rbegin(),
                            selectableComponents.rbegin() + availableDisconnectedComponents,
                            ConnectionSize{}, ConnectionSizeAccumulator<false>());
    } else {
      disconnectedSize
          = std::accumulate(selectableComponents.rbegin(), disconnectedNodeEnd,
                            ConnectionSize{}, ConnectionSizeAccumulator<false>());
    }
  } else if(getComponentSize(connectedSize) == 0 && m_connectedComponentSizes.front() == 0) {
    RIMACS_ASSERT(false && "Never  use this, at least one component should be made available by "
                        "the selection of a new extendable node");
  }
  return add(connectedSize, disconnectedSize);
}


template<bool weighted>
bool SubgraphMapping<weighted>::canBeDisconnectedExpanded(
    const VertexMappingType& vm)
{
  // do not test whether the actual component size is positive.
  // Doing so would mark the mapping as invalid for longer than actual needed
  m_currentExpansionIsValid
      = getNodeComp().mappingIsValid(m_nodeMappings.end() - m_connectedComponentSizes.front(),
                                     m_nodeMappings.end());
  if(!m_currentExpansionIsValid) {
    return false;
  }
  if(!m_use_cc_feature) {
    return m_config.isDisconnected()
           && m_config.componentSize() <= m_connectedComponentSizes.front()
           && m_config.componentSize() <= getComponentSize(vm.getNofExtendableNodes())
           && m_connectedComponentSizes.size() < m_config.nofComponents();
  }
  if(m_config.isDisconnected()
     && m_config.componentSize() <= m_connectedComponentSizes.front()
     && m_connectedComponentSizes.size() < m_config.nofComponents()) {
    auto extendableNodes = countExtendableNodes();
    if(!getComponentSize(extendableNodes)) {
      return false;
    }
    return isExtendable(vm);
  }
  return false;
}


template<bool weighted>
bool SubgraphMapping<weighted>::isExtendable(
    const VertexMappingType& vm) const
{
  if(m_mcs.empty() && m_equivalentMcs.empty()) {
    return true;
  }
  double max_weight{0.0};
  MappingIndex max_size{0};
  if(m_nodeMappings.empty()
     || (m_config.isDisconnected()
         && m_connectedComponentSizes.size() - static_cast<MappingIndex>(m_connectedComponentSizes.front() == 0)
            < m_config.nofComponents()
         && m_connectedComponentSizes.front() >= m_config.componentSize())) {
    auto extendableNodes = vm.getNofExtendableNodes();
    max_size = m_nodeMappings.size() + getComponentSize(extendableNodes);
    max_weight = !weighted ? max_size :  currentSize() + getComponentWeight(extendableNodes);
  } else {
    auto extendableNodes = vm.getNofConnectedExtendableNodes();
    max_size = m_nodeMappings.size() + getComponentSize(extendableNodes);
    max_weight = !weighted ? max_size :  currentSize() + getComponentWeight(extendableNodes);
  }
  double bestWeight = 0.0;
  if(!m_mcs.empty()) {
    bestWeight = m_mcs.front().weight();
  } else if(!m_equivalentMcs.empty()) {
    bestWeight = m_equivalentMcs.begin()->weight();
  }
  return m_config.getMinimumResultSize() <= max_size
         && (!weighted || m_config.getMinimumResultWeight() <= max_weight)
         && (m_config.resultType() == Config::ResultType::Maximum
             ? max_weight > bestWeight
             : max_weight >= bestWeight);
}


template<bool weighted>
bool SubgraphMapping<weighted>::addResult()
{
  RIMACS_ASSERT(m_mcs.empty() || m_config.isDisconnected()
                || weighted || m_nodeMappings.size() >= m_mcs.front().size());
  if(!m_currentExpansionIsValid
     || (m_connectedComponentSizes.size() > 1
         && m_connectedComponentSizes.front() < m_config.componentSize())
     || m_connectedComponentSizes.front() == 0) {
    return false;
  }
  bool handleEquivalentResults = m_queryEquivalenceClasses && m_targetEquivalenceClasses;

  double realWeight
      = getNodeComp().estimateFinalWeight(m_nodeMappings.begin(), m_nodeMappings.end(),
                                          currentSize());
  RIMACS_ASSERT(realWeight <= currentSize()
                && "Final weight must not be greater than the sum of the node mapping weights");
  if(m_nodeMappings.size() < m_config.getMinimumResultSize()
     || (weighted && realWeight < m_config.getMinimumResultWeight())) {
    return true;
  }

  if(m_mcs.empty() && m_equivalentMcs.empty()) {
  } else {
    const MCSResult& best = handleEquivalentResults ? *m_equivalentMcs.begin()
                                                    : m_mcs.front();
    if(m_config.resultType() == Config::ResultType::Maximum
       && (best.weight() < realWeight
           || (best.weight() == realWeight &&
               m_connectedComponentSizes.size() < best.getComponentSizes().size()))) {
      m_mcs.clear();
      m_equivalentMcs.clear();
    } else if(m_config.resultType() == Config::ResultType::AllMaximum
              && best.weight() <= realWeight) {

      if((handleEquivalentResults ? m_equivalentMcs.begin()->weight()
                                  : m_mcs.front().weight()) < realWeight) {
        m_mcs.clear();
        m_equivalentMcs.clear();
        m_resultCount = 0;
      }
    } else {
      return true;
    }
  }
  if(handleEquivalentResults) {
    const auto& compare = m_equivalentMcs.key_comp();
    MCSResult res(m_nodeMappings, m_connectedComponentSizes, realWeight,
                  m_initialMappingSize, compare.lessCompareCast());
    auto iteratorPair = m_equivalentMcs.equal_range(res);
    ptrdiff_t dist = std::distance(iteratorPair.first, iteratorPair.second);
    ++m_resultCount;
    if(dist < m_config.equivalentResults()) {
      m_equivalentMcs.insert(std::move(res));
    }
    if(LOGGING_IS_ENABLED(m_info, Base::LogLevel::Debug) && m_resultCount % 100 == 0) {
      LOG_DEBUG(m_info, "Enumerated " << m_resultCount << "equivalent results.");
    }
  } else {
    m_mcs.emplace_back(m_nodeMappings, m_connectedComponentSizes,
                       realWeight, m_initialMappingSize);
  }
  return true;
}

inline void recount_component_sizes(
    std::vector<SubgraphMapping<true>::ConnectionSize>& connectedSizes,
    const std::vector<ComponentIndex>& componentIndices,
    const VertexMapping<true>& vm,
    bool isDisconnected)
{
  for(auto it = connectedSizes.begin() + 1, last = connectedSizes.end(); it != last; ++it) {
    if(it->first != 0) {
      it->second = 0.0;
    }
  }
  for(MappingIndex i = 0, last = componentIndices.size(); i < last; ++i) {
    const ComponentIndex& component = componentIndices[i];
    SubgraphMapping<true>::ConnectionSize& cs = connectedSizes.at(std::abs(component) - 1);
    if(component == SubgraphMapping<true>::InitiallyUnmappableNode) {
      continue; // -1 is the marker for unmappable nodes
    }
    const auto& compNodes = vm.getCompatibleNodes().at(i);
    RIMACS_ASSERT(!compNodes.empty() || (!isDisconnected && cs.first < 0) || cs.first == 0);
    if(cs.first > 0) {
      cs.second -= compNodes.empty() ? 0.0 : compNodes.front().second;
    } else if(cs.first < 0 && isDisconnected) {
      cs.second += compNodes.empty() ? 0.0 : compNodes.front().second;
    }
  }
}

inline void recount_component_sizes(
    std::vector<SubgraphMapping<false>::ConnectionSize>& /*connectedSizes*/,
    const std::vector<ComponentIndex>& /*componentIndices*/,
    const VertexMapping<false>& /*vm*/,
    bool /*isDisconnected*/){}

inline void increase_weight(
    std::vector<double>& weights,
    const MappingTuple& t)
{
  // Weight is stored at last position
  weights.push_back(weights.back() + t.m_weight);
}

inline void increase_weight(
    std::vector<double>&,
    const MappingPair&)
{}


template<bool weighted>
void SubgraphMapping<weighted>::recountDisconnectedComponents(
    const VertexMappingType& vm)
{
  if(weighted && m_use_cc_feature) {
    recount_component_sizes(m_connectedSizes, m_connectedComponentIds, vm, m_config.isDisconnected());
  }
}


template<bool weighted>
typename SubgraphMapping<weighted>::VertexMappingType SubgraphMapping<weighted>::addMatchUpdateVM(
    const VertexMappingType& vmm,
    const typename VertexMappingType::RealMappingTuple& nodeMappings,
    std::vector<MappingIndex>& newUnmappableNodes)
{
  m_nodeStatus.at(nodeMappings.m_from) = m_nodeMappings.size();
  m_nodeMappings.push_back(nodeMappings);
  ++m_counter;
  ++m_connectedComponentSizes.front();
  VertexMappingType result = VertexMappingType(vmm, m_query, m_target,
                                               m_compatibility.second, nodeMappings,
                                               m_nodeStatus, newUnmappableNodes);
  for(const MappingIndex& unmappableNode : newUnmappableNodes | boost::adaptors::reversed) {
    updateConnectedComponents(unmappableNode, result);
  }

  if(m_connectedComponentSizes.front() == 1) {
    result.reinitializeAdjacentNodes(m_query, nodeMappings, m_nodeStatus);
    identifyDisconnectedComponents();
  }
  if(m_use_cc_feature) {
    if(weighted) {
      recount_component_sizes(m_connectedSizes, m_connectedComponentIds, result,
                              m_config.isDisconnected());
    }
    result.setExtendableNodes(countExtendableNodes());
  }
  increase_weight(m_currentWeight, nodeMappings);
  m_currentExpansionIsValid = true;
  return result;
}


template<>
inline typename SubgraphMapping<true>::VertexMappingType SubgraphMapping<true>::handleInitialMapping(
    const VertexMappingType& vm,
    const MappingPair& initialMapping)
{
  ++m_initialMappingSize;
  std::vector<MappingIndex> locallyUnmappableNodes;
  double weight = getNodeComp().nodesAreCompatible(initialMapping)
                ? getNodeComp().getWeight(initialMapping)
                : 1.0;
  return addMatchUpdateVM(vm, MappingTuple(initialMapping.m_from, initialMapping.m_to, weight),
                          locallyUnmappableNodes);
}

template<>
inline typename SubgraphMapping<false>::VertexMappingType SubgraphMapping<false>::handleInitialMapping(
    const VertexMappingType& vm,
    const MappingPair& initialMapping)
{
  ++m_initialMappingSize;
  std::vector<MappingIndex> locallyUnmappableNodes;
  return addMatchUpdateVM(vm, initialMapping, locallyUnmappableNodes);
}


namespace {

template<bool weighted>
struct connected_component_size_filter {

  bool operator()(MappingIndex node) const
  {
    return getComponentSize(m_componentSizes[std::abs(m_connectedComponentIds[node]) - 1]) == 0;
  }

  const std::vector<ComponentIndex>& m_connectedComponentIds;
  const std::vector<typename SubgraphMapping<weighted>::ConnectionSize>& m_componentSizes;
};

} // namespace


template<bool weighted>
bool SubgraphMapping<weighted>::prepareDisconnectionUpdateVM(
    VertexMappingType& vm)
{
  m_connectedComponentSizes.insert(m_connectedComponentSizes.begin(), 0);
  vm.prepareDisconnected(m_query);
  vm.filterAdjacentNodes(connected_component_size_filter<weighted>{m_connectedComponentIds,
                                                                   m_connectedSizes});
  return getComponentSize(vm.getNofExtendableNodes()) >= m_config.componentSize();
}


template<bool weighted>
void SubgraphMapping<weighted>::resetDisconnection()
{
  RIMACS_ASSERT(m_connectedComponentSizes.front() == 0);
  m_connectedComponentSizes.erase(m_connectedComponentSizes.begin());
  auto setToNegativeValue = [](ComponentIndex& l) {l = -std::abs(l);};
  auto setComponentToNegativeValue = [](ConnectionSize& l) {
    getComponentSize(l) = -std::abs(getComponentSize(l));
    if(weighted) {
      getComponentWeight(l) = -std::fabs(getComponentWeight(l));
    }};
  std::for_each(m_connectedComponentIds.begin(), m_connectedComponentIds.end(), setToNegativeValue);
  std::for_each(m_connectedSizes.begin(),
                m_connectedSizes.end(), setComponentToNegativeValue);
}


namespace {

struct unmappable_node_policy {

  MappingIndex getNextComponentId(
      MappingIndex) const
  {
    return ++m_maxComponentIndex;
  }

  bool hasValidComponentId(
      MappingIndex neighbour,
      ComponentIndex /*nextComponentId*/) const
  {
    return m_componentIds[neighbour] == m_componentId;
  }

  bool neighbourHasValidComponentId(
      MappingIndex neighbour) const
  {
    return m_componentIds[neighbour] == m_componentId;
  }

  constexpr bool push_back_node() const
  {
    return true;
  }


  Index& m_maxComponentIndex;
  std::vector<ComponentIndex>& m_componentIds;
  ComponentIndex m_componentId;
};


struct mark_node_unmappable_policy {

  ComponentIndex getNextComponentId(
      Index neighbour) const
  {
    return m_componentIds[neighbour];
  }

  bool hasValidComponentId(
      Index neighbour,
      ComponentIndex nextComponentId) const
  {
    return m_componentIds[neighbour] == nextComponentId;
  }

  bool neighbourHasValidComponentId(
      MappingIndex neighbour) const
  {
    return m_componentId < m_componentIds[neighbour];
  }

  constexpr bool push_back_node() const
  {
    return false;
  }

  std::vector<ComponentIndex>& m_componentIds;
  ComponentIndex m_componentId;
};


} // namespace


template<bool weighted>
template<typename ConnectedComponentDFSPolicy>
inline void  SubgraphMapping<weighted>::run_dfs(
    const AdjacencyFunctor& graph,
    std::vector<MappingIndex>& nodeStack,
    VertexMapping<weighted>& vm,
    std::vector<bool>& visited,
    const ConnectedComponentDFSPolicy& policy,
    ConnectionSize& componentSize,
    const MappingIndex& nextComponentId,
    bool& hasConnectionToMappedNodes)
{
  MappingIndex node = nodeStack.back();
  nodeStack.pop_back();
  m_connectedComponentIds.at(node) = nextComponentId;
  if(weighted) {
    getComponentWeight(componentSize)
        -= getComponentWeight(vm.getCompatibleNodes().at(node).empty()
                              ? std::pair<MappingIndex, double>()
                              : vm.getCompatibleNodes().at(node).front());
  }
  ++getComponentSize(componentSize);
  for(Index neighbour : graph.nodeGetNeighbours(node)) {
    if(!visited.at(neighbour) &&
       policy.hasValidComponentId(neighbour, nextComponentId)) {
      nodeStack.push_back(neighbour);
      visited.at(neighbour) = true;
    } else if(m_nodeStatus.at(neighbour)
              <= m_nodeMappings.size()) {
      hasConnectionToMappedNodes = true;
    }
  }
}


template<bool weighted>
inline void SubgraphMapping<weighted>::invalidate_connected_component(
    VertexMapping<weighted>& vm,
    const ComponentIndex& nextComponentId)
{
  // Invalidate components
  for(MappingIndex i = 0, last = m_connectedComponentIds.size(); i < last; ++i) {
    if(m_connectedComponentIds.at(i) == nextComponentId) {
      m_connectedComponentIds.at(i) = -nextComponentId;
      m_nodeStatus.at(i) = static_cast<MappingIndex>(NodeMappingStatus::Unmappable);
      // In case of disconnected mappings, removing a node can decrease the size of an disconnected
      // component below the threshold. Then it is considered unmappable. If this component
      // is the currently extended one, then there are nodes within this component, which
      // are considered adjacent for mapping extension. This has to be avoided.
      vm.setNodeUnmappable(i, m_connectedComponentSizes.front() == 0);
    }
  }
  getComponentSize(m_connectedSizes.at(nextComponentId - 1)) = 0;
  if(weighted) {
    getComponentWeight(m_connectedSizes, nextComponentId - 1)
        = -getComponentWeight(m_connectedSizes, nextComponentId - 1);
  }
}


template<bool weighted>
inline void SubgraphMapping<weighted>::update_connected_component_size(
    const ConnectionSize& componentSize,
    const MappingIndex& adjacentNode,
    const ComponentIndex& nextComponentId)
{
  getComponentSize(m_connectedSizes, nextComponentId - 1) = -getComponentSize(componentSize);
  if(weighted) {
    getComponentWeight(m_connectedSizes, nextComponentId - 1)
        = -getComponentWeight(m_connectedSizes, nextComponentId - 1);
  }
  // Ensure, that components are correctly marked as disconnected, if this is a disconnected expansion
  if(getComponentSize(componentSize) == 1) {
    RIMACS_ASSERT(m_connectedComponentIds.at(adjacentNode) > 0);
    m_connectedComponentIds.at(adjacentNode) = -std::abs(m_connectedComponentIds.at(adjacentNode));
  } else {
    for(MappingIndex i = 0, last = m_connectedComponentIds.size(); i < last; ++i) {
      if(m_connectedComponentIds.at(i) == nextComponentId) {
        m_connectedComponentIds.at(i) = -nextComponentId;
      }
    }
  }
}


template<bool weighted>
template<typename ConnectedComponentDFSPolicy>
inline void SubgraphMapping<weighted>::run_connected_component_dfs_for_node(
    const AdjacencyFunctor& graph,
    VertexMapping<weighted>& vm,
    std::vector<MappingIndex>& nodeStack,
    std::vector<bool>& visited,
    const ConnectedComponentDFSPolicy& policy,
    MappingIndex adjacentNode,
    ComponentIndex componentId)
{
  // Ensure, that no unvisited neighbour has a higher component id, which might affect the restored state
  RIMACS_ASSERT(visited.at(adjacentNode)
                || std::abs(m_connectedComponentIds.at(adjacentNode)) <= std::abs(componentId)
                || componentId == InitialComponentId
                                /* special value 1 as it is the initial component id for all nodes.
                                 * For nodes that are unmappable unmappable at all, the assumption
                                 * that no unvisited neighbour has a higher componentId might be
                                 * violated.
                                 * Those nodes are assigned -1 as componentId
                                 * This doesn't matter since they are never restored. */
             || (std::is_same<ConnectedComponentDFSPolicy, mark_node_unmappable_policy>::value));
  if(visited.at(adjacentNode) || !policy.neighbourHasValidComponentId(adjacentNode)) {
    return;
  }
  nodeStack.push_back(adjacentNode);
  visited.at(adjacentNode) = true;
  ConnectionSize componentSize{};
  bool hasConnectionToMappedNodes = false;
  ComponentIndex nextComponentId = policy.getNextComponentId(adjacentNode);
  while(!nodeStack.empty()) {
    run_dfs(graph, nodeStack, vm, visited, policy, componentSize,
            nextComponentId, hasConnectionToMappedNodes);
  }
  if(policy.push_back_node()) {
    m_connectedSizes.push_back(componentSize);
  }
  if(!hasConnectionToMappedNodes) {
    // Here we have to ensure that we do not invalidate disconnected components directly after
    // a disconnected expansion. In this case, there is no connected component, thus the
    // number of connected components criterion needs to be relaxed.
    if(m_connectedComponentSizes.size() - static_cast<MappingIndex>(m_connectedComponentSizes.front() == 0)
       >= m_config.nofComponents()
       || (getComponentSize(m_connectedSizes.at(nextComponentId - 1))
           < static_cast<ComponentIndex>(m_config.componentSize())
           && m_connectedComponentSizes.front() + m_connectedComponentSizes.size() > 1
           /* prevent failures ignoring a too small component,
            * if there is no mapping of sufficient size*/)) {
      invalidate_connected_component(vm, nextComponentId);
    } else { // only if it is no new disconnected extension
      update_connected_component_size(componentSize, adjacentNode, nextComponentId);
    }
  }
}

template<bool weighted>
template<typename ConnectedComponentDFSPolicy>
void SubgraphMapping<weighted>::connected_component_dfs(
    ConnectedComponentDFSPolicy& policy,
    VertexMappingType& vm,
    MappingIndex unmappableNode)
{
  const AdjacencyFunctor& graph = m_query;
  std::vector<bool> visited(graph.getNofNodes(), false);
  visited[unmappableNode] = true;
  const ComponentIndex componentId = std::abs(m_connectedComponentIds.at(unmappableNode));
  std::vector<MappingIndex> nodeStack;
  m_connectedComponentIds.at(unmappableNode) = -componentId;
  getComponentSize(m_connectedSizes, componentId - 1) = 0;
  if(weighted) {
    getComponentWeight(m_connectedSizes, componentId - 1)
        = -getComponentWeight(m_connectedSizes, componentId - 1);
  }
  for(MappingIndex adjacentNode : graph.nodeGetNeighbours(unmappableNode)) {
    run_connected_component_dfs_for_node(graph, vm, nodeStack, visited,
                                         policy, adjacentNode, componentId);
  }
}


template<bool weighted>
void SubgraphMapping<weighted>::updateConnectedComponents(
    MappingIndex unmappableNode,
    VertexMappingType& vm,
    bool isMappedNode)
{
  if(!m_use_cc_feature) {
    if(m_nodeStatus.at(unmappableNode) == static_cast<MappingIndex>(NodeMappingStatus::Unmappable)) {
      m_nodeStatus.at(unmappableNode) = static_cast<MappingIndex>(NodeMappingStatus::Mappable);
    }
    // Connected component feature is disabled.
    // -> return without doing anything
    return;
  }
  if(isMappedNode) {
    m_nodeStatus[unmappableNode] = 0;
  }
  unmappable_node_policy policy{m_maxConnectedComponentIndex,
                                m_connectedComponentIds,
                                m_connectedComponentIds[unmappableNode]};
  // handle the special case for initial disconnected component detection
  // Those nodes have an initial component id of 1  an all adjacent nodes initially are in component 2
  policy.m_componentId += static_cast<ComponentIndex>(policy.m_componentId == 1);
  connected_component_dfs(policy, vm, unmappableNode);
}


template<bool weighted>
void SubgraphMapping<weighted>::setNodeUnmappable(
    MappingIndex idx,
    VertexMappingType& vm)
{
  vm.finishedNode(idx);
  m_nodeStatus.at(idx) = static_cast<MappingIndex>(NodeMappingStatus::Unmappable);
  if(m_use_cc_feature) {
    mark_node_unmappable_policy policy{m_connectedComponentIds,
                                       std::abs(m_connectedComponentIds[idx])};
    RIMACS_ASSERT(m_connectedComponentIds.at(idx) < 0);
    connected_component_dfs(policy, vm, idx);
  }
}


template<bool weighted>
void SubgraphMapping<weighted>::restoreConnectedComponents(
    MappingIndex unmappableNode)
{
  if(!m_use_cc_feature) {
    return;
  }
  RIMACS_ASSERT(m_connectedComponentIds.at(unmappableNode) < 0);
  std::vector<MappingIndex>& nodeStatus = m_nodeStatus ;
  const AdjacencyFunctor& graph = m_query;
  const ComponentIndex componentId = -m_connectedComponentIds.at(unmappableNode);
  bool isAdjacentToMapping = false;
  for(MappingIndex neighbour : graph.nodeGetNeighbours(unmappableNode)) {
    if(m_connectedComponentIds.at(neighbour) > componentId) {
      getComponentSize(m_connectedSizes, m_connectedComponentIds.at(neighbour) - 1)
          = std::numeric_limits<ComponentIndex>::max();
      isAdjacentToMapping = true;
    } else if(m_connectedComponentIds.at(neighbour) < -componentId) {
      getComponentSize(m_connectedSizes, -m_connectedComponentIds.at(neighbour) - 1)
          = std::numeric_limits<ComponentIndex>::max();
    } else {
      isAdjacentToMapping |= nodeStatus[neighbour] < m_nodeMappings.size();
    }
  }
  ComponentIndex restoredComponentId = isAdjacentToMapping ? componentId : -componentId;
  ComponentIndex& componentSize = getComponentSize(m_connectedSizes, componentId - 1);
  for(auto it = m_connectedComponentIds.begin(), last = m_connectedComponentIds.end()
      ; it != last
      ; ++it) {
    if(*it == -componentId
       || (getComponentSize(m_connectedSizes, std::abs(*it) - 1)
           == std::numeric_limits<ComponentIndex>::max())) {
      *it = restoredComponentId;
      ++componentSize;
      nodeStatus.at(std::distance(m_connectedComponentIds.begin(), it))
          = static_cast<MappingIndex>(NodeMappingStatus::Mappable);
    }
  }
  if(!isAdjacentToMapping && componentSize < static_cast<ComponentIndex>(m_config.componentSize())
     && (m_nodeMappings.size() >= m_config.componentSize()
         || (!m_mcs.empty() && m_mcs.front().size() > m_config.componentSize())
         || (!m_equivalentMcs.empty() && m_equivalentMcs.begin()->size()
             > m_config.componentSize()))) {
    componentSize = 0;
  } else if(restoredComponentId == -componentId) {
    componentSize = -componentSize;
  }
  while(getComponentSize(m_connectedSizes.back()) == std::numeric_limits<ComponentIndex>::max()) {
    m_connectedSizes.pop_back();
    --m_maxConnectedComponentIndex;
  }
  RIMACS_ASSERT(std::find_if(m_connectedSizes.begin(),
                             m_connectedSizes.end(), [](const auto& cs) {
                return getComponentSize(cs) == std::numeric_limits<ComponentIndex>::max();
              })
                == m_connectedSizes.end());
  ComponentIndex maxElem = *std::max_element(m_connectedComponentIds.begin(), m_connectedComponentIds.end());
  ComponentIndex minElem = *std::min_element(m_connectedComponentIds.begin(), m_connectedComponentIds.end());
  RIMACS_ASSERT(std::max(maxElem, std::abs(minElem))
                == static_cast<ComponentIndex>(m_connectedSizes.size()));
}


namespace {

template<typename T, typename V>
void wrap_emplace_front(
    std::vector<T>& vec,
    T&& t,
    V&&)
{
  vec.emplace(vec.begin(), std::forward<T>(t));
}
template<typename T, typename V>
void wrap_emplace_front(
    std::vector<std::pair<T, V>>& vec,
    T&& t,
    V&& v)
{
  vec.emplace(vec.begin(), std::forward<T>(t), std::forward<V>(v));
}

} // namespace


template<bool weighted>
void SubgraphMapping<weighted>::markUnmappableNodes(
    const VertexMappingType& vm)
{
  RIMACS_ASSERT(m_maxConnectedComponentIndex == MappingIndex{1});
  if(!m_use_cc_feature) {
    return;
  }
  ++m_maxConnectedComponentIndex;
  std::vector<MappingIndex> unmappableNodes;
  for(MappingIndex i = 0, last = vm.getCompatibleNodes().size(); i < last; ++i) {
    if(vm.getCompatibleNodes()[i].empty()) {
      unmappableNodes.push_back(i);
      m_connectedComponentIds.at(i) = 1;
    } else {
      m_connectedComponentIds.at(i) = 2;
    }
  }
  getComponentSize(m_connectedSizes, 0) -= unmappableNodes.size();
  wrap_emplace_front(m_connectedSizes, static_cast<ComponentIndex>(unmappableNodes.size()), 0.0);
  for(const MappingIndex& unmappableNode : unmappableNodes) {
    // since the node is already unmappable for the VertexMapping, it won't be modified
    // The update of connected components can handle the case, that those unmappable nodes are
    // assumed to be part,or no part of the actual mapping.
    updateConnectedComponents(unmappableNode, const_cast<VertexMappingType&>(vm), false);
    m_nodeStatus[unmappableNode] = static_cast<MappingIndex>(NodeMappingStatus::Unmappable);
  }
  getComponentSize(m_connectedSizes, 0) = -unmappableNodes.size();
  if(unmappableNodes.empty()) {
    if(weighted) {
      getComponentWeight(m_connectedSizes.back()) = getComponentWeight(vm.getNofExtendableNodes());
    }
  } else {
    // hack to remove the "virtual" component 2 which was used for the mappable nodes
    getComponentSize(m_connectedSizes, 1) = 0;
    getComponentWeight(m_connectedSizes, 1) = 0.0;
  }
}


template<bool weighted>
void SubgraphMapping<weighted>::identifyDisconnectedComponents()
{
  RIMACS_ASSERT(m_connectedComponentSizes.front() == 1);
  if(getComponentSize(m_connectedSizes, 0) == 0) {
    // There are no previously unmappable nodes, thus there are no disconnected components.
    return;
  }
  MappingIndex currentNode = m_nodeMappings.back().m_from;
  SortedVector<MappingIndex> adjacentComponents;
  for(const MappingIndex& neighbour : m_query.get().nodeGetNeighbours(currentNode)) {
    ComponentIndex id = m_connectedComponentIds.at(neighbour);
    if(id > 0) {
      adjacentComponents.insert_unique(id);
    }
  }
  SortedVector<MappingIndex> disconnectedComponentIds;
  for(ComponentIndex& id : m_connectedComponentIds) {
    if(id > 0 && !adjacentComponents.contains(id)) {
      disconnectedComponentIds.insert_unique(id);
      id = -id;
    }
  }
  for(const MappingIndex& id : disconnectedComponentIds) {
    getComponentWeight(m_connectedSizes, id - 1) = -getComponentWeight(m_connectedSizes, id - 1);
    if(weighted) {
      getComponentWeight(m_connectedSizes, id - 1) = -getComponentWeight(m_connectedSizes, id - 1);
    }
  }
}


} // namespace RIMACS

