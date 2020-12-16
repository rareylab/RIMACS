
#pragma once

#include <algorithm>
#include <numeric>
#include <type_traits>
#include <utility>
#include <cmath>

#include <RIMACS/ASSERT.hpp>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Structs.hpp>
#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Algorithm/CompareUtils.hpp>

#include <RIMACS/Algorithm/VertexMapping.hpp>
#include <RIMACS/Algorithm/SubgraphMapping.hpp>
#include <RIMACS/Algorithm/WeightTypeTraits.hpp>

namespace RIMACS {

template<bool weighted>
VertexMapping<weighted>::VertexMapping(
    const VertexMapping& other,
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const EdgeCompatibilityFunctor& edgeComp,
    const MappingPair& newNode,
    std::vector<MappingIndex>& mappingStatus,
    std::vector<MappingIndex>& newUnmappableNodes)
  : m_compatibleNodes(other.m_compatibleNodes.size())
  , m_adjacentNodes(other.m_adjacentNodes)
  , m_hasAdjacentNodes(other.m_hasAdjacentNodes)
{
  doCopyConstruct(other, query, target, edgeComp, newNode, mappingStatus, newUnmappableNodes);
}


template<bool weighted>
void VertexMapping<weighted>::doCopyConstruct(
    const VertexMapping& other,
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const EdgeCompatibilityFunctor& edgeComp,
    const MappingPair& newNode,
    std::vector<MappingIndex>& mappingStatus,
    std::vector<MappingIndex>& newUnmappableNodes)
{
  this->updateCompatibleNodes(other, query, target, edgeComp, mappingStatus, newNode);
  RIMACS_ASSERT(newUnmappableNodes.empty());
  newUnmappableNodes = updateExtendableNodesAndStatus(mappingStatus);
  updateAdjacentNodes(query, newNode, mappingStatus);
}


template<bool weighted>
bool VertexMapping<weighted>::getNextMatch(
    MappingPair& matchIndices,
    RealMappingTuple& mappedNodes)
{
  mappedNodes.m_from = matchIndices.m_from;
  {
    ++matchIndices.m_to;
    if(matchIndices.m_to == m_compatibleNodes.at(mappedNodes.m_from).size()) {
      matchIndices.m_to = 0;
      mappedNodes.m_to = getCompatibleNode(m_compatibleNodes.at(mappedNodes.m_from), 0);
      return false;
    }
  }
  mappedNodes.m_to = getCompatibleNode(m_compatibleNodes.at(mappedNodes.m_from),
                                       matchIndices.m_to);
  if(weighted) {
    set_weight(mappedNodes, -getCompatibleNodeWeight(m_compatibleNodes.at(mappedNodes.m_from),
                                                     matchIndices.m_to));
  }
  return true;
}


template<bool weighted>
MappingIndex VertexMapping<weighted>::getNextNode()
{
  if(m_adjacentNodes.empty()) {
    return NODE_END;
  }
  MappingIndex result = m_adjacentNodes.back().m_node;
  m_adjacentNodes.pop_back();
  return result;
}


template<bool weighted>
void VertexMapping<weighted>::finishedNode(
    MappingIndex node)
{
  // if the node has no mapping partner in any of the graphs, then is wasn't counted as extendable
  if(!m_compatibleNodes.at(node).empty()) {
    if(weighted) {
      getCompatibleNodeWeight(m_nofExtendableNodes)
        += getCompatibleNodeWeight(m_compatibleNodes[node], 0);
    }
    --getCompatibleNode(m_nofExtendableNodes);
  }
  m_compatibleNodes.at(node).clear();
  m_hasAdjacentNodes = !m_adjacentNodes.empty();
}

template<bool weighted>
void VertexMapping<weighted>::setNodeUnmappable(
    MappingIndex node,
    bool ensureNodeIsNotAdjacent)
{
  if(!m_compatibleNodes.at(node).empty()) {
    m_compatibleNodes[node].clear();
    --getCompatibleNode(m_nofExtendableNodes);
    if(ensureNodeIsNotAdjacent) {
      auto it = std::find_if(m_adjacentNodes.begin(), m_adjacentNodes.end(),
                             [node](const auto& n) {
        return n.m_node == node;
      });
      if(it != m_adjacentNodes.end()) {
        m_adjacentNodes.erase(it);
      }
    }
  }
}


template<bool weighted>
void VertexMapping<weighted>::initAdjacentNodes(
    const AdjacencyFunctor& query)
{
  m_adjacentNodes.reserve(query.getNofNodes());
  for(MappingIndex i = 0; i < query.getNofNodes(); ++i) {
    MappingIndex count = m_compatibleNodes.at(i).size();
    if(count) {
      if(weighted) {
        m_adjacentNodes.emplace_back(i, count, 0, getCompatibleNodeWeight(m_compatibleNodes[i], 0),
            -query.getNodeOrderSuggestion(i));
      } else {
        m_adjacentNodes.emplace_back(i, count, 0, -query.getNodeOrderSuggestion(i),
            -getCompatibleNodeWeight(m_compatibleNodes[i], 0));
      }
    }
  }
  std::sort(m_adjacentNodes.begin(), m_adjacentNodes.end(), std::greater<>());
  if(weighted) {
    getCompatibleNodeWeight(m_nofExtendableNodes)
        = -std::accumulate(m_adjacentNodes.begin(), m_adjacentNodes.end(), 0.0,
                           [](const double& currentValue, const AdjacentNode& t) {
            return t.m_weight + currentValue;
          });
  }
  getCompatibleNode(m_nofExtendableNodes) = m_adjacentNodes.size();
  m_hasAdjacentNodes = static_cast<bool>(getCompatibleNode(m_nofExtendableNodes));
}


inline void emplace_sort_element(
    std::vector<std::tuple<double, int, int, double, MappingIndex>>& sortNodeBuffer,
    const double& weight,
    const int& degreeDelta,
    const int& extendedNodeDegreeDelta,
    const double estimatedSimilarity,
    const MappingIndex node)
{
  // Enable disable extended node degree here
  sortNodeBuffer.emplace_back(weight, degreeDelta, extendedNodeDegreeDelta, estimatedSimilarity, node);
}


inline void emplace_sort_element(
    std::vector<std::tuple<int, int, double, MappingIndex>>& sortNodeBuffer,
    const double&,
    const int& degreeDelta,
    const int& extendedNodeDegreeDelta,
    const double& estimatedSimilarity,
    const MappingIndex node)
{
  // Enable disable extended node degree here
  sortNodeBuffer.emplace_back(degreeDelta, extendedNodeDegreeDelta, estimatedSimilarity, node);
}

inline ComponentIndex get_node_neighbour_degree(
    const AdjacencyFunctor& graph,
    MappingIndex node)
{
  ComponentIndex result = 0;
  for(const MappingIndex& neighbour : graph.nodeGetNeighbours(node)) {
    result += graph.nodeGetNofEdges(neighbour);
  }
  return result;
}


template<bool weighted>
void VertexMapping<weighted>::initCompatibleNodes(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp)
{
  using SortNodeBufferType
      = std::vector<std::conditional_t<weighted, std::tuple<double, int, int, double, MappingIndex>,
                                                 std::tuple<int, int, double, MappingIndex>>>;
  std::vector<ComponentIndex> g2NodeNeighbourDegrees;
  g2NodeNeighbourDegrees.reserve(target.getNofNodes());
  for(MappingIndex i = 0; i < target.getNofNodes(); ++i) {
    g2NodeNeighbourDegrees.push_back(get_node_neighbour_degree(target, i));
  }
  std::vector<double> g2NodeRelativePositions;
  g2NodeRelativePositions.reserve(target.getNofNodes());
  for(MappingIndex i = 0; i < target.getNofNodes(); ++i) {
    g2NodeRelativePositions.push_back(1 - static_cast<double>(target.getNodeOrderSuggestion(i))
                                                              / (target.getMaxNodeOrderSuggestion() + 0.1));
  }
  SortNodeBufferType sortNodeBuffer;
  sortNodeBuffer.reserve(query.getNofNodes());
  for(MappingIndex idx1 = 0; idx1 < query.getNofNodes(); ++idx1) {
    m_compatibleNodes.at(idx1).reserve(target.getNofNodes() / 2);
    ComponentIndex g0NodeDegree = query.nodeGetNofEdges(idx1);
    int neighbourDegree = get_node_neighbour_degree(query, idx1);
    double rel_node_pos = 1 - static_cast<double>(query.getNodeOrderSuggestion(idx1))
                                                  / (query.getMaxNodeOrderSuggestion() + 0.1);
    for(MappingIndex idx2 = 0; idx2 < target.getNofNodes(); ++idx2) {
      if(nodeComp.nodesAreCompatible(idx1, idx2)) {
        double sim = nodeComp.estimateNodeSimilarity(idx1, idx2);
        double pos_sim = 1 - std::fabs(rel_node_pos - g2NodeRelativePositions[idx2]);
        if(weighted) {
          emplace_sort_element(sortNodeBuffer,
              -nodeComp.getWeight(idx1, idx2),
              std::abs(g0NodeDegree - static_cast<ComponentIndex>(target.nodeGetNofEdges(idx2))),
              std::abs(neighbourDegree - g2NodeNeighbourDegrees[idx2]), -sim * pos_sim, idx2);
        } else {
          emplace_sort_element(sortNodeBuffer,
              0.0,
              std::abs(g0NodeDegree - static_cast<ComponentIndex>(target.nodeGetNofEdges(idx2))),
              std::abs(neighbourDegree - g2NodeNeighbourDegrees[idx2]), -sim * pos_sim, idx2);
        }
      }
    }
    std::sort(sortNodeBuffer.begin(), sortNodeBuffer.end());
    std::transform(
        sortNodeBuffer.begin(), sortNodeBuffer.end(), std::back_inserter(m_compatibleNodes[idx1]),
        [](const auto& t) {
          constexpr size_t SortNodeTypeNodeIdx = std::tuple_size<typename SortNodeBufferType::value_type>::value - 1;
          constexpr size_t NodeSimilarityEstIdx = SortNodeTypeNodeIdx - 1;
          return std::make_pair(std::get<SortNodeTypeNodeIdx>(t), std::get<weighted ? 0 : NodeSimilarityEstIdx>(t));
          });
    sortNodeBuffer.clear();
  }
}


template<bool weighted>
void VertexMapping<weighted>::initCompatibleNodes(
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const NodeCompatibilityFunctor& nodeComp,
    const std::vector<std::vector<MappingIndex>>& hint)
{
  RIMACS_ASSERT(hint.size() == query.getNofNodes());

  MappingIndex currentNode = 0;
  std::vector<std::pair<MappingIndex, double>> compNodes;
  m_compatibleNodes.reserve(query.getNofNodes());

  auto unweightedInsertionLambda = [&nodeComp,&currentNode,&compNodes](const MappingIndex& hintNode) {
    if(nodeComp.nodesAreCompatible(currentNode, hintNode)) {
      compNodes.emplace_back(hintNode, nodeComp.estimateNodeSimilarity(currentNode, hintNode));
    }
  };
  auto weightedInsertionLambda = [&nodeComp,&currentNode,&compNodes](const MappingIndex& hintNode) {
    if(nodeComp.nodesAreCompatible(currentNode, hintNode)) {
      compNodes.emplace_back(hintNode, -nodeComp.getWeight(currentNode, hintNode));
    }
  };
  for(MappingIndex last = query.getNofNodes(); currentNode < last; ++currentNode) {
    const std::vector<MappingIndex>& currentHint = hint[currentNode];
    compNodes.clear();
    compNodes.reserve(currentHint.size());
    if(weighted) {
      std::for_each(currentHint.begin(), currentHint.end(), weightedInsertionLambda);
    } else {
      std::for_each(currentHint.begin(), currentHint.end(), unweightedInsertionLambda);
    }
    m_compatibleNodes.push_back(std::move(compNodes));
  }
}


template<bool weighted>
std::vector<MappingIndex> VertexMapping<weighted>::updateExtendableNodesAndStatus(
    std::vector<MappingIndex>& status)
{
  std::vector<MappingIndex> newUnmappableNodes;
  m_nofExtendableNodes = CompatibleNodesType{};
  for(MappingIndex i = 0; i < status.size(); ++i) {
    if(status.at(i) == static_cast<MappingIndex>(NodeMappingStatus::Mappable)) {
      if(!m_compatibleNodes.at(i).empty()) {
        ++getCompatibleNode(m_nofExtendableNodes);
        if(weighted) {
          getCompatibleNodeWeight(m_nofExtendableNodes)
              -= getCompatibleNodeWeight(m_compatibleNodes[i], 0);
        }
      } else {
        newUnmappableNodes.push_back(i);
        status.at(i) = static_cast<MappingIndex>(NodeMappingStatus::Unmappable);
      }
    }
  }
  m_hasAdjacentNodes = static_cast<bool>(getCompatibleNode(m_nofExtendableNodes));
  return newUnmappableNodes;
}


template<bool weighted>
void VertexMapping<weighted>::updateAdjacentNodes(
    const AdjacencyFunctor& graph,
    const MappingPair& nodeMappings,
    const std::vector<MappingIndex>& mappingStatus)
{
  for(MappingIndex otherNode : graph.nodeGetNeighbours(nodeMappings.m_from)) {
    // insert adjacent nodes if they are new
    if(mappingStatus.at(otherNode) == static_cast<MappingIndex>(NodeMappingStatus::Mappable)) {
      m_adjacentNodes.emplace_back(MappingIndex{otherNode}, MappingIndex{0} , 1u, double{0} , 0.0);
    }
  }
  for(AdjacentNode& adjacentNode : m_adjacentNodes) {
    adjacentNode.m_nofPartners = m_compatibleNodes.at(adjacentNode.m_node).size();
    if(weighted) {
      adjacentNode.m_weight
          = m_compatibleNodes.at(adjacentNode.m_node).empty()
            ? 0.0
            : getCompatibleNodeWeight(m_compatibleNodes.at(adjacentNode.m_node), 0);
    }
  }
  // Sort solely by node, in order to identify and merge duplicates
  std::sort(m_adjacentNodes.begin(),
            m_adjacentNodes.end(), [] (const AdjacentNode& lhs, const AdjacentNode& rhs) {
      return lhs.m_node < rhs.m_node;
    });
  bool foundRing = false;
  auto adjacentEnd = std::unique(m_adjacentNodes.begin(),
                                 m_adjacentNodes.end(),
                                 [&foundRing](AdjacentNode& lhs, AdjacentNode& rhs) {
      if(lhs.m_node == rhs.m_node && lhs.m_nofPartners == rhs.m_nofPartners) {
        lhs.m_nofMappedNeighbours += rhs.m_nofMappedNeighbours;
        foundRing = true;
        return true;
      }
      return false;
  });
  m_adjacentNodes.erase(std::remove_if(m_adjacentNodes.begin(),
                                       adjacentEnd,
                                       [](const AdjacentNode& n) {
                                         return n.m_nofPartners == 0;
                                       }),
                        m_adjacentNodes.end());
  //  Sort and restore the desired node order
  std::sort(m_adjacentNodes.begin(),
            m_adjacentNodes.end(), std::greater<>());

  m_hasAdjacentNodes = !m_adjacentNodes.empty();
}



template<bool weighted>
void VertexMapping<weighted>::updateCompatibleNodes(
    const VertexMapping& other,
    const AdjacencyFunctor& query,
    const AdjacencyFunctor& target,
    const EdgeCompatibilityFunctor& edgeComp,
    const std::vector<MappingIndex>& status,
    const MappingPair& newNode)
{
  for(MappingIndex nodeIdx = 0, last = m_compatibleNodes.size() ; nodeIdx < last; ++nodeIdx) {
    if(status.at(nodeIdx) == static_cast<MappingIndex>(NodeMappingStatus::Mappable))
    {
      MappingIndex edgeIdx = query.getEdgeId(newNode.m_from, nodeIdx);
      bool adjacentInQuery = edgeIdx != AdjacencyFunctor::NO_EDGE;

      for(const StoredCompatibleNodesType& compatibleNode : other.m_compatibleNodes.at(nodeIdx)) {
        if(newNode.m_to != getCompatibleNode(compatibleNode)
           && target.adjacent(newNode.m_to, getCompatibleNode(compatibleNode)) == adjacentInQuery
           && (!adjacentInQuery
               || edgeComp.edgesAreCompatible(
                      edgeIdx,
                      target.getEdgeId(newNode.m_to, getCompatibleNode(compatibleNode))))) {
          m_compatibleNodes.at(nodeIdx).push_back(compatibleNode);
        }
      }
    }
  }
}


namespace {


template<typename CompatibleNodesType>
static std::enable_if_t<std::is_same<CompatibleNodesType, MappingIndex>::value, const MappingIndex&>
expand_with_weight(const MappingIndex& compNode)
{
  return compNode;
}
template<typename CompatibleNodesType>
static std::enable_if_t<std::is_same<CompatibleNodesType, std::pair<MappingIndex, double>>::value,
                        std::pair<MappingIndex, double>>
expand_with_weight(const MappingIndex& compNode)
{
  return std::make_pair(compNode, 1.0);
}

} // namespace

template<bool weighted>
void VertexMapping<weighted>::insertInitialMapping(
    const MappingPair& initialMapping)
{
  if(m_compatibleNodes.at(initialMapping.m_from).empty()) {
    m_compatibleNodes[initialMapping.m_from].push_back(
        expand_with_weight<StoredCompatibleNodesType>(initialMapping.m_to));
  }
}


} // namespace RIMACS
