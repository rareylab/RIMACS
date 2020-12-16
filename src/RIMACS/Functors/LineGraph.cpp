
#include <algorithm>
#include <unordered_map>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <RIMACS/ASSERT.hpp>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Config.hpp>
#include <RIMACS/Functors/SwappedMapping.hpp>
#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Algorithm/CompareUtils.hpp>
#include <RIMACS/Utils/SortedVector.hpp>

#include <RIMACS/Functors/LineGraph.hpp>

namespace RIMACS {

// At least one method of the functor should be implemented outside the header file
LineGraphEdgeCompatibilityFunctor::~LineGraphEdgeCompatibilityFunctor() = default;
LineGraph::~LineGraph() = default;

bool BasicLineGraph::adjacent(
    MappingIndex fromNode,
    MappingIndex toNode) const
{
  auto it = find_edge(fromNode, toNode);
  return it != m_edges.end() && it->m_from == fromNode && it->m_to == toNode;
}


MappingIndex BasicLineGraph::getEdgeId(
    MappingIndex nodeFromId,
    MappingIndex nodeToId) const
{
  auto it = find_edge(nodeFromId, nodeToId);
  return it != m_edges.end() && it->m_from == nodeFromId && it->m_to == nodeToId
         ? static_cast<MappingIndex>(std::distance(m_edges.begin(), it))
         : NO_EDGE;
}


const MappingIndex& BasicLineGraph::getMappedNode(
    MappingIndex nodeFromId,
    MappingIndex nodeToId) const
{
  return m_edges.at(getEdgeId(nodeFromId, nodeToId)).m_commonNode;
}


const MappingIndex& BasicLineGraph::getLineNodeFromNode(
    MappingIndex node) const
{
  return m_nodes.at(node).m_fromNode;
}


const MappingIndex& BasicLineGraph::getLineNodeToNode(
    MappingIndex node) const
{
  return m_nodes.at(node).m_toNode;
}


MappingIndex BasicLineGraph::edgeGetFromId(
    MappingIndex edgeId) const
{
  return m_edges.at(edgeId).m_from;
}


MappingIndex BasicLineGraph::edgeGetToId(
    MappingIndex edgeId) const
{
  return m_edges.at(edgeId).m_to;
}


MappingIndex BasicLineGraph::nodeGetNofEdges(
    MappingIndex nodeIdx) const
{
  return m_nodeEdges.at(nodeIdx).size();
 }


std::vector<MappingIndex> BasicLineGraph::nodeGetEdges(
    MappingIndex nodeIdx) const
{
  std::vector<MappingIndex> result;
  result.reserve(m_nodeEdges.at(nodeIdx).size());
  std::transform(m_nodeEdges[nodeIdx].begin(), m_nodeEdges[nodeIdx].end(),
                 std::back_inserter(result), [] (const node_edge_info& e) { return e.m_edgeId; });
  return result;
}


std::vector<MappingIndex> BasicLineGraph::nodeGetNeighbours(
    MappingIndex nodeIdx) const
{
  std::vector<MappingIndex> result;
  result.reserve(m_nodeEdges.at(nodeIdx).size());
  std::transform(m_nodeEdges[nodeIdx].begin(), m_nodeEdges[nodeIdx].end(),
                 std::back_inserter(result), [] (const node_edge_info& e) { return e.m_to; });
  return result;
}


MappingIndex BasicLineGraph::getNofNodes() const
{
  return m_nodes.size();
}


MappingIndex BasicLineGraph::getNofEdges() const
{
  return m_edges.size();
}


inline typename std::vector<BasicLineGraph::line_graph_edge>::const_iterator BasicLineGraph::find_edge(
    MappingIndex& nodeFrom,
    MappingIndex& nodeTo) const
{
  if(nodeTo < nodeFrom) {
    std::swap(nodeFrom, nodeTo);
  }
  auto it = std::lower_bound(m_edges.begin(), m_edges.end(), MappingPair(nodeFrom, nodeTo),
                             [] (const line_graph_edge& edge,
                                 const MappingPair& searched) {
      return searched.m_from != edge.m_from ? edge.m_from < searched.m_from
                                            : edge.m_to < searched.m_to;
    });
  return it;
}


void BasicLineGraph::init(
    const AdjacencyFunctor& baseFunctor)
{
  m_nodes.reserve(baseFunctor.getNofEdges());
  for(MappingIndex i = 0, last = baseFunctor.getNofEdges(); i < last; ++i) {
    MappingIndex from = baseFunctor.edgeGetFromId(i);
    MappingIndex to = baseFunctor.edgeGetToId(i);
    if(to < from) {
      std::swap(from, to);
    }
    m_nodes.emplace_back(i, from, to);
  }
  m_nodeEdges.resize(m_nodes.size());
  std::vector<MappingIndex> edgeCounts(m_nodes.size(), 0);
  for(MappingIndex i = 0, nodeEnd = baseFunctor.getNofNodes(); i < nodeEnd; ++i) {
    std::vector<MappingIndex> edges = baseFunctor.nodeGetEdges(i);
    for(MappingIndex edgeI = 0, last = edges.size(), lastIdx = last ? last - 1 : 0; edgeI < lastIdx; ++edgeI) {
      for(MappingIndex edgeJ = edgeI + 1; edgeJ < last; ++edgeJ) {
        m_edges.emplace_back(std::min(edges[edgeI], edges[edgeJ]),
                             std::max(edges[edgeI], edges[edgeJ]), i);
        ++edgeCounts.at(edgeI);
        ++edgeCounts.at(edgeJ);
      }
    }
  }
  std::sort(m_edges.begin(), m_edges.end());
  for(MappingIndex i = 0, last = m_nodes.size(); i < last; ++i) {
    m_nodeEdges[i].reserve(edgeCounts[i]);
  }
  for(MappingIndex i = 0, last = m_edges.size(); i < last; ++i) {
    const line_graph_edge& e = m_edges[i];
    m_nodeEdges.at(e.m_from).emplace_back(i, e.m_to);
    m_nodeEdges.at(e.m_to).emplace_back(i, e.m_from);
  }
}


bool LineGraphNodeCompatibilityFunctor::mappingIsValid(
    typename std::vector<MappingPair>::const_iterator beginMapping,
    typename std::vector<MappingPair>::const_iterator endMapping) const
{
  return isValidImpl(beginMapping, endMapping);
}


bool LineGraphNodeCompatibilityFunctor::mappingIsValid(
    SwappedMappingIterator beginMapping,
    SwappedMappingIterator endMapping) const
{
  return isValidImpl(beginMapping, endMapping);
}


namespace {

template<typename Iterator>
bool prepare_mapped_nodes(
    Iterator beginMapping,
    Iterator endMapping,
    const BasicLineGraph& query,
    const BasicLineGraph& target,
    SortedVector<MappingPair>& mappedNodes,
    const SortedVector<MappingIndex>& queryEdges,
    const SortedVector<MappingIndex>& availableTargetNodes)
{
  // The mapping is not necessarily sorted, we always have to sort is for the following index tricks
  SortedVector<MappingPair> sortedMapping;
  sortedMapping.reserve(std::distance(beginMapping, endMapping));
  std::transform(beginMapping, endMapping, std::back_inserter(sortedMapping),
                 [](auto& p) { return MappingPair(p.m_from, p.m_to); });
  for(Iterator begin = beginMapping; begin != endMapping; ++begin) {
    for(MappingIndex neighbour : query.nodeGetNeighbours(begin->m_from)) {
      if(queryEdges.count(neighbour)) {
        MappingIndex mappedNeighbourEdgeId
            = std::lower_bound(sortedMapping.begin(), sortedMapping.end(),
                               MappingPair(neighbour, std::numeric_limits<MappingIndex>::max()))[-1].m_to;
        mappedNodes.emplace_unique(query.getMappedNode(begin->m_from, neighbour),
                                   target.getMappedNode(begin->m_to, mappedNeighbourEdgeId));
        if(mappedNodes.size() > availableTargetNodes.size()) {
          return false;
        }
      }
    }
  }
  return true;
}

bool evaluate_mapped_nodes_and_incident_edges(
    const SortedVector<MappingPair>& mappedNodes,
    SortedVector<MappingIndex>& availableTargetNodes,
    const SortedVector<MappingIndex>& queryNodes,
    const SortedVector<MappingIndex>& targetNodes)
{
  MappingIndex lastNode = std::numeric_limits<MappingIndex>::max();
  for(const MappingPair& mp : mappedNodes) {
    // ensure that each node of both graphs is mapped at most once
    // and count the number of incident nodes. The numbers must be identical
    if(lastNode == mp.m_from || ! availableTargetNodes.count(mp.m_to)
       || queryNodes.contains(mp.m_from) != targetNodes.contains(mp.m_to)) {
      return false;
    }
    availableTargetNodes.erase(mp.m_to);
    lastNode = mp.m_from;
  }
  // The loop above only tests for nodes that are definitely mapped, because at least two
  // incident edges are mapped. Nodes with only one incident edge are not accounted for
  return true;
}

} // namespace

template<typename Iterator>
inline bool LineGraphNodeCompatibilityFunctor::isValidImpl(
    Iterator beginMapping,
    Iterator endMapping) const
{
  // Perform triangle triod (delta Y) exchange test
  switch(std::distance(beginMapping, endMapping)) {
    case(3): {
      SortedVector<MappingIndex> queryNodes(4);
      SortedVector<MappingIndex> targetNodes(4);
      for(; beginMapping != endMapping; ++beginMapping) {
        queryNodes.insert_unique(m_query.m_nodes[beginMapping->m_from].m_fromNode);
        queryNodes.insert_unique(m_query.m_nodes[beginMapping->m_from].m_toNode);
        targetNodes.insert_unique(m_target.m_nodes[beginMapping->m_to].m_fromNode);
        targetNodes.insert_unique(m_target.m_nodes[beginMapping->m_to].m_toNode);
      }
      return queryNodes.size() == targetNodes.size();
    } case(4):
    case(5):
    case(6): {
      // perform symmetry tests for extended triangles, diamonds and tetrahedrons
      // e.g. C1CC1C, C12CC1C2, C12C3C1C23

      /*
       * There are some cases on small graphs, the cases listed above, where the line graph
       * contains more symmetry than the original graph. In such a situation the restored node
       * mapping could be invalid. Therefore we handle some special cases described in
       * Nicholson, V.; Tsai, C.-C.; Johnson, M.; Naim, M. In Graph Theory and Topology in Chemistry.
       * A Subgraph Isomorphism Theorem for Molecular Graphs.; Elsevier: Amsterdam 1987; pp 226–230
       * (Available in Chemistry library on the campus)
       */

      SortedVector<MappingIndex> queryNodes(14);
      SortedVector<MappingIndex> queryEdges(6);
      SortedVector<MappingIndex> targetNodes(14);
      SortedVector<MappingPair> mappedNodes(7);
      // extract the incident nodes for each edge
      // This is used such that nodes that could be mapped based on the current result,
      // have the same number of incident mapped edges
      for(Iterator begin = beginMapping; begin != endMapping; ++begin) {
        queryNodes.insert(m_query.m_nodes[begin->m_from].m_fromNode);
        queryNodes.insert(m_query.m_nodes[begin->m_from].m_toNode);
        queryEdges.insert_unique(begin->m_from);
        targetNodes.insert(m_target.m_nodes[begin->m_to].m_fromNode);
        targetNodes.insert(m_target.m_nodes[begin->m_to].m_toNode);
      }
      // extract unique nodes from target, they are used as upper bound for mapping sizes
      SortedVector<MappingIndex> availableTargetNodes(targetNodes);
      availableTargetNodes.make_unique();
      // ensure that the number of mapped nodes of query is the same as from target
      if(SortedVector<MappingIndex>(queryNodes).make_unique().size()
         != availableTargetNodes.size()) {
        return false;
      } else if(static_cast<std::ptrdiff_t>(availableTargetNodes.size())
                > std::distance(beginMapping, endMapping)) {
        // All of the considered cases have in common, that there are more edges than nodes,
        // the mappings contain cycles
        return true;
      }

      if(!prepare_mapped_nodes(beginMapping, endMapping, m_query, m_target, mappedNodes, queryEdges, availableTargetNodes)) {
        return false;
      }
      if(!evaluate_mapped_nodes_and_incident_edges(mappedNodes, availableTargetNodes, queryNodes, targetNodes)) {
        return false;
      }
    } default:
      return true;
  }
}


namespace {

std::vector<MappingPair> convert_component(
    const std::vector<MappingIndex>& component,
    const BasicLineGraph& query,
    const BasicLineGraph& target,
    const std::unordered_map<MappingIndex, MappingIndex>& mapping)
{
  RIMACS_ASSERT(component.size() > 1);
  std::unordered_map<MappingIndex, MappingIndex> res;
  for(const MappingIndex& node : component) {
    auto anyResultIt = res.end();
    for(const MappingIndex& neighbour : query.nodeGetNeighbours(node)) {
      auto mappingIt = mapping.find(neighbour);
      if(mappingIt == mapping.end()) {
        continue;
      }
      anyResultIt = res.emplace(query.getMappedNode(node, neighbour),
                                target.getMappedNode(mapping.at(node), mapping.at(neighbour))).first;
    }
    RIMACS_ASSERT(anyResultIt != res.end());
    auto mappedIt = res.find(query.getLineNodeFromNode(node));
    if(mappedIt == res.end()) {
      const MappingIndex& otherNode = mapping.at(node);
      const MappingIndex& otherLineFrom = target.getLineNodeFromNode(otherNode);
      const MappingIndex& otherLineTo = target.getLineNodeToNode(otherNode);
      res.emplace(query.getLineNodeFromNode(node),
                  otherLineFrom == anyResultIt->second ? otherLineTo : otherLineFrom);
    }
    mappedIt = res.find(query.getLineNodeToNode(node));
    if(mappedIt == res.end()) {
      const MappingIndex& otherNode = mapping.at(node);
      const MappingIndex& otherLineFrom = target.getLineNodeFromNode(otherNode);
      const MappingIndex& otherLineTo = target.getLineNodeToNode(otherNode);
      res.emplace(query.getLineNodeToNode(node),
                  otherLineFrom == anyResultIt->second ? otherLineTo : otherLineFrom);
    }
  }
  std::vector<MappingPair> result;
  result.reserve(res.size());
  std::transform(res.begin(), res.end(), std::back_inserter(result), [](const auto& p) {
      return MappingPair(p.first, p.second);
    });
  std::sort(result.begin(), result.end());
  return result;
}

std::vector<MappingPair> convert_component(
    const MappingPair& mapping,
    const BasicLineGraph& query,
    const BasicLineGraph& target,
    const NodeCompatibilityFunctor& functor)
{
  std::vector<MappingPair> result;
  result.reserve(2);
  RIMACS_ASSERT(!dynamic_cast<const LineGraphNodeCompatibilityFunctor*>(&functor)
                && "Node compatibility functor must be of the base graphs, not of the line graph");
  const MappingIndex qLineFrom = query.getLineNodeFromNode(mapping.m_from);
  const MappingIndex qLineTo = query.getLineNodeToNode(mapping.m_from);
  const MappingIndex tLineFrom = target.getLineNodeFromNode(mapping.m_to);
  const MappingIndex tLineTo = target.getLineNodeToNode(mapping.m_to);
  if(functor(qLineFrom, tLineFrom) && functor(qLineTo, tLineTo)) {
    result.emplace_back(qLineFrom, tLineFrom);
    result.emplace_back(qLineTo, tLineTo);
  }
  if(functor(qLineFrom, tLineTo) && functor(qLineTo, tLineFrom)) {
    result.emplace_back(qLineFrom, tLineTo);
    result.emplace_back(qLineTo, tLineFrom);
  }
  return result;
}


std::unordered_map<MappingIndex, MappingIndex> to_map(
    const std::vector<MappingPair>& mapping)
{
  std::unordered_map<MappingIndex, MappingIndex> mappingMap;
  mappingMap.reserve(mapping.size());
  for(const MappingPair& mp : mapping) {
    mappingMap.emplace(mp.m_from, mp.m_to);
  }
  return mappingMap;
}

} // namespace


std::vector<std::pair<std::vector<MappingPair>, std::vector<MappingIndex>>> EdgeMappingConverter::operator()(
    const std::vector<MappingPair>& mapping,
    const std::vector<MappingIndex>& cumulatedComponentSizes,
    bool expandLoneBondMappings) const
{
  std::vector<std::pair<std::vector<MappingPair>, std::vector<MappingIndex>>> resultVector;
  std::vector<std::tuple<MappingIndex, MappingPair, MappingPair>> alternativeResults;

  std::vector<MappingPair> result;
  result.reserve(mapping.size() + cumulatedComponentSizes.size());
  std::vector<MappingIndex> cumComponentSizes;
  cumComponentSizes.reserve(cumulatedComponentSizes.size());
  std::unordered_map<MappingIndex, MappingIndex> mappingMap = to_map(mapping);

  MappingIndex first = 0;
  for(auto it = cumulatedComponentSizes.begin(); it != cumulatedComponentSizes.end(); ++it) {
    std::vector<MappingIndex> component;
    component.reserve(*it - first);
    std::transform(std::next(mapping.begin(), first), std::next(mapping.begin(), *it),
                   std::back_inserter(component), mapping_from());
    std::vector<MappingPair> componentMapping;
    if(component.size() == 1) {
      componentMapping = convert_component(mapping[first], m_query, m_target, m_functor);
      // Handle the following case, where a bond in mapped, but the mapping is valid in both directions
      // In this case, the convert component method returns two results of size two. But they are
      // concatenated, thus the result size is four.
      if(componentMapping.size() == 4) {
        // Store the second result variant in the alternative results
        alternativeResults.emplace_back(std::distance(cumulatedComponentSizes.begin(), it),
                                        componentMapping[2], componentMapping[3]);
        // keep the first result variant in the mapping
        componentMapping.erase(componentMapping.begin() + 2, componentMapping.end());
      }
    } else {
      componentMapping = convert_component(component, m_query, m_target, mappingMap);
    }
    result.insert(result.end(), componentMapping.begin(), componentMapping.end());
    cumComponentSizes.push_back(result.size());
    first = *it;
  }
  resultVector.emplace_back(std::move(result), std::move(cumComponentSizes));

  if(expandLoneBondMappings) {
    resultVector.resize(1 << alternativeResults.size(), resultVector.front());
    MappingIndex increment = 1;
    MappingIndex finalSize = resultVector.size();
    for(const auto& alternativeResultTuple : alternativeResults) {
      for(MappingIndex outerIdx = 0; outerIdx < finalSize; outerIdx += increment << 1) {
        for(MappingIndex idx = 0; idx < increment; ++idx) {
          std::pair<std::vector<MappingPair>, std::vector<MappingIndex>>& singleResult
              = resultVector[outerIdx + idx];
          MappingIndex componentBegin = std::get<0>(alternativeResultTuple)
                               ? singleResult.second.at(std::get<0>(alternativeResultTuple) - 1)
                               : 0;
          singleResult.first.at(componentBegin) = std::get<1>(alternativeResultTuple);
          singleResult.first.at(componentBegin + 1) = std::get<2>(alternativeResultTuple);
        }
      }
      increment <<= 1;
    }
  }

  return  resultVector;
}

double PreferMaxNumberOfMappedNodesLineGraphNodeCompatibilityFunctor::estimateFinalWeight(
    typename std::vector<MappingPair>::const_iterator beginMapping,
    typename std::vector<MappingPair>::const_iterator endMapping,
    double currentWeight) const
{
  return estimateFinalWeight_impl(beginMapping, endMapping,
                                  currentWeight);
}

double PreferMaxNumberOfMappedNodesLineGraphNodeCompatibilityFunctor::estimateFinalWeight(
    SwappedMappingIterator beginMapping,
    SwappedMappingIterator endMapping,
    double currentWeight) const
{
  return estimateFinalWeight_impl(beginMapping, endMapping,
                                  currentWeight);
}

template<typename Iterator>
double PreferMaxNumberOfMappedNodesLineGraphNodeCompatibilityFunctor::estimateFinalWeight_impl(
    Iterator beginMapping,
    Iterator endMapping,
    double currentWeight) const
{
  boost::dynamic_bitset<> mappedNodes(this->m_query.getNofNodes() + 1);
  size_t nofMappedEdges = std::distance(beginMapping, endMapping);
  for(; beginMapping != endMapping; ++beginMapping) {
    mappedNodes.set(this->m_query.getLineNodeFromNode(beginMapping->m_from));
    mappedNodes.set(this->m_query.getLineNodeToNode(beginMapping->m_from));
  }
  return currentWeight - calculatePenaltyValue(mappedNodes.count(), nofMappedEdges,
                                               calculateNofComponents(nofMappedEdges));
}

double PreferMaxNumberOfMappedNodesLineGraphNodeCompatibilityFunctor::initializePenaltyValue(
    MappingIndex lineGraphNofNodes)
{
  size_t upper_bound = 1;
  while(upper_bound < lineGraphNofNodes) {
    upper_bound <<= 1;
  }
  return 1.0 / static_cast<double>(upper_bound);
}

inline MappingIndex PreferMaxNumberOfMappedNodesLineGraphNodeCompatibilityFunctor::calculateNofComponents(
    MappingIndex mappedEdges) const
{
  // The upper bound estimation of the number of connected components 'm_query.getNofNodes() - mappedEdges + 1'
  // is only valid, because the graph is a linegraph. On a generic graph, there could be an arbitrary number
  // of connected components, if only a single node is not mapped. (E.g. in a star)
  return std::max(MappingIndex{1},
                  std::min(std::min(m_config.nofComponents(), m_query.getNofNodes() - mappedEdges + 1),
                           mappedEdges / m_config.componentSize()));
}


} // namespace RIMACS
