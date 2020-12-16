
#pragma once

#include <type_traits>
#include <algorithm>

#include <RIMACS/ASSERT.hpp>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Structs.hpp>
#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Algorithm/CompareUtils.hpp>

namespace RIMACS {


template<bool weighted>
class VertexMapping {
public:

  using CompatibleNodesType = typename std::conditional<weighted, std::pair<MappingIndex, double>,
                                                                  MappingIndex>::type;
  using StoredCompatibleNodesType = std::pair<MappingIndex, double>;
  using RealMappingTuple = typename std::conditional_t<weighted, MappingTuple, MappingPair>;

  enum : MappingIndex { NODE_END = std::numeric_limits<MappingIndex>::max() };

  VertexMapping(
      const VertexMapping& other,
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const EdgeCompatibilityFunctor& edgeComp,
      const MappingPair& newNode,
      std::vector<MappingIndex>& mappingStatus,
      std::vector<MappingIndex>& newUnmappableNodes);

  VertexMapping(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp)
    : m_compatibleNodes(query.getNofNodes())
  {
    initCompatibleNodes(query, target, nodeComp);
    initAdjacentNodes(query);
  }

  VertexMapping(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const MappingPair& initialMapping)
    : VertexMapping(query, target, nodeComp)
  {
    insertInitialMapping(initialMapping);
  }

  VertexMapping(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const std::vector<std::vector<MappingIndex>>& hint)
  {
    initCompatibleNodes(query, target, nodeComp, hint);
    initAdjacentNodes(query);
  }

  VertexMapping(VertexMapping&&) noexcept = default;
  VertexMapping& operator=(VertexMapping&&) noexcept = default;
  //VertexMapping(const VertexMapping&) = delete;
  VertexMapping(const VertexMapping&) = default;
  VertexMapping& operator=(const VertexMapping&) = delete;

  static MappingPair initIndexTuple(
      MappingIndex g0Node)
  {
    // Produce an overflow on the last index in order to start with an all zero tuple after the
    // first index increment.
    return {g0Node, NODE_END};
  }

  /**
   * @brief getNextNode
   * @return The index of the next node or NODE_END
   */
  MappingIndex getNextNode();

  template<typename Policy>
  void filterAdjacentNodes(
      const Policy& policy) {
    m_adjacentNodes.erase(std::remove_if(m_adjacentNodes.begin(), m_adjacentNodes.end(),
                                         [&policy](const AdjacentNode& adjacencyTuple) {
                                           return policy(adjacencyTuple.m_node);
                                         }), m_adjacentNodes.end());
    m_hasAdjacentNodes = !m_adjacentNodes.empty();
  }

  virtual ~VertexMapping() = default;

public:


  void prepareDisconnected(
      const AdjacencyFunctor& query)
  {
    RIMACS_ASSERT(m_adjacentNodes.empty());
    initAdjacentNodes(query);
  }


  bool getNextMatch(
      MappingPair& matchIndices,
      RealMappingTuple& mappedNodes);


  const CompatibleNodesType& getNofExtendableNodes() const
  {
    return m_nofExtendableNodes;
  }


  CompatibleNodesType getNofConnectedExtendableNodes() const
  {
    auto extendable = getNofExtendableNodes();
    // we can't use the container of adjacent nodes since it's modified during mapping iteration
    // thus for the last extension, the current mapping could be considered as not connected extendable
    return m_hasAdjacentNodes ? extendable : CompatibleNodesType{};
  }

  void setExtendableNodes(
      const CompatibleNodesType& nofNodes)
  {
    RIMACS_ASSERT(nofNodes <= m_nofExtendableNodes);
    m_nofExtendableNodes = nofNodes;
  }

  void reinitializeAdjacentNodes(
      const AdjacencyFunctor& graph,
      const MappingPair& newNode,
      const std::vector<MappingIndex>& mappingStatus)
  {
    m_adjacentNodes.clear();
    updateAdjacentNodes(graph, newNode, mappingStatus);
  }

  bool nodeIsMappable(MappingIndex node) const
  {
    return !m_compatibleNodes.at(node).empty();
  }

  void finishedNode(
      MappingIndex node);


  const std::vector<std::vector<StoredCompatibleNodesType>>& getCompatibleNodes() const
  {
    return m_compatibleNodes;
  }

  void setNodeUnmappable(
      MappingIndex node,
      bool ensureNodeIsNotAdjacent);


  friend class ::VertexMappingTester;
  friend class ::SubgraphMappingTester;


  std::vector<MappingIndex> updateExtendableNodesAndStatus(
      std::vector<MappingIndex>& status);


  void updateCompatibleNodes(
      const VertexMapping& other,
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const EdgeCompatibilityFunctor& edgeComp,
      const std::vector<MappingIndex>& status,
      const MappingPair& newNode);


  void updateAdjacentNodes(
      const AdjacencyFunctor& graph,
      const MappingPair& nodeMappings,
      const std::vector<MappingIndex>& mappingStatus);

protected:

  void insertInitialMapping(
      const MappingPair& initialMapping);

  void initAdjacentNodes(
      const AdjacencyFunctor& query);


  void initCompatibleNodes(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp);

  void initCompatibleNodes(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const std::vector<std::vector<MappingIndex>>& hint);

  void doCopyConstruct(
      const VertexMapping& other,
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const EdgeCompatibilityFunctor& edgeComp,
      const MappingPair& newNode,
      std::vector<MappingIndex>& mappingStatus,
      std::vector<MappingIndex>& newUnmappableNodes);

  template<typename DerivedMapping,
      typename=std::enable_if_t<std::is_base_of<VertexMapping, DerivedMapping>::value
                                && !std::is_same<VertexMapping, DerivedMapping>::value>>
  VertexMapping(
      const DerivedMapping& other);

private:

  // for each node of the first graph, there is an array for each other graph containing the list
  // of compatible nodes
  std::vector<std::vector<StoredCompatibleNodesType>> m_compatibleNodes;
  std::vector<AdjacentNode> m_adjacentNodes;
  CompatibleNodesType m_nofExtendableNodes;
  bool m_hasAdjacentNodes = false;
};


} // namespace RIMACS

#include <RIMACS/Algorithm/VertexMappingTemplate.hpp>
