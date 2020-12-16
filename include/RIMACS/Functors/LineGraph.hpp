
#pragma once

#include <RIMACS/Forward.hpp>
#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Functors/GraphCachingFunctor.hpp>

namespace RIMACS {

class LineGraphNodeCompatibilityFunctor;
class LineGraphEdgeCompatibilityFunctor;

/**
 * @brief The LineGraph class
 *  The LineGraph is necessary to compute a noninduced edge based MCS on two graphs.
 *  Each edge of the original graph is represented by a node in the line graph.
 *  Each edge in the line graph corresponds to a node, but nodes with a degree greater than two
 *  are represented by multiple edges in the line graph and nodes with a degree of one or less
 *  are not represented by any edge in the line graph.
 */
class BasicLineGraph : public AdjacencyFunctor {
public:

  explicit BasicLineGraph(
      const AdjacencyFunctor& base)
    : m_baseFunctor(base)
    , m_maxOrderSuggestion((m_baseFunctor.getMaxNodeOrderSuggestion()
                            + m_baseFunctor.getMaxNodeOrderSuggestion()) / 2)
  {
    init(base);
  }

  bool adjacent(
      MappingIndex fromNode,
      MappingIndex toNode) const override;

  MappingIndex getEdgeId(
      MappingIndex nodeFromId,
      MappingIndex nodeToId) const override;

  const MappingIndex& getMappedNode(
    MappingIndex nodeFromId,
    MappingIndex nodeToId) const;

  const MappingIndex& getLineNodeFromNode(
      MappingIndex node) const;

  const MappingIndex& getLineNodeToNode(
      MappingIndex node) const;

  MappingIndex edgeGetFromId(
      MappingIndex edgeId) const override;

  MappingIndex edgeGetToId(
      MappingIndex edgeId) const override;

  MappingIndex nodeGetNofEdges(
      MappingIndex nodeIdx) const override;

  std::vector<MappingIndex> nodeGetEdges(
      MappingIndex nodeIdx) const override;

  std::vector<MappingIndex> nodeGetNeighbours(
      MappingIndex nodeIdx) const override;

  MappingIndex getNofNodes() const override;

  MappingIndex getNofEdges() const override;


  ComponentIndex getNodeOrderSuggestion(
      MappingIndex node) const override
  {
    RIMACS_ASSERT(false && "Do not use this function directly! Use a LineGraph instead");
    return ComponentIndex{-1};
  }

  ComponentIndex getMaxNodeOrderSuggestion() const override
  {
    RIMACS_ASSERT(false && "Do not use this function directly! Use a LineGraph instead");
    return m_maxOrderSuggestion;
  }

private:

  friend class LineGraphNodeCompatibilityFunctor;
  friend class LineGraphEdgeCompatibilityFunctor;

  /**
   * @brief The line_graph_edge struct
   *  Internal struct representing an edge in the line graph
   *  storing basic edge information, the from and to linegraph node indices as well as
   *  the (real) node between the edges of the original graph.
   */
  struct line_graph_edge {

    line_graph_edge(line_graph_edge&&) noexcept = default;
    line_graph_edge(const line_graph_edge&) noexcept = default;
    line_graph_edge& operator=(line_graph_edge&&) noexcept = default;
    line_graph_edge& operator=(const line_graph_edge&) = default;

    line_graph_edge(
        MappingIndex from,
        MappingIndex to,
        MappingIndex commonNode)
      : m_from(from)
      , m_to(to)
      , m_commonNode(commonNode)
    {}

    bool operator<(
        const line_graph_edge& other) const
    {
      if(m_from != other.m_from) {
        return m_from < other.m_from;
      }
      if(m_to != other.m_to) {
        return m_to < other.m_to;
      }
      return m_commonNode < other.m_commonNode;
    }

    MappingIndex m_from;
    MappingIndex m_to;
    MappingIndex m_commonNode;
  };


  /**
   * @brief The line_graph_node struct
   *  Internal struct representing a node of the linegraph, thus an edge in the (real) graph.
   *  It stores the index of the represented edge in the (real) graph as well as the edges incident
   *  node indices.
   */
  struct line_graph_node {

    line_graph_node(line_graph_node&&) noexcept = default;
    line_graph_node(const line_graph_node&) noexcept = default;
    line_graph_node& operator=(line_graph_node&&) noexcept = default;
    line_graph_node& operator=(const line_graph_node&) = default;

    line_graph_node(
        MappingIndex edgeId,
        MappingIndex fromNode,
        MappingIndex toNode)
      : m_baseEdgeId(edgeId)
      , m_fromNode(fromNode)
      , m_toNode(toNode)
    {}
    MappingIndex m_baseEdgeId;
    MappingIndex m_fromNode;
    MappingIndex m_toNode;
  };


  inline typename std::vector<line_graph_edge>::const_iterator find_edge(
      MappingIndex& nodeFrom,
      MappingIndex& nodeTo) const;

  void init(
      const AdjacencyFunctor& baseFunctor);

  struct node_edge_info{
    node_edge_info(
        MappingIndex id,
        MappingIndex to)
      : m_edgeId(id)
      , m_to(to)
    {}
    MappingIndex m_edgeId;
    MappingIndex m_to;
  };

  std::vector<line_graph_node> m_nodes;
  std::vector<line_graph_edge> m_edges;
  std::vector<std::vector<node_edge_info>> m_nodeEdges;
  const AdjacencyFunctor& m_baseFunctor;
  ComponentIndex m_maxOrderSuggestion;
};


class LineGraph : public GraphCachingAdjacencyFunctor<BasicLineGraph> {
public:

  LineGraph(
      const AdjacencyFunctor& base)
    : GraphCachingAdjacencyFunctor(base)
  {}

  ~LineGraph();

  operator const BasicLineGraph&() const
  {
    return this->getBaseFunctor();
  }

};


/**
 * @brief The LineGraphNodeCompatibilityFunctor class
 *  Overload for node compatibility in line graphs.
 *  Nodes in a line graph are compatible, if the edges and their incident nodes are compatible.
 */
class LineGraphNodeCompatibilityFunctor : public NodeCompatibilityFunctor {
public:

  LineGraphNodeCompatibilityFunctor(
      const BasicLineGraph& query,
      const BasicLineGraph& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp)
    : m_query(query)
    , m_target(target)
    , m_nodeCompatibility(nodeComp)
    , m_edgeCompatibility(edgeComp)
  {}

  bool nodesAreCompatible(
      MappingIndex queryNode,
      MappingIndex targetNode) const override
{
    const BasicLineGraph::line_graph_node& qn = m_query.m_nodes.at(queryNode);
    const BasicLineGraph::line_graph_node& tn = m_target.m_nodes.at(targetNode);
    return m_edgeCompatibility(qn.m_baseEdgeId, tn.m_baseEdgeId)
           && ((m_nodeCompatibility(qn.m_fromNode, tn.m_fromNode)
                && m_nodeCompatibility(qn.m_toNode, tn.m_toNode))
               || (m_nodeCompatibility(qn.m_fromNode, tn.m_toNode)
                   && m_nodeCompatibility(qn.m_toNode, tn.m_fromNode)));
  }

  /**
   * @brief mappingIsValid
   *  Overload for the mapping is valid method. It is necessary to ensure that no Delta Y exchange
   *  occurs during the line graph mapping.
   */
  bool mappingIsValid(
      typename std::vector<MappingPair>::const_iterator beginMapping,
      typename std::vector<MappingPair>::const_iterator endMapping) const override;

  /**
   * @brief mappingIsValid
   *  Overload for the other mappingIsValid method as well.
   *  Overloading only one method can lead to undefined matching behaviour.
   */
  bool mappingIsValid(
      SwappedMappingIterator beginMapping,
      SwappedMappingIterator endMapping) const override;

  double estimateNodeSimilarity(
      MappingIndex queryNode,
      MappingIndex targetNode) const override
  {
    const BasicLineGraph::line_graph_node& qn = m_query.m_nodes.at(queryNode);
    const BasicLineGraph::line_graph_node& tn = m_target.m_nodes.at(targetNode);
    return std::max((m_nodeCompatibility.estimateNodeSimilarity(qn.m_fromNode, tn.m_fromNode)
                     + m_nodeCompatibility.estimateNodeSimilarity(qn.m_toNode, tn.m_toNode)),
                    (m_nodeCompatibility.estimateNodeSimilarity(qn.m_fromNode, tn.m_toNode)
                     + m_nodeCompatibility.estimateNodeSimilarity(qn.m_toNode, tn.m_fromNode)))
           / 2;
  }

protected:

  /**
   * @brief isValidImpl
   *  Implementation of the validity test for a given mapping. Since there are
   *  two possible iterator types, it is advantageous to combine both variants in
   *  a shared template based implementation. Both types under the iterator provide
   *  a 'first' and a 'second' member, corresponding to the indices of the mapped
   *  elements in the first, respective second graph as expected by this functor.
   * @tparam Iterator Iterator template corresponding to the vector iterator
             or the SwappedMappingIterator
   * @param begin begin iterator
   * @param end end iterator
   * @return false, if the range [begin, end) describes a delta- Y exchange, true otherwise
   */
  template<typename Iterator>
  inline bool isValidImpl(
      Iterator begin,
      Iterator end) const;

  const BasicLineGraph& m_query;
  const BasicLineGraph& m_target;
  const NodeCompatibilityFunctor& m_nodeCompatibility;
  const EdgeCompatibilityFunctor& m_edgeCompatibility;
};

class PreferMaxNumberOfMappedNodesLineGraphNodeCompatibilityFunctor
    : public LineGraphNodeCompatibilityFunctor {
public:

  PreferMaxNumberOfMappedNodesLineGraphNodeCompatibilityFunctor(
      const BasicLineGraph& query,
      const BasicLineGraph& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const Config& config)
    : LineGraphNodeCompatibilityFunctor(query, target, nodeComp, edgeComp)
    , m_config(config)
    , m_unmappedNodePenaltyValue(initializePenaltyValue(query.getNofNodes()))
  {}

  double estimateFinalWeight(
      typename std::vector<MappingPair>::const_iterator beginMapping,
      typename std::vector<MappingPair>::const_iterator endMapping,
      double currentWeight) const override;

  double estimateFinalWeight(
      SwappedMappingIterator beginMapping,
      SwappedMappingIterator endMapping,
      double currentWeight) const override;

private:

  double calculatePenaltyValue(
      size_t mappedRealNodes,
      size_t mappedEdges,
      size_t nofComponents) const
  {
    return std::min(m_unmappedNodePenaltyValue
                    * static_cast<double>(mappedEdges + nofComponents - mappedRealNodes),
                    1 - m_unmappedNodePenaltyValue);
  }

  inline MappingIndex calculateNofComponents(
      MappingIndex mappedEdges) const;

  template<typename Iterator>
  double estimateFinalWeight_impl(
      Iterator beginMapping,
      Iterator endMapping,
      double currentWeight) const;

  static double initializePenaltyValue(
      MappingIndex lineGraphNofNodes);

  const Config& m_config;
  double m_unmappedNodePenaltyValue;
};

class LineGraphEdgeCompatibilityFunctor : public EdgeCompatibilityFunctor {
public:

  LineGraphEdgeCompatibilityFunctor(
      const BasicLineGraph& query,
      const BasicLineGraph& target,
      const NodeCompatibilityFunctor& nodeComp)
    : m_query(query)
    , m_target(target)
    , m_nodeCompatibility(nodeComp)
    , m_cache(2 * query.getNofEdges() * target.getNofEdges(), false)
    , m_queryNofEdges(query.getNofEdges())
  {}

  ~LineGraphEdgeCompatibilityFunctor() override;

  bool edgesAreCompatible(
      MappingIndex queryEdge,
      MappingIndex targetEdge) const override
  {
    auto edge_cache_it = m_cache.begin() + 2 * (m_queryNofEdges * targetEdge + queryEdge);
    if(*edge_cache_it) {
      return edge_cache_it[1];
    }
    *edge_cache_it = true;
    const BasicLineGraph::line_graph_edge& q = m_query.m_edges.at(queryEdge);
    const BasicLineGraph::line_graph_edge& t = m_target.m_edges.at(targetEdge);
    return (edge_cache_it[1] = m_nodeCompatibility(q.m_commonNode, t.m_commonNode));
  }

private:
  const BasicLineGraph& m_query;
  const BasicLineGraph& m_target;
  const NodeCompatibilityFunctor& m_nodeCompatibility;
  mutable std::vector<bool> m_cache;
  MappingIndex m_queryNofEdges;
};


/**
 * @brief The EdgeMappingConverter class
 *  Converts an edge mapping to a node based mapping
 */
class EdgeMappingConverter {
public:

  EdgeMappingConverter(
      const BasicLineGraph& query,
      const BasicLineGraph& target,
      const NodeCompatibilityFunctor& functor)
    : m_query(query)
    , m_target(target)
    , m_functor(functor)
  {}

  /**
   * @brief operator ()
   *  Converts a mappings of edges indices to a node index mapping
   *  If there are components consisting of exactly one node, they can be expanded to
   *  node mappings in both directions.
   * @param mapping The mapping to be converted
   * @param cumulatedComponentSizes Components of the node mapping
   * @param expandLoneBondMappings Flag whether mapping consisting of exactly one edge but both node
   *        are wise versa compatible should be expanded in both directions
   * @return The vector of pair of edge mappings and updated cumulated component sizes.
   */
  std::vector<std::pair<std::vector<MappingPair>, std::vector<MappingIndex>>> operator()(
      const std::vector<MappingPair>& mapping,
      const std::vector<MappingIndex>& cumulatedComponentSizes,
      bool expandLoneBondMappings) const;

  /**
   * @brief operator ()
   *  Converts a mapping of edges indices to a node index mapping
   * @param mapping The mapping to be converted
   * @param cumulatedComponentSizes Components of the node mapping
   * @return The edge mapping and the updated cumulated component sizes.
   */
  std::pair<std::vector<MappingPair>, std::vector<MappingIndex>> operator()(
      const std::vector<MappingPair>& mapping,
      const std::vector<MappingIndex>& cumulatedComponentSizes) const
  {
    return operator()(mapping, cumulatedComponentSizes, false).front();
  }

private:
  const BasicLineGraph& m_query;
  const BasicLineGraph& m_target;
  const NodeCompatibilityFunctor& m_functor;
};

} // namespace RIMACS
