
#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>
#include <RIMACS/ASSERT.hpp>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Structs.hpp>
#include <RIMACS/Functors/SwappedMapping.hpp>


namespace RIMACS {

class CachingRequiredMinimalAdjacencyFunctor {
public:
  // Return value, if there is no edge between two tested atoms
  enum : MappingIndex { NO_EDGE = std::numeric_limits<MappingIndex>::max() };

  virtual ~CachingRequiredMinimalAdjacencyFunctor();

  /**
   * @brief edgeGetFromId
   * @param edgeId Internal index of an edge
   * @return Internal index of the first incident node of the edge
   */
  virtual MappingIndex edgeGetFromId(
      MappingIndex edgeId) const = 0;

  /**
   * @brief edgeGetToId
   * @param edgeId Internal index of an edge
   * @return Internal index of the second incident node of the edge
   */
  virtual MappingIndex edgeGetToId(
      MappingIndex edgeId) const = 0;

  /**
   * @brief getNofNodes
   * @return Number of nodes in the graph
   *         (Last internal node index + 1)
   */
  virtual MappingIndex getNofNodes() const = 0;

  /**
   * @brief getNofEdges
   * @return number of edges in the graph
   *         (Last internal edge index + 1)
   */
  virtual MappingIndex getNofEdges() const = 0;

};

/**
 * @brief The AdjacencyFunctor class
 *  The adjacency functor wraps the implementation of an actual graph for the MCS calculation.
 *  It has to provide indices for all nodes and edges and represent consistent index based
 *  adjacency.
 *  All internal indices of nodes and edges have to be continues starting at 0.
 *  Graphs are assumed to be undirected.
 */
class AdjacencyFunctor : public CachingRequiredMinimalAdjacencyFunctor {
public:

  ~AdjacencyFunctor() override;

  /**
   * @brief adjacent
   * @param node1 Internal index of first node
   * @param node2 Internal index of second node
   * @return true, if the nodes at the indices node1 and node2 are adjacent
   */
  virtual bool adjacent(
      MappingIndex node1,
      MappingIndex node2) const = 0;

  /**
   * @brief operator ()
   * Forward for adjacent
   */
  bool operator()(
      MappingIndex node1,
      MappingIndex node2) const
  {
    return adjacent(node1, node2);
  }

  /**
   * @brief nodeGetEdges
   * @param nodeIdx internal index of the node
   * @return Vector containing the internal indices of all incident edges
   */
  virtual std::vector<MappingIndex> nodeGetEdges(
      MappingIndex nodeIdx) const = 0;

  /**
   * @brief nodeGetNofEdges
   * @param nodeIdx Internal index of the node
   * @return degree of the specified node
   */
  virtual MappingIndex nodeGetNofEdges(
      MappingIndex nodeIdx) const
  {
    return static_cast<MappingIndex>(nodeGetEdges(nodeIdx).size());
  }

  /**
   * @brief getEdgeId
   *  Get the internal index of an edge connecting two nodes
   * @param nodeFromId Internal index of the first node
   * @param nodeToId Internal index of the second node
   * @return The internal index of the edge or NO_EDGE
   */
  virtual MappingIndex getEdgeId(
      MappingIndex nodeFromId,
      MappingIndex nodeToId) const = 0;

  /**
   * @brief nodeGetNeighbours
   * @param nodeIdx Internal node index
   * @return Vector containing the internal indices of all adjacent nodes
   */
  virtual std::vector<MappingIndex> nodeGetNeighbours(
      MappingIndex nodeIdx) const = 0;

  virtual ComponentIndex getNodeOrderSuggestion(MappingIndex) const
  {
    return 0;
  }

  virtual ComponentIndex getMaxNodeOrderSuggestion() const
  {
    return std::numeric_limits<ComponentIndex>::max();
  }

protected:

  std::vector<ComponentIndex> getNodeOrderSuggestionVector(
      double partial) const;
};


/**
 * @brief The NodeCompatibilityFunctor class
 *  A wrapper class providing information about node compatibility for the graphs used in the MCS
 *  calculation.
 */
class NodeCompatibilityFunctor {
public:

  virtual ~NodeCompatibilityFunctor();

  /**
   * @brief nodesAreCompatible
   * @param fstGraphNode Internal index of the node in the first graph
   * @param secGraphNode Internal index of the node in the second graph
   * @return true, if nodes are allowed to be mapped onto each other, false otherwise
   */
  virtual bool nodesAreCompatible(
      MappingIndex fstGraphNode,
      MappingIndex secGraphNode) const = 0;

  bool nodesAreCompatible(
      const MappingPair& graphNodes) const
  {
    return nodesAreCompatible(graphNodes.m_from, graphNodes.m_to);
  }

  /**
   * @brief operator ()
   * Wrapper for nodesAreCompatible
   */
  bool operator()(
      MappingIndex fstGraphNode,
      MappingIndex secGraphNode) const
  {
    return nodesAreCompatible(fstGraphNode, secGraphNode);
  }

  bool operator()(
      const MappingPair& graphNodes) const
  {
    return nodesAreCompatible(graphNodes.m_from, graphNodes.m_to);
  }

  /**
   * @brief getWeight
   *  Provide a weight for a given node mapping.
   *  If no weight there is no overload, node will be weighted 1.0 if they are compatible
   *  and 0.0 otherwise. The algorithm ensures that nodes are compatible if the weight is requested.
   *  Weights MUST be positive values. Negative values and a weight of 0.0 are used internally.
   *  Those weights may lead to incorrect results.
   * @param fstGraphNode Internal index of the node in the first graph
   * @param secGraphNode Internal index of the node in the second graph
   * @return A weight for the mapping.
   *         (This should only be overloaded if weights are not 1 and 0)
   */
  virtual double getWeight(
      MappingIndex fstGraphNode,
      MappingIndex secGraphNode) const
  {
    return nodesAreCompatible(fstGraphNode, secGraphNode) ? 1.0 : 0.0;
  }

  double getWeight(
      const MappingPair& nodes) const
  {
    return getWeight(nodes.m_from, nodes.m_to);
  }

  /**
   * Make an estimation for internal node ordering.
   * This does not affect the computed mappings, it's about to improve the state space traversal
   * strategies. If no overload is provided, the node mapping weight is used as approximation
   * @param node1
   * @param node2
   * @return
   */
  virtual double estimateNodeSimilarity(
      MappingIndex fstGraphNode,
      MappingIndex secGraphNode) const
  {
    return getWeight(fstGraphNode, secGraphNode);
  }

  /**
   * @brief mappingIsValid
   *  Test whether a connected component of a mapping yields a valid mapping
   *  This might be affected due to the delta-Y-exchange on line graphs
   *  or other context specific constraints
   *  IMPORTANT: If you want to overload this function, do not forget to overload
   *             the other mappingIsValid method as well. If the graph order was swapped
   *             then the mappingIsValid method accepting SwappedPairIterators is called.
   *             It is recommended to implement a template based mappingIsValid_impl method
   *             that just accepts iterators. This method then can access the ->first and
   *             ->second members of the iterator. This works for both variants.
   * If you are looking for an example, consider the LineGraphNodeCompatibilityFunctor
   * @param beginMapping Iterator to the first mapping index pair
   * @param endMapping Corresponding end iterator
   * @return true is mapping is valid, false otherwise
   */
  virtual bool mappingIsValid(
      typename std::vector<MappingPair>::const_iterator /*beginMapping*/,
      typename std::vector<MappingPair>::const_iterator /*endMapping*/) const
  {
    return true;
  }

  /**
   * @brief mappingIsValid
   * Support for mappingIsValid if the order of the molecules was swapped.
   * Consider the documentation of the base variant of mappingIsValid using the vector iterators.
   */
  virtual bool mappingIsValid(
      SwappedMappingIterator /*beginMapping*/,
      SwappedMappingIterator /*endMapping*/) const
  {
    return true;
  }

  /**
   * @brief estimateFinalWeight given the mapping, components and the sum of the weighted
   *  node mapping, this function allows to update the final weight estimation as an
   *  customisation point
   *
   *  IMPORTANT: If you want to overload this function, do not forget to overload
   *             the other mappingIsValid method as well. If the graph order was swapped
   *             then the mappingIsValid method accepting SwappedPairIterators is called.
   *             It is recommended to implement a template based mappingIsValid_impl method
   *             that just accepts iterators. This method then can access the ->first and
   *             ->second members of the iterator. This works for both variants.
   * @param beginMapping the begin iterator to the mapping, the graphs are access
   *                     using the fields first and second
   * @param endMapping the corresponding end iterator for the mapping
   * @param currentWeight The weight estimation based on the weight for the node mappings.
   *                      The returned weight MUST NOT be greater.
   * @return The updated weight value. It this function is not overridden, the estimation is returned
   */
  virtual double estimateFinalWeight(
      typename std::vector<MappingPair>::const_iterator /*beginMapping*/,
      typename std::vector<MappingPair>::const_iterator /*endMapping*/,
      double currentWeight) const
  {
    return currentWeight;
  }

  /**
   * @brief estimateFinalWeight
   * Support for estimateFinalWeight if the order of the molecules was swapped.
   * Consider the documentation of the base variant of estimateFinalWeight using the vector iterators.
   */
  virtual double estimateFinalWeight(
      SwappedMappingIterator /*beginMapping*/,
      SwappedMappingIterator /*endMapping*/,
      double currentWeight) const
  {
    return currentWeight;
  }
};


/**
 * @brief The EdgeCompatibilityFunctor class
 *  A wrapper class providing information about edge compatibility for the graphs used in the MCS
 *  calculation.
 */
class EdgeCompatibilityFunctor {
public:
  virtual ~EdgeCompatibilityFunctor();

  /**
   * @brief edgesAreCompatible
   * @param fstGraphNode Internal index of the node in the first graph
   * @param secGraphNode Internal index of the node in the second graph
   * @return true if edges are allowed to be mapped onto each other, false otherwise
   *         (compatibility of incident nodes is ensured by the algorithm)
   */
  virtual bool edgesAreCompatible(
      MappingIndex fstGraphEdge,
      MappingIndex secGraphEdge) const = 0;

  /**
   * @brief operator ()
   * Wrapper for edgesAreCompatible
   */
  bool operator()(
      MappingIndex fstGraphEdge,
      MappingIndex secGraphEdge) const
  {
    return edgesAreCompatible(fstGraphEdge, secGraphEdge);
  }
};

} // namespace RIMACS
