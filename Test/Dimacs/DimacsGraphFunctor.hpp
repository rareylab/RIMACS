
#pragma once

#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Functors/GraphCachingFunctor.hpp>

#include "DimacsGraph.hpp"

namespace Dimacs {

class BasicDimacsAdjacencyFunctor : public RIMACS::CachingRequiredMinimalAdjacencyFunctor {
public:

  ~BasicDimacsAdjacencyFunctor() override;

  BasicDimacsAdjacencyFunctor(
      const DimacsGraph& graph)
    : m_graph(graph)
  {}

  BasicDimacsAdjacencyFunctor(
      DimacsGraph&& graph)
    : BasicDimacsAdjacencyFunctor(std::make_unique<DimacsGraph>(std::move(graph)))
  {}

  BasicDimacsAdjacencyFunctor(
      const std::string& file)
    : BasicDimacsAdjacencyFunctor(std::make_unique<DimacsGraph>(file))
  {}

  RIMACS::MappingIndex edgeGetFromId(
      RIMACS::MappingIndex edgeId) const override
  {
    return m_graph.m_edges.at(edgeId).m_from;
  }

  RIMACS::MappingIndex edgeGetToId(
      RIMACS::MappingIndex edgeId) const override
  {
    return m_graph.m_edges.at(edgeId).m_to;
  }

  RIMACS::MappingIndex getNofNodes() const override
  {
    return m_graph.m_nodeLabels.size();
  }

  RIMACS::MappingIndex getNofEdges() const override
  {
    return m_graph.m_edges.size();
  }

  operator const DimacsGraph&() const
  {
    return m_graph;
  }

private:

  BasicDimacsAdjacencyFunctor(std::unique_ptr<DimacsGraph> ownedGraph)
    : m_graph(*ownedGraph)
    , m_ownedGraph(ownedGraph.get())
  {
    // manual transfer of ownership with explicit release of the graph, in order to guarantee
    // that the reference member is initialized with a valid reference and not with a nullptr
    ownedGraph.release();
  }

  const DimacsGraph& m_graph;
  std::unique_ptr<DimacsGraph> m_ownedGraph;
};

class DimacsAdjacencyFunctor : public RIMACS::GraphCachingAdjacencyFunctor<BasicDimacsAdjacencyFunctor> {
public:

  ~DimacsAdjacencyFunctor() override;

  explicit DimacsAdjacencyFunctor(
      const DimacsGraph& graph)
    : GraphCachingAdjacencyFunctor(graph)
  {}

  explicit DimacsAdjacencyFunctor(
      DimacsGraph&& graph)
    : GraphCachingAdjacencyFunctor(std::move(graph))
  {}

  operator const DimacsGraph&() const
  {
    return this->m_realInstance;
  }
};


class DimacsNodeCompatibilityFunctor : public RIMACS::NodeCompatibilityFunctor {
public:

  DimacsNodeCompatibilityFunctor(
      const DimacsAdjacencyFunctor& query,
      const DimacsAdjacencyFunctor& target)
    : m_queryGraph(query)
    , m_targetGraph(target)
  {}

  ~DimacsNodeCompatibilityFunctor() override;

  bool nodesAreCompatible(
      RIMACS::MappingIndex queryNode,
      RIMACS::MappingIndex targetNode) const override
  {
    const DimacsGraph& query = m_queryGraph;
    const DimacsGraph& target = m_targetGraph;
    return query.m_nodeLabels.at(queryNode) == target.m_nodeLabels.at(targetNode);
  }

  double estimateNodeSimilarity(
      RIMACS::MappingIndex queryNode,
      RIMACS::MappingIndex targetNode) const override;

protected:
  const DimacsAdjacencyFunctor& m_queryGraph;
  const DimacsAdjacencyFunctor& m_targetGraph;
};


class DimacsEdgeCompatibilityFunctor : public RIMACS::EdgeCompatibilityFunctor {
public:

  ~DimacsEdgeCompatibilityFunctor() override;

  DimacsEdgeCompatibilityFunctor(
      const DimacsGraph& query,
      const DimacsGraph& target)
    : m_queryGraph(query)
    , m_targetGraph(target)
  {}

  bool edgesAreCompatible(
      RIMACS::MappingIndex queryEdge,
      RIMACS::MappingIndex targetEdge) const override
  {
    return m_queryGraph.m_edges.at(queryEdge).m_label == m_targetGraph.m_edges.at(targetEdge).m_label;
  }

private:
  const DimacsGraph& m_queryGraph;
  const DimacsGraph& m_targetGraph;
};

} // namespace Dimacs
