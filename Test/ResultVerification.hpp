
#pragma once

#include <unordered_set>

#include <gtest/gtest.h>


#include <RIMACS/Forward.hpp>
#include <RIMACS/MCSResult.hpp>
#include <unordered_map>

#include "Dimacs/DimacsGraph.hpp"
#include "Dimacs/DimacsGraphFunctor.hpp"

extern std::string path;

namespace {

bool connected(
    const Dimacs::DimacsGraph& graph,
    RIMACS::MappingIndex from,
    RIMACS::MappingIndex to)
{
  return std::find_if(graph.m_edges.begin(), graph.m_edges.end(),
                      [from, to](const auto& e) {return (e.m_from == from && e.m_to == to)
                                                     || (e.m_from == to && e.m_to == from);})
         != graph.m_edges.end();
}

void verify_connected_components(
    const std::vector<std::unordered_set<RIMACS::MappingIndex>>& connectivity,
    std::vector<RIMACS::MappingIndex> connectedComponents)
{
  if(connectivity.size() < 2 || connectedComponents.empty()) {
    return;
  }
  connectedComponents.insert(connectedComponents.begin(), RIMACS::MappingIndex{0});
  for(auto it = connectedComponents.begin(), next = std::next(it), last = connectedComponents.end()
      ; next != last
      ; ++it, ++next) {
    std::set<RIMACS::MappingIndex> lastComponent;
    std::set<RIMACS::MappingIndex> component = {connectivity.at(*it).begin(), connectivity.at(*it).end()};
    while(component != lastComponent) {
      lastComponent = component;
      std::for_each(lastComponent.begin(), lastComponent.end(),
                   [&connectivity, &component] (RIMACS::MappingIndex neighbour) {
                     component.insert(connectivity.at(neighbour).begin(), connectivity.at(neighbour).end());
      });
    }
    component.insert(*it);
    bool identical = component.size() == *next - *it;
    EXPECT_EQ(component.size(), *next - *it);
    EXPECT_EQ(*component.begin(), *it);
    EXPECT_EQ(*component.rbegin(), (*next) - 1);
    EXPECT_TRUE(identical); // debug var
  }
}

int get_common_node_idx(
    const Dimacs::DimacsEdge& e1,
    const Dimacs::DimacsEdge& e2)
{
  if(e1 == e2) {return -1;}
  if(e1.m_from == e2.m_from || e1.m_from == e2.m_to) {
    return e1.m_from;
  }
  if(e1.m_to == e2.m_from || e1.m_to == e2.m_to) {
    return e1.m_to;
  }
  return -1;
}

int get_common_node_idx(
    const Dimacs::DimacsGraph& graph,
    int e1,
    int e2)
{
  return get_common_node_idx(graph.m_edges.at(e1), graph.m_edges.at(e2));
}

int get_other_node_idx(
    const Dimacs::DimacsEdge& e1,
    int node_idx)
{
  if(e1.m_from == node_idx) {
    return e1.m_to;
  }
  if(e1.m_to == node_idx) {
    return e1.m_from;
  }
  return -1;
}

template<bool induced=true>
void verify_result(
    const Dimacs::DimacsAdjacencyFunctor& first,
    const Dimacs::DimacsAdjacencyFunctor& second,
    const Dimacs::DimacsNodeCompatibilityFunctor& nodeComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& edgeComp,
    const std::vector<RIMACS::MappingPair>& mapping,
    const std::vector<RIMACS::MappingIndex>& connectedComponents = std::vector<RIMACS::MappingIndex>())
{
  const Dimacs::DimacsGraph& q = first;
  const Dimacs::DimacsGraph& t = second;
  std::vector<std::unordered_set<RIMACS::MappingIndex>> connectivity(mapping.size());
  for(size_t i = 0; i < mapping.size(); ++i) {
    const RIMACS::MappingPair& mp = mapping[i];
    EXPECT_TRUE(nodeComp.nodesAreCompatible(mp.m_from, mp.m_to));
    for(size_t j = 0; j < mapping.size(); ++j) {
      const RIMACS::MappingPair& toPair = mapping[j];
      RIMACS::MappingIndex queryEdgeIdx = first.getEdgeId(mp.m_from, toPair.m_from);
      RIMACS::MappingIndex targetEdgeIdx = second.getEdgeId(mp.m_to, toPair.m_to);
      ASSERT_EQ(queryEdgeIdx != std::numeric_limits<RIMACS::MappingIndex>::max() && induced,
                targetEdgeIdx != std::numeric_limits<RIMACS::MappingIndex>::max() && induced);
      if(queryEdgeIdx != std::numeric_limits<RIMACS::MappingIndex>::max()
         && targetEdgeIdx != std::numeric_limits<RIMACS::MappingIndex>::max()) {
        ASSERT_TRUE(edgeComp(queryEdgeIdx, targetEdgeIdx));
        connectivity.at(i).insert(j);
      }
    }
  }
  verify_connected_components(connectivity, connectedComponents);
}


template<bool induced=true>
void verify_result(
    const Dimacs::DimacsAdjacencyFunctor& first,
    const Dimacs::DimacsAdjacencyFunctor& second,
    const Dimacs::DimacsNodeCompatibilityFunctor& nodeComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& edgeComp,
    const RIMACS::MCSResult& result)
{
  verify_result<induced>(first, second, nodeComp, edgeComp, result.getMappings(), result.getComponentSizes());
}


inline size_t get_number_of_incident_bonds(
    const std::unordered_map<RIMACS::MappingIndex, size_t>& mappedAtoms,
    const RIMACS::MappingIndex& a)
{
  return mappedAtoms.count(a) ? mappedAtoms.at(a) : 1;
}

void verify_mces_result(
    const Dimacs::DimacsAdjacencyFunctor& first,
    const Dimacs::DimacsAdjacencyFunctor& second,
    const Dimacs::DimacsNodeCompatibilityFunctor& nodeComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& edgeComp,
    const std::vector<RIMACS::MappingPair>& mapping,
    const std::vector<RIMACS::MappingIndex>& connectedComponents = std::vector<RIMACS::MappingIndex>())
{
  const Dimacs::DimacsGraph& q = first;
  const Dimacs::DimacsGraph& t = second;
  std::vector<std::unordered_set<RIMACS::MappingIndex>> connectivity(mapping.size());
  std::unordered_set<RIMACS::MappingIndex> mappedQueryEdges, mappedTargetEdges;
  std::unordered_map<RIMACS::MappingIndex, size_t> mappedQueryNodes, mappedTargetNodes;
  std::set<RIMACS::MappingPair> mappedAtomPairs;
  for(size_t i = 0; i < mapping.size(); ++i) {
    const RIMACS::MappingPair& mp = mapping[i];
    EXPECT_TRUE(edgeComp.edgesAreCompatible(mp.m_from, mp.m_to));
    const Dimacs::DimacsEdge& queryEdge = q.m_edges.at(mp.m_from);
    const Dimacs::DimacsEdge& targetEdge = t.m_edges.at(mp.m_to);
    mappedQueryEdges.emplace(mp.m_from);
    mappedTargetEdges.emplace(mp.m_to);
    mappedQueryNodes.emplace(queryEdge.m_from, 0); mappedQueryNodes.emplace(queryEdge.m_to, 0);
    mappedTargetNodes.emplace(targetEdge.m_from, 0); mappedTargetNodes.emplace(targetEdge.m_to, 0);
    bool nodesCompatibleFrom_QMapsFrom_T = nodeComp(queryEdge.m_from, targetEdge.m_from)
                                           && nodeComp(queryEdge.m_to, targetEdge.m_to);
    bool nodesCompatibleFrom_QMapsTo_T = nodeComp(queryEdge.m_from, targetEdge.m_to)
                                         && nodeComp(queryEdge.m_to, targetEdge.m_from);
    EXPECT_TRUE(nodesCompatibleFrom_QMapsFrom_T
                || nodesCompatibleFrom_QMapsTo_T);
    int addedQueryAtom = -1, addedTargetAtom = -1;
    for(size_t j = 0; j < mapping.size(); ++j) {
      const RIMACS::MappingPair& toPair = mapping[j];
      int commonQueryNode = get_common_node_idx(q, mp.m_from, toPair.m_from);
      int commonTargetNode = get_common_node_idx(t, mp.m_to, toPair.m_to);
      ASSERT_EQ(commonQueryNode != -1, commonTargetNode != -1);
      if(commonQueryNode != -1) {
        EXPECT_TRUE(nodeComp(commonQueryNode, commonTargetNode));
        int otherQueryNode = get_other_node_idx(queryEdge, commonQueryNode);
        int otherTargetNode = get_other_node_idx(targetEdge, commonTargetNode);
        EXPECT_TRUE(nodeComp(otherQueryNode, otherTargetNode));
        connectivity.at(i).insert(j);
        addedQueryAtom = commonQueryNode, addedTargetAtom = commonTargetNode;
      }
    }
    // store all valid and enumerated atom mappings
    if(addedQueryAtom != -1) {
      mappedAtomPairs.emplace(addedQueryAtom, addedTargetAtom);
      mappedAtomPairs.emplace(get_other_node_idx(queryEdge, addedQueryAtom), get_other_node_idx(targetEdge, addedTargetAtom));
    } else {
      if(nodesCompatibleFrom_QMapsFrom_T) {
        mappedAtomPairs.emplace(queryEdge.m_from, targetEdge.m_from);
      }
      if(nodesCompatibleFrom_QMapsTo_T) {
        mappedAtomPairs.emplace(queryEdge.m_from, targetEdge.m_to);
      }
    }
  }
  for(std::pair<const RIMACS::MappingIndex, size_t>& ap : mappedQueryNodes) {
    auto edges = first.nodeGetEdges(ap.first);
    ap.second =
        std::count_if(edges.begin(), edges.end(),
                      [&mappedQueryEdges](const RIMACS::MappingIndex& e) {return mappedQueryEdges.count(e);});
  }
  for(std::pair<const RIMACS::MappingIndex, size_t>& ap : mappedTargetNodes) {
    auto edges = second.nodeGetEdges(ap.first);
    ap.second =
        std::count_if(edges.begin(), edges.end(),
                      [&mappedTargetEdges](const RIMACS::MappingIndex& e) {return mappedTargetEdges.count(e);});
  }
  for(const RIMACS::MappingPair& ap : mappedAtomPairs) {
    size_t incidentBondsInQueryMoleculeForAtom = get_number_of_incident_bonds(mappedQueryNodes, ap.m_from);
    size_t incidentBondsInTargetMoleculeForAtom = get_number_of_incident_bonds(mappedTargetNodes, ap.m_to);
    EXPECT_EQ(incidentBondsInQueryMoleculeForAtom, incidentBondsInTargetMoleculeForAtom);
    EXPECT_TRUE(nodeComp(ap.m_from, ap.m_to));
  }
  verify_connected_components(connectivity, connectedComponents);
}


void verify_mces_result(
    const Dimacs::DimacsAdjacencyFunctor& first,
    const Dimacs::DimacsAdjacencyFunctor& second,
    const Dimacs::DimacsNodeCompatibilityFunctor& nodeComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& edgeComp,
    const RIMACS::MCSResult& result)
{
  verify_mces_result(first, second, nodeComp, edgeComp, result.getMappings(), result.getComponentSizes());
}

} // namespace
