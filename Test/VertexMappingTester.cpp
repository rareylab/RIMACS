
#include <unordered_set>

#include <gtest/gtest.h>

#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Algorithm/VertexMapping.hpp>

#include "InitializationUtils.hpp"
#include "Dimacs/DimacsGraphFunctor.hpp"


class VertexMappingTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph>>
{
public:
  template<typename VertexMapping>
  static auto& compatible_nodes(VertexMapping&& vm)
  {
    return vm.m_compatibleNodes;
  }
  template<typename VertexMapping>
  static auto& adjacent_nodes(VertexMapping&& vm)
  {
    return vm.m_adjacentNodes;
  }
  template<typename VertexMapping>
  static auto& has_adjacent_nodes(VertexMapping&& vm)
  {
    return vm.m_hasAdjacentNodes;
  }
};

namespace {

template<typename VertexMapping>
auto& adjacent_nodes(
    VertexMapping&& vm)
{
  return VertexMappingTester::adjacent_nodes(vm);
}
template<typename VertexMapping>
auto& has_adjacent_nodes(
    VertexMapping&& vm)
{
  return VertexMappingTester::has_adjacent_nodes(vm);
}
template<typename VertexMapping>
auto& compatible_nodes(
    VertexMapping&& vm)
{
  return VertexMappingTester::compatible_nodes(vm);
}

} // namespace

INSTANTIATE_TEST_SUITE_P(GraphVariations, VertexMappingTester,
    Data::generic_disconnected_test_data());

template<bool weighted>
void test_initialization_tpl(
    const Dimacs::DimacsAdjacencyFunctor& query,
    const Dimacs::DimacsAdjacencyFunctor& target,
    const Dimacs::DimacsNodeCompatibilityFunctor& atomComp)
{
  RIMACS::VertexMapping<weighted> vm(query, target, atomComp);

  EXPECT_TRUE(has_adjacent_nodes(vm));
  EXPECT_EQ(adjacent_nodes(vm).size(), query.getNofNodes());
  EXPECT_EQ(compatible_nodes(vm).size(), query.getNofNodes());
  EXPECT_EQ(RIMACS::getCompatibleNode(vm.getNofExtendableNodes()), query.getNofNodes());
  EXPECT_EQ(RIMACS::getCompatibleNode(vm.getNofConnectedExtendableNodes()), query.getNofNodes());

  for(RIMACS::MappingIndex i = 0; i < query.getNofNodes(); ++i) {
    auto weightedTransform = std::bind(&std::pair<RIMACS::MappingIndex, double>::first,
                                       std::placeholders::_1);
    std::unordered_set<RIMACS::MappingIndex> compatibleNodes;
    std::transform(compatible_nodes(vm).at(i).begin(), compatible_nodes(vm).at(i).end(),
                   std::inserter(compatibleNodes, compatibleNodes.end()),
                   weightedTransform);
    for(RIMACS::MappingIndex otherNode = 0; otherNode < target.getNofNodes(); ++otherNode) {
      EXPECT_EQ(static_cast<bool>(compatibleNodes.count(otherNode)),
               atomComp.nodesAreCompatible(i, otherNode));
    }
  }
}

TEST_P(VertexMappingTester, testInitialization)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor q(query);
  Dimacs::DimacsAdjacencyFunctor t(target);
  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(q, t);

  test_initialization_tpl<true>(q, t, nodeComp);
  test_initialization_tpl<false>(q, t, nodeComp);
}


template<bool weighted>
void test_initialized_access_tpl(
    const Dimacs::DimacsAdjacencyFunctor& query,
    const Dimacs::DimacsAdjacencyFunctor& target,
    const Dimacs::DimacsNodeCompatibilityFunctor& atomComp)
{
  RIMACS::VertexMapping<weighted> vm(query, target, atomComp);
  const RIMACS::VertexMapping<weighted> orgVm(query, target, atomComp);

  std::unordered_set<RIMACS::MappingIndex> usedNodes;
  RIMACS::MappingIndex nextNode = vm.getNextNode();
  do {
    EXPECT_TRUE(usedNodes.insert(nextNode).second);
    const auto& compNodes = compatible_nodes(orgVm).at(nextNode);
    RIMACS::MappingPair indices = vm.initIndexTuple(nextNode);
    RIMACS::MappingTupleType<weighted> mapping;
    size_t idx = 0;
    while(vm.getNextMatch(indices, mapping)) {
      EXPECT_TRUE(idx < compNodes.size());
      EXPECT_EQ(indices.m_from, nextNode);
      EXPECT_EQ(indices.m_to , idx);
      EXPECT_EQ(mapping.m_from, nextNode);
      EXPECT_EQ(mapping.m_to, RIMACS::getCompatibleNode(compNodes.at(idx)));
      ++idx;
    }
    EXPECT_EQ(indices.m_to, RIMACS::MappingIndex{0});
    nextNode = vm.getNextNode();
  } while(nextNode != std::numeric_limits<RIMACS::MappingIndex>::max());
  EXPECT_EQ(usedNodes.size(), query.getNofNodes());
  for(RIMACS::MappingIndex i = 0; i < query.getNofNodes(); ++i) {
    EXPECT_TRUE(static_cast<bool>(usedNodes.count(i)));
  }
}


TEST_P(VertexMappingTester, initializedAccess)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor q(query);
  Dimacs::DimacsAdjacencyFunctor t(target);
  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(q, t);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(q, t);

  test_initialized_access_tpl<true>(q, t, nodeComp);
  test_initialized_access_tpl<false>(q, t, nodeComp);
}


template<bool weighted>
void test_single_update_tpl(
    const Dimacs::DimacsAdjacencyFunctor& query,
    const Dimacs::DimacsAdjacencyFunctor& target,
    const Dimacs::DimacsNodeCompatibilityFunctor& atomComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& bondComp)
{
  RIMACS::VertexMapping<weighted> vm(query, target, atomComp);
  const RIMACS::VertexMapping<weighted> orgVm(query, target, atomComp);
  std::vector<RIMACS::MappingIndex> status(query.getNofNodes(), Mappable);

  RIMACS::MappingIndex selectedNode = adjacent_nodes(vm).at(adjacent_nodes(vm).size() / 2).m_node;
  EXPECT_FALSE(compatible_nodes(vm).at(selectedNode).empty());
  RIMACS::MappingIndex selectedCompatibleNode = get_compatible_node_idx(compatible_nodes(vm), selectedNode);
  adjacent_nodes(vm) = {{selectedNode,  1, 0, 0.0, 0.0}};
  status.at(selectedNode) = 0;

  EXPECT_EQ(vm.getNextNode(), selectedNode);
  EXPECT_EQ(vm.getNextNode(), std::numeric_limits<RIMACS::MappingIndex>::max());

  std::vector<RIMACS::MappingIndex> newUnmappableNodes;
  RIMACS::VertexMapping<weighted> updatedMapping(vm, query, target, bondComp,
                                                 RIMACS::MappingPair(selectedNode, selectedCompatibleNode),
                                                 status, newUnmappableNodes);
  EXPECT_EQ(adjacent_nodes(updatedMapping).size(), query.nodeGetEdges(selectedNode).size());
  std::unordered_set<RIMACS::MappingIndex> adjNodes;
  for(const RIMACS::MappingIndex& incidentEdge : query.nodeGetEdges(selectedNode)) {
    RIMACS::MappingIndex adjacentNode = query.edgeGetFromId(incidentEdge) != selectedNode
                                   ? query.edgeGetFromId(incidentEdge) : query.edgeGetToId(incidentEdge);
    adjNodes.insert(adjacentNode);
    auto it = std::find_if(adjacent_nodes(updatedMapping).begin(),
                           adjacent_nodes(updatedMapping).end(),
                           [adjacentNode] (const RIMACS::AdjacentNode& adn) {
                             return adn.m_node == adjacentNode;
                           });
    EXPECT_TRUE(it != adjacent_nodes(updatedMapping).end() && "Each adjacent node must be present");
  }
  for(RIMACS::MappingIndex nodeId = 0; nodeId < query.getNofNodes(); ++nodeId) {
    const auto& originalCompNodes = compatible_nodes(orgVm).at(nodeId);
    const auto& updatedCompNodes = compatible_nodes(updatedMapping).at(nodeId);
    if(adjNodes.count(nodeId)) {
      auto weightedTransform = std::bind(&std::pair<RIMACS::MappingIndex, double>::first,
                                         std::placeholders::_1);
      std::unordered_set<RIMACS::MappingIndex> compNodes;
      std::transform(updatedCompNodes.begin(), updatedCompNodes.end(),
                     std::inserter(compNodes, compNodes.end()), weightedTransform);
      for(const auto& otherCompNode : originalCompNodes) {
        RIMACS::MappingIndex otherCompNodeIdx = RIMACS::getCompatibleNode(otherCompNode) ;
        // After update, only nodes fulfilling all adjacency relations are compatible
        bool adjacent = target.adjacent(selectedCompatibleNode, otherCompNodeIdx);
        EXPECT_EQ(static_cast<bool>(compNodes.count(otherCompNodeIdx)), adjacent);
      }
    } else if(nodeId == selectedNode) {
      // already mapped nodes have no compatible nodes
      EXPECT_TRUE(compatible_nodes(updatedMapping).at(selectedNode).empty());
    } else {
      for(auto upCompIt = updatedCompNodes.begin(), orgCompIt = originalCompNodes.begin()
          ; upCompIt != updatedCompNodes.end()
          ; ++orgCompIt) {
        EXPECT_TRUE(orgCompIt != originalCompNodes.end()
                    && "Each compatible node of the updated mapping must be part of the original");
        if(*upCompIt == *orgCompIt) {
          EXPECT_FALSE(target.adjacent(selectedCompatibleNode, RIMACS::getCompatibleNode(*orgCompIt)));
          ++upCompIt;
        } else {
          EXPECT_TRUE(RIMACS::getCompatibleNode(*orgCompIt) == selectedCompatibleNode
                  || target.adjacent(selectedCompatibleNode, RIMACS::getCompatibleNode(*orgCompIt)));
        }
      }
    }
  }
}


TEST_P(VertexMappingTester, singleUpdate)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor q(query);
  Dimacs::DimacsAdjacencyFunctor t(target);
  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(q, t);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(q, t);

  test_single_update_tpl<true>(q, t, nodeComp, edgeComp);
  test_single_update_tpl<false>(q, t, nodeComp, edgeComp);
}



template<bool weighted>
void test_disconnected_initialization_tpl(
    const Dimacs::DimacsAdjacencyFunctor& query,
    const Dimacs::DimacsAdjacencyFunctor& target,
    const Dimacs::DimacsNodeCompatibilityFunctor& atomComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& bondComp)
{
  RIMACS::VertexMapping<weighted> vm(query, target, atomComp);
  const RIMACS::VertexMapping<weighted> orgVm(query, target, atomComp);
  std::vector<RIMACS::MappingIndex> status(query.getNofNodes(), Mappable);

  RIMACS::MappingIndex selectedNode = adjacent_nodes(vm).at(adjacent_nodes(vm).size() / 2).m_node;
  status.at(selectedNode) = 0;
  EXPECT_FALSE(compatible_nodes(vm).at(selectedNode).empty());

  RIMACS::MappingIndex selectedCompatibleNode = get_compatible_node_idx(compatible_nodes(vm), selectedNode);
  adjacent_nodes(vm) = {{selectedNode, 1, 0, 0.0, 0.0}};

  vm.finishedNode(vm.getNextNode());
  EXPECT_TRUE(compatible_nodes(vm).at(selectedNode).empty());
  EXPECT_EQ(vm.getNextNode(), std::numeric_limits<RIMACS::MappingIndex>::max());

  // Update mapping
  std::vector<RIMACS::MappingIndex> newUnmappableNodes;
  vm = RIMACS::VertexMapping<weighted>(vm, query, target, bondComp,
                                       RIMACS::MappingPair(selectedNode, selectedCompatibleNode), status,
                                       newUnmappableNodes);
  for(RIMACS::MappingIndex nextNode = vm.getNextNode()
      ; nextNode != std::numeric_limits<RIMACS::MappingIndex>::max()
      ; nextNode = vm.getNextNode()) {
    vm.finishedNode(nextNode);
  }
  vm.prepareDisconnected(query);

  std::unordered_set<RIMACS::MappingIndex> usedNodes;
  RIMACS::MappingIndex nextNode = vm.getNextNode();
  do {
    EXPECT_TRUE(usedNodes.insert(nextNode).second);
    const auto& originalCompNodes = compatible_nodes(orgVm).at(nextNode);
    const auto& updatedCompNodes = compatible_nodes(vm).at(nextNode);
    for(auto orgCompIt = originalCompNodes.begin(), upCompIt = updatedCompNodes.begin()
        ; upCompIt != updatedCompNodes.end(); ++orgCompIt) {
      // Each node of the updated VM must be processed before the original vector might be at end.
      EXPECT_TRUE(orgCompIt != originalCompNodes.end());
      if(*orgCompIt == *upCompIt) {
        EXPECT_FALSE(target.adjacent(selectedCompatibleNode, RIMACS::getCompatibleNode(*orgCompIt)));
        ++upCompIt;
      } else {
        EXPECT_TRUE(RIMACS::getCompatibleNode(*orgCompIt) == selectedCompatibleNode
                || target.adjacent(selectedCompatibleNode, RIMACS::getCompatibleNode(*orgCompIt)));
      }
    }
    nextNode = vm.getNextNode();
  } while(nextNode != std::numeric_limits<RIMACS::MappingIndex>::max());
  EXPECT_TRUE(usedNodes.size() <= query.getNofNodes() - query.nodeGetEdges(selectedNode).size() - 1);
  for(RIMACS::MappingIndex i = 0; i < query.getNofNodes(); ++i) {
    if(!usedNodes.count(i)) {
      EXPECT_TRUE(compatible_nodes(vm).at(i).empty()
                  && "Processed nodes do not need compatible node entries");
    }
  }
}

TEST_P(VertexMappingTester, disconnectedInitialization)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor q(query);
  Dimacs::DimacsAdjacencyFunctor t(target);
  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(q, t);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(q, t);

  test_disconnected_initialization_tpl<true>(q, t, nodeComp, edgeComp);
  test_disconnected_initialization_tpl<false>(q, t, nodeComp, edgeComp);
}
