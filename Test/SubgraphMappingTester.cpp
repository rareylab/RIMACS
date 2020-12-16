
#include <gtest/gtest.h>

#include <RIMACS/Algorithm/SubgraphMapping.hpp>

#include "Dimacs/DimacsGraph.hpp"
#include "Dimacs/DimacsGraphFunctor.hpp"

#include <RIMACS/Logger.hpp>

#include "InitializationUtils.hpp"


class SubgraphMappingTester
     : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph>> {
public:
  template<typename SubgraphMapping>
  static auto& node_mappings(SubgraphMapping&& mapping)
  {
    return mapping.m_nodeMappings;
  }
  template<typename SubgraphMapping>
  static auto& node_status(SubgraphMapping&& mapping)
  {
    return mapping.m_nodeStatus;
  }
  template<typename SubgraphMapping>
  static auto& connected_component_sizes(SubgraphMapping&& mapping)
  {
    return mapping.m_connectedComponentSizes;
  }

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
};


namespace {

template<typename SubgraphMapping>
auto& node_mappings(
    SubgraphMapping&& vm)
{
  return SubgraphMappingTester::node_mappings(vm);
}
template<typename SubgraphMapping>
auto& node_status(
    SubgraphMapping&& vm)
{
  return SubgraphMappingTester::node_status(vm);
}
template<typename SubgraphMapping>
auto& connected_component_sizes(
    SubgraphMapping&& vm)
{
  return SubgraphMappingTester::connected_component_sizes(vm);
}
template<typename VertexMapping>
auto& adjacent_nodes(
    VertexMapping&& vm)
{
  return SubgraphMappingTester::adjacent_nodes(vm);
}
template<typename VertexMapping>
auto& compatible_nodes(
    VertexMapping&& vm)
{
  return SubgraphMappingTester::compatible_nodes(vm);
}

} // namespace


INSTANTIATE_TEST_SUITE_P(GraphVariations, SubgraphMappingTester,
                          Data::generic_test_data());


template<bool weighted>
void test_initialization_tpl(
    const Dimacs::DimacsAdjacencyFunctor& query,
    const Dimacs::DimacsAdjacencyFunctor& target,
    const Dimacs::DimacsNodeCompatibilityFunctor& atomComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& bondComp,
    const RIMACS::Config& config)
{
  RIMACS::SubgraphMapping<weighted> mapping(query, target, atomComp, bondComp, config);

  EXPECT_EQ(mapping.getCount(), size_t{0});
  EXPECT_TRUE(mapping.moveMCSResults().empty());
  EXPECT_FALSE(mapping.recursionLimitReached());
  EXPECT_TRUE(node_mappings(mapping).empty());
  EXPECT_EQ(node_status(mapping).size(), query.getNofNodes());
  EXPECT_EQ(std::count(node_status(mapping).begin(), node_status(mapping).end(), Mappable),
           query.getNofNodes());
  EXPECT_EQ(connected_component_sizes(mapping).size(), size_t{1});
  EXPECT_EQ(connected_component_sizes(mapping).front(), size_t{0});
}

TEST_P(SubgraphMappingTester,testInitialization)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor q(query);
  Dimacs::DimacsAdjacencyFunctor t(target);
  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(q, t);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(q, t);
  const auto& config = RIMACS::getDefaultConfig();

  test_initialization_tpl<true>(q, t, nodeComp, edgeComp, config);
  test_initialization_tpl<false>(q, t, nodeComp, edgeComp, config);
}


template<bool weighted>
void test_extend_mapping_tpl(
    const Dimacs::DimacsAdjacencyFunctor& query,
    const Dimacs::DimacsAdjacencyFunctor& target,
    const Dimacs::DimacsNodeCompatibilityFunctor& atomComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& bondComp,
    const RIMACS::Config& config)
{
  RIMACS::SubgraphMapping<weighted> mapping(query, target, atomComp, bondComp, config);
  RIMACS::VertexMapping<weighted> vm(query, target, atomComp);

  RIMACS::MappingIndex selectedNode = adjacent_nodes(vm).at(adjacent_nodes(vm).size() / 2).m_node;
  EXPECT_FALSE(compatible_nodes(vm).at(selectedNode).empty());
  const RIMACS::MappingIndex selectedCompatibleNode = get_compatible_node_idx(compatible_nodes(vm), selectedNode);

  std::vector<RIMACS::MappingIndex> newUnmappableNodes;
  vm = mapping.addMatchUpdateVM(
        vm,
        RIMACS::MappingTuple(selectedNode, selectedCompatibleNode, 0.0),
        newUnmappableNodes);

  EXPECT_TRUE(node_status(mapping).at(selectedNode) != Mappable);
  EXPECT_EQ(std::count(node_status(mapping).begin(), node_status(mapping).end(), Mappable),
           query.getNofNodes() - 1);
  EXPECT_EQ(connected_component_sizes(mapping).size(), size_t{1});
  EXPECT_EQ(connected_component_sizes(mapping).front(), size_t{1});
  EXPECT_EQ(node_mappings(mapping).size(), size_t{1});
  EXPECT_EQ(node_mappings(mapping).front(),
            RIMACS::MappingPair(selectedNode, selectedCompatibleNode));
  EXPECT_EQ(mapping.getCount(), size_t{1});
  EXPECT_TRUE(mapping.isExtendable(vm));
  EXPECT_FALSE(mapping.canBeDisconnectedExpanded(vm));
}


TEST_P(SubgraphMappingTester, extendMapping)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor q(query);
  Dimacs::DimacsAdjacencyFunctor t(target);
  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(q, t);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(q, t);
  const auto& config = RIMACS::getDefaultConfig();

  test_extend_mapping_tpl<true>(q, t, nodeComp, edgeComp, config);
  test_extend_mapping_tpl<false>(q, t, nodeComp, edgeComp, config);
}


template<bool weighted>
void test_remove_mapping_tpl(
    const Dimacs::DimacsAdjacencyFunctor& query,
    const Dimacs::DimacsAdjacencyFunctor& target,
    const Dimacs::DimacsNodeCompatibilityFunctor& atomComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& bondComp,
    const RIMACS::Config& config)
{
  RIMACS::SubgraphMapping<weighted> mapping(query, target, atomComp, bondComp, config);
  RIMACS::VertexMapping<weighted> vm(query, target, atomComp);

  RIMACS::MappingIndex selectedNode = adjacent_nodes(vm).at(adjacent_nodes(vm).size() / 2).m_node;
  EXPECT_FALSE(compatible_nodes(vm).at(selectedNode).empty());
  const RIMACS::MappingIndex selectedCompatibleNode = get_compatible_node_idx(compatible_nodes(vm), selectedNode);


  std::vector<RIMACS::MappingIndex> newUnmappableNodes;
  mapping.addMatchUpdateVM(
        vm,
        RIMACS::MappingTuple(selectedNode, selectedCompatibleNode, 0.0),
        newUnmappableNodes);

  EXPECT_TRUE(mapping.isExtendable(vm));
  EXPECT_FALSE(mapping.canBeDisconnectedExpanded(vm));

  mapping.removeMapping(newUnmappableNodes);
  // node is mappable again, after mapping is removed
  EXPECT_EQ(node_status(mapping).at(selectedNode), Mappable);
  EXPECT_EQ(std::count(node_status(mapping).begin(), node_status(mapping).end(), Mappable),
            query.getNofNodes());
  EXPECT_EQ(connected_component_sizes(mapping).size(), size_t{1});
  EXPECT_EQ(connected_component_sizes(mapping).front(), size_t{0});
  EXPECT_TRUE(node_mappings(mapping).empty());
  EXPECT_EQ(mapping.getCount(), size_t{1});
  EXPECT_TRUE(mapping.isExtendable(vm));
  EXPECT_FALSE(mapping.canBeDisconnectedExpanded(vm));
}

TEST_P(SubgraphMappingTester, removeMapping)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor q(query);
  Dimacs::DimacsAdjacencyFunctor t(target);
  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(q, t);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(q, t);
  const auto& config = RIMACS::getDefaultConfig();

  test_remove_mapping_tpl<true>(q, t, nodeComp, edgeComp, config);
  test_remove_mapping_tpl<false>(q, t, nodeComp, edgeComp, config);
}


template<bool weighted>
void test_disconnected_extension_tpl(
    const Dimacs::DimacsAdjacencyFunctor& query,
    const Dimacs::DimacsAdjacencyFunctor& target,
    const Dimacs::DimacsNodeCompatibilityFunctor& atomComp,
    const Dimacs::DimacsEdgeCompatibilityFunctor& bondComp,
    const RIMACS::Config& config)
{
  RIMACS::SubgraphMapping<weighted> mapping(query, target, atomComp, bondComp, config);
  RIMACS::VertexMapping<weighted> vm(query, target, atomComp);

  RIMACS::MappingIndex selectedNode = adjacent_nodes(vm).at(adjacent_nodes(vm).size() / 2).m_node;
  EXPECT_FALSE(compatible_nodes(vm).at(selectedNode).empty());
  const RIMACS::MappingIndex selectedCompatibleNode = get_compatible_node_idx(compatible_nodes(vm), selectedNode);

  std::vector<RIMACS::MappingIndex> locallyUnmappableNodes;
  vm = mapping.addMatchUpdateVM(
        vm,
        RIMACS::MappingTuple(selectedNode, selectedCompatibleNode, 0.0),
        locallyUnmappableNodes);
  for(const RIMACS::AdjacentNode& adn : adjacent_nodes(vm)) {
    vm.finishedNode(adn.m_node);
  }
  adjacent_nodes(vm).clear();

  EXPECT_TRUE(mapping.canBeDisconnectedExpanded(vm));

  mapping.prepareDisconnectionUpdateVM(vm);

  EXPECT_EQ(connected_component_sizes(mapping).size(), size_t{2});
  EXPECT_EQ(connected_component_sizes(mapping).front(), size_t{0});
  EXPECT_EQ(connected_component_sizes(mapping).back(), size_t{1});

  RIMACS::MappingPair indices = vm.initIndexTuple(vm.getNextNode());
  RIMACS::MappingTupleType<weighted> nodeMapping;
  EXPECT_TRUE(vm.getNextMatch(indices, nodeMapping));

  std::vector<RIMACS::MappingIndex> newUnmappableNodes;
  mapping.addMatchUpdateVM(vm, nodeMapping, newUnmappableNodes);

  EXPECT_EQ(node_mappings(mapping).size(), size_t{2});
  EXPECT_EQ(connected_component_sizes(mapping).size(), size_t{2});
  EXPECT_EQ(connected_component_sizes(mapping).front(), size_t{1});
  EXPECT_EQ(connected_component_sizes(mapping).back(), size_t{1});

  mapping.removeMapping(newUnmappableNodes);

  EXPECT_EQ(node_mappings(mapping).size(), size_t{1});
  EXPECT_EQ(connected_component_sizes(mapping).size(), size_t{2});
  EXPECT_EQ(connected_component_sizes(mapping).front(), size_t{0});
  EXPECT_EQ(connected_component_sizes(mapping).back(), size_t{1});

  mapping.resetDisconnection();

  EXPECT_EQ(connected_component_sizes(mapping).size(), size_t{1});
  EXPECT_EQ(connected_component_sizes(mapping).front(), size_t{1});
}


TEST_P(SubgraphMappingTester, extendDisconnected)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor q(query);
  Dimacs::DimacsAdjacencyFunctor t(target);
  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(q, t);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(q, t);
  const auto& config = RIMACS::getDisconnectedConfig(RIMACS::Config::ResultType::Maximum, 2, 1);


  test_disconnected_extension_tpl<true>(q, t, nodeComp, edgeComp, config);
  test_disconnected_extension_tpl<false>(q, t, nodeComp, edgeComp, config);
}
