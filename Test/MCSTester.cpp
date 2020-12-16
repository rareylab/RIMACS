
#include <numeric>

#include <gtest/gtest.h>

#include <RIMACS/MCSRunner.hpp>
#include <RIMACS/Algorithm/CompareUtils.hpp>

#include "ResultVerification.hpp"
#include "InitializationUtils.hpp"


class SimpleMCSTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph, size_t, size_t>> {};

TEST_P(SimpleMCSTester, connectedMCS)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;
  size_t size, nofResults;

  std::tie(query, target, size, nofResults) = GetParam();

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  RIMACS::Config conf = RIMACS::getConnectedConfig(RIMACS::Config::ResultType::AllMaximum);
  auto result = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, conf);

  EXPECT_EQ(result.front().size(), size);
  EXPECT_EQ(result.size(), nofResults);

  std::vector<RIMACS::MappingIndex> expectedComponent;
  expectedComponent.push_back(size);
  for(const RIMACS::MCSResult& res : result) {
    verify_result(queryAdj, targetAdj, nodeComp, edgeComp, res.getMappings(),expectedComponent);
  }

  // redo the test for the opposite molecule order
  Dimacs::DimacsNodeCompatibilityFunctor reverseNodeComp(targetAdj, queryAdj);
  Dimacs::DimacsEdgeCompatibilityFunctor reverseEdgeComp(targetAdj, queryAdj);
  result = RIMACS::MCSRunner::runGenericMCS(targetAdj, queryAdj, reverseNodeComp, reverseEdgeComp, conf);

  EXPECT_EQ(result.front().size(), size);
  EXPECT_EQ(result.size(), nofResults);

  for(const RIMACS::MCSResult& res : result) {
    verify_result(targetAdj, queryAdj, reverseNodeComp, reverseEdgeComp, res.getMappings(),expectedComponent);
  }
}

INSTANTIATE_TEST_SUITE_P(Graphs, SimpleMCSTester,
    testing::Values(std::make_tuple(Data::SimpleThreeMemberChain, Data::SimpleThreeMemberChain,
                                    size_t{3}, size_t{2}),

                    std::make_tuple(Data::SimpleThreeMemberChain, Data::IsoThreeMemberChain,
                                    size_t{3}, size_t{2}),

                    std::make_tuple(Data::MixedThreeMemberChain, Data::MixedThreeMemberChain,
                                    size_t{3}, size_t{2}),

                    std::make_tuple(Data::get_simple_chain(6, 6, 8), Data::get_simple_chain(6, 6, 7, 8),
                                    size_t{2}, size_t{2}),

                    std::make_tuple(Data::get_graph("dMCSdemo_1"), Data::get_graph("dMCSdemo_2"),
                                    size_t{10}, size_t{72}),

                    std::make_tuple(Data::get_graph("caffeine"), Data::get_graph("theobromine"),
                                    size_t{13}, size_t{1}),

                    std::make_tuple(Data::get_graph("TestCases_1_f"), Data::get_graph("TestCases_1_t"),//
                                    size_t{14}, size_t{1}),

                    std::make_tuple(Data::get_graph("TestCases_2_f"), Data::get_graph("TestCases_2_t"),
                                    size_t{7}, size_t{6})));


class SimpleDMCSTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph, size_t, size_t>> {};

TEST_P(SimpleDMCSTester, defaultDisconnectedMCS)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;
  size_t size, nofComponents;

  std::tie(query, target, size, nofComponents) = GetParam();

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  RIMACS::Config conf = RIMACS::getDefaultDisconnectedConfig();

  auto result = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, conf);

  ASSERT_EQ(result.size(), 1);
  EXPECT_EQ(result.front().size(), size);
  EXPECT_EQ(result.front().getComponentSizes().size(), nofComponents);

  for(const RIMACS::MCSResult& res : result) {
    verify_result(queryAdj, targetAdj, nodeComp, edgeComp, res.getMappings(),res.getComponentSizes());
  }

  // redo the test for the opposite molecule order
  Dimacs::DimacsNodeCompatibilityFunctor reverseNodeComp(targetAdj, queryAdj);
  Dimacs::DimacsEdgeCompatibilityFunctor reverseEdgeComp(targetAdj, queryAdj);
  result = RIMACS::MCSRunner::runGenericMCS(targetAdj, queryAdj, reverseNodeComp, reverseEdgeComp, conf);

  ASSERT_EQ(result.size(), 1);
  EXPECT_EQ(result.front().size(), size);
  EXPECT_EQ(result.front().getComponentSizes().size(), nofComponents);

  for(const RIMACS::MCSResult& res : result) {
    verify_result(targetAdj, queryAdj, reverseNodeComp, reverseEdgeComp, res.getMappings(), res.getComponentSizes());
  }
}

INSTANTIATE_TEST_SUITE_P(Graphs, SimpleDMCSTester,
    testing::Values(std::make_tuple(Data::SimpleThreeMemberChain, Data::SimpleThreeMemberChain,
                                    size_t{3}, size_t{1}),

                    std::make_tuple(Data::get_simple_chain(6, 8, 6), Data::get_simple_chain(6, 7, 6),
                                    size_t{1}, size_t{1}),

                    std::make_tuple(Data::get_simple_chain(6, 6, 8, 6), Data::get_simple_chain(6, 7, 6, 6),
                                    size_t{2}, size_t{1}),

                    std::make_tuple(Data::get_simple_chain(6, 6, 8, 6, 6), Data::get_simple_chain(6, 6, 7, 6, 6),
                                    size_t{2}, size_t{1}),

                    std::make_tuple(Data::get_simple_chain(6, 6, 6, 8, 6, 6), Data::get_simple_chain(6, 6, 7, 6, 6, 6),
                                    size_t{3}, size_t{1}),

                    std::make_tuple(Data::get_simple_chain(6, 6, 6, 8, 6, 6, 6), Data::get_simple_chain(6, 6, 6, 7, 6, 6, 6),
                                    size_t{6}, size_t{2}),

                    std::make_tuple(Data::get_simple_chain(6, 16, 6, 8, 6, 15, 6), Data::get_simple_chain(6, 15, 6, 7, 6, 16, 6),
                                    size_t{6}, size_t{2}),

                    std::make_tuple(Data::get_graph("dMCSdemo_1"), Data::get_graph("dMCSdemo_2"),
                                    size_t{11}, size_t{3}),

                    std::make_tuple(Data::get_graph("dMCSdemo_1"), Data::get_graph("dMCSdemo_3"),
                                    size_t{11}, size_t{3}),

                    std::make_tuple(Data::get_graph("dMCSdemo_2"), Data::get_graph("dMCSdemo_1"),
                                    size_t{11}, size_t{3}),

                    std::make_tuple(Data::get_graph("dMCSdemo_3"), Data::get_graph("dMCSdemo_1"),
                                    size_t{11}, size_t{3}),

                    std::make_tuple(Data::get_graph("caffeine"), Data::get_graph("theobromine"),
                                    size_t{13}, size_t{1}),

                    std::make_tuple(Data::get_graph("graphs_1"), Data::get_graph("graphs_2"),
                                    size_t{8}, size_t{1}),

                    std::make_tuple(Data::get_graph("graphs_3"), Data::get_graph("graphs_4"),
                                    size_t{9}, size_t{3})));



namespace {

int ring_dfs(
    const Dimacs::DimacsAdjacencyFunctor& graph,
    std::vector<int>& discovery,
    std::vector<bool>& ringMember,
    RIMACS::MappingIndex current,
    RIMACS::MappingIndex parent,
    size_t& time)
{
  discovery[current] = time;
  int lowpoint = std::numeric_limits<int>::max();
  for(RIMACS::MappingIndex neighbour : graph.nodeGetNeighbours(current)) {
    if(neighbour == parent) {
      continue;
    } else if(discovery[neighbour] > 0) {
      lowpoint = std::min(lowpoint, discovery[neighbour]);
    } else {
      lowpoint = std::min(ring_dfs(graph, discovery, ringMember, neighbour, current, ++time), lowpoint);
    }
  }
  ringMember[current] = lowpoint <= discovery[current];
  return lowpoint;
}

class ring_weighter: public Dimacs::DimacsNodeCompatibilityFunctor {
public:
  ring_weighter(
      const Dimacs::DimacsAdjacencyFunctor& query,
      const Dimacs::DimacsAdjacencyFunctor& target)
    : DimacsNodeCompatibilityFunctor(query, target)
  {
    m_queryRingMember.resize(query.getNofNodes());
    m_targetRingMember.resize(target.getNofNodes());
    std::vector<int> discovery(std::max(query.getNofNodes(), target.getNofNodes()), -1);
    size_t time = 1;
    ring_dfs(query, discovery, m_queryRingMember, 0, std::numeric_limits<RIMACS::MappingIndex>::max(), time);
    discovery.clear();
    discovery.resize(target.getNofNodes(), -1);
    ring_dfs(target, discovery, m_targetRingMember, 0, std::numeric_limits<RIMACS::MappingIndex>::max(), time);
  }

  double getWeight(
      RIMACS::MappingIndex queryNode,
      RIMACS::MappingIndex targetNode) const override
  {
    RIMACS_ASSERT(this->nodesAreCompatible(queryNode, targetNode));
    return (this->m_queryRingMember.at(queryNode)
            && this->m_targetRingMember.at(targetNode)) ? 10.0 : 1.0;
  }

private:
  std::vector<bool> m_queryRingMember;
  std::vector<bool> m_targetRingMember;
};

class label_weighter: public Dimacs::DimacsNodeCompatibilityFunctor {
public:
  label_weighter(
      const Dimacs::DimacsAdjacencyFunctor& query,
      const Dimacs::DimacsAdjacencyFunctor& target)
    : DimacsNodeCompatibilityFunctor(query, target)
  {}

  double getWeight(
      RIMACS::MappingIndex queryNode,
      RIMACS::MappingIndex targetNode) const override
  {
    RIMACS_ASSERT(this->nodesAreCompatible(queryNode, targetNode));
    return static_cast<Dimacs::DimacsGraph>(this->m_queryGraph).m_nodeLabels.at(queryNode);
  }
};

} // namespace

class SimpleWeightedMCSTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph, size_t, size_t, size_t>> {};

TEST_P(SimpleWeightedMCSTester, weightedMCS)
{

  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;
  size_t disconnectedParam, size, labelledSize;

  std::tie(query, target, disconnectedParam, size, labelledSize) = GetParam();

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  label_weighter weightedNodeComp(queryAdj, targetAdj);


  {
    ring_weighter weightedNodeComp(queryAdj, targetAdj);
    auto res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, weightedNodeComp, edgeComp,
                                             RIMACS::getDefaultConnectedConfig(), true);
    ASSERT_EQ(res.size(), size_t{1});
    EXPECT_EQ(res.front().size(), size);

    res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, weightedNodeComp, edgeComp,
                                        RIMACS::getDefaultConnectedConfig(),
                                        true, nullptr, nullptr, true);
    ASSERT_EQ(res.size(), size_t{1});
    EXPECT_EQ(res.front().size(), size);
    verify_result(queryAdj, targetAdj, weightedNodeComp, edgeComp, res.front());
  }
  {
    RIMACS::Config config
        = RIMACS::getDisconnectedConfig(RIMACS::Config::ResultType::Maximum,
                                     disconnectedParam, disconnectedParam);

    auto res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, weightedNodeComp, edgeComp, config, true);
    ASSERT_EQ(res.size(), size_t{1});
    EXPECT_EQ(res.front().size(), labelledSize);

    res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, weightedNodeComp, edgeComp, config,
                                        true, nullptr, nullptr, true);
    ASSERT_EQ(res.size(), size_t{1});
    EXPECT_EQ(res.front().size(), labelledSize);
    verify_result(queryAdj, targetAdj, weightedNodeComp, edgeComp, res.front());
  }
}

INSTANTIATE_TEST_SUITE_P(Graphs, SimpleWeightedMCSTester,
    testing::Values(std::make_tuple(Data::get_graph("ring_and_chain"),
                                    Data::get_graph("chain_and_ring"),
                                    size_t{2}, size_t{3}, size_t{9}),

                    std::make_tuple(Data::get_graph("ring_and_chains_1"),
                                    Data::get_graph("chains_and_ring_1"),
                                    size_t{2}, size_t{3}, size_t{8}),

                    std::make_tuple(Data::get_graph("ring_and_chains_2"),
                                    Data::get_graph("chains_and_ring_2"),
                                    size_t{2}, size_t{3}, size_t{4}),

                    std::make_tuple(Data::get_graph("ring_and_chains_3"),
                                    Data::get_graph("chains_and_ring_3"),
                                    size_t{42}, size_t{3}, size_t{4}),

                    std::make_tuple(Data::get_graph("ring_and_chains_4"),
                                    Data::get_graph("chains_and_ring_4"),
                                    size_t{42}, size_t{4}, size_t{3})
));

class AcceptanceTester : public testing::TestWithParam<std::tuple<std::string, std::string, size_t, size_t, size_t>> {};

TEST_P(AcceptanceTester, diverseTests)
{
  std::string queryFile, targetFile;
  size_t nofConnectedResults, connected_size, disconnected_size;

  std::tie(queryFile, targetFile, nofConnectedResults, connected_size, disconnected_size) = GetParam();

  Dimacs::DimacsGraph query = Data::get_graph(queryFile);
  Dimacs::DimacsGraph target = Data::get_graph(targetFile);

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  RIMACS::MCSResults res;
  res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp,
                                      RIMACS::getConnectedConfig(RIMACS::Config::ResultType::AllMaximum));

  ASSERT_GE(res.size(), 1);
  EXPECT_EQ(res.size(), nofConnectedResults);
  for(const RIMACS::MCSResult& result : res) {
    ASSERT_EQ(result.getComponentSizes().size(), 1);
    EXPECT_EQ(result.size(), connected_size);
    EXPECT_EQ(result.getComponentSizes().front(), connected_size);
    verify_result(queryAdj, targetAdj, nodeComp, edgeComp, result.getMappings(), result.getComponentSizes());
  }

  res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp,
                                      RIMACS::getDisconnectedConfig(RIMACS::Config::ResultType::Maximum));

  ASSERT_EQ(res.size(), 1);
  for(const RIMACS::MCSResult& result : res) {
    EXPECT_EQ(result.size(), disconnected_size);
    ASSERT_GE(result.getComponentSizes().size(), 1);
    EXPECT_LE(result.getComponentSizes().size(), RIMACS::Config::Defaults::NofComponents);
    verify_result(queryAdj, targetAdj, nodeComp, edgeComp, result.getMappings(), result.getComponentSizes());
  }
}

INSTANTIATE_TEST_SUITE_P(Graphs, AcceptanceTester,
                         testing::Values(std::make_tuple("graphs_5", "graphs_6", size_t{8}, size_t{7}, size_t{11}),
                    std::make_tuple("graphs_7", "graphs_6", size_t{3}, size_t{10}, size_t{14}),
                    std::make_tuple("graphs_8", "graphs_6", size_t{5}, size_t{9}, size_t{13}),
                    std::make_tuple("graphs_9", "caffeine", size_t{4}, size_t{2}, size_t{2}),
                    std::make_tuple("graphs_10", "graphs_11", size_t{6}, size_t{7}, size_t{12})));


INSTANTIATE_TEST_SUITE_P(OnTestCases, AcceptanceTester,
                         testing::Values(std::make_tuple("TestCases_3_f", "TestCases_3_t", size_t{2}, size_t{7}, size_t{15}),
                    std::make_tuple("TestCases_4_f", "TestCases_4_t", size_t{6}, size_t{7}, size_t{15}),
                    std::make_tuple("TestCases_5_f", "TestCases_5_t", size_t{2}, size_t{6}, size_t{10}),
                    std::make_tuple("TestCases_6_f", "TestCases_6_t", size_t{2}, size_t{6}, size_t{12}),
                    std::make_tuple("TestCases_7_f", "TestCases_7_t", size_t{2}, size_t{12}, size_t{20}),
                    std::make_tuple("TestCases_8_f", "TestCases_8_t", size_t{2}, size_t{6}, size_t{11}),
                    std::make_tuple("TestCases_9_f", "TestCases_9_t", size_t{10}, size_t{7}, size_t{17}),
                    std::make_tuple("TestCases_10_f", "TestCases_10_t", size_t{12}, size_t{6}, size_t{9}),
                    std::make_tuple("TestCases_11_f", "TestCases_11_t", size_t{1}, size_t{9}, size_t{18}),
                    std::make_tuple("TestCases_12_f", "TestCases_12_t", size_t{18}, size_t{8}, size_t{17})));


TEST(MCS_Tester, equivalenceClassDemo)
{
  Dimacs::DimacsGraph g = Data::get_simple_cycle(6, 6, 6, 6, 6, 6);
  INITIALIZE_DEFAULT_NAMED_FUNCTORS(g, g);

  RIMACS::Config conf = RIMACS::getConnectedConfig(RIMACS::Config::ResultType::AllMaximum);
  auto res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, conf);
  EXPECT_EQ(res.size(), 12);

  std::vector<RIMACS::MappingIndex> equivalence_classes(6, 0);
  RIMACS::Config equivalence_config = RIMACS::getConnectedConfig(RIMACS::Config::ResultType::AllMaximum, 3);
  res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, equivalence_config, false,
                                        &equivalence_classes, &equivalence_classes);

  EXPECT_EQ(res.size(), 3);

  res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, equivalence_config, false,
                                        &equivalence_classes, nullptr);

  EXPECT_EQ(res.size(), 3);

  res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, equivalence_config, false,
                                        nullptr, &equivalence_classes);

  EXPECT_EQ(res.size(), 3);
}


class EquivalenceClassTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph,
                                               std::vector<RIMACS::MappingIndex>, std::vector<RIMACS::MappingIndex>,
                                               int, int, int>> {};

TEST_P(EquivalenceClassTester, testClasses)
{
  Dimacs::DimacsGraph query, target;
  std::vector<RIMACS::MappingIndex> queryEquivalenceClasses, targetEquivalenceClasses;
  int maxNofEquivalentResults, nofEquivalentResults, size;

  std::tie(query, target,
           queryEquivalenceClasses, targetEquivalenceClasses,
           maxNofEquivalentResults, nofEquivalentResults, size) = GetParam();

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  ASSERT_EQ(query.m_nodeLabels.size(), queryEquivalenceClasses.size());
  ASSERT_EQ(target.m_nodeLabels.size(), targetEquivalenceClasses.size());

  RIMACS::Config conf
      = RIMACS::getConnectedConfig(RIMACS::Config::ResultType::AllMaximum, maxNofEquivalentResults);
  auto res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, conf, false,
                                             &queryEquivalenceClasses, &targetEquivalenceClasses);

  EXPECT_EQ(res.size(), nofEquivalentResults);
  for(const RIMACS::MCSResult& result : res) {
    EXPECT_EQ(result.size(), size);
    verify_result(queryAdj, targetAdj, nodeComp, edgeComp, result);
  }
}

namespace {

using EQclass = std::vector<RIMACS::MappingIndex>;

INSTANTIATE_TEST_SUITE_P(Graphs, EquivalenceClassTester,
    testing::Values(std::make_tuple(Data::SimpleThreeMemberChain, Data::SimpleFourMemberChain,
                                    EQclass({0, 1, 0}), EQclass({0, 1, 1, 0}), 1, 1, 3),
                    std::make_tuple(Data::get_graph("triangle"), Data::get_graph("triod"),
                                    EQclass({0, 0, 0}), EQclass({0, 1, 0, 0}), 1, 1, 2),
                    std::make_tuple(Data::get_graph("triangle"), Data::get_graph("linked_triangle"),
                                    EQclass({0, 0, 0}), EQclass({0, 1, 2, 2}), 1, 1, 3),
                    std::make_tuple(Data::get_graph("triangle"), Data::get_graph("linked_triangle"),
                                    EQclass({0, 0, 0}), EQclass({0, 1, 2, 2}), 2, 2, 3),
                    // add an artificial marker to the equivalence classes of the triangle
                    std::make_tuple(Data::get_graph("triangle"), Data::get_graph("linked_triangle"),
                                    EQclass({1, 0, 0}), EQclass({0, 1, 2, 2}), 2, 4, 3),

                    // Three different results that are not equivalent
                    std::make_tuple(Data::get_graph("linked_triod"), Data::get_graph("diamond"),
                                    EQclass({0, 1, 2, 3, 3}), EQclass({0, 1, 0, 1}), 1, 3, 3),

                    // Two results per equivalence class on the results
                    std::make_tuple(Data::get_graph("linked_triod"), Data::get_graph("diamond"),
                                    EQclass({0, 1, 2, 3, 3}), EQclass({0, 1, 0, 1}), 2, 6, 3)));

} // namespace


class MappingMCSTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph,
                                    std::vector<RIMACS::MappingPair>,
                                    size_t, size_t>> {};

TEST_P(MappingMCSTester, testInitialMapping)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;
  std::vector<RIMACS::MappingPair> initialMapping;
  size_t resultSize, singleMappingSize;

  std::tie(query, target, initialMapping, resultSize, singleMappingSize) = GetParam();

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  RIMACS::MCSResults res
      = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, initialMapping);

  ASSERT_EQ(res.size(), 1);
  EXPECT_EQ(res.front().size(), resultSize);
  EXPECT_EQ(res.front().weight(), static_cast<double>(resultSize));

  verify_result(queryAdj, targetAdj, nodeComp, edgeComp, res.front().getMappings());

  RIMACS::MCSResults singleMappingResult
      = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, initialMapping.front());

  ASSERT_EQ(singleMappingResult.size(), 1);
  EXPECT_EQ(singleMappingResult.front().size(), singleMappingSize);
  EXPECT_EQ(singleMappingResult.front().weight(), static_cast<double>(singleMappingSize));

  verify_result(queryAdj, targetAdj, nodeComp, edgeComp, singleMappingResult.front().getMappings());
}

TEST_P(MappingMCSTester, testWeightedInitialMapping)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;
  std::vector<RIMACS::MappingPair> initialMapping;
  size_t resultSize, singleMappingSize;

  std::tie(query, target, initialMapping, resultSize, singleMappingSize) = GetParam();

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  const RIMACS::Config& conf = RIMACS::getDefaultConnectedConfig();

  RIMACS::MCSResults res
      = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, initialMapping, conf, true);

  ASSERT_EQ(res.size(), 1);
  EXPECT_EQ(res.front().size(), resultSize);
  EXPECT_EQ(res.front().weight(), static_cast<double>(resultSize));

  verify_result(queryAdj, targetAdj, nodeComp, edgeComp, res.front().getMappings());

  RIMACS::MCSResults singleMappingResult
      = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, initialMapping.front(), conf, true);

  ASSERT_EQ(singleMappingResult.size(), 1);
  EXPECT_EQ(singleMappingResult.front().size(), singleMappingSize);
  EXPECT_EQ(singleMappingResult.front().weight(), static_cast<double>(singleMappingSize));

  verify_result(queryAdj, targetAdj, nodeComp, edgeComp, singleMappingResult.front().getMappings());
}

namespace {

using IniMapping = std::vector<RIMACS::MappingPair>;

INSTANTIATE_TEST_SUITE_P(Graphs, MappingMCSTester,
    testing::Values(std::make_tuple(Data::get_simple_chain(6, 8, 7, 16, 6, 8, 7),
                                    Data::get_simple_chain(6, 8, 7, 15, 6, 8, 7),
                                    IniMapping({{0,0}, {4, 4}}), size_t{6}, size_t{3}),

                    std::make_tuple(Data::get_graph("triangle"), Data::get_graph("triangle"),
                                    IniMapping({{1,2}}), size_t{3}, size_t{3}),

                    std::make_tuple(Data::get_simple_chain(6, 8, 6, 6, 6),
                                    Data::get_simple_chain(6, 6, 6, 8, 6),
                                    IniMapping({{1,3}, {0, 2}}), size_t{3}, size_t{5})
                                    ));

} // namespace


class HintedMappingTester
    : public testing::TestWithParam<std::tuple<std::string, std::string,
                                               std::vector<std::vector<RIMACS::MappingIndex>>,
                                               size_t, size_t>> {};

TEST_P(HintedMappingTester, testEffectOfHint)
{
  std::string queryFile, targetFile;
  std::vector<std::vector<RIMACS::MappingIndex>> hint;
  size_t normalSize, hintedSize;

  std::tie(queryFile, targetFile, hint, normalSize, hintedSize) = GetParam();

  Dimacs::DimacsGraph query = Data::get_graph(queryFile);
  Dimacs::DimacsGraph target = Data::get_graph(targetFile);

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target)

  auto result = RIMACS::MCSRunner::runGenericMCSIncludingCount(queryAdj, targetAdj, nodeComp, edgeComp);

  ASSERT_FALSE(result.first.empty());
  EXPECT_EQ(result.first.front().size(), normalSize);

  result = RIMACS::MCSRunner::runHintedMCS(queryAdj, targetAdj, nodeComp, edgeComp,
                                          hint, RIMACS::getDefaultConnectedConfig());

  ASSERT_FALSE(result.first.empty());
  EXPECT_EQ(result.first.front().size(), hintedSize);
}

namespace {

  using Hint = std::vector<std::vector<RIMACS::MappingIndex>>;


INSTANTIATE_TEST_SUITE_P(Graphs, HintedMappingTester,
    testing::Values(std::make_tuple("chain_and_ring", "chain_and_ring",
                                    Hint({{9},{8},{7},{6},{5},{4},{3},{2},{1},{0}}), 10, 2),
                    std::make_tuple("chains_and_ring_3", "chains_and_ring_1",
                                    Hint({{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13}}), 6, 4),
                    std::make_tuple("chains_and_ring_4", "chains_and_ring_2",
                                    Hint({{10},{9},{8},{7},{6},{5},{4},{3},{2},{1},{0}}), 6, 6)));

} // namespace

class HintedCorrectnessTester
    : public testing::TestWithParam<std::tuple<std::string, std::string, size_t>> {};

TEST_P(HintedCorrectnessTester, testHintDoesNotGenerateCompatibility)
{
  std::string queryFile, targetFile;
  size_t size;

  std::tie(queryFile, targetFile, size) = GetParam();

  Dimacs::DimacsGraph query = Data::get_graph(queryFile);
  Dimacs::DimacsGraph target = Data::get_graph(targetFile);

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target)

  auto result = RIMACS::MCSRunner::runGenericMCSIncludingCount(queryAdj, targetAdj, nodeComp, edgeComp);

  ASSERT_FALSE(result.first.empty());
  EXPECT_EQ(result.first.front().size(), size);

  std::vector<RIMACS::MappingIndex> singleHint(targetAdj.getNofNodes());
  std::iota(singleHint.begin(), singleHint.end(), 0);
  std::vector<std::vector<RIMACS::MappingIndex>> hint(queryAdj.getNofNodes(), singleHint);

  result = RIMACS::MCSRunner::runHintedMCS(queryAdj, targetAdj, nodeComp, edgeComp,
                                          hint, RIMACS::getDefaultConnectedConfig());

  ASSERT_FALSE(result.first.empty());
  EXPECT_EQ(result.first.front().size(), size);
}

INSTANTIATE_TEST_SUITE_P(Graphs, HintedCorrectnessTester,
    testing::Values(std::make_tuple("graphs_1", "graphs_2", 8),
                    std::make_tuple("graphs_3", "graphs_4", 3),
                    std::make_tuple("graphs_5", "graphs_6", 7),
                    std::make_tuple("graphs_7", "graphs_8", 16),
                    std::make_tuple("graphs_9", "graphs_10", 3)));
