
#include <gtest/gtest.h>

#include "RIMACS/MCSRunner.hpp"

#include "Dimacs/DimacsGraphFunctor.hpp"
#include "InitializationUtils.hpp"

class MinimumResultSizeTester
    : public testing::TestWithParam<std::tuple<std::string, std::string, int, size_t, size_t>> {};

TEST_P(MinimumResultSizeTester, testMinSize)
{
  std::string queryGraphFile;
  std::string targetGraphFile;
  int nofComponents;
  size_t minSize;
  size_t resultSize;

  std::tie(queryGraphFile, targetGraphFile, nofComponents, minSize, resultSize) = GetParam();

  Dimacs::DimacsGraph query = Data::get_graph(queryGraphFile);
  Dimacs::DimacsGraph target = Data::get_graph(targetGraphFile);

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  RIMACS::Config conf = RIMACS::getConfig(RIMACS::Config::ResultType::Maximum, nofComponents);
  conf.setMinimumSize(minSize);
  auto result = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp, conf);

  if(minSize > resultSize) {
    EXPECT_EQ(result.size(), size_t{0});
  } else {
    ASSERT_EQ(result.size(), size_t{1});
    EXPECT_EQ(result.front().size(), resultSize);
  }
}


INSTANTIATE_TEST_SUITE_P(MinSizeGraphs, MinimumResultSizeTester,
     testing::Values(std::make_tuple("triangle", "triangle", 1, size_t{3}, size_t{3}),
                     std::make_tuple("triangle", "triangle", 1, size_t{4}, size_t{3}),
                     std::make_tuple("chains_and_ring_13", "chains_and_ring_14", 1, size_t{9}, size_t{9}),
                     std::make_tuple("chains_and_ring_13", "chains_and_ring_14", 1, size_t{10}, size_t{9}),
                     std::make_tuple("chains_and_ring_13", "chains_and_ring_14", 2, size_t{10}, size_t{10})));

namespace {

class minimum_node_position_weighter : public Dimacs::DimacsNodeCompatibilityFunctor {
public:
  minimum_node_position_weighter(
      const Dimacs::DimacsAdjacencyFunctor& query,
      const Dimacs::DimacsAdjacencyFunctor& target)
    : DimacsNodeCompatibilityFunctor(query, target)
  {}

  double getWeight(
      RIMACS::MappingIndex queryGraphNode,
      RIMACS::MappingIndex targetGraphNode) const override
  {
    return std::min(queryGraphNode, targetGraphNode) + 0.125;
  }
};

} // namespace

class MinimumResultWeightTester
    : public testing::TestWithParam<std::tuple<std::string, std::string, int, double, double>> {};

TEST_P(MinimumResultWeightTester, testMinWeight)
{
  std::string queryGraphFile;
  std::string targetGraphFile;
  int nofComponents;
  double minWeight;
  double resultWeight;

  std::tie(queryGraphFile, targetGraphFile, nofComponents, minWeight, resultWeight) = GetParam();

  Dimacs::DimacsGraph query = Data::get_graph(queryGraphFile);
  Dimacs::DimacsGraph target = Data::get_graph(targetGraphFile);

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);
  minimum_node_position_weighter weightedNodeComp(queryAdj, targetAdj);

  RIMACS::Config conf = RIMACS::getConfig(RIMACS::Config::ResultType::Maximum, nofComponents);
  conf.setMinimumWeight(minWeight);
  auto result = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, weightedNodeComp, edgeComp, conf, true);

  if(minWeight > resultWeight) {
    EXPECT_EQ(result.size(), size_t{0});
  } else {
    ASSERT_EQ(result.size(), size_t{1});
    EXPECT_EQ(result.front().weight(), resultWeight);
  }
}


INSTANTIATE_TEST_SUITE_P(MinWeightGraphs, MinimumResultWeightTester,
       testing::Values(std::make_tuple("triangle", "triangle", 1, 3.375, 3.375),
                       std::make_tuple("triangle", "triangle", 1, 3.376, 3.375),
                       std::make_tuple("chains_and_ring_13", "chains_and_ring_14", 1, 35.0, 35.0),
                       std::make_tuple("chains_and_ring_13", "chains_and_ring_14", 1, 35.125, 35.0),
                       std::make_tuple("chains_and_ring_13", "chains_and_ring_14", 2, 51.5, 51.25),
                       std::make_tuple("chains_and_ring_13", "chains_and_ring_14", 2, 51.25, 51.25)));

//void MinimumResultConstraintTester::testMinimumSizeMCS_data()
//{
//  QTest::newRow("Small components") << "C" << "C" << 1 << size_t{2} << size_t{1};
//  QTest::newRow("Different Substitution c f") << "NCCc1ccc(cc1)CCO" << "NCCc1cccc(c1)CCO" << 1 << size_t{10} << size_t{9};
//  QTest::newRow("Different Substitution c s") << "NCCc1ccc(cc1)CCO" << "NCCc1cccc(c1)CCO" << 1 << size_t{9} << size_t{9};
//  QTest::newRow("Different Substitution d2 s") << "NCCc1ccc(cc1)CCO" << "NCCc1cccc(c1)CCO" << 2 << size_t{10} << size_t{10};
//  QTest::newRow("Different Substitution d2 f") << "NCCc1ccc(cc1)CCO" << "NCCc1cccc(c1)CCO" << 2 << size_t{11} << size_t{10};
//}

//void MinimumResultConstraintTester::testMinimumWeightMCS_data()
//{
//  QTest::newRow("Small components, no result") << "C" << "C" << 1 << 0.25 << 0.125;
//  QTest::newRow("Small components, success") << "C" << "C" << 1 << 0.0 << 0.125;
//  QTest::newRow("Different Substitution c s") << "c1(ccc(cc1)CCO)CCN" << "c1(cccc(c1)CCN)CCO" << 1 << 35.0 << 35.0;
//  QTest::newRow("Different Substitution c f") << "c1(ccc(cc1)CCO)CCN" << "c1(cccc(c1)CCN)CCO" << 1 << 35.125 << 35.0;
//  QTest::newRow("Different Substitution d2 f") << "c1(ccc(cc1)CCO)CCN" << "c1(cccc(c1)CCN)CCO" << 2 << 51.5 << 51.25;
//  QTest::newRow("Different Substitution s2 s") << "c1(ccc(cc1)CCO)CCN" << "c1(cccc(c1)CCN)CCO" << 2 << 51.0 << 51.25;
//}

