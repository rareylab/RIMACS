
#include <gtest/gtest.h>

#include <RIMACS/MCSRunner.hpp>

#include "InitializationUtils.hpp"

#include "ResultVerification.hpp"

using EqClass = std::vector<RIMACS::MappingIndex>;

class EquivalenceClassTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph,
                                               EqClass, EqClass, unsigned, unsigned>> {};

TEST_P(EquivalenceClassTester, equivalenceClassMatch)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;
  EqClass queryClasses;
  EqClass targetClasses;
  unsigned nofResults;
  unsigned resultSize;

  std::tie(query, target, queryClasses, targetClasses, nofResults, resultSize) = GetParam();

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  const EqClass* queryClassPtr = queryClasses.empty() ? nullptr : &queryClasses;
  const EqClass* targetClassPtr = targetClasses.empty() ? nullptr : &targetClasses;

  auto config = RIMACS::getConnectedConfig(RIMACS::Config::ResultType::AllMaximum, 1);
  RIMACS::MCSResults res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj, nodeComp, edgeComp,
                                                            queryClassPtr, targetClassPtr, config);

  EXPECT_EQ(res.size(), nofResults);
  for(const RIMACS::MCSResult& singleRes : res) {
    EXPECT_EQ(singleRes.size(), resultSize);
    EXPECT_EQ(singleRes.getComponentSizes().size(), 1);
    verify_result(queryAdj, targetAdj, nodeComp, edgeComp, singleRes);
  }
}

INSTANTIATE_TEST_SUITE_P(SimplePairs, EquivalenceClassTester,
    testing::Values(std::make_tuple(Data::get_simple_chain(1, 1), Data::get_simple_chain(1, 1),
                                    EqClass(), EqClass(), 2, 2),
                    std::make_tuple(Data::get_simple_chain(1, 1), Data::get_simple_chain(1, 1),
                                    EqClass({0, 0}), EqClass(), 1, 2),
                    std::make_tuple(Data::get_simple_chain(1, 1), Data::get_simple_chain(1, 1),
                                    EqClass(), EqClass({0, 0}), 1, 2),
                    std::make_tuple(Data::get_simple_chain(1, 1), Data::get_simple_chain(1, 1),
                                    EqClass({0, 0}), EqClass({0, 0}), 1, 2)));


INSTANTIATE_TEST_SUITE_P(SixMemberRing, EquivalenceClassTester,
    testing::Values(std::make_tuple(Data::get_simple_cycle(1, 1, 1, 1, 1, 1),
                                    Data::get_simple_cycle(1, 1, 1, 1, 1, 2),
                                    EqClass(), EqClass(), 12, 5),
                    std::make_tuple(Data::get_simple_cycle(1, 1, 1, 1, 1, 1),
                                    Data::get_simple_cycle(1, 1, 1, 1, 1, 2),
                                    EqClass({0, 0, 0, 0, 0, 0}), EqClass(), 1, 5),
                    std::make_tuple(Data::get_simple_cycle(1, 1, 1, 1, 1, 1),
                                    Data::get_simple_cycle(1, 1, 1, 1, 1, 2),
                                    EqClass(), EqClass({0, 1, 2, 1, 0, 3}), 6, 5),
                    std::make_tuple(Data::get_simple_cycle(1, 1, 1, 1, 1, 1),
                                    Data::get_simple_cycle(1, 1, 1, 1, 1, 2),
                                    EqClass({0, 0, 0, 0, 0, 0}), EqClass({0, 1, 2, 1, 0, 3}), 1, 5)));


INSTANTIATE_TEST_SUITE_P(SimpleTree, EquivalenceClassTester,
    testing::Values(std::make_tuple(Data::get_graph("tree_1"),
                                    Data::get_graph("tree_2"),
                                    EqClass(), EqClass(), 72, 9),
                    std::make_tuple(Data::get_graph("tree_1"),
                                    Data::get_graph("tree_2"),
                                    EqClass({0, 1, 1, 1, 2, 3, 4, 5, 6, 7, 4, 5}), EqClass(), 6, 9),
                    std::make_tuple(Data::get_graph("tree_1"),
                                    Data::get_graph("tree_2"),
                                    EqClass(), EqClass({0, 1, 1, 1, 2, 0, 1, 1, 1}), 1, 9),
                    std::make_tuple(Data::get_graph("tree_1"),
                                    Data::get_graph("tree_2"),
                                    EqClass({0, 1, 1, 1, 2, 3, 4, 5, 6, 7, 4, 5}),
                                    EqClass({0, 1, 1, 1, 2, 0, 1, 1, 1}), 1, 9)));
