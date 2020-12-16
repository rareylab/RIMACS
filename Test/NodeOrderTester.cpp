
#include <gtest/gtest.h>

#include "InitializationUtils.hpp"

#include <RIMACS/Functors/LineGraph.hpp>

class NodeOrderTester
    : public testing::TestWithParam<std::pair<Dimacs::DimacsGraph,
                                              std::vector<RIMACS::MappingIndex>>> {};

TEST_P(NodeOrderTester, testSuggestion)
{
  Dimacs::DimacsGraph graph;
  std::vector<RIMACS::MappingIndex> nodeOrderSuggestion;
  std::tie(graph, nodeOrderSuggestion) = GetParam();

  Dimacs::DimacsAdjacencyFunctor adj(graph);
  std::vector<RIMACS::MappingIndex> actualSuggestion;
  actualSuggestion.reserve(adj.getNofNodes());
  for(size_t i = 0, last = adj.getNofNodes(); i < last; ++i) {
    actualSuggestion.push_back(adj.getNodeOrderSuggestion(i));
  }
  EXPECT_EQ(actualSuggestion, nodeOrderSuggestion);
}

namespace {
using Order = std::vector<RIMACS::MappingIndex>;
INSTANTIATE_TEST_SUITE_P(Graphs, NodeOrderTester,
   testing::Values(std::make_pair(Data::get_simple_chain(1, 1),
                                  Order({0, 0})),
                    std::make_pair(Data::get_simple_chain(1, 1, 1, 1, 1, 1),
                                  Order({0, 1, 2, 2, 1, 0})),
                    std::make_pair(Data::get_simple_chain(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                                  Order({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0})),
                   std::make_pair(Data::get_graph("tree"),
                                  Order({4, 3, 2, 1, 0, 1, 0, 2, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 2, 1, 0, 1, 0})),
                   std::make_pair(Data::get_graph("large_tree"),
                                  Order({4, 3, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0})),
                   std::make_pair(Data::get_graph("deep_tree"),
                                  Order({7, 6, 5, 4, 3, 2, 1, 0, 2, 1, 0, 4, 3, 2, 1, 0, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0, 2, 1, 0, 4, 3, 2, 1, 0, 2, 1, 0})),
                   std::make_pair(Data::get_graph("caffeine"),
                                  Order({0, 1, 1, 1, 2, 2, 2, 1, 2, 1, 0, 2, 1, 1})),
                   std::make_pair(Data::get_graph("remdesivir"),
                                  Order({0, 1, 2, 1, 0, 3, 4, 5, 4, 6, 5, 7, 8, 7, 9, 8, 7, 6, 5, 5, 6, 4, 3, 4, 3, 2, 2, 3, 2, 1, 0, 1, 0, 4, 5, 7, 6, 5, 4, 3, 4, 5}))));

} // namespace



class EdgeOrderTester : public testing::TestWithParam<std::pair<Dimacs::DimacsGraph, std::vector<RIMACS::MappingIndex>>> {};

TEST_P(EdgeOrderTester, testSuggestion)
{
  Dimacs::DimacsGraph graph;
  std::vector<RIMACS::MappingIndex> edgeOrderSuggestion;
  std::tie(graph, edgeOrderSuggestion) = GetParam();

  Dimacs::DimacsAdjacencyFunctor adj(graph);
  RIMACS::LineGraph line(adj);
  std::vector<RIMACS::MappingIndex> actualSuggestion;
  actualSuggestion.reserve(line.getNofNodes());
  for(size_t i = 0, last = line.getNofNodes(); i < last; ++i) {
    actualSuggestion.push_back(line.getNodeOrderSuggestion(i));
  }
  EXPECT_EQ(actualSuggestion, edgeOrderSuggestion);
}

namespace {
using Order = std::vector<RIMACS::MappingIndex>;

INSTANTIATE_TEST_SUITE_P(Graphs, EdgeOrderTester,
    testing::Values(std::make_pair(Data::get_simple_chain(1, 1),
                                   Order({0})),
                    std::make_pair(Data::get_simple_chain(1, 1, 1, 1, 1, 1),
                                   Order({1, 0, 1, 0, 1})),
                    std::make_pair(Data::get_simple_chain(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                                   Order({3, 2, 1, 0, 1, 2, 3, 4, 5, 6, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3})),
                    std::make_pair(Data::get_graph("tree"),
                                   Order({2, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 2, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1})),
                    std::make_pair(Data::get_graph("large_tree"),
                                   Order({2, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 2, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 2, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1})),
                    std::make_pair(Data::get_graph("deep_tree"),
                                   Order({4, 3, 2, 1, 0, 1, 2, 0, 1, 2, 2, 1, 0, 1, 2, 0, 1, 2, 4, 3, 2, 1, 0, 1, 2, 0, 1, 2, 2, 1, 0, 1, 2, 0, 1, 2})),
                    std::make_pair(Data::get_graph("caffeine"),
                                   Order({1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1})),
                    std::make_pair(Data::get_graph("remdesivir"),
                                   Order({2, 1, 1, 2, 0, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6, 5, 4, 3, 2, 3, 4, 2, 1, 2, 1, 0, 1, 0, 1, 0, 1, 2, 2, 1, 2, 2, 3, 5, 4, 3, 2, 1, 1, 2, 3}))));
} // namespace
