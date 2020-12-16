
#include <gtest/gtest.h>

#include <RIMACS/Functors/GraphCachingFunctor.hpp>
#include <RIMACS/Utils/SortedVector.hpp>

#include "InitializationUtils.hpp"
#include "ResultVerification.hpp"



class FunctorTester : public testing::TestWithParam<Dimacs::DimacsGraph> {};

TEST_P(FunctorTester, basicAdjacency)
{
  Dimacs::DimacsGraph graph = GetParam();

  Dimacs::BasicDimacsAdjacencyFunctor baseFunctor(graph);

  EXPECT_EQ(baseFunctor.getNofNodes(), graph.m_nodeLabels.size());
  EXPECT_EQ(baseFunctor.getNofEdges(), graph.m_edges.size());
  for(size_t i = 0, last = graph.m_edges.size(); i < last; ++i) {
    const Dimacs::DimacsEdge& e = graph.m_edges[i];
    EXPECT_EQ(baseFunctor.edgeGetFromId(i), e.m_from);
    EXPECT_EQ(baseFunctor.edgeGetToId(i), e.m_to);
  }
}

TEST_P(FunctorTester, basicAdjacencyForCachingFunctor)
{
  Dimacs::DimacsGraph graph = GetParam();

  Dimacs::DimacsAdjacencyFunctor functor(graph);
  const RIMACS::CachingRequiredMinimalAdjacencyFunctor& baseFunctor = functor;

  EXPECT_EQ(baseFunctor.getNofNodes(), graph.m_nodeLabels.size());
  EXPECT_EQ(baseFunctor.getNofEdges(), graph.m_edges.size());
  for(size_t i = 0, last = graph.m_edges.size(); i < last; ++i) {
    const Dimacs::DimacsEdge& e = graph.m_edges[i];
    EXPECT_EQ(baseFunctor.edgeGetFromId(i), e.m_from);
    EXPECT_EQ(baseFunctor.edgeGetToId(i), e.m_to);
  }
}

TEST_P(FunctorTester, cachingAdjacency)
{
  Dimacs::DimacsGraph graph = GetParam();

  Dimacs::DimacsAdjacencyFunctor functor(graph);

  EXPECT_EQ(functor.getNofNodes(), graph.m_nodeLabels.size());
  EXPECT_EQ(functor.getNofEdges(), graph.m_edges.size());

  std::vector<RIMACS::SortedVector<RIMACS::MappingIndex>> neighbours(graph.m_nodeLabels.size());
  for(size_t from = 0, last = graph.m_nodeLabels.size(); from < last; ++from) {
    for(size_t to = 0; to < last; ++to) {
      EXPECT_EQ(functor.adjacent(from, to), connected(graph, from, to));
      if(connected(graph, from, to)) {
        neighbours.at(from).insert_unique(to);
        neighbours.at(to).insert_unique(from);
        EXPECT_NE(functor.getEdgeId(from, to), RIMACS::AdjacencyFunctor::NO_EDGE);
        EXPECT_NE(functor.getEdgeId(to, from), RIMACS::AdjacencyFunctor::NO_EDGE);
      } else {
        EXPECT_EQ(functor.getEdgeId(from, to), RIMACS::AdjacencyFunctor::NO_EDGE);
        EXPECT_EQ(functor.getEdgeId(to, from), RIMACS::AdjacencyFunctor::NO_EDGE);
      }
    }
  }
  std::vector<RIMACS::SortedVector<RIMACS::MappingIndex>> incidentEdges(graph.m_nodeLabels.size());
  for(size_t i = 0, last = graph.m_edges.size(); i < last; ++i) {
    const Dimacs::DimacsEdge& e = graph.m_edges[i];
    incidentEdges.at(e.m_from).insert_unique(i);
    incidentEdges.at(e.m_to).insert_unique(i);
    EXPECT_EQ(functor.getEdgeId(e.m_from, e.m_to), i);
    EXPECT_EQ(functor.getEdgeId(e.m_to, e.m_from), i);
  }

  for(size_t i = 0, last = graph.m_nodeLabels.size(); i < last; ++i) {
    const std::vector<RIMACS::MappingIndex>& expectedNeighbours = neighbours[i];
    std::vector<RIMACS::MappingIndex> neighbours = functor.nodeGetNeighbours(i);
    std::sort(neighbours.begin(), neighbours.end());
    EXPECT_EQ(neighbours, expectedNeighbours);

    const std::vector<RIMACS::MappingIndex>& expectedEdgeIds = incidentEdges[i];
    std::vector<RIMACS::MappingIndex> incidentEdges = functor.nodeGetEdges(i);
    std::sort(incidentEdges.begin(), incidentEdges.end());
    EXPECT_EQ(incidentEdges, expectedEdgeIds);
    EXPECT_EQ(functor.nodeGetNofEdges(i), expectedEdgeIds.size());
  }
}

INSTANTIATE_TEST_SUITE_P(Graphs, FunctorTester,
    testing::Values(Data::get_graph("caffeine"),
                    Data::get_graph("chain_and_ring"),
                    Data::get_graph("chains_and_ring_1"),
                    Data::get_graph("chains_and_ring_2"),
                    Data::get_graph("chains_and_ring_3"),
                    Data::get_graph("chains_and_ring_4"),
                    Data::get_graph("diamond"),
                    Data::get_graph("dMCSdemo_1"),
                    Data::get_graph("dMCSdemo_2"),
                    Data::get_graph("dMCSdemo_3"),
                    Data::get_graph("graphs_10"),
                    Data::get_graph("graphs_11"),
                    Data::get_graph("graphs_1"),
                    Data::get_graph("graphs_2"),
                    Data::get_graph("graphs_3"),
                    Data::get_graph("graphs_4"),
                    Data::get_graph("graphs_5"),
                    Data::get_graph("graphs_6"),
                    Data::get_graph("graphs_7"),
                    Data::get_graph("graphs_8"),
                    Data::get_graph("graphs_9"),
                    Data::get_graph("TestCases_10_f"),
                    Data::get_graph("TestCases_10_t"),
                    Data::get_graph("TestCases_11_f"),
                    Data::get_graph("TestCases_11_t"),
                    Data::get_graph("TestCases_12_f"),
                    Data::get_graph("TestCases_12_t"),
                    Data::get_graph("TestCases_1_f"),
                    Data::get_graph("TestCases_1_t"),
                    Data::get_graph("TestCases_2_f"),
                    Data::get_graph("TestCases_2_t"),
                    Data::get_graph("TestCases_3_f"),
                    Data::get_graph("TestCases_3_t"),
                    Data::get_graph("TestCases_4_f"),
                    Data::get_graph("TestCases_4_t"),
                    Data::get_graph("TestCases_5_f"),
                    Data::get_graph("TestCases_5_t"),
                    Data::get_graph("TestCases_6_f"),
                    Data::get_graph("TestCases_6_t"),
                    Data::get_graph("TestCases_7_f"),
                    Data::get_graph("TestCases_7_t"),
                    Data::get_graph("TestCases_8_f"),
                    Data::get_graph("TestCases_8_t"),
                    Data::get_graph("TestCases_9_f"),
                    Data::get_graph("TestCases_9_t"),
                    Data::get_graph("linked_diamond"),
                    Data::get_graph("linked_triangle"),
                    Data::get_graph("linked_triod"),
                    Data::get_graph("marked_tetrahedron"),
                    Data::get_graph("marked_triangle"),
                    Data::get_graph("marked_triod"),
                    Data::get_graph("ring_and_chain"),
                    Data::get_graph("ring_and_chains_1"),
                    Data::get_graph("ring_and_chains_2"),
                    Data::get_graph("ring_and_chains_3"),
                    Data::get_graph("ring_and_chains_4"),
                    Data::get_graph("tetrahedron"),
                    Data::get_graph("theobromine"),
                    Data::get_graph("triangle"),
                    Data::get_graph("triangle_triod"),
                    Data::get_graph("triod")));



