
#include <gtest/gtest.h>


#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Functors/LineGraph.hpp>
#include <RIMACS/MCSRunner.hpp>

#include "Dimacs/DimacsGraph.hpp"
#include "Dimacs/DimacsGraphFunctor.hpp"

#include "InitializationUtils.hpp"
#include "ResultVerification.hpp"


class LineGraphTester : public testing::TestWithParam<Dimacs::DimacsGraph> {};

TEST_P(LineGraphTester, testAdjacency)
{
  Dimacs::DimacsGraph graph = GetParam();

  RIMACS::LineGraph line{Dimacs::DimacsAdjacencyFunctor(graph)};

  EXPECT_EQ(line.getNofNodes(), graph.m_edges.size());

  for(size_t i = 0; i < line.getNofNodes(); ++i) {
    for(size_t j : line.nodeGetNeighbours(i)) {
      EXPECT_NE(get_common_node_idx(graph, i, j), -1);
    }
  }
  for(size_t i = 0, last = graph.m_edges.size(), lastI = last - 1; i < lastI ; ++i) {
    EXPECT_EQ(line.getEdgeId(i, i), std::numeric_limits<RIMACS::MappingIndex>::max());
    for(size_t j = i + 1; j < last; ++j) {
      bool connectedInLineGraph = line.getEdgeId(i, j) < line.getNofEdges();
      bool connectedInMolecule
          = get_common_node_idx(graph, i, j) != -1;
      EXPECT_EQ(connectedInLineGraph, connectedInMolecule);
      connectedInLineGraph = line.getEdgeId(j, i) < line.getNofEdges();
      EXPECT_EQ(connectedInLineGraph, connectedInMolecule);
    }
  }
}


INSTANTIATE_TEST_SUITE_P(graphs, LineGraphTester,
    testing::Values(Data::get_graph("caffeine"),
                    Data::get_graph("triangle"),
                    Data::get_graph("triod"),
                    Data::SimpleFourMemberChain));

namespace {

class LoggingNodeCompatibilityFunctor : public RIMACS::NodeCompatibilityFunctor {
public:

  explicit LoggingNodeCompatibilityFunctor(
      const RIMACS::NodeCompatibilityFunctor& base)
    : m_baseFunctor(base)
  {}

  std::vector<RIMACS::MappingPair>&& moveLog()
  {
    return std::move(m_log);
  }

private:
  bool nodesAreCompatible(RIMACS::MappingIndex from, RIMACS::MappingIndex to) const override
  {
    m_log.emplace_back(from, to);
    return m_baseFunctor.nodesAreCompatible(from, to);
  }

  const RIMACS::NodeCompatibilityFunctor& m_baseFunctor;
  mutable std::vector<RIMACS::MappingPair> m_log;
};

class LoggingEdgeCompatibilityFunctor : public RIMACS::EdgeCompatibilityFunctor {
public:

  explicit LoggingEdgeCompatibilityFunctor(
      const RIMACS::EdgeCompatibilityFunctor& base)
    : m_baseFunctor(base)
  {}

  std::vector<RIMACS::MappingPair>&& moveLog()
  {
    return std::move(m_log);
  }

private:
  bool edgesAreCompatible(RIMACS::MappingIndex from, RIMACS::MappingIndex to) const override
  {
    m_log.emplace_back(from, to);
    return m_baseFunctor.edgesAreCompatible(from, to);
  }

  const RIMACS::EdgeCompatibilityFunctor& m_baseFunctor;
  mutable std::vector<RIMACS::MappingPair> m_log;
};


} // namespace


class LineGraphCompatibilityTester
    : public testing::TestWithParam<std::pair<Dimacs::DimacsGraph, Dimacs::DimacsGraph>> {};

TEST_P(LineGraphCompatibilityTester, testLineGraphAtomCompatibility)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor queryAdj(query);
  Dimacs::DimacsAdjacencyFunctor targetAdj(target);
  RIMACS::LineGraph queryLine{queryAdj};
  RIMACS::LineGraph targetLine{targetAdj};

  Dimacs::DimacsNodeCompatibilityFunctor atomComp(queryAdj, targetAdj);
  Dimacs::DimacsEdgeCompatibilityFunctor bondComp(queryAdj, targetAdj);
  LoggingNodeCompatibilityFunctor nodeComp(atomComp);
  LoggingEdgeCompatibilityFunctor edgeComp(bondComp);

  RIMACS::LineGraphNodeCompatibilityFunctor lineNodeComp{queryLine, targetLine, nodeComp, edgeComp};

  ASSERT_EQ(query.m_edges.size(), queryLine.getNofNodes());
  ASSERT_EQ(target.m_edges.size(), targetLine.getNofNodes());

  for(RIMACS::MappingIndex i = 0, lastI = queryLine.getNofNodes(); i < lastI; ++i) {
    const Dimacs::DimacsEdge& queryEdge = query.m_edges.at(i);
    ASSERT_LT(queryEdge.m_from, queryEdge.m_to);
    for(RIMACS::MappingIndex j = 0, lastJ = targetLine.getNofNodes(); j < lastJ; ++j) {
      const Dimacs::DimacsEdge& targetEdge = target.m_edges.at(j);
      ASSERT_LT(targetEdge.m_from, targetEdge.m_to);
      bool compatible = lineNodeComp(i, j);
      std::vector<RIMACS::MappingPair> nodeComparisons = std::move(nodeComp.moveLog());
      std::vector<RIMACS::MappingPair> edgeComparisons = std::move(edgeComp.moveLog());
      bool bondsCompatible = bondComp(i, j);
      bool atomsCompatibleInOrder
          = (atomComp(queryEdge.m_from, targetEdge.m_from)
             && atomComp(queryEdge.m_to, targetEdge.m_to));
      bool atomsCompatibleInReverseOrder
          = (atomComp(queryEdge.m_to, targetEdge.m_from)
             && atomComp(queryEdge.m_from, targetEdge.m_to));
      bool atomsCompatible
          = atomsCompatibleInOrder || atomsCompatibleInReverseOrder;
      EXPECT_EQ(compatible, bondsCompatible && atomsCompatible);
      if(compatible) {
        EXPECT_EQ(edgeComparisons.size(), size_t{1});
        EXPECT_EQ(edgeComparisons.front(), RIMACS::MappingPair(i, j));
        if(atomsCompatibleInOrder) {
          EXPECT_NE(std::find(nodeComparisons.begin(), nodeComparisons.end(),
                            RIMACS::MappingPair(queryEdge.m_from,
                                               targetEdge.m_from)),
                    nodeComparisons.end());
          EXPECT_NE(std::find(nodeComparisons.begin(), nodeComparisons.end(),
                            RIMACS::MappingPair(queryEdge.m_to,
                                               targetEdge.m_to)),
                    nodeComparisons.end());
        } else if(atomsCompatibleInReverseOrder) {
          EXPECT_NE(std::find(nodeComparisons.begin(), nodeComparisons.end(),
                            RIMACS::MappingPair(queryEdge.m_from,
                                               targetEdge.m_to)),
                   nodeComparisons.end());
          EXPECT_NE(std::find(nodeComparisons.begin(), nodeComparisons.end(),
                            RIMACS::MappingPair(queryEdge.m_to,
                                               targetEdge.m_from)),
                    nodeComparisons.end());
        }
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Graphs, LineGraphCompatibilityTester,
    testing::Values(std::make_pair(Data::SimpleFourMemberChain, Data::SimpleFourMemberChain),
                    std::make_pair(Data::get_graph("triangle"), Data::get_graph("triangle")),
                    std::make_pair(Data::get_graph("triod"), Data::get_graph("triod")),
                    std::make_pair(Data::get_graph("triangle"), Data::get_graph("triod")),
                    std::make_pair(Data::get_graph("triod"), Data::get_graph("triangle")),
                    std::make_pair(Data::get_graph("caffeine"), Data::get_graph("theobromine"))));


TEST_P(LineGraphCompatibilityTester, testLineGraphBondCompatibility)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;

  std::tie(query, target) = GetParam();

  Dimacs::DimacsAdjacencyFunctor queryAdj(query);
  Dimacs::DimacsAdjacencyFunctor targetAdj(target);
  RIMACS::LineGraph queryLine{queryAdj};
  RIMACS::LineGraph targetLine{targetAdj};

  Dimacs::DimacsNodeCompatibilityFunctor atomComp(queryAdj, targetAdj);
  Dimacs::DimacsEdgeCompatibilityFunctor bondComp(queryAdj, targetAdj);
  LoggingNodeCompatibilityFunctor nodeComp(atomComp);
  LoggingEdgeCompatibilityFunctor edgeComp(bondComp);

  RIMACS::LineGraphEdgeCompatibilityFunctor lineEdgeComp{queryLine, targetLine, nodeComp};

  for(RIMACS::MappingIndex i = 0, lastI = queryLine.getNofEdges(); i < lastI; ++i) {
    const Dimacs::DimacsEdge& qEdgeFrom = query.m_edges.at(queryLine.edgeGetFromId(i));
    const Dimacs::DimacsEdge& qEdgeTo   = query.m_edges.at(queryLine.edgeGetToId(i));
    int commonQueryNode = get_common_node_idx(qEdgeFrom, qEdgeTo);
    for(RIMACS::MappingIndex j = 0, lastJ = targetLine.getNofEdges(); j < lastJ; ++j) {
      const Dimacs::DimacsEdge& tEdgeFrom = target.m_edges.at(targetLine.edgeGetFromId(j));
      const Dimacs::DimacsEdge& tEdgeTo   = target.m_edges.at(targetLine.edgeGetToId(j));
      int commonTargetNode = get_common_node_idx(tEdgeFrom, tEdgeTo);
      bool compatible = lineEdgeComp(i, j);
      std::vector<RIMACS::MappingPair> nodeComparisons = std::move(nodeComp.moveLog());
      bool bondsAdjacent = commonQueryNode != -1 && commonTargetNode != -1;
      bool atomsCompatible = bondsAdjacent
                             && atomComp(commonQueryNode, commonTargetNode);
      EXPECT_EQ(compatible, atomsCompatible);
      EXPECT_EQ(nodeComparisons.size(), size_t{1});
      EXPECT_EQ(nodeComparisons.front(), (RIMACS::MappingPair(commonQueryNode, commonTargetNode)));
    }
  }
}


class BasicMCESTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph, size_t>> {};

TEST_P(BasicMCESTester, testMoleculeMCES)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;
  size_t size;

  std::tie(query, target, size) = GetParam();

  Dimacs::DimacsAdjacencyFunctor queryAdj(query);
  Dimacs::DimacsAdjacencyFunctor targetAdj(target);
  RIMACS::LineGraph queryLine{queryAdj};
  RIMACS::LineGraph targetLine{targetAdj};

  Dimacs::DimacsNodeCompatibilityFunctor atomComp(queryAdj, targetAdj);
  Dimacs::DimacsEdgeCompatibilityFunctor bondComp(queryAdj, targetAdj);
  RIMACS::LineGraphNodeCompatibilityFunctor nodeComp(queryLine, targetLine, atomComp, bondComp);
  RIMACS::LineGraphEdgeCompatibilityFunctor edgeComp(queryLine, targetLine, atomComp);

  Dimacs::DimacsNodeCompatibilityFunctor reverseAtomComp(targetAdj, queryAdj);
  Dimacs::DimacsEdgeCompatibilityFunctor reverseBondComp(targetAdj, queryAdj);
  RIMACS::LineGraphNodeCompatibilityFunctor reverseNodeComp(targetLine, queryLine, reverseAtomComp, reverseBondComp);
  RIMACS::LineGraphEdgeCompatibilityFunctor reverseEdgeComp(targetLine, queryLine, reverseAtomComp);

  auto result = RIMACS::MCSRunner::runGenericMCSIncludingCount(queryLine, targetLine, nodeComp, edgeComp);
  EXPECT_EQ(result.first.size(), size_t{1});
  EXPECT_EQ(result.first.front().size(), size);
  verify_mces_result(queryAdj, targetAdj, atomComp, bondComp,result.first.front());

  // redo the test with different molecule order
  result = RIMACS::MCSRunner::runGenericMCSIncludingCount(targetLine, queryLine, reverseNodeComp, reverseEdgeComp);
  EXPECT_EQ(result.first.size(), size_t{1});
  EXPECT_EQ(result.first.front().size(), size);
  verify_mces_result(targetAdj, queryAdj, reverseAtomComp, reverseBondComp, result.first.front());
}



INSTANTIATE_TEST_SUITE_P(Graphs, BasicMCESTester,
    testing::Values(std::make_tuple(Data::SimpleFourMemberChain, Data::SimpleFourMemberChain, 3),
                    std::make_tuple(Data::get_graph("triod"), Data::get_graph("triangle"), 2),
                    std::make_tuple(Data::get_graph("triangle"), Data::get_graph("triod"), 2),
                    std::make_tuple(Data::get_graph("linked_triangle"), Data::get_graph("triod"), 3),
                    std::make_tuple(Data::get_graph("linked_triod"), Data::get_graph("triangle"), 2),
                    std::make_tuple(Data::get_graph("triangle"), Data::get_graph("linked_triod"), 2),
                    std::make_tuple(Data::get_graph("triod"), Data::get_graph("linked_triangle"), 3),
                    std::make_tuple(Data::get_graph("linked_triangle"), Data::get_graph("linked_triod"), 3),
                    std::make_tuple(Data::get_graph("linked_triod"), Data::get_graph("linked_triangle"), 3),
                    std::make_tuple(Data::get_graph("caffeine"), Data::get_graph("theobromine"), 14)));


class EdgeMappingConversionTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph,
                                    std::vector<RIMACS::MappingPair>, std::vector<RIMACS::MappingIndex>,
                                    std::vector<std::vector<RIMACS::MappingPair>>,
                                    std::vector<std::vector<RIMACS::MappingIndex>> > > {};

TEST_P(EdgeMappingConversionTester, testMappingConversion)
{
  using MappingPair = RIMACS::MappingPair;
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;
  std::vector<MappingPair> mapping;
  std::vector<RIMACS::MappingIndex> components;
  std::vector<std::vector<MappingPair>> resultMappings;
  std::vector<std::vector<RIMACS::MappingIndex>> resultComponents;

  std::tie(query, target, mapping, components, resultMappings, resultComponents) = GetParam();

  ASSERT_EQ(resultMappings.size(), resultComponents.size());


  Dimacs::DimacsAdjacencyFunctor queryAdj(query);
  Dimacs::DimacsAdjacencyFunctor targetAdj(target);
  Dimacs::DimacsNodeCompatibilityFunctor atomComp(queryAdj, targetAdj);
  RIMACS::LineGraph queryLine(queryAdj);
  RIMACS::LineGraph targetLine(targetAdj);

  RIMACS::EdgeMappingConverter conv(queryLine, targetLine, atomComp);

  auto convertedResults = conv(mapping, components, true);
  ASSERT_EQ(convertedResults.size(), resultMappings.size());
  for(size_t i = 0; i < convertedResults.size(); ++i) {
    std::vector<MappingPair> actualMapping;
    std::vector<RIMACS::MappingIndex> actualComponents;
    std::tie(actualMapping, actualComponents) = convertedResults.at(i);
    EXPECT_EQ(actualComponents, resultComponents.at(i));
    EXPECT_EQ(actualMapping, resultMappings.at(i));
  }
}

namespace {
  using MappingPair = RIMACS::MappingPair;
  using Mappings = std::vector<MappingPair>;
  using MappingsVector = std::vector<Mappings>;
  using Component = std::vector<RIMACS::MappingIndex>;
  using Components = std::vector<Component>;

INSTANTIATE_TEST_SUITE_P(Graphs, EdgeMappingConversionTester,
    testing::Values(std::make_tuple(Data::get_graph("triangle"), Data::get_graph("triangle"),
                                    Mappings({{0, 1}, {1, 0}, {2, 2}}), Component({3}),
                                    MappingsVector({{{0, 0}, {1, 2}, {2, 1}}}), Components({{3}})),
                    std::make_tuple(Data::get_graph("triod"), Data::get_graph("triod"),
                                    Mappings({{0, 0}, {1, 1}, {2, 2}}), Component({3}),
                                    MappingsVector({{{0, 0}, {1, 1}, {2, 2}, {3, 3}}}), Components({{4}})),
                    std::make_tuple(Data::get_graph("triangle"), Data::SimpleThreeMemberChain,
                                    Mappings({{0, 0}, {1, 1}}), Component({2}),
                                    MappingsVector({{{0, 1}, {1, 0}, {2, 2}}}), Components({{3}})),
                    std::make_tuple(Data::get_graph("triod"), Data::get_graph("triod"),
                                    Mappings({{0, 0}, {1, 1}, {2, 2}}), Component({3}),
                                    MappingsVector({{{0, 0}, {1, 1}, {2, 2}, {3, 3}}}), Components({{4}}))));




INSTANTIATE_TEST_SUITE_P(Components, EdgeMappingConversionTester,
    testing::Values(std::make_tuple(Data::get_graph("triangle_triod"), Data::get_graph("triangle_triod"),
                                    Mappings({{0,0}, {1, 1}, {2, 2}, {5, 5}, {6, 6}, {7, 7}}), Component({3, 6}),
                                    MappingsVector({{{0,0}, {1, 1}, {2, 2}, {4, 4}, {5, 5}, {6, 6}, {7, 7}}}),
                                    Components({{3, 7}}))));

INSTANTIATE_TEST_SUITE_P(LoneEdges, EdgeMappingConversionTester,
    testing::Values(std::make_tuple(Data::get_simple_chain(7, 6, 6),
                                    Data::get_simple_chain(8, 6, 6),
                                    Mappings({{1, 1}}), Component({1}),
                                    MappingsVector({{{1, 2}, {2, 1}}, {{1, 1}, {2, 2}}}),
                                    Components({{2}, {2}})),
                    std::make_tuple(Data::get_simple_chain(16, 7, 6, 6),
                                    Data::get_simple_chain(7, 16, 8, 6, 6),
                                    Mappings({{0, 0}, {2, 3}}), Component({1, 2}),
                                    MappingsVector({{{0, 1}, {1, 0}, {2, 4}, {3, 3}},
                                                    {{0, 1}, {1, 0}, {2, 3}, {3, 4}}}),
                                    Components({{2, 4}, {2, 4}}))));
} // namespace


class MCESNodeResultTester
    : public testing::TestWithParam<std::tuple<Dimacs::DimacsGraph, Dimacs::DimacsGraph,
                                               RIMACS::Config::ResultType, size_t, size_t,
                                               std::vector<RIMACS::MappingIndex> > > {};

TEST_P(MCESNodeResultTester, testMoleculeMCESNodeResult)
{
  Dimacs::DimacsGraph query;
  Dimacs::DimacsGraph target;
  RIMACS::Config::ResultType res;
  size_t nComp;
  size_t compS;
  std::vector<RIMACS::MappingIndex> sizes;

  std::tie(query, target, res, nComp, compS, sizes) = GetParam();

  Dimacs::DimacsAdjacencyFunctor queryAdj(query);
  Dimacs::DimacsAdjacencyFunctor targetAdj(target);

  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(queryAdj, targetAdj);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(queryAdj, targetAdj);

  RIMACS::MCSResults nodeResult;
  RIMACS::MCSResults edgeResult;

  RIMACS::Config conf = RIMACS::getConfig(res, nComp, compS);
  std::tie(nodeResult, std::ignore, std::ignore)
      = RIMACS::MCSRunner::runGenericNodeMappingMCES(queryAdj, targetAdj, nodeComp, edgeComp, conf);
  edgeResult = RIMACS::MCSRunner::runGenericMCES(queryAdj, targetAdj, nodeComp, edgeComp, conf);

  EXPECT_GE(nodeResult.size(), edgeResult.size());

  std::vector<RIMACS::MappingIndex> actualResultSizes;
  std::transform(nodeResult.begin(), nodeResult.end(), std::back_inserter(actualResultSizes),
                 [](const RIMACS::MCSResult& r){
    return r.size();
  });
  for(size_t i = 0, last = std::max(nodeResult.size(), edgeResult.size()); i < last; ++i) {
    const RIMACS::MCSResult& atomRes = nodeResult[std::min(i, nodeResult.size() - 1)];
    const RIMACS::MCSResult& bondRes = edgeResult[std::min(i, edgeResult.size() - 1)];
    verify_mces_result(queryAdj, targetAdj, nodeComp, edgeComp, bondRes);
    verify_result<false>(queryAdj, targetAdj, nodeComp, edgeComp, atomRes);
  }
  EXPECT_EQ(actualResultSizes, sizes);
}

namespace {
using MappingSizes = std::vector<RIMACS::MappingIndex>;

INSTANTIATE_TEST_SUITE_P(Graphs, MCESNodeResultTester,
    testing::Values(std::make_tuple(Data::get_graph("linked_triangle"), Data::get_graph("linked_triangle"),
                                    RIMACS::Config::ResultType::AllMaximum,
                                    size_t{10}, size_t{10},
                                    MappingSizes({4, 4})),

                    std::make_tuple(Data::get_graph("triangle"), Data::get_graph("triangle"),
                                    RIMACS::Config::ResultType::Maximum,
                                    size_t{3}, size_t{1},
                                    MappingSizes({3})),

                    std::make_tuple(Data::get_graph("marked_triangle"), Data::get_graph("marked_triod"),
                                    RIMACS::Config::ResultType::AllMaximum,
                                    size_t{3}, size_t{1},
                                    MappingSizes({3, 3, 3, 3, 3, 3}))
                                         ));

INSTANTIATE_TEST_SUITE_P(LoneBonds, MCESNodeResultTester,
    testing::Values(std::make_tuple(Data::get_simple_chain(6, 6, 7, 6, 6),
                                    Data::get_simple_chain(6, 6, 8, 6, 6),
                                    RIMACS::Config::ResultType::AllMaximum,
                                    size_t{3}, size_t{1},
                                    MappingSizes({4, 4, 4, 4, 4, 4, 4, 4}))));


INSTANTIATE_TEST_SUITE_P(EdgeSymmetry, MCESNodeResultTester,
    testing::Values(std::make_tuple(Data::get_graph("linked_triangle"), Data::get_graph("linked_triangle"),
                                    RIMACS::Config::ResultType::AllMaximum,
                                    size_t{10}, size_t{10},
                                    MappingSizes({4, 4})),

                    std::make_tuple(Data::get_graph("diamond"), Data::get_graph("diamond"),
                                    RIMACS::Config::ResultType::AllMaximum,
                                    size_t{10}, size_t{10},
                                    MappingSizes({4, 4, 4, 4})),

                    std::make_tuple(Data::get_graph("marked_tetrahedron"), Data::get_graph("marked_tetrahedron"),
                                    RIMACS::Config::ResultType::AllMaximum,
                                    size_t{10}, size_t{10},
                                    MappingSizes({4, 4, 4, 4, 4, 4})),

                    std::make_tuple(Data::get_graph("tetrahedron"), Data::get_graph("tetrahedron"),
                                    RIMACS::Config::ResultType::AllMaximum,
                                    size_t{10}, size_t{10},
                                    MappingSizes({4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4})),

                    std::make_tuple(Data::get_graph("linked_diamond"), Data::get_graph("diamond"),
                                    RIMACS::Config::ResultType::AllMaximum,
                                    size_t{10}, size_t{10},
                                    MappingSizes({4, 4, 4, 4})),

                    std::make_tuple(Data::get_graph("diamond"), Data::get_graph("linked_diamond"),
                                    RIMACS::Config::ResultType::AllMaximum,
                                    size_t{10}, size_t{10},
                                    MappingSizes({4, 4, 4, 4}))
                                    ));
} // namespace

namespace {

std::vector<RIMACS::MappingIndex> get_intra_edges(
    const Dimacs::DimacsGraph& graph,
    const std::vector<RIMACS::MappingIndex>& nodes)
{
  std::vector<RIMACS::MappingIndex> result;
  std::vector<bool> mappedNodes(graph.m_nodeLabels.size(), false);
  for(const auto& node : nodes) {
    mappedNodes.at(node) = true;
  }
  for(size_t i = 0, last = graph.m_edges.size(); i < last; ++i) {
    if(mappedNodes.at(graph.m_edges[i].m_from)
       && mappedNodes.at(graph.m_edges[i].m_to)) {
      result.push_back(i);
    }
  }
  return result;
}

std::vector<RIMACS::MappingIndex> to_node_vector(
    const Dimacs::DimacsGraph& graph,
    const std::vector<RIMACS::MappingIndex>& edges)
{
  std::vector<RIMACS::MappingIndex> result;
  std::vector<bool> nodesIncidentToMappedEdges(graph.m_nodeLabels.size(), false);
  for(const auto& edge : edges) {
    nodesIncidentToMappedEdges.at(graph.m_edges.at(edge).m_from) = true;
    nodesIncidentToMappedEdges.at(graph.m_edges.at(edge).m_to) = true;
  }
  for(size_t i = 0, last = graph.m_nodeLabels.size(); i < last; ++i) {
    if(nodesIncidentToMappedEdges.at(i)) {
      result.push_back(i);
    }
  }
  return result;
}

} // namespace

class MaximumNumberOfMappedNodesMCESTester
    : public testing::TestWithParam<std::tuple<std::string, std::string, int, size_t, size_t>> {};

TEST_P(MaximumNumberOfMappedNodesMCESTester, testMaximizeNumberOfMappedNodes)
{
  std::string queryGraphFile;
  std::string targetGraphFile;
  int nofComponents;
  size_t nodeSize;
  size_t edgeSize;

  std::tie(queryGraphFile, targetGraphFile, nofComponents, nodeSize, edgeSize) = GetParam();

  Dimacs::DimacsGraph query = Data::get_graph(queryGraphFile);
  Dimacs::DimacsGraph target = Data::get_graph(targetGraphFile);

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(query, target);

  const RIMACS::Config conf = RIMACS::getConfig(RIMACS::Config::ResultType::Maximum, nofComponents);

  auto nodeResult = RIMACS::MCSRunner::runGenericNodeMappingMCES(queryAdj, targetAdj,
                                                                nodeComp, edgeComp, conf);

  auto edgeResult = RIMACS::MCSRunner::runGenericMCESIncludingCount(queryAdj, targetAdj,
                                                                   nodeComp, edgeComp, conf);

  ASSERT_EQ(std::get<0>(nodeResult).size(), size_t{1});
  ASSERT_EQ(edgeResult.first.size(), size_t{1});

  const RIMACS::MCSResult nodeMappingResult = std::get<0>(nodeResult).front();
  const RIMACS::MCSResult edgeMappingResult = edgeResult.first.front();
  EXPECT_EQ(nodeMappingResult.size(), nodeSize);
  EXPECT_EQ(edgeMappingResult.size(), edgeSize);

  size_t nofMappedBondsGraph1 = 0;
  size_t nofMappedBondsGraph2 = 0;
  for(size_t compIdx = 0, currentIdx = 0, last = nodeMappingResult.getComponentSizes().size()
      ; compIdx < last
      ; ++compIdx) {
    // Extract the mapped atoms per component in order to guarantee that no unmapped bond is assumed as mapped
    std::vector<RIMACS::MappingIndex> mappedQueryNodes;
    std::vector<RIMACS::MappingIndex> mappedTargetNodes;
    for(size_t endIdx = nodeMappingResult.getComponentSizes()[compIdx]; currentIdx < endIdx; ++currentIdx) {
      mappedQueryNodes.push_back(nodeMappingResult.getMappings().at(currentIdx).m_from);
      mappedTargetNodes.push_back(nodeMappingResult.getMappings().at(currentIdx).m_to);
    }
    nofMappedBondsGraph1 += get_intra_edges(query, mappedQueryNodes).size();
    nofMappedBondsGraph2 += get_intra_edges(target, mappedTargetNodes).size();
  }
  EXPECT_EQ(nofMappedBondsGraph1, edgeSize);
  EXPECT_EQ(nofMappedBondsGraph2, edgeSize);

  std::vector<RIMACS::MappingIndex> mappedQueryEdges;
  std::vector<RIMACS::MappingIndex> mappedTargetEdges;
  for(const RIMACS::MappingPair& edgeMapping: edgeMappingResult.getMappings()) {
    mappedQueryEdges.push_back(edgeMapping.m_from);
    mappedTargetEdges.push_back(edgeMapping.m_to);
  }
  EXPECT_LT(to_node_vector(query, mappedQueryEdges).size(), nodeSize
               /*&& "The unguided noninduced MCS accidentally found the optimum node mapping result. "
                  "Please fix the example of the test"*/);
  EXPECT_LT(to_node_vector(target, mappedTargetEdges).size(), nodeSize);
  EXPECT_LT(edgeResult.second, std::get<1>(nodeResult)
              /*&& "Searching for the maximum number of atoms while keeping the bond mapping "
                 "maximal should be more expensive than just looking for a bond maximum"*/);
}

INSTANTIATE_TEST_SUITE_P(Graphs, MaximumNumberOfMappedNodesMCESTester,
    testing::Values(std::make_tuple("chains_and_ring_5", "chains_and_ring_6", 1, size_t{4}, size_t{3}),
                    std::make_tuple("chains_and_ring_7", "chains_and_ring_8", 1, size_t{4}, size_t{3}),
                    std::make_tuple("chains_and_ring_9", "chains_and_ring_10", 1, size_t{4}, size_t{3}),
                    std::make_tuple("chains_and_ring_12", "chains_and_ring_11", 1, size_t{4}, size_t{3}),
                    std::make_tuple("chains_and_ring_12", "chains_and_ring_11", 2, size_t{8}, size_t{6}),
                    std::make_tuple("chains_and_ring_12", "chains_and_ring_11", 3, size_t{11}, size_t{9})));


class DeltaYTester
    : public testing::TestWithParam<std::tuple<std::string, unsigned, unsigned, unsigned>> {};

TEST_P(DeltaYTester, testAllVariants)
{
  std::string graphFile;
  unsigned nofResults;
  unsigned nodeResultSize;
  unsigned edgeResultSize;

  std::tie(graphFile, nofResults, nodeResultSize, edgeResultSize) = GetParam();
  Dimacs::DimacsGraph graph = Data::get_graph(graphFile);

  INITIALIZE_DEFAULT_NAMED_FUNCTORS(graph, graph);

  RIMACS::Config conf = RIMACS::getConnectedConfig(RIMACS::Config::ResultType::AllMaximum);
  RIMACS::MCSResults res = RIMACS::MCSRunner::runGenericMCS(queryAdj, targetAdj,
                                                            nodeComp, edgeComp, conf);

  EXPECT_EQ(res.size(), nofResults);
  for(const RIMACS::MCSResult& singleRes : res) {
    EXPECT_EQ(singleRes.size(), nodeResultSize);
    verify_result(queryAdj, targetAdj, nodeComp, edgeComp, singleRes);
  }
  RIMACS::MCSResults edgeRes = RIMACS::MCSRunner::runGenericMCES(queryAdj, targetAdj,
                                                                 nodeComp, edgeComp, conf);

  EXPECT_EQ(edgeRes.size(), nofResults);
  for(const RIMACS::MCSResult& singleRes : edgeRes) {
    EXPECT_EQ(singleRes.size(), edgeResultSize);
    verify_mces_result(queryAdj, targetAdj, nodeComp, edgeComp, singleRes);
  }

  RIMACS::LineGraph queryLine(queryAdj);
  RIMACS::LineGraph targetLine(targetAdj);
  RIMACS::LineGraphNodeCompatibilityFunctor lineNodeComp(queryLine, targetLine, nodeComp, edgeComp);
  RIMACS::LineGraphEdgeCompatibilityFunctor lineEdgeComp(queryLine, targetLine, nodeComp);

  RIMACS::MCSResults realEdgeRes = RIMACS::MCSRunner::runGenericMCS(queryLine, targetLine,
                                                                    lineNodeComp, lineEdgeComp, conf);

  EXPECT_EQ(edgeRes.size(), nofResults);
  for(const RIMACS::MCSResult& singleRes : realEdgeRes) {
    EXPECT_EQ(singleRes.size(), edgeResultSize);
    verify_mces_result(queryAdj, targetAdj, nodeComp, edgeComp, singleRes);
  }
}

INSTANTIATE_TEST_SUITE_P(Examples, DeltaYTester,
    testing::Values(std::make_tuple("linked_triangle", 2, 4, 4),
                    std::make_tuple("diamond", 4, 4, 5),
                    std::make_tuple("tetrahedron", 24, 4, 6)));
