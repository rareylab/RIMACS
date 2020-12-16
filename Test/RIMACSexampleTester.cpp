
#include <gtest/gtest.h>

#include <RIMACS/MCSRunner.hpp>

#include "Dimacs/DimacsGraph.hpp"
#include "Dimacs/DimacsGraphFunctor.hpp"
#include "InitializationUtils.hpp"

TEST(MCSexample, testSimpleMCS)
{
  const Dimacs::DimacsAdjacencyFunctor query(Data::SimpleCarboxylicAcidLike);
  const Dimacs::DimacsAdjacencyFunctor target(Data::SimpleAmideBondLike);

  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(query, target);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(query, target);

  RIMACS::MCSResults res = RIMACS::MCSRunner::runGenericMCS(query, target,
                                                          nodeComp, edgeComp,
                                                          RIMACS::getDefaultConnectedConfig());

  // Maximum MCS (default one) has at most one result
  EXPECT_EQ(res.size(), size_t{1});
  ASSERT_GE(res.size(), size_t{1});

  // Retrieve the actual atom mapping
  const std::vector<RIMACS::MappingPair>& nodeMapping = res.front().getMappings();
  std::vector<RIMACS::MappingPair> refNodeMapping({{0, 0},
                                                {1, 1},
                                                {2, 2}});
  EXPECT_EQ(nodeMapping, refNodeMapping);

  // Test the component sizes of the mapping
  EXPECT_EQ(res.front().getComponentSizes(), std::vector<RIMACS::MappingIndex>(1, 3));
}


TEST(MCSexample, testDisconnectedMCS)
{
  const Dimacs::DimacsAdjacencyFunctor query(Data::TwiceCarboxylicAcidLike);
  const Dimacs::DimacsAdjacencyFunctor target(Data::TwiceAmideBondLike);

  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(query, target);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(query, target);

  RIMACS::MCSResults res = RIMACS::MCSRunner::runGenericMCS(query, target, nodeComp, edgeComp,
                                                      RIMACS::getDefaultDisconnectedConfig());

  // Maximum MCS (default one) has at most one result
  EXPECT_EQ(res.size(), size_t{1});

  // Retrieve the actual atom mapping
  const std::vector<RIMACS::MappingPair>& nodeMapping = res.front().getMappings();
  std::vector<RIMACS::MappingPair> refNodeMapping({{0, 0},
                                                  {1, 1},
                                                  {2, 2},
                                                  {4, 4},
                                                  {5, 5},
                                                  {6, 6}});
  EXPECT_EQ(nodeMapping, refNodeMapping);

  // Test the component sizes of the mapping
  EXPECT_EQ(res.front().getComponentSizes(), std::vector<RIMACS::MappingIndex>({3, 6}));

  // Set a different number of allowed component (2) and a different minimum size for additional
  // connected components (5)
  res = RIMACS::MCSRunner::runGenericMCS(query, target, nodeComp, edgeComp,
                                        RIMACS::getDisconnectedConfig(RIMACS::Config::ResultType::Maximum, 2, 5));

  // Maximum MCS (default one) has at most one result
  EXPECT_EQ(res.size(), size_t{1});
  ASSERT_GE(res.size(), size_t{1});

  // Now there is only the connected mapping, since the disconnected extension does not have 5 atoms
  EXPECT_EQ(res.front().getComponentSizes(), std::vector<RIMACS::MappingIndex>(1, 3));
}


TEST(MCSexample, testAllMaximumMCS)
{
  const Dimacs::DimacsAdjacencyFunctor query(Data::LocalizedBenzene);
  const Dimacs::DimacsAdjacencyFunctor target(Data::LocalizedPhenol);

  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(query, target);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(query, target);

  RIMACS::MCSResults res
      = RIMACS::MCSRunner::runGenericMCS(query, target, nodeComp, edgeComp,
                                        RIMACS::getConnectedConfig(RIMACS::Config::ResultType::AllMaximum));

  // All maximum MCS may have several results
  EXPECT_EQ(res.size(), size_t{6});
  ASSERT_GE(res.size(), size_t{1});

  // iterate over all base MCS results
  for(auto it = res.begin(), last = res.end(); it != last; ++it) {
    // Test for the size
    EXPECT_EQ(it->size(), size_t{6});
    // Test for the component sizes of the mapping
    EXPECT_EQ(it->getComponentSizes(), std::vector<RIMACS::MappingIndex>(1, 6));
  }

  std::vector<RIMACS::MappingIndex> queryEquivalenceClasses({0, 0, 0, 0, 0, 0});
  std::vector<RIMACS::MappingIndex> targetEquivalenceClasses({0, 0, 0, 0, 0, 0, 1});

  // Allow some equivalent results (3)
  RIMACS::Config config = RIMACS::getConnectedConfig(RIMACS::Config::ResultType::AllMaximum, 3);

  // If we provide the config, then we can only choose between atom MCS (induced) and edge MCS (noninduced)
  res = RIMACS::MCSRunner::runGenericMCS(query, target, nodeComp, edgeComp, config,
                                      false, &queryEquivalenceClasses, &targetEquivalenceClasses);

  // All maximum MCS may have several results
  EXPECT_EQ(res.size(), size_t{3});
  ASSERT_GE(res.size(), size_t{1});

  // iterate over all base MCS results
  for(auto it = res.begin(), last = res.end(); it != last; ++it) {
    // Test for the size
    EXPECT_EQ(it->size(), size_t{6});
    // Test for the component sizes of the mapping
    EXPECT_EQ(it->getComponentSizes(), std::vector<RIMACS::MappingIndex>(1, 6));
  }
}


TEST(MCSexample, testMCES)
{
  const Dimacs::DimacsAdjacencyFunctor query(Data::FourMemberExampleCycle);
  const Dimacs::DimacsAdjacencyFunctor target(Data::FourMemberExampleChain);

  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(query, target);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(query, target);

  RIMACS::MCSResults res = RIMACS::MCSRunner::runGenericMCES(query, target,
                                                           nodeComp, edgeComp,
                                                           RIMACS::getDefaultConfig());

  // Maximum MCS (default one) has at most one result
  EXPECT_EQ(res.size(), size_t{1});
  ASSERT_GE(res.size(), size_t{1});

  // Retrieve the actual bond mapping
  const std::vector<RIMACS::MappingPair>& edgeMapping = res.front().getMappings();
  std::vector<RIMACS::MappingPair> refEdgeMapping({{0, 0},
                                                  {1, 1},
                                                  {2, 2},
                                                  {3, 3}});
  EXPECT_EQ(edgeMapping, refEdgeMapping);

  // Test the component sizes of the mapping
  EXPECT_EQ(res.front().getComponentSizes(),
            std::vector<RIMACS::MappingIndex>(size_t{1}, RIMACS::MappingIndex{4}));
}


TEST(MCSexampleTester, testNodeMappingMCES)
{
  const Dimacs::DimacsAdjacencyFunctor query(Data::get_simple_cycle(8, 6, 6, 7));
  const Dimacs::DimacsAdjacencyFunctor target(Data::get_simple_chain(8, 6, 6, 7));

  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(query, target);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(query, target);
  // Run a noninduced MCS and convert the result to an atom mapping
  RIMACS::MCSResults res;
  double edgeBasedSimilarity;
  std::tie(res, std::ignore, edgeBasedSimilarity)
      = RIMACS::MCSRunner::runGenericNodeMappingMCES(query, target,
                                                    nodeComp, edgeComp,
                                                    RIMACS::getDefaultConfig());

  // The underlying edge mapping does not cover all edges, but the node mapping does
  EXPECT_LT(edgeBasedSimilarity, 1);
  // Maximum MCS (default one) has at most one result
  EXPECT_EQ(res.size(), size_t{1});
  ASSERT_GE(res.size(), 1);

  RIMACS::MCSResult& result = res.front();

  // Retrieve the actual atom mapping
  const std::vector<RIMACS::MappingPair>& nodeMapping = result.getMappings();
  std::vector<RIMACS::MappingPair> refNodeMapping({{0, 0}, {1, 1}, {2, 2}, {3, 3}});
  EXPECT_EQ(nodeMapping, refNodeMapping);

  // Test the component sizes of the mapping
  EXPECT_EQ(result.getComponentSizes(),
            std::vector<RIMACS::MappingIndex>(size_t{1}, RIMACS::MappingIndex{4}));
}

TEST(MCSexampleTester, testAllMaximumNodeMappingMCES)
{
  const Dimacs::DimacsAdjacencyFunctor query(Data::get_graph("MCES_size_demo_query"));
  const Dimacs::DimacsAdjacencyFunctor target(Data::get_graph("MCES_size_demo_target"));

  Dimacs::DimacsNodeCompatibilityFunctor nodeComp(query, target);
  Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(query, target);
  // Run a noninduced MCS and convert the result to an atom mapping
  RIMACS::MCSResults res;
  double edgeBasedSimilarity;
  // Run a disconnected MCS and convert the result to an atom mapping
  std::tie(res, std::ignore, edgeBasedSimilarity)
      = RIMACS::MCSRunner::runGenericNodeMappingMCES(query, target,
                                                    nodeComp, edgeComp,
                                                    RIMACS::getDisconnectedConfig(RIMACS::Config::ResultType::AllMaximum, 10, 2));

  // The graphs allow for one connected result and 16 permutations of a result with 4 connected components
  ASSERT_GE(res.size(), size_t{1});
  EXPECT_EQ(res.size(), size_t{17});

  // The first result has only one component
  EXPECT_EQ(res.front().getComponentSizes(), std::vector<RIMACS::MappingIndex>(1, 12));
  for(int i = 1 ; i < res.size() ; ++i) {
    // The remaining result consist of 4 components
    EXPECT_EQ(res.at(i).getComponentSizes(), std::vector<RIMACS::MappingIndex>({3, 6, 9, 12}));
  }

  // Search for the maximum result for the same pair of graphs
  std::tie(res, std::ignore, edgeBasedSimilarity)
      = RIMACS::MCSRunner::runGenericNodeMappingMCES(query, target,
                                                    nodeComp, edgeComp,
                                                    RIMACS::getDisconnectedConfig(RIMACS::Config::ResultType::Maximum, 10, 2));
  ASSERT_EQ(res.size(), size_t{1});
  // The result will always be the connected result
  EXPECT_EQ(res.front().getComponentSizes(), std::vector<RIMACS::MappingIndex>(1, 12));
}

