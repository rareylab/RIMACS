
#pragma once

#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Config.hpp>
#include <RIMACS/Algorithm/WeightTypeTraits.hpp>

#include "Dimacs/DimacsGraphFunctor.hpp"

extern std::string path;
extern char delimiter;

#define INITIALIZE_DEFAULT_NAMED_FUNCTORS(GRAPH1, GRAPH2) \
    Dimacs::DimacsAdjacencyFunctor queryAdj(GRAPH1); \
    Dimacs::DimacsAdjacencyFunctor targetAdj(GRAPH2); \
    Dimacs::DimacsNodeCompatibilityFunctor nodeComp(queryAdj, targetAdj); \
    Dimacs::DimacsEdgeCompatibilityFunctor edgeComp(queryAdj, targetAdj);

namespace Data {
using Nodes = std::vector<int>;
using Edges = std::vector<Dimacs::DimacsEdge>;
using Graph = Dimacs::DimacsGraph;


inline std::string get_path(
    const std::string& name)
{
  std::stringstream str;
  str << path << "Examples" << delimiter << name << ".dimacs";
  return str.str();
}

inline Dimacs::DimacsGraph get_graph(
    const std::string& name)
{
  Dimacs::DimacsGraph ret(get_path(name));
  if(ret.name.empty()) {
    ret.name = name;
  }
  return ret;
}

template<typename... IntTypes>
inline Graph get_simple_chain(
    int fstNodeLabel,
    IntTypes... remainingNodeLabels)
{
  Nodes labels({fstNodeLabel, remainingNodeLabels...});
  std::vector<std::tuple<int, int, int>> edges;
  edges.reserve(labels.size() - 1);
  for(size_t i = 1; i < labels.size(); ++i) {
    edges.emplace_back(static_cast<int>(i - 1), static_cast<int>(i), 1);
  }
  return Graph(labels, edges);
}

template<typename... IntTypes>
inline Graph get_simple_cycle(
    int fstNodeLabel,
    IntTypes... remainingNodeLabels)
{
  Nodes labels({fstNodeLabel, remainingNodeLabels...});
  std::vector<std::tuple<int, int, int>> edges;
  edges.reserve(labels.size() - 1);
  for(size_t i = 1; i < labels.size(); ++i) {
    edges.emplace_back(static_cast<int>(i - 1), static_cast<int>(i), 1);
  }
  edges.insert(edges.begin() + 1, {0, static_cast<int>(labels.size() - 1), 1});
  return Graph(labels, edges);
}

static const Graph IsoThreeMemberChain = Graph{Nodes{6, 6, 6}, Edges{Dimacs::DimacsEdge(0, 1), Dimacs::DimacsEdge(0, 2)}};
static const Graph SimpleThreeMemberChain = get_simple_chain(6, 6, 6);
static const Graph SimpleFourMemberChain = get_simple_chain(6, 6, 6, 6);

static const Graph MixedFourMemberChain = get_simple_chain(6, 8, 6, 6);
static const Graph MixedSixMemberChain = get_simple_chain(6, 6, 8, 6, 6, 6);
static const Graph MixedThreeMemberChain = get_simple_chain(6, 8, 6);

static const Graph SimpleCarboxylicAcidLike = Graph{Nodes{6, 6, 8, 8}, Edges{{0, 1, 1}, {1, 2, 2}, {1, 3, 1}}};
static const Graph SimpleAmideBondLike = Graph{Nodes{6, 6, 8, 7}, Edges{{0, 1, 1}, {1, 2, 2}, {1, 3, 1}}};

static const Graph TwiceCarboxylicAcidLike = Graph{Nodes{6, 6, 8, 8, 6, 6, 8, 8}, Edges{{0, 1, 1}, {1, 2, 2}, {1, 3, 1}, {3, 4, 1}, {4, 5, 1}, {5, 6, 2}, {5, 7, 1}}};
static const Graph TwiceAmideBondLike = Graph{Nodes{6, 6, 8, 7, 6, 6, 8, 7}, Edges{{0, 1, 1}, {1, 2, 2}, {1, 3, 1}, {3, 4, 1}, {4, 5, 1}, {5, 6, 2}, {5, 7, 1}}};

static const Graph LocalizedBenzene = Graph{Nodes{6, 6, 6, 6, 6, 6}, Edges{{0, 1, 2}, {1, 2, 1}, {2, 3, 2}, {3, 4, 1}, {4, 5, 2}, {5, 0, 1}}};
static const Graph LocalizedPhenol = Graph{Nodes{6, 6, 6, 6, 6, 6, 8}, Edges{{0, 1, 2}, {1, 2, 1}, {2, 3, 2}, {3, 4, 1}, {4, 5, 2}, {5, 6, 1}, {5, 0, 1}}};

static const Graph FourMemberExampleCycle = get_simple_cycle(8, 6, 6, 7);
static const Graph FourMemberExampleChain = get_simple_cycle(8, 6, 6, 7);

inline auto generic_test_data()
{
  return testing::Combine(testing::Values(SimpleFourMemberChain, MixedFourMemberChain, SimpleThreeMemberChain, MixedSixMemberChain),
                          testing::Values(SimpleFourMemberChain, MixedFourMemberChain, SimpleThreeMemberChain, MixedSixMemberChain));
}


inline auto generic_disconnected_test_data()
{
  return testing::Values(std::make_tuple(SimpleFourMemberChain, SimpleThreeMemberChain),
                         std::make_tuple(MixedFourMemberChain, MixedSixMemberChain));
}

inline auto example_test_data()
{
  return testing::Values(std::make_tuple(SimpleCarboxylicAcidLike, SimpleAmideBondLike));
}

} // namespace Data

namespace {

constexpr size_t Mappable = static_cast<size_t>(RIMACS::NodeMappingStatus::Mappable);
constexpr size_t Unmappable = static_cast<size_t>(RIMACS::NodeMappingStatus::Unmappable);


template<typename CompatibleNodesType>
size_t get_compatible_node_idx(
    const std::vector<std::vector<CompatibleNodesType>>& compNodes,
    size_t node)
{
  size_t nofCompatibleNodes = compNodes.at(node).size();
  const auto& selectedCompNodeWeightPair = compNodes.at(node).at(nofCompatibleNodes / 2);
  return RIMACS::getCompatibleNode(selectedCompNodeWeightPair);
}


} // namespace
