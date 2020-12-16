
#pragma once

#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace Dimacs {

struct DimacsEdge {
  DimacsEdge(
      int from,
      int to,
      int label = 0)
    : m_from(from)
    , m_to(to)
    , m_label(label)
  {}

  int m_from = -1;
  int m_to = -1;
  int m_label = 0;
};

inline bool operator==(
    const DimacsEdge& e1,
    const DimacsEdge& e2)
{
  return e1.m_from  == e2.m_from
      && e1.m_to    == e2.m_to
      && e1.m_label == e2.m_label;
}

struct DimacsGraph {
  DimacsGraph(
      const std::string& file);

  DimacsGraph(
      const std::vector<int>& nodeLabels,
      const std::vector<DimacsEdge>& edges);

  DimacsGraph(
      const std::vector<int>& nodeLabels,
      const std::vector<std::tuple<int, int,int>>& edges);

  DimacsGraph() = default;

  const char* get_name() const;

  std::vector<int> m_nodeLabels;
  std::vector<DimacsEdge> m_edges;
  std::string name;
};

std::ostream& operator<<(
    std::ostream& os,
    const DimacsGraph& g);

} // namespace Dimacs
