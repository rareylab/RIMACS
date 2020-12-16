
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <tuple>

#include <RIMACS/ASSERT.hpp>

#include "DimacsGraph.hpp"

namespace Dimacs {

std::ostream& operator<<(
    std::ostream& os,
    const DimacsGraph& g)
{
  os << "Dimacs Graph named '" << g.name << "' with " << g.m_nodeLabels.size() << " nodes and " << g.m_edges.size() << " edges." << std::endl;
  for(const int nodeLabel : g.m_nodeLabels) {
    os << nodeLabel << '\t';
  }
  os << std::endl;
  for(const DimacsEdge& e : g.m_edges) {
    os << e.m_from << "," << e.m_to << ":" << e.m_label << '\t';
  }
  os << std::endl;
  return os;
}

DimacsGraph::DimacsGraph(
    const std::string& file)
{
  std::ifstream in(file);
  std::string line;
  bool parsed = false;

  while(std::getline(in, line)) {
    if(line.empty() || line[0] == 'c') {
      continue;
    }
    std::stringstream lineBuffer;
    lineBuffer << line;
    std::string type;
    lineBuffer >> type;
    switch(type[0]) {
      case('p'): {
        std::string problem;
        int nofEdges;
        int nofNodes;
        lineBuffer >> problem >> nofNodes >> nofEdges;
        m_nodeLabels.reserve(nofNodes);
        m_edges.reserve(nofEdges);
        RIMACS_ASSERT(nofNodes > 0);
        RIMACS_ASSERT(nofEdges > 0);
        break;
      } case(static_cast<int>('n')): {
        int label;
        lineBuffer >> label;
        m_nodeLabels.push_back(label);
        RIMACS_ASSERT(label > 0);
        parsed = true;
        break;
      } case(static_cast<int>('e')): {
        int from = 0, to = 0, label = 0;
        lineBuffer >> from >> to >> label;
        if(to < from) {
          std::swap(from, to);
        }
        m_edges.push_back(DimacsEdge{from - 1, to - 1, label});
        RIMACS_ASSERT(0 < from);
        RIMACS_ASSERT(to <= static_cast<int>(m_nodeLabels.size()));
        RIMACS_ASSERT(0 < label);
        parsed = true;
        break;
      }
      default: {
        std::cout << line << std::endl;
        RIMACS_ASSERT(false && "Unknown line in dimacs file, abort");
      }
    }
  }
  RIMACS_ASSERT(parsed && "The file does not contain any graph");
}

DimacsGraph::DimacsGraph(
    const std::vector<int>& nodes,
    const std::vector<DimacsEdge>& edges)
  : m_nodeLabels(nodes)
  , m_edges(edges)
{
  std::for_each(m_edges.begin(), m_edges.end(), [](auto& e) {e.m_label = std::max(1, e.m_label);});
}

DimacsGraph::DimacsGraph(
    const std::vector<int>& nodes,
    const std::vector<std::tuple<int, int, int>>& edges)
    : m_nodeLabels(nodes)
{
  m_edges.reserve(edges.size());
  std::transform(edges.begin(), edges.end(), std::back_inserter(m_edges),
      [](const std::tuple<int, int, int>& t) {
    return DimacsEdge(std::get<0>(t), std::get<1>(t), std::max(1, std::get<2>(t)));
  });
}

// If an IDE has problems with the representation of native std::string
const char* DimacsGraph::get_name() const {return name.c_str();}

} // namespace Dimacs
