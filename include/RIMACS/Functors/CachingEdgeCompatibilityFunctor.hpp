
#pragma once

#include <RIMACS/Functors/Functors.hpp>

namespace RIMACS {

class CachingEdgeCompatibilityFunctor : public EdgeCompatibilityFunctor {
public:

  CachingEdgeCompatibilityFunctor(
      const EdgeCompatibilityFunctor& base,
      MappingIndex g1NofEdges,
      MappingIndex g2NofEdges)
    : m_base(base)
    , m_cache(2 * g1NofEdges * g2NofEdges, false)
    , m_g1NofEdges(g1NofEdges)
  {}

  bool edgesAreCompatible(
      MappingIndex edge1,
      MappingIndex edge2) const override
  {
    auto edge_cache_it = m_cache.begin() + 2 * (m_g1NofEdges * edge2 + edge1);
    RIMACS_ASSERT(edge_cache_it < m_cache.end());
    if(*edge_cache_it) {
      return edge_cache_it[1];
    }
    *edge_cache_it = true;
    return (edge_cache_it[1] = m_base(edge1, edge2));
  }

private:
  const EdgeCompatibilityFunctor& m_base;
  mutable std::vector<bool> m_cache;
  MappingIndex m_g1NofEdges;
};

} // namespace RIMACS
