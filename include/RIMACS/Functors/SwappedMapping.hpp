
#pragma once

#include <RIMACS/Forward.hpp>
#include <RIMACS/Structs.hpp>

namespace RIMACS {


/**
 * Helper class for the mappingIsValidFunction of the
 * NodeCompatibilityFunctor.
 * Since the smaller graph is always used as query, the
 * algorithm may swap the order of the graphs.
 * The order of the index arguments for the CompatibilityFunctors
 * is handled by a wrapping functor, but the final test,
 * whether a mapping is valid would need manual swapping of
 * all mapping pairs. This cast-based class just changes the
 * interpretation of the mapped index pairs.
 */
class swapped_mapping {
public:
  MappingIndex m_to;
  MappingIndex m_from;

  static swapped_mapping* wrapMapping(
      MappingPair* p)
  {
    checkTemplateParameters();
    return reinterpret_cast<swapped_mapping*>(p);
  }

  static const swapped_mapping* wrapMapping(
      const MappingPair* p)
  {
    checkTemplateParameters();
    return reinterpret_cast<const swapped_mapping*>(p);
  }

  static swapped_mapping& wrapMapping(
      MappingPair& p)
  {
    return *wrapMapping(&p);
  }

  static const swapped_mapping& wrapMapping(
      const MappingPair& p)
  {
    return *wrapMapping(&p);
  }

private:
  template<typename... Args>
  explicit swapped_mapping(Args&&...)
  {}

  static void checkTemplateParameters() {
    static_assert(sizeof(swapped_mapping) == sizeof(MappingPair), "Pair sizes must be identical");
    static_assert(sizeof(swapped_mapping) == 2 * sizeof(MappingIndex), "Pair size must be twice the member size");
  }
};

using SwappedMappingIterator = typename std::vector<swapped_mapping>::const_iterator;
inline SwappedMappingIterator do_swap(
    const typename std::vector<MappingPair>::const_iterator& it)
{
  return reinterpret_cast<const SwappedMappingIterator&>(it);
}

template<typename T>
bool operator<(
    const swapped_mapping& sp1,
    const swapped_mapping& sp2)
{
  return sp1.m_from != sp2.m_from
              ? sp1.m_from < sp2.m_from
              : sp1.m_to < sp2.m_to;
}

} // namespace RIMACS
