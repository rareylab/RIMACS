
#pragma once

#include <algorithm>

#include <RIMACS/ASSERT.hpp>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Structs.hpp>
#include <RIMACS/ResultCompare.hpp>


namespace RIMACS {

/**
 * @brief The MCSResult class
 */
class MCSResult {
public:

  using iterator = typename std::vector<MappingPair>::const_iterator;

  MCSResult(const MCSResult&) = default;
  MCSResult(MCSResult&&) noexcept = default;
  MCSResult& operator=(const MCSResult&) = default;
  MCSResult& operator=(MCSResult&&) noexcept = default;

  MCSResult() noexcept = default;


  MCSResult(
      const std::vector<MappingPair>& nodeMapping,
      const std::vector<MappingIndex>& componentSizes,
      double weight,
      MappingIndex initialMappingSize,
      const std::vector<MappingIndex>* queryEquivalenceClasses = nullptr,
      const std::vector<MappingIndex>* targetEquivalenceClasses = nullptr)
    : MCSResult(nodeMapping, componentSizes, weight)
  {
    if(queryEquivalenceClasses && targetEquivalenceClasses) {
      sortResults(initialMappingSize, Compare::equivalence_class_less{queryEquivalenceClasses,
                                                                      targetEquivalenceClasses});
    } else {
      sortResults(initialMappingSize);
    }
  }


  MCSResult(
      std::vector<MappingPair>&& nodeMapping,
      std::vector<MappingIndex>&& componentSizes,
      double weight)
    : m_nodeMapping(std::move(nodeMapping))
    , m_cumulatedComponentSizes(std::move(componentSizes))
    , m_weight(weight)
  {}


  template<typename Compare>
  MCSResult(
      const std::vector<MappingPair>& nodeMapping,
      const std::vector<MappingIndex>& componentSizes,
      double weight,
      MappingIndex initialMappingSize,
      const Compare& compare)
    : MCSResult(nodeMapping, componentSizes, weight)
  {
    sortResults(initialMappingSize, compare);
  }

  /**
   * @brief valid
   * @return true, if there is any result, false otherwise
   */
  bool valid() const
  {
    return !m_cumulatedComponentSizes.empty();
  }


  /**
   * @brief size
   * @return number of mapped nodes or edges, depending on the performed algorithm
   */
  Index size() const
  {
    return m_nodeMapping.size();
  }


  /**
   * @brief weight
   * @return returns the weight of the MCS mapping.
   */
  const double& weight() const
  {
    return m_weight;
  }


  /**
   * @brief isDisconnected
   * @return true, if the mapping consists of more than one connected component, false otherwise
   */
  bool isDisconnected() const
  {
    return m_cumulatedComponentSizes.size() > 1;
  }


  /**
   * @brief begin
   * @return Iterator to the first mapping index pair
   */
  iterator begin() const
  {
    return m_nodeMapping.begin();
  }

  /**
   * @brief end
   * @return End iterator of the mapping index pairs
   */
  iterator end() const
  {
    return m_nodeMapping.end();
  }


  /**
   * @brief getMappings
   * @return The underlying mapping container
   */
  const std::vector<MappingPair>& getMappings() const
  {
    return m_nodeMapping;
  }


  /**
   * @brief getComponentSizes
   * @return The number of connected components of the mapping
   */
  const std::vector<MappingIndex>& getComponentSizes() const
  {
    return m_cumulatedComponentSizes;
  }


  /**
   * @brief operator <
   *  Compare one mapping to another for mapping unification
   *  The operator < is based on the number of components, their size and the node mappings
   * @param other
   * @return true, if the mapping is considered smaller, false otherwise
   */
  bool operator<(
      const MCSResult& other) const
  {
    return operatorLess(other);
  }


  /**
   * @brief operatorLess
   *  Template based implementation of the operator<
   *  It accepts the additional template parameters for an equality and less compare
   *  functor that are used to compare the node (or edge) index mappings.
   * @tparam EqualCompare operator== like struct for MappingPair instances
   * @tparam LessCompare operator< like struct for MappingPair instances
   * @param other
   * @param equal Instance of the equality compare operator
   * @param less Instance of the less compare operator
   * @return true, if the mapping is considered smaller, false otherwise
   */
  template<typename EqualCompare = std::equal_to<MappingPair>,
           typename LessCompare = std::less<MappingPair>>
  inline bool operatorLess(
      const MCSResult& other,
      const EqualCompare& equal = EqualCompare(),
      const LessCompare& less = LessCompare()) const;


  /**
   * @brief operator ==
   * @param other
   * @return true if the mappings are considered equal, false otherwise
   */
  bool operator==(
      const MCSResult& other) const
  {
    return operatorEqual(other);
  }


  /**
   * @brief operatorEqual
   *  Template based implementation of the operator==
   *  It accepts the additional template parameter for an equality
   *  functor that is used to compare the node (or edge) index mappings.
   * @tparam Compare operator== like struct for MappingPair instances
   * @param other
   * @param compare Instance of the equality compare operator
   * @return true, if the mapping is considered smaller, false otherwise
   */
  template<typename Compare = std::equal_to<MappingPair>>
  inline bool operatorEqual(
      const MCSResult& other,
      const Compare& compare = Compare()) const;

  /**
   * @brief operatorEqual
   *  Implementation of the operator== with additional provided equivalence class lists
   *  Mappings (1,2) and (3,4) are considered equal, if the equivalence classes for node 1 and 3 of
   *  the first graph are equal and the equivalence classes for the nodes 2 and 4 of the second graph
   *  are equal.
   * @param other
   * @param queryEquivalenceClasses Equivalence classes for nodes of the first graph
   * @param targetEquivalenceClasses Equivalence classes for nodes of the second graph
   * @return true, if the mappings are considered equal using the equivalence class mapping,
   *         false otherwise
   */
  inline bool operatorEqual(
      const MCSResult& other,
      const std::vector<MappingIndex>& queryEquivalenceClasses,
      const std::vector<MappingIndex>& targetEquivalenceClasses) const;

  /**
   * @brief orderDisconnected
   *  Reorder the internal mapping, such that the connected components are ordered by
   *  increasing size. If the components have equal size, then the mapping indices are
   *  compared to determine an unique order.
   */
  inline void orderDisconnected();


  /**
   * @brief doSwap
   *  Swap all mapping pairs. Then sort them.
   * @tparam Compare operator== like struct for MappingPair instances
   * @param initialMappingSize Number of nodes that are not sorted, because the mapping was
   *        provided as a constraint
   * @param compare Instance of the equality compare operator
   */
  template<typename Compare = std::less<>>
  void doSwap(
      const MappingIndex& initialMappingSize = MappingIndex{0},
      const Compare& compare = Compare())
  {
    for(MappingPair& mp : m_nodeMapping) {
      std::swap(mp.m_from, mp.m_to);
    }
    m_cumulatedComponentSizes.insert(m_cumulatedComponentSizes.begin(), initialMappingSize);
    sortResults<false, Compare>(initialMappingSize, compare);
  }

  /**
   * Provide a mutable reference to the internal mapping.
   * Used for adaptors of the RIMACS library to other contexts
   */
  std::vector<MappingPair>& getMutableMapping()
  {
    return m_nodeMapping;
  }

private:
  friend class MCSRunner;
  template<bool weighted>
  friend class SubgraphMapping;


  // Private constructor setting all values
  MCSResult(
      const std::vector<MappingPair>& nodeMapping,
      const std::vector<MappingIndex>& componentSizes,
      double weight)
    : m_nodeMapping(nodeMapping)
    , m_weight(weight)
  {
    m_cumulatedComponentSizes.reserve(componentSizes.size() + 1);
    m_cumulatedComponentSizes.push_back(0);
    m_cumulatedComponentSizes.insert(m_cumulatedComponentSizes.end(),
                                     componentSizes.rbegin(), componentSizes.rend());
  }


  /**
   * @brief sortResults
   *  Sort the node index mapping per connected component.
   * @tparam doCumulate Template flag indicating whether the component sizes member represents
   *         actual sizes or already the accumulated indices that can be used for sorting.
   * @tparam Compare operator== like struct for MappingPair instances
   * @param initialMappingSize Number of nodes that are not sorted, because the mapping was
   *        provided as a constraint
   * @param compare Instance of the equality compare operator
   */
  template<bool doCumulate=true, typename Compare=std::less<MappingPair>>
  inline void sortResults(
      const MappingIndex& initialMappingSize,
      const Compare& compare = Compare());


  std::vector<MappingPair> m_nodeMapping;
  std::vector<MappingIndex> m_cumulatedComponentSizes;
  double m_weight = std::numeric_limits<double>::quiet_NaN();
};


template<typename EqualCompare,
         typename LessCompare>
inline bool MCSResult::operatorLess(
    const MCSResult& other,
    const EqualCompare& equal,
    const LessCompare& less) const
{
  if(m_cumulatedComponentSizes.size() != other.m_cumulatedComponentSizes.size()) {
    return m_cumulatedComponentSizes.size() < other.m_cumulatedComponentSizes.size();
  }
  for(auto it = m_cumulatedComponentSizes.begin(),
      oIt = other.m_cumulatedComponentSizes.begin(),
      last = m_cumulatedComponentSizes.end()
      ; it != last
      ; ++it, ++oIt) {
    if(*it != *oIt) {
      return *it < *oIt;
    }
  }
  if(m_weight != other.m_weight) {
    return m_weight < other.m_weight;
  }
  RIMACS_ASSERT(other.m_nodeMapping.size() == m_nodeMapping.size());
  for(auto it = m_nodeMapping.begin(), oIt = other.m_nodeMapping.begin(), last = m_nodeMapping.end()
      ; it != last
      ; ++it, ++oIt) {
    if(!equal(*it, *oIt)) {
      // Return false, if this is an equality compare
      return std::is_same<LessCompare, EqualCompare>::value ? false : less(*it, *oIt);
    }
  }
  // Return true, if this is an equality compare
  return std::is_same<LessCompare, EqualCompare>::value;
}


inline bool MCSResult::operatorEqual(
    const MCSResult& other,
    const std::vector<MappingIndex>& queryEquivalenceClasses,
    const std::vector<MappingIndex>& targetEquivalenceClasses) const
{
  return m_cumulatedComponentSizes == other.m_cumulatedComponentSizes
         && m_weight == other.m_weight
         && std::equal(m_nodeMapping.begin(), m_nodeMapping.end(), other.m_nodeMapping.begin(),
                       Compare::equivalence_class_equal{&queryEquivalenceClasses, &targetEquivalenceClasses});
}


template<typename Compare>
inline bool MCSResult::operatorEqual(
    const MCSResult& other,
    const Compare& compare) const
{
  return m_cumulatedComponentSizes == other.m_cumulatedComponentSizes
         && m_weight == other.m_weight
         && std::equal(m_nodeMapping.begin(), m_nodeMapping.end(), other.m_nodeMapping.begin(),
                       compare);
}


inline void MCSResult::orderDisconnected()
{
  std::vector<std::vector<MappingPair>> components;
  m_cumulatedComponentSizes.insert(m_cumulatedComponentSizes.begin(), MappingIndex{0});
  for(auto it = m_cumulatedComponentSizes.begin(),
      next = std::next(m_cumulatedComponentSizes.begin()),
      last = m_cumulatedComponentSizes.end()
      ; next != last
      ; ++next, ++it) {
    components.emplace_back(m_nodeMapping.begin() + *it, m_nodeMapping.begin() + *next);
  }
  std::sort(components.begin(), components.end(), [](auto& v1, auto& v2) {
    for(MappingIndex i = 0, last = std::min(v1.size(), v2.size()); i < last; ++i) {
      if(v1[i] != v2[i]) {
        return v1[i] < v2[i];
      }
    }
    return v1.size() < v2.size();
  });
  m_cumulatedComponentSizes.clear();
  m_nodeMapping.clear();
  m_cumulatedComponentSizes = {0};
  for(const auto& c : components) {
    m_cumulatedComponentSizes.push_back(m_cumulatedComponentSizes.back() + c.size());
    m_nodeMapping.insert(m_nodeMapping.end(), c.begin(), c.end());
  }
  m_cumulatedComponentSizes.erase(m_cumulatedComponentSizes.begin());
}

namespace {

template<typename Compare>
inline void sort_components(
    std::vector<std::vector<MappingPair>>& components,
    const Compare& compare)
{
  auto componentCompare = [&compare] (const std::vector<MappingPair>& lhs,
                                      const std::vector<MappingPair>& rhs) {
      auto lIt = lhs.begin(), lEnd = lhs.end(), rIt = rhs.begin(), rEnd = rhs.end();
      for(; lIt != lEnd && rIt != rEnd; ++lIt, ++rIt) {
        if(compare(*lIt, *rIt)) {
          return true;
        } else if(compare(*rIt, *lIt)) {
          return false;
        }
      }
      return lIt != lEnd && rIt == rEnd;
    };
  std::sort(components.begin(), components.end(), componentCompare);
}

template<typename Compare>
void ensure_components_are_sorted(
    const Compare& compare,
    std::vector<MappingPair>& nodeMapping,
    std::vector<MappingIndex>& componentSizes,
    MappingIndex initialMappingSize)
{
  if(componentSizes.size() > 2u + std::min(initialMappingSize, MappingIndex{1})) {
    std::vector<std::vector<MappingPair>> components;
    size_t nofComponentsToSort = componentSizes.size() - std::min(initialMappingSize, MappingIndex{1}) - 1;
    components.reserve(nofComponentsToSort);
    for(size_t i = 0; i < nofComponentsToSort; ++i) {
      components.emplace_back(std::next(nodeMapping.begin(), *std::next(componentSizes.rbegin())),
                              std::next(nodeMapping.begin(), *componentSizes.rbegin()));
      componentSizes.pop_back();
      nodeMapping.erase(std::next(nodeMapping.begin(), *componentSizes.rbegin()),
                        nodeMapping.end());
    }
    sort_components(components, compare);
    for(const std::vector<MappingPair>& component : components) {
      nodeMapping.insert(nodeMapping.end(), component.begin(), component.end());
      componentSizes.push_back(nodeMapping.size());
    }
  }
}

} // namespace

template<bool doCumulate, typename Compare>
inline void MCSResult::sortResults(
    const MappingIndex& initialMappingSize,
    const Compare& compare)
{
  // If there is an initial mapping, then the first nodes are kept as part of the initial mapping
  if(doCumulate && initialMappingSize) {
    m_cumulatedComponentSizes.front() += initialMappingSize;
    m_cumulatedComponentSizes.at(1) -= initialMappingSize;
  }
  for(auto it = m_cumulatedComponentSizes.begin() + 1, last = m_cumulatedComponentSizes.end()
      ; it != last
      ; ++it) {
    if(doCumulate) {
      *it += it[-1];
    }
    std::sort(std::next(m_nodeMapping.begin(), it[-1]), std::next(m_nodeMapping.begin(), *it), compare);
  }
  ensure_components_are_sorted(compare, m_nodeMapping, m_cumulatedComponentSizes, initialMappingSize);
  m_cumulatedComponentSizes.erase(m_cumulatedComponentSizes.begin());
}


} // namespace RIMACS

#include <RIMACS/ResultCompareTemplate.hpp>
