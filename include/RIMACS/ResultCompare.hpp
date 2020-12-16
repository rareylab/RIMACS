
#pragma once

#include <functional>
#include <RIMACS/Forward.hpp>

namespace RIMACS {

namespace Compare {

template<typename Derived, typename Derived2>
struct derived_wrapper : public Derived, public Derived2 {

  derived_wrapper(const derived_wrapper&) = default;
  derived_wrapper(derived_wrapper&&) noexcept = default;
  derived_wrapper& operator=(const derived_wrapper&) = default;
  derived_wrapper& operator=(derived_wrapper&&) noexcept = default;


  template<typename... Args>
  derived_wrapper(
      Args&& ... args)
    : Derived(std::forward<Args>(args)...)
    , Derived2(std::forward<Args>(args)...)
  {}
};


template<typename Derived>
struct derived_wrapper<Derived, Derived> : public Derived {

  derived_wrapper(const derived_wrapper&) = default;
  derived_wrapper(derived_wrapper&&) noexcept = default;
  derived_wrapper& operator=(const derived_wrapper&) = default;
  derived_wrapper& operator=(derived_wrapper&&) noexcept = default;


  template<typename... Args>
  explicit derived_wrapper(
      Args&& ... args)
      : Derived(std::forward<Args>(args)...) {}

};


template<typename Compare, typename EqualCompare = Compare>
struct equivalence_class_compare;


template<template<typename> class Compare,
    template<typename> class EqualCompare>
struct equivalence_class_compare<Compare<MappingIndex>, EqualCompare<MappingIndex>>
    : public derived_wrapper<Compare<MappingIndex>, EqualCompare<MappingIndex>> {

  using EQ = EqualCompare<MappingIndex>;
  using LESS = Compare<MappingIndex>;

  equivalence_class_compare(const equivalence_class_compare&) = default;
  equivalence_class_compare(equivalence_class_compare&&) noexcept = default;
  equivalence_class_compare& operator=(const equivalence_class_compare&) = default;
  equivalence_class_compare& operator=(equivalence_class_compare&&) noexcept = default;


  equivalence_class_compare(
      const std::vector<MappingIndex>* queryEq,
      const std::vector<MappingIndex>* targetEq)
    : m_queryEquivalenceClasses(queryEq)
    , m_targetEquivalenceClasses(targetEq)
  {}


  bool valid() const
  {
    return m_queryEquivalenceClasses && m_targetEquivalenceClasses;
  }


  bool operator()(
      const MappingPair& lhs,
      const MappingPair& rhs) const;

  const std::vector<MappingIndex>* m_queryEquivalenceClasses;
  const std::vector<MappingIndex>* m_targetEquivalenceClasses;
};


using equivalence_class_equal = equivalence_class_compare<std::equal_to<MappingIndex>>;
using equivalence_class_less = equivalence_class_compare<std::less<MappingIndex>, std::equal_to<MappingIndex>>;


template<typename EqualCompare = equivalence_class_equal, typename LessCompare = equivalence_class_less>
struct custom_compare_MCSResult_compare
    : derived_wrapper<EqualCompare, LessCompare> {

  using BaseClass = derived_wrapper<EqualCompare, LessCompare>;


  custom_compare_MCSResult_compare(const custom_compare_MCSResult_compare&) = default;
  custom_compare_MCSResult_compare(custom_compare_MCSResult_compare&&) noexcept = default;
  custom_compare_MCSResult_compare& operator=(const custom_compare_MCSResult_compare&) = default;
  custom_compare_MCSResult_compare& operator=(custom_compare_MCSResult_compare&&) noexcept = default;


  template<typename... Args>
  custom_compare_MCSResult_compare(
      Args&& ... args)
    : BaseClass(std::forward<Args>(args)...)
  {}


  const EqualCompare& eqCompareCast() const
  {
    return *this;
  }


  const LessCompare& lessCompareCast() const
  {
    return *this;
  }


  bool operator()(
      const MCSResult& m1,
      const MCSResult& m2) const;
};


} } // RIMACS::Compare


#include <RIMACS/MCSResult.hpp>
