
#pragma once

#include <RIMACS/Forward.hpp>
#include <RIMACS/ResultCompare.hpp>

namespace RIMACS {
namespace Compare {


template<typename EqualCompare, typename LessCompare>
bool custom_compare_MCSResult_compare<EqualCompare, LessCompare>::operator()(
    const MCSResult& m1,
    const MCSResult& m2) const
{
  return m1.operatorLess(m2, static_cast<const EqualCompare&>(*this),
                         static_cast<const LessCompare&>(*this));
}


template<template<typename> class Compare, template<typename> class EqualCompare>
bool equivalence_class_compare<Compare<MappingIndex>, EqualCompare<MappingIndex>>::operator()(
    const MappingPair& lhs,
    const MappingPair& rhs) const
{
  if(!EQ::operator()(m_queryEquivalenceClasses->at(lhs.m_from),
                     m_queryEquivalenceClasses->at(rhs.m_from))) {
    return LESS::operator()(m_queryEquivalenceClasses->at(lhs.m_from),
                            m_queryEquivalenceClasses->at(rhs.m_from));
  }
  return LESS::operator()(m_targetEquivalenceClasses->at(lhs.m_to),
                          m_targetEquivalenceClasses->at(rhs.m_to));
}

} } // RIMACS::Compare
