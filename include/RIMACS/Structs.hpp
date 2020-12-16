
#pragma once

#include "Forward.hpp"

namespace RIMACS {

struct MappingPair {

  MappingPair() = default;

  MappingPair(
      MappingIndex from,
      MappingIndex to)
    : m_from(from)
    , m_to(to)
  {}

  bool operator==(
      const MappingPair& other) const
  {
    return m_from == other.m_from
           && m_to == other.m_to;
  }

  bool operator!=(
      const MappingPair& other) const
  {
    return !(*this == other);
  }

  bool operator<(
      const MappingPair& other) const
  {
    return m_from != other.m_from ? m_from < other.m_from
                                  : m_to < other.m_to;
  }

  bool operator>(
      const MappingPair& other) const
  {
    return other < *this;
  }

  MappingIndex m_from = s_invalidId;
  MappingIndex m_to = s_invalidId;
private:
  static constexpr unsigned s_invalidId = std::numeric_limits<MappingIndex>::max();
};

struct MappingTuple : MappingPair {

  MappingTuple() = default;

  MappingTuple(
      MappingIndex from,
      MappingIndex to,
      double weight)
    : MappingPair(from, to)
    , m_weight(weight)
  {}

  MappingTuple(
      MappingPair base)
    : MappingPair(base)
    , m_weight(1.0)
  {}

  bool operator==(
      const MappingTuple& other) const
  {
    return MappingPair::operator==(other)
           && m_weight == other.m_weight;
  }

  bool operator!=(
      const MappingTuple& other) const
  {
    return !(*this == other);
  }

  bool operator<(
      const MappingTuple& other) const
  {
    return MappingPair::operator!=(other) ? MappingPair::operator<(other)
                                          : m_weight < other.m_weight;
  }

  bool operator>(
      const MappingTuple& other) const
  {
    return other < *this;
  }

  double m_weight = s_noWeight;
private:
  static constexpr double s_noWeight = std::numeric_limits<double>::signaling_NaN();
};

namespace {

inline void set_weight(
    const MappingPair&,
    double)
{}

inline void set_weight(
    MappingTuple& mt,
    const double& weight)
{
  mt.m_weight = weight;
}

} // namespace

struct AdjacentNode {

  AdjacentNode(
      MappingIndex node,
      MappingIndex nofPartners,
      MappingIndex nofMappedNeighbours,
      double weight = s_noWeight,
      double initialOrderWeight = s_noWeight)
    : m_node(node),
      m_nofPartners(nofPartners),
      m_nofMappedNeighbours(nofMappedNeighbours),
      m_weight(weight),
      m_initialOrderWeight(initialOrderWeight)
  {}

  bool operator==(
      const AdjacentNode& other) const
  {
    return m_node == other.m_node
           && m_nofPartners == other.m_nofPartners
           && m_nofMappedNeighbours == other.m_nofMappedNeighbours
           && m_weight == other.m_weight
           && m_initialOrderWeight == other.m_initialOrderWeight;
  }

  bool operator<(
      const AdjacentNode& other) const
  {
    if(m_nofMappedNeighbours != other.m_nofMappedNeighbours) {
      return m_nofMappedNeighbours > other.m_nofMappedNeighbours;
    }
    if(m_nofPartners != other.m_nofPartners) {
      return m_nofPartners > other.m_nofPartners;
    }
    if(m_weight != other.m_weight) {
      return m_weight < other.m_weight;
    }
    if(m_initialOrderWeight != other.m_initialOrderWeight) {
      return m_initialOrderWeight < other.m_initialOrderWeight;
    }
    return m_node < other.m_node;
  }

  bool operator>(
      const AdjacentNode& other) const
  {
    return other < *this;
  }

  MappingIndex m_node;
  MappingIndex m_nofPartners;
  MappingIndex m_nofMappedNeighbours;
  double m_weight = s_noWeight;
  double m_initialOrderWeight = s_noWeight;
private:
  static constexpr double s_noWeight = std::numeric_limits<double>::signaling_NaN();
};

} // namespace RIMACS
