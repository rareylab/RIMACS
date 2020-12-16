
#pragma once

#include <functional>
#include <RIMACS/Forward.hpp>

namespace RIMACS {

template<bool connected>
struct ConnectionSizeAccumulator {
  ComponentIndex operator()(
      const ComponentIndex& initial,
      const ComponentIndex& value)
  {
    return connected ? initial + value : initial - value;
  }

  std::pair<ComponentIndex, double> operator()(
      const std::pair<ComponentIndex, double>& initial,
      const std::pair<ComponentIndex, double>& value)
  {
    return connected ? std::pair<ComponentIndex, double>(initial.first + value.first,
                                                         initial.second + value.second)
                     : std::pair<ComponentIndex, double>(initial.first - value.first,
                                                         initial.second - value.second);
  }
};

struct WeightedNodeCompare {
  bool operator()(
      const std::pair<ComponentIndex, double>& p1,
      const std::pair<ComponentIndex, double>& p2)
  {
    if(p1.second != p2.second) {
      return p2.second < p1.second;
    }
    return p2.first < p1.first;
  }

  bool operator()(
      const ComponentIndex& i,
      const ComponentIndex& j)
  {
    return j < i;
  }
};


struct mapping_identifier {
  template<typename MappingType>
  void check_template_parameter(
      const MappingType& value) const
  {
    static_assert(std::is_same<decltype(value.m_from),
                               decltype(value.m_to)>::value
                  && std::is_same<decltype(value.m_to),
                                  decltype(value.m_to)>::value
                  && sizeof(MappingType) == sizeof(decltype(value.m_from)) + sizeof(decltype(value.m_to)),
                  "Mapping Type argument must have m_from and m_to members");
  }
};

struct mapping_from : private mapping_identifier {
  template<typename MappingType>
  MappingIndex operator()(
      const MappingType& value) const
  {
    check_template_parameter(value);
    return value.m_from;
  }
};


struct mapping_from_compare : private mapping_identifier {
  template<typename MappingType>
  bool operator()(
      const MappingType& lhs,
      const MappingType& rhs) const
  {
    check_template_parameter(lhs);
    return lhs.m_from < rhs.m_from;
  }
};

struct mapping_to_compare : private mapping_identifier {
  template<typename MappingType>
  bool operator()(
      const MappingType& lhs,
      const MappingType& rhs) const
  {
    check_template_parameter(lhs);
    return lhs.m_to < rhs.m_to;
  }
};

} // namespace RIMACS
