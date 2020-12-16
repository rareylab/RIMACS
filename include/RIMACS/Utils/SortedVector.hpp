
#pragma once

#include <algorithm>
#include <vector>
#include <type_traits>

namespace RIMACS {


template<typename T, typename Compare=std::less<void>>
class SortedVector: private Compare {
public:

  using iterator = typename std::vector<T>::const_iterator;
  using mutable_iterator = typename std::vector<T>::iterator;
  using value_type = typename std::vector<T>::value_type;
  using BaseClass = Compare;

  SortedVector(
      const Compare& comp = Compare())
    : BaseClass(comp)
  {}

  SortedVector(
      size_t reserved_size)
  {
    m_container.reserve(reserved_size);
  }

  constexpr operator const std::vector<value_type>& () const
  {
    return m_container;
  }

  const T& at(
      size_t idx) const
  {
    return m_container.at(idx);
  }

  const T& operator[](
      size_t idx) const
  {
    return m_container[idx];
  }

  iterator begin() const
  {
    return m_container.begin();
  }

  iterator end() const
  {
    return m_container.end();
  }

  void push_back(
      const value_type& value)
  {
    if(empty() || !key_comp()(value, m_container.back())) {
      m_container.push_back(value);
    } else {
      insert(value);
    }
  }

  void erase(
      iterator it)
  {
    m_container.erase(it);
  }

  void erase(
      const value_type& value)
  {
    iterator pos = lower_bound(value);
    if(pos != end() && !key_comp()(value, *pos)) {
      m_container.erase(pos);
    }
  }

  void reserve(
      size_t n)
  {
    m_container.reserve(n);
  }

  bool empty() const
  {
    return m_container.empty();
  }

  void clear()
  {
    m_container.clear();
  }

  SortedVector& make_unique()
  {
    const Compare& comp = key_comp();
    m_container.erase(std::unique(m_container.begin(),
                                  m_container.end(), [&comp](const auto& a, const auto& b) {return !comp(a, b);}),
                      m_container.end());
    return *this;
  }

  iterator insert(
      const T& value)
  {
    mutable_iterator pos = upper_bound(value);
    pos = m_container.insert(pos, value);
    return pos;
  }

  iterator insert(
      T&& value)
  {
    mutable_iterator pos = upper_bound(value);
    pos = m_container.insert(pos, std::move(value));
    return pos;
  }

  template<typename... Args>
  iterator emplace(
      Args&&... args)
  {
    value_type value(std::forward<Args>(args)...);
    mutable_iterator pos = upper_bound(value);
    pos = m_container.insert(pos, std::move(value));
    return pos;
  }


  iterator insert_unique(
      const T& value)
  {
    mutable_iterator pos = lower_bound(value);
    if(pos == m_container.end() || Compare::operator()(value, *pos)) {
      pos = m_container.insert(pos, value);
    }
    return pos;
  }

  iterator insert_unique(
      T&& value)
  {
    mutable_iterator pos = lower_bound(value);
    if(pos == m_container.end() || Compare::operator()(value, *pos)) {
      pos = m_container.insert(pos, std::move(value));
    }
    return pos;
  }

  template<typename... Args>
  iterator emplace_unique(
      Args&&... args)
  {
    value_type value(std::forward<Args>(args)...);
    mutable_iterator pos = lower_bound(value);
    if(pos == m_container.end() || Compare::operator()(value, *pos)) {
      pos = m_container.insert(pos, std::move(value));
    }
    return pos;
  }

  size_t size() const
  {
    return m_container.size();
  }

  bool contains(
      const T& value) const
  {
    return std::binary_search(begin(), end(), value, key_comp());
  }

  size_t count(
      const T& value) const
  {
    auto range_iterators = equal_range(value);
    return std::distance(range_iterators.first, range_iterators.second);
  }

  std::pair<iterator, iterator> equal_range(
      const T& value) const
  {
    return std::equal_range(begin(), end(), value, key_comp());
  }

  const Compare& key_comp() const
  {
    return *this;
  }

  std::vector<T> moveOut()
  {
    return std::move(m_container);
  }

private:

  mutable_iterator lower_bound(
      const T& value)
  {
    return std::lower_bound(m_container.begin(), m_container.end(), value, key_comp());
  }

  mutable_iterator upper_bound(
      const T& value)
  {
    return std::upper_bound(m_container.begin(), m_container.end(), value, key_comp());
  }

  std::vector<T> m_container;
};

} // namespace RIMACS
