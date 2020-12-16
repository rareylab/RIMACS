
#pragma once

#include <type_traits>
#include <RIMACS/Forward.hpp>

namespace RIMACS {

/*
 * In this file, type handling for the different types according to the weighted and unweighted
 * MCS calculation is implemented.
 * Whereas in the weighted case it is necessary to carry the component weight along, the
 * unweighted calculations don't need this kind of information and are more efficient in memory
 * and sometimes in the amount of computation since the weights are not updated since they are
 * not present.
 */


template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, const T&> getCompatibleNode(
    const std::pair<T, double>& compNode)
{
  return compNode.first;
}
template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, T&> getCompatibleNode(
    std::pair<T, double>& compNode)
{
  return compNode.first;
}
template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, const T&> getCompatibleNode(
    const T& compNode)
{
  return compNode;
}
template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, T&> getCompatibleNode(
    T& compNode)
{
  return compNode;
}

template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, T&> getCompatibleNode(
    std::vector<T>& compNodes,
    MappingIndex i)
{
  return compNodes[i];
}
template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, const T&> getCompatibleNode(
    const std::vector<T>& compNodes,
    MappingIndex i)
{
  return compNodes[i];
}

template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, T&> getCompatibleNode(
    std::vector<std::pair<T, double>>& compNodes,
    MappingIndex i)
{
  return compNodes[i].first;
}
template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, const T&> getCompatibleNode(
    const std::vector<std::pair<T, double>>& compNodes,
    MappingIndex i)
{
  return compNodes[i].first;
}


template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, double&> getCompatibleNodeWeight(
    std::pair<T, double>& compNode)
{
  return compNode.second;
}
template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, const double&> getCompatibleNodeWeight(
    const std::pair<T, double>& compNode)
{
  return compNode.second;
}

template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, const double&> getCompatibleNodeWeight(
    const std::vector<std::pair<T, double>>& compNodes,
    MappingIndex i)
{
  return compNodes[i].second;
}
template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, double&> getCompatibleNodeWeight(
    std::vector<std::pair<T, double>>& compNodes,
    MappingIndex i)
{
  return compNodes[i].second;
}


template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, double> getCompatibleNodeWeight(
    const T& /*compNode*/)
{
  RIMACS_ASSERT(false && "Never call this weight wrapper on unweighted mappings");
  return 0.0;
}
template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, double&> getCompatibleNodeWeight(
    T& compNode)
{
  RIMACS_ASSERT(false && "Never call this weight wrapper on unweighted mappings");
  return reinterpret_cast<double&>(compNode);
}

template<typename T,
         typename=std::enable_if_t<
           sizeof(decltype(getCompatibleNode(*static_cast<std::remove_reference_t<T>*>(nullptr)))) != 0>>
inline auto getComponentSize(T&& value) -> decltype(getCompatibleNode(value))
{
  return getCompatibleNode(value);
}

template<typename T,
         typename=std::enable_if_t<
           sizeof(decltype(getCompatibleNode(*static_cast<std::remove_reference_t<T>*>(nullptr), MappingIndex{0})))
           != 0>>
inline auto getComponentSize(
    T&& value,
    MappingIndex idx) -> decltype(getCompatibleNode(value, idx))
{
  return getCompatibleNode(value, idx);
}


template<typename T,
         typename=std::enable_if_t<
           sizeof(decltype(getCompatibleNodeWeight(*static_cast<std::remove_reference_t<T>*>(nullptr))))
           != 0 && !std::is_integral<T>::value>>
inline auto getComponentWeight(
    T&& value) -> decltype(getCompatibleNodeWeight(value))
{
  return getCompatibleNodeWeight(value);
}


template<typename T>
inline std::enable_if_t<std::is_integral<T>::value, T&> getComponentWeight(
    T& size)
{
  return size;
}

template<typename T,
         typename=std::enable_if_t<
           sizeof(decltype(getCompatibleNodeWeight(*static_cast<std::remove_reference_t<T>*>(nullptr), MappingIndex{0})))
           != 0>>
inline auto getComponentWeight(
    T&& value,
    MappingIndex idx) -> decltype(getCompatibleNodeWeight(value, idx))
{
  return getCompatibleNodeWeight(value, idx);
}

inline ComponentIndex& getComponentWeight(
    std::vector<ComponentIndex>& sizes,
    MappingIndex idx)
{
  return sizes[idx];
}

inline const ComponentIndex& getComponentWeight(
    const std::vector<ComponentIndex>& sizes,
    MappingIndex idx)
{
  return sizes[idx];
}


} // namespace RIMACS
