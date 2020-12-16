
#pragma once

#include <RIMACS/Forward.hpp>
#include <RIMACS/Functors/Functors.hpp>


namespace RIMACS {

struct SparseMatrixCache;

/**
 * @brief The CachingAdjacencyFunctor class
 *  This functor provides efficient access to the adjacency of a graph.
 *  It is used since most graphs in the Naomi library, e.g. Molecule or SmartsGraph
 *  are designed to be expressive storing nodes (e.g. Atoms) and edges (e.g. Bonds) explicitly
 *  and referring onto each other using pointers.
 *  Those indirections can be erased storing the adjacency in a matrix.
 */
class CachingAdjacencyFunctor : public AdjacencyFunctor {
public:
  ~CachingAdjacencyFunctor() override;

  // GraphCachingAdjacencyFunctor has to be unmovable, since it might hold a reference
  // to its implementation located in a derived object. Moving would invalidate the reference
  // If derived classes implement copy or assignment operators, they do handle this problem
  CachingAdjacencyFunctor(CachingAdjacencyFunctor&&) = delete;
  CachingAdjacencyFunctor(const CachingAdjacencyFunctor&) = delete;
  CachingAdjacencyFunctor& operator=(CachingAdjacencyFunctor&&) = delete;
  CachingAdjacencyFunctor& operator=(const CachingAdjacencyFunctor&) = delete;

  CachingAdjacencyFunctor(
      CachingAdjacencyFunctor&& other,
      const CachingRequiredMinimalAdjacencyFunctor& newBaseRef) noexcept;

  bool adjacent(
      MappingIndex node1,
      MappingIndex node2) const override;

  MappingIndex getEdgeId(
      MappingIndex nodeFromId,
      MappingIndex nodeToId) const override;

  MappingIndex edgeGetFromId(
      MappingIndex edgeId) const override;

  MappingIndex edgeGetToId(
      MappingIndex edgeId) const override;

  std::vector<MappingIndex> nodeGetEdges(
      MappingIndex nodeIdx) const override;

  std::vector<MappingIndex> nodeGetNeighbours(
      MappingIndex nodeIdx) const override;

  MappingIndex getNofNodes() const override;

  MappingIndex getNofEdges() const override;

  MappingIndex nodeGetNofEdges(
      MappingIndex nodeIdx) const override;

  ComponentIndex getNodeOrderSuggestion(
      MappingIndex nodeIdx) const override
  {
    return m_nodeOrder[nodeIdx];
  }

  ComponentIndex getMaxNodeOrderSuggestion() const override
  {
    return m_nodeOrder.back();
  }


protected:
  explicit CachingAdjacencyFunctor(
      const CachingRequiredMinimalAdjacencyFunctor& realInstance);

  void doInitialize();

private:

  const CachingRequiredMinimalAdjacencyFunctor& m_instance;
  std::unique_ptr<SparseMatrixCache> m_cache;
  std::vector<ComponentIndex> m_nodeOrder;
};


/**
 * @brief The GraphCachingAdjacencyFunctor class
 *
 * The GraphCachingAdjacencyFunctor variants derive from the CachingAdjacencyFunctor
 * and wrap the template describing the actually wrapped functor.
 * Instances of GraphCachingAdjacencyFunctor can store an actual instance of the template
 * parameter Functor as well as shared_ptr or reference_wrappers to an existing functor.
 * @tparam Functor The real instance of the graph functor.
 */
template<typename Functor=std::reference_wrapper<const AdjacencyFunctor>>
class GraphCachingAdjacencyFunctor : public CachingAdjacencyFunctor {
public:
  static_assert(std::is_base_of<CachingRequiredMinimalAdjacencyFunctor, Functor>::value,
                "Cached functor must be derived of AdjacencyFunctor");

  /* cppcheck does not detect that m_realInstance is initialized in the constructor */
  /* cppcheck-suppress uninitMemberVar */
  explicit GraphCachingAdjacencyFunctor(
      const Functor& f)
    : CachingAdjacencyFunctor(m_realInstance)
    , m_realInstance(f)
  {
    this->doInitialize();
  }

  /* cppcheck-suppress uninitMemberVar */
  explicit GraphCachingAdjacencyFunctor(
      Functor&& f)
    : CachingAdjacencyFunctor(m_realInstance)
    , m_realInstance(std::move(f))
  {
    this->doInitialize();
  }

  /* cppcheck-suppress uninitMemberVar */
  GraphCachingAdjacencyFunctor(
      GraphCachingAdjacencyFunctor&& other) noexcept
    : CachingAdjacencyFunctor(std::move(static_cast<CachingAdjacencyFunctor&>(other)),
                              static_cast<const CachingRequiredMinimalAdjacencyFunctor&>(m_realInstance))
    , m_realInstance(std::move(other.m_realInstance))
  {}

  template<typename... Args,
           typename = typename std::enable_if<std::is_constructible<Functor, Args...>::value>::type>
  /* cppcheck-suppress uninitMemberVar */
  explicit GraphCachingAdjacencyFunctor(
      Args&&... args)
    : CachingAdjacencyFunctor(static_cast<const CachingRequiredMinimalAdjacencyFunctor&>(m_realInstance))
    , m_realInstance(std::forward<Args>(args)...)
  {
    static_assert(std::is_constructible<Functor, Args...>::value,
                  "Real functor has to be constructable with given arguments");
    this->doInitialize();
  }

  const Functor& getBaseFunctor() const
  {
    return m_realInstance;
  }

protected:
  Functor m_realInstance;
};


template<typename Functor>
class GraphCachingAdjacencyFunctor<std::shared_ptr<Functor>> : public CachingAdjacencyFunctor {
public:
  static_assert(std::is_base_of<CachingRequiredMinimalAdjacencyFunctor, Functor>::value,
                "Cached functor must be derived of AdjacencyFunctor");

  /* cppcheck-suppress uninitMemberVar */
  explicit GraphCachingAdjacencyFunctor(
      const std::shared_ptr<Functor>& f)
    : CachingAdjacencyFunctor(*f)
    , m_realInstance(f)
  {
    this->doInitialize();
  }

  /* cppcheck-suppress uninitMemberVar */
  explicit GraphCachingAdjacencyFunctor(std::shared_ptr<Functor>&& f)
      : CachingAdjacencyFunctor(*f)
  {
    m_realInstance = std::move(f);
    this->doInitialize();
  }

private:
  std::shared_ptr<Functor> m_realInstance;
};


template<typename Functor>
class GraphCachingAdjacencyFunctor<std::unique_ptr<Functor>> : public CachingAdjacencyFunctor {
public:
  static_assert(std::is_base_of<CachingRequiredMinimalAdjacencyFunctor, Functor>::value,
                "Cached functor must be derived of AdjacencyFunctor");

  /* cppcheck-suppress uninitMemberVar */
  explicit GraphCachingAdjacencyFunctor(
      std::unique_ptr<Functor>&& f)
    : CachingAdjacencyFunctor(*f)
  {
    m_realInstance = std::move(f);
    this->doInitialize();
  }

private:
  std::unique_ptr<Functor> m_realInstance;
};


template<typename Functor>
class GraphCachingAdjacencyFunctor<std::reference_wrapper<Functor>> : public CachingAdjacencyFunctor {
public:

  explicit GraphCachingAdjacencyFunctor(
      const Functor& f)
    : CachingAdjacencyFunctor(f)
    , m_realInstance(f)
  {
    this->doInitialize();
    static_assert(std::is_base_of<CachingRequiredMinimalAdjacencyFunctor, Functor>::value,
                  "Cached functor must be derived of AdjacencyFunctor");
  }

private:
  const Functor& m_realInstance;
};


/**
 * Several helpers creating GraphCachingAdjacencyFunctors provided some construction arguments
 */
template<typename Functor, typename... Args>
std::shared_ptr<GraphCachingAdjacencyFunctor<Functor>> makeSharedGraphCachingFunctor(
    Args&&... args)
{
  return std::make_shared<GraphCachingAdjacencyFunctor<Functor>>(std::forward<Args>(args)...);
}

template<typename Functor, typename... Args>
std::unique_ptr<GraphCachingAdjacencyFunctor<Functor>> makeUniqueGraphCachingFunctor(
    Args&&... args)
{
  return std::make_unique<GraphCachingAdjacencyFunctor<Functor>>(std::forward<Args>(args)...);
}

template<typename Functor>
std::shared_ptr<GraphCachingAdjacencyFunctor<Functor>> makeSharedReferenceGraphCachingFunctor(
    const Functor& f)
{
  return std::make_shared<GraphCachingAdjacencyFunctor<std::reference_wrapper<const Functor>>>(f);
}

template<typename Functor>
std::unique_ptr<GraphCachingAdjacencyFunctor<std::reference_wrapper<const Functor>>>
makeUniqueReferenceGraphCachingFunctor(
    const Functor& f)
{
  return std::make_unique<GraphCachingAdjacencyFunctor<std::reference_wrapper<const Functor>>>(f);
}

} // namespace RIMACS
