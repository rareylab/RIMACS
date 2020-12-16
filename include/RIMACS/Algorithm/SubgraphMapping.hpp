
#pragma once

#include <numeric>
#include <set>
#include <RIMACS/Logger.hpp>

#include <RIMACS/ASSERT.hpp>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Functors/Functors.hpp>
#include <RIMACS/Algorithm/VertexMapping.hpp>
#include <RIMACS/Config.hpp>
#include <RIMACS/MCSResult.hpp>
#include <RIMACS/Algorithm/WeightTypeTraits.hpp>
#include <RIMACS/Utils/SortedVector.hpp>


namespace RIMACS {


template<bool weighted>
class SubgraphMapping {
public:
  using ConnectionSize = typename std::conditional<weighted, std::pair<ComponentIndex, double>,
                                                             ComponentIndex>::type;
  using VertexMappingType = VertexMapping<weighted>;

  static constexpr ComponentIndex InitialComponentId = 1;
  static constexpr ComponentIndex InitiallyUnmappableNode = -InitialComponentId;

  SubgraphMapping(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const CompatibilityFunctors& functors,
      const Config& config,
      const std::vector<MappingIndex>* queryEquivalenceClasses = nullptr,
      const std::vector<MappingIndex>* targetEquivalenceClasses = nullptr,
      Base::UserLog* info = nullptr)
    : SubgraphMapping(query, target, functors.first, functors.second, config,
                      queryEquivalenceClasses, targetEquivalenceClasses, info)
  {}

  SubgraphMapping(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const Config& config,
      const std::vector<MappingIndex>* queryEquivalenceClasses = nullptr,
      const std::vector<MappingIndex>* targetEquivalenceClasses = nullptr,
      Base::UserLog* info = nullptr);

  SubgraphMapping(SubgraphMapping&&) noexcept = default;

  virtual ~SubgraphMapping() = default;

  void enableConnectedComponentTrackingFeature(
      bool enabled)
  {
    m_use_cc_feature = enabled;
  }

  void setUserLog(
      Base::UserLog* info)
  {
    m_info = info;
  }

  Base::UserLog* getUserLog() const
  {
    return m_info;
  }

  const Config& getConfig() const
  {
    return m_config;
  }

  const AdjacencyFunctor& getQueryGraph() const
  {
    return m_query;
  }

  const AdjacencyFunctor& getTargetGraph() const
  {
    return m_target;
  }

  const NodeCompatibilityFunctor& getNodeComp() const
  {
    return m_compatibility.first;
  }

  const EdgeCompatibilityFunctor& getEdgeComp() const
  {
    return m_compatibility.second;
  }

  void removeMapping(
      const std::vector<MappingIndex>& locallyUnmappableNodes);


  VertexMappingType addMatchUpdateVM(
      const VertexMappingType& vmm,
      const typename VertexMappingType::RealMappingTuple& nodeMappings,
      std::vector<MappingIndex>& newUnmappableNodes);


  VertexMappingType handleInitialMapping(
      const VertexMappingType& vmm,
      const MappingPair& initialMapping);


  bool prepareDisconnectionUpdateVM(
      VertexMappingType& vm);

  void recountDisconnectedComponents(
      const VertexMappingType& vm);

  void resetDisconnection();


  const Index& getCount() const
  {
    return m_counter;
  }


  bool recursionLimitReached() const
  {
    return m_counter > m_config.recursionLimit();
  }


  bool isExtendable(
      const VertexMappingType& vm) const;


  bool canBeDisconnectedExpanded(
      const VertexMappingType& vm);

  /**
   * Convert the current node mappings into a result and store them in the respective
   * containers. This function takes care not to override existing results with higher
   * weights or more mapped nodes.
   * @return false if the mapping size constrains are violated, true otherwise
   */
  virtual bool addResult();

  void markUnmappableNodes(
      const VertexMappingType& vm);

  std::vector<MCSResult>&& moveMCSResults();


  void setNodeMappable(
      MappingIndex idx)
  {
    m_nodeStatus.at(idx) = static_cast<MappingIndex>(NodeMappingStatus::Mappable);
  }


  void setNodeUnmappable(
      MappingIndex idx,
      VertexMappingType& vm);

  void updateConnectedComponents(
      MappingIndex unmappableNode,
      VertexMappingType& vm,
      bool isMappedNode=false);

  void restoreConnectedComponents(
      MappingIndex unmappableNode);

  bool nodeIsMappable(
      MappingIndex node) const
  {
    return m_nodeStatus[node] == static_cast<MappingIndex>(NodeMappingStatus::Mappable);
  }


private:

  friend class ::SubgraphMappingTester;

  template<bool is_weighted>
  std::enable_if_t<is_weighted, const double&> current_size_impl() const
  {
    return m_currentWeight.back();
  }
  template<bool is_weighted>
  std::enable_if_t<!is_weighted, MappingIndex> current_size_impl() const
  {
    return m_nodeMappings.size();
  }
  auto currentSize() const -> decltype(current_size_impl<weighted>())
  {
    return current_size_impl<weighted>();
  }

  template<typename ConnectedComponentDFSPolicy>
  void connected_component_dfs(
      ConnectedComponentDFSPolicy& policy,
      VertexMappingType& vm,
      MappingIndex unmappableNode);

  template<typename ConnectedComponentDFSPolicy>
  inline void run_connected_component_dfs_for_node(
      const AdjacencyFunctor& graph,
      VertexMapping<weighted>& vm,
      std::vector<MappingIndex>& nodeStack,
      std::vector<bool>& visited,
      const ConnectedComponentDFSPolicy& policy,
      MappingIndex adjacentNode,
      ComponentIndex componentId);

  template<typename ConnectedComponentDFSPolicy>
  inline void  run_dfs(
      const AdjacencyFunctor& graph,
      std::vector<MappingIndex>& nodeStack,
      VertexMapping<weighted>& vm,
      std::vector<bool>& visited,
      const ConnectedComponentDFSPolicy& policy,
      ConnectionSize& componentSize,
      const MappingIndex& nextComponentId,
      bool& hasConnectionToMappedNodes);

  inline void invalidate_connected_component(
      VertexMapping<weighted>& vm,
      const ComponentIndex& nextComponentId);

  inline void update_connected_component_size(
      const ConnectionSize& componentSize,
      const MappingIndex& adjacentNode,
      const ComponentIndex& nextComponentId) ;

  ConnectionSize countExtendableNodes() const;

  void identifyDisconnectedComponents();

  static bool getConnectedComponentFeatureSuggestion(
      const Config& conf)
  {
    return conf.componentSize() < 3 && conf.nofComponents() > 6;
  }


private:

  std::vector<MappingPair> m_nodeMappings;
  std::reference_wrapper<const AdjacencyFunctor> m_query;
  std::reference_wrapper<const AdjacencyFunctor> m_target;
  CompatibilityFunctors m_compatibility;
  std::vector<MappingIndex> m_nodeStatus;
  std::vector<MCSResult> m_mcs;
  SortedVector<MCSResult, Compare::custom_compare_MCSResult_compare<>> m_equivalentMcs;
  std::vector<MappingIndex> m_connectedComponentSizes;
  Config m_config;
  std::vector<double> m_currentWeight;
  Index m_counter = 0;
  Index m_resultCount = 0;
  MappingIndex m_initialMappingSize = Index{0};
  const std::vector<MappingIndex>* m_queryEquivalenceClasses = nullptr;
  const std::vector<MappingIndex>* m_targetEquivalenceClasses = nullptr;
  std::vector<ComponentIndex> m_connectedComponentIds;
  std::vector<ConnectionSize> m_connectedSizes;
  Index m_maxConnectedComponentIndex;
  bool m_use_cc_feature = false;
  mutable bool m_currentExpansionIsValid = true;
  Base::UserLog* m_info = nullptr;
};

} // namespace RIMACS

#include <RIMACS/Algorithm/SubgraphMappingTemplate.hpp>
