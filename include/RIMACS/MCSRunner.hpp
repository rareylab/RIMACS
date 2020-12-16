
#pragma once

#include <type_traits>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Config.hpp>
#include <RIMACS/MCSResult.hpp>

#include <RIMACS/Functors/LineGraph.hpp>


namespace RIMACS {

/**
 * @brief The MCSRunner class
 *  This MCSRunner provides methods for an induced MCS named runGenericMCSIncludingCount or
 *  runGenericMCS and methods for the disconnected MCS named runGenericNodeMappingMCES,
 *  runGenericMCES, or runGenericMCESIncludingCount.
 *  The generic (induced) MCS methods differ in return types. Some provide the results and
 *  the total number of backtracking steps, the 'IncludingCount' variants.
 *  Beyond that, the initial mapping parameter differs. Otherwise all induces methods are the same.
 *  The initial mapping can be an vector of mapping pairs, a single mapping pair or it can be omitted.
 *
 *  For the noninduced MCES, initial mappings and weights are unsupported.
 *  The 'NodeMapping' variant transforms the edge based mapping into a node mapping.
 *  The return type difference of the 'IncludingCount' variant is the same as for the (induced) MCS.
 *
 *  In general the input graphs are named query and target, although this assignment is
 *  arbitrarily set. This shall help to understand whats going on and how the mapping is calculated.
 */
class MCSRunner {
public:

  /**
   * @brief runGenericMCSIncludingCount
   *  Compute a MCS between two graphs.
   * @param query Functor describing the query graph
   * @param target Functor describing the target graph
   * @param nodeComp Functor for index based node compatibility
   * @param edgeComp Functor for index based edge compatibility
   * @param initialMapping List of initially mapped nodes
   * @param config Config for the connectedness of the result, result type and runtime cutoffs
   * @param weighted flag whether weights provided by the NodeCompatibility functors should be used
   * @param queryEquivalenceClasses (optional) pointer to a vector of equivalence class indices
   *        for the nodes of the first graph (there must be exactly one entry per node)
   * @param targetEquivalenceClasses (optional) pointer to a vector of equivalence class indices
   *        for the nodes of the second graph (there must be exactly one entry per node)
   * @param ccFeature Flag whether an internal optimization tracking available connected component
   *        in the query graph should be used
   * @param assumeInitialMappingAsCompatible A flag which can skip the test for node compatibility
   *        if exactly one initial node mapping is provided
   * @return the list of @ref MCSResult instances describing the computed mappings and the total
   *         number of iterations
   */
  static std::pair<std::vector<MCSResult>, Index> runGenericMCSIncludingCount(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const std::vector<MappingPair>& initialMapping,
      const Config& config = getDefaultConnectedConfig(),
      bool weighted = false,
      const EquivalenceClasses* queryEquivalenceClasses = nullptr,
      const EquivalenceClasses* targetEquivalenceClasses = nullptr,
      bool ccFeature = true,
      bool assumeInitialMappingAsCompatible = false);


  /**
   * @brief runGenericMCS
   *  Same as runGenericMCSIncludingCount but with different return type.
   * @return the list of @ref MCSResult instances describing the computed mappings
   */
  static std::vector<MCSResult> runGenericMCS(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const std::vector<MappingPair>& initialMapping,
      const Config& config = getDefaultConnectedConfig(),
      bool weighted = false,
      const EquivalenceClasses* queryEquivalenceClasses = nullptr,
      const EquivalenceClasses* targetEquivalenceClasses = nullptr,
      bool ccFeature = true,
      bool assumeInitialMappingAsCompatible = false)
  {
    return std::move(runGenericMCSIncludingCount(query, target, nodeComp, edgeComp, initialMapping, config,
                                                 weighted, queryEquivalenceClasses,
                                                 targetEquivalenceClasses, ccFeature,
                                                 assumeInitialMappingAsCompatible).first);
  }


  /**
   * @brief runGenericMCS
   *  Same as runGenericMCSIncludingCount but with different return type.
   * @return the list of @ref MCSResult instances describing the computed mappings
   */
  static std::vector<MCSResult> runGenericMCS(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const EquivalenceClasses* queryEquivalenceClasses,
      const EquivalenceClasses* targetEquivalenceClasses,
      const Config& config = getDefaultConnectedConfig(),
      bool ccFeature = true)
  {
    return runGenericMCS(query, target, nodeComp, edgeComp, std::vector<MappingPair>{}, config,
                         false, queryEquivalenceClasses, targetEquivalenceClasses, ccFeature, false);
  }


  /**
   * @brief runGenericMCSIncludingCount
   *
   * @param query Functor describing the query graph
   * @param target Functor describing the target graph
   * @param nodeComp Functor for index based node compatibility
   * @param edgeComp Functor for index based edge compatibility
   * @param initialMapping A pair describing the initial node mapping
   * @param config Config for the connectedness of the result, result type and runtime cutoffs
   * @param weighted flag whether weights provided by the NodeCompatibility functors should be used
   * @param queryEquivalenceClasses (optional) pointer to a vector of equivalence class indices
   *        for the nodes of the first graph (there must be exactly one entry per node)
   * @param targetEquivalenceClasses (optional) pointer to a vector of equivalence class indices
   *        for the nodes of the second graph (there must be exactly one entry per node)
   * @param ccFeature Flag whether an internal optimization tracking available connected component
   *        in the query graph should be used
   * @param assumeInitialMappingAsCompatible A flag whether node compatibility should be
   *        ignored for the initial mapping
   * @return the list of @ref MCSResult instances describing the computed mappings and the total
   *         number of iterations
   */
  static std::pair<std::vector<MCSResult>, Index> runGenericMCSIncludingCount(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const MappingPair& initialMapping,
      const Config& config = getDefaultConnectedConfig(),
      bool weighted = false,
      const EquivalenceClasses* queryEquivalenceClasses = nullptr,
      const EquivalenceClasses* targetEquivalenceClasses = nullptr,
      bool ccFeature = true,
      bool assumeInitialMappingAsCompatible = false);


  /**
   * @brief runGenericMCS
   *  Same as runGenericMCSIncludingCount but with different return type.
   * @return the list of @ref MCSResult instances describing the computed mappings
   */
  static std::vector<MCSResult> runGenericMCS(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const MappingPair& initialMapping,
      const Config& config = getDefaultConnectedConfig(),
      bool weighted = false,
      const EquivalenceClasses* queryEquivalenceClasses = nullptr,
      const EquivalenceClasses* targetEquivalenceClasses = nullptr,
      bool ccFeature = true,
      bool assumeInitialMappingAsCompatible = false)
  {
    return std::move(runGenericMCSIncludingCount(query, target, nodeComp, edgeComp, initialMapping, config,
                                                 weighted, queryEquivalenceClasses,
                                                 targetEquivalenceClasses, ccFeature,
                                                 assumeInitialMappingAsCompatible).first);
  }


  /**
   * @brief runGenericMCSIncludingCount
   *
   * @param query Functor describing the query graph
   * @param target Functor describing the target graph
   * @param nodeComp Functor for index based node compatibility
   * @param edgeComp Functor for index based edge compatibility
   * @param config Config for the connectedness of the result, result type and runtime cutoffs
   * @param weighted flag whether weights provided by the NodeCompatibility functors should be used
   * @param queryEquivalenceClasses (optional) pointer to a vector of equivalence class indices
   *        for the nodes of the first graph (there must be exactly one entry per node)
   * @param targetEquivalenceClasses (optional) pointer to a vector of equivalence class indices
   *        for the nodes of the second graph (there must be exactly one entry per node)
   * @param ccFeature Flag whether an internal optimization tracking available connected component
   *        in the query graph should be used
   * @return the list of @ref MCSResult instances describing the computed mappings and the total
   *         number of iterations
   */
  static std::pair<std::vector<MCSResult>, Index> runGenericMCSIncludingCount(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const Config& config = getDefaultConnectedConfig(),
      bool weighted = false,
      const EquivalenceClasses* queryEquivalenceClasses = nullptr,
      const EquivalenceClasses* targetEquivalenceClasses = nullptr,
      bool ccFeature = true);

  /**
   * @brief runGenericMCS
   *  Same as runGenericMCSIncludingCount but with different return type.
   * @return the list of @ref MCSResult instances describing the computed mappings
   */
  static std::vector<MCSResult> runGenericMCS(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const Config& config = getDefaultConnectedConfig(),
      bool weighted = false,
      const std::vector<MappingIndex>* queryEquivalenceClasses = nullptr,
      const std::vector<MappingIndex>* targetEquivalenceClasses = nullptr,
      bool ccFeature = true)
  {
    return std::move(runGenericMCSIncludingCount(query, target, nodeComp, edgeComp, config,
                                                 weighted, queryEquivalenceClasses,
                                                 targetEquivalenceClasses, ccFeature).first);
  }


  /**
   * @brief runGenericMCESIncludingCount
   *
   * @param query Functor describing the query graph
   * @param target Functor describing the target graph
   * @param nodeComp Functor for index based node compatibility
   * @param edgeComp Functor for index based edge compatibility
   * @param config Config for the connectedness of the result, result type and runtime cutoffs
   * @param ccFeature Flag whether an internal optimization tracking available connected component
   *        in the query graph should be used
   * @return the list of @ref MCSResult instances describing the computed edge mappings and the total
   *         number of iterations
   */
  static std::pair<std::vector<MCSResult>, Index> runGenericMCESIncludingCount(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const Config& config = getDefaultConnectedConfig(),
      bool ccFeature = true);


  /**
   * @brief runGenericMCES
   *  Same as runGenericMCESIncludingCount but with different return type.
   * @return the list of @ref MCSResult instances describing the computed edge mappings
   * @return
   */
  static std::vector<MCSResult> runGenericMCES(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const Config& config = getDefaultConnectedConfig(),
      bool ccFeature = true)
  {
    return runGenericMCESIncludingCount(query, target, nodeComp, edgeComp, config,
                                        ccFeature).first;
  }


  /**
   * @brief runGenericNodeMappingMCES
   *  Compute a noninduced MCES between two graphs, but provide a node mapping as result.
   *
   * @param query Functor describing the query graph
   * @param target Functor describing the target graph
   * @param nodeComp Functor for index based node compatibility
   * @param edgeComp Functor for index based edge compatibility
   * @param config Config for the connectedness of the result, result type and runtime cutoffs
   * @param ccFeature Flag whether an internal optimization tracking available connected component
   *        in the query graph should be used
   * @return the list of @ref MCSResult instances describing the node mappings and the total
   *         number of iterations
   *         The node mappings are the generated transforming edge mappings to node mappings
   */
  static std::tuple<std::vector<MCSResult>, Index, double> runGenericNodeMappingMCES(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const Config& config = getDefaultConnectedConfig(),
      bool ccFeature = true);


  /**
   * @brief runHintedMCS run the MCS algorithm providing a hint for the node compatibility.
   *                     The hint is used to restrict nodeCompatibility, but the algorithm
   *                     will also ensure node compatibility as provided by the functor
   * @param query Functor describing the query graph
   * @param target Functor describing the target graph
   * @param nodeComp Functor for index based node compatibility
   * @param edgeComp Functor for index based edge compatibility
   * @param config Config for the connectedness of the result, result type and runtime cutoffs
   * @param compatibilityHint The hint for the node compatibility. Only nodes that are assigned
   *                          compatible by this hint are potentially mapped
   * @param config Config for the connectedness of the result, result type and runtime cutoffs
   * @param initialMapping List of initially mapped nodes
   * @param weighted flag whether weights provided by the NodeCompatibility functors should be used
   * @param queryEquivalenceClasses (optional) pointer to a vector of equivalence class indices
   *        for the nodes of the first graph (there must be exactly one entry per node)
   * @param targetEquivalenceClasses (optional) pointer to a vector of equivalence class indices
   *        for the nodes of the second graph (there must be exactly one entry per node)
   * @return the list of @ref MCSResult instances describing the node mappings and the total
   *         number of iterations
   *         The node mappings are the generated transforming edge mappings to node mappings
   */
  static std::pair<MCSResults, Index> runHintedMCS(
      const AdjacencyFunctor& query,
      const AdjacencyFunctor& target,
      const NodeCompatibilityFunctor& nodeComp,
      const EdgeCompatibilityFunctor& edgeComp,
      const CompatibilityHint& compatibilityHint,
      const Config& config = getDefaultConfig(),
      const std::vector<MappingPair>& initialMapping = std::vector<MappingPair>(),
      bool weighted = false,
      const EquivalenceClasses* queryEquivalenceClasses = nullptr,
      const EquivalenceClasses* targetEquivalenceClasses = nullptr);
};

} // namespace RIMACS
