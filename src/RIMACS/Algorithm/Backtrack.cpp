
#include <vector>
#include <limits>

#include <RIMACS/ASSERT.hpp>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Config.hpp>
#include <RIMACS/Algorithm/SubgraphMapping.hpp>
#include <RIMACS/Algorithm/VertexMapping.hpp>

#include <RIMACS/Algorithm/Backtrack.hpp>


namespace RIMACS {

template<bool weighted>
bool backtrack_impl(
    SubgraphMapping<weighted>& mapping,
    VertexMapping<weighted>& vm)
{
  if(LOGGING_IS_ENABLED(mapping.getUserLog(), Base::LogLevel::Debug)
     && mapping.getCount() && !mapping.getCount() % 1000000 ) {
    LOG_DEBUG(mapping.getUserLog(), "Backtrack is at iteration " << mapping.getCount()
              << " of at most " << mapping.getConfig().recursionLimit()
              << "(" << 100 * mapping.getCount() / mapping.getConfig().recursionLimit() << ")");
  }
  if(mapping.recursionLimitReached()) {
    // calculation timed out
    LOG_WARNING(mapping.getUserLog(), "Calculation exceeded the limit of "
                << mapping.getConfig().recursionLimit() << " recursions");
    return false;
  }

  bool extended = false;
  // get next node and remove it from the adjacent nodes
  MappingIndex g0node = vm.getNextNode();
  if(g0node == VertexMapping<weighted>::NODE_END) {
    /* if a subgraph can not be extended the following */
    /* operations are depending on the current backtracking mode */
    if(mapping.canBeDisconnectedExpanded(vm)) {
      if(mapping.prepareDisconnectionUpdateVM(vm)) {
        extended = backtrack_impl(mapping, vm);
      }
      mapping.resetDisconnection();
    }
    if(!extended) {
      extended = mapping.addResult();
    }
    return extended;
  }
  mapping.updateConnectedComponents(g0node, vm, true);
  MappingPair mappingIndices = vm.initIndexTuple(g0node);
  MappingTupleType<weighted> newNodeMapping;
  bool extendable = !weighted;
  do {
    /* matched all nodes with g0node? */
    if(!vm.getNextMatch(mappingIndices, newNodeMapping)) {
      break;
    }
    // nodes that become unmappable after the extension with the new mapping
    std::vector<MappingIndex> newUnmappableNodes;
    // match g0node with a compatible node
    VertexMapping<weighted> updatedMap = mapping.addMatchUpdateVM(vm, newNodeMapping,
                                                                  newUnmappableNodes);

    if(mapping.isExtendable(updatedMap)) {
      extended |= backtrack_impl(mapping, updatedMap);
      mapping.removeMapping(newUnmappableNodes);
      if(!weighted && extended) {
        // in the unweighted case, extendability of the mapping can only change after successful result extension
        extendable = mapping.isExtendable(vm);
      }
    } else {
      mapping.removeMapping(newUnmappableNodes);
    }
    if(weighted) {
      // In the weighted case, each extension can potentially influence the upper bounds
      mapping.recountDisconnectedComponents(vm);
      extendable = mapping.isExtendable(vm);
    }
  } while(extendable);
  // Backtrack leaving g0node unmapped
  mapping.setNodeUnmappable(g0node, vm);
  if(mapping.isExtendable(vm)) {
    extended |= backtrack_impl(mapping, vm);
  }
  mapping.setNodeMappable(g0node);
  mapping.restoreConnectedComponents(g0node);
  return extended;
}


bool backtrack(
    SubgraphMapping<false>& mapping,
    VertexMapping<false>& vm)
{
  return backtrack_impl(mapping, vm);
}


bool backtrack(
    SubgraphMapping<true>& mapping,
    VertexMapping<true>& vm)
{
  return backtrack_impl(mapping, vm);
}


} // namespace RIMACS
