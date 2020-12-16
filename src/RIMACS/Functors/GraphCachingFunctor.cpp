
#include <Eigen/SparseCore>

#include <RIMACS/Forward.hpp>
#include <RIMACS/Functors/LineGraph.hpp>

#include <RIMACS/Functors/GraphCachingFunctor.hpp>


namespace RIMACS {

// Experimental evaluation for optimal node order suggestion as performed by
// Schmidt et.al.: Disconnected Maximum Common Substructures Under Constraints
// https://dx.doi.org/10.1021/acs.jcim.0c00741
static constexpr double NODE_ORDER_RELATIVE_POSITION = 1.0;
static constexpr double EDGE_ORDER_RELATIVE_POSITION = 0.75;

// Con- Destructors have to be defined in cpp because the SparseMatrixCache struct is undefined in the
// hpp. But the unique_ptr con- destructors require, that the struct is defined.
CachingAdjacencyFunctor::CachingAdjacencyFunctor(
    const CachingRequiredMinimalAdjacencyFunctor& realInstance)
  : m_instance(realInstance)
{}

CachingAdjacencyFunctor::~CachingAdjacencyFunctor() = default;
CachingAdjacencyFunctor::CachingAdjacencyFunctor(
    CachingAdjacencyFunctor&& other,
    const CachingRequiredMinimalAdjacencyFunctor& newBaseRef) noexcept
  : m_instance(newBaseRef)
  , m_cache(std::move(other.m_cache))
  , m_nodeOrder(std::move(other.m_nodeOrder))
{}

struct SparseMatrixCache {
  SparseMatrixCache(
      ComponentIndex nofNodes)
    : m_matrix(nofNodes, nofNodes)
  {}

  Eigen::SparseMatrix<MappingIndex> m_matrix;
};


void CachingAdjacencyFunctor::doInitialize()
{
  m_cache = std::make_unique<SparseMatrixCache>(m_instance.getNofNodes());
  auto& mat = m_cache->m_matrix;
  std::vector<Eigen::Triplet<unsigned>> edges;
  edges.reserve(2 * m_instance.getNofEdges());
  for(MappingIndex i = 0, last = m_instance.getNofEdges(); i < last; ++i) {
    MappingIndex from = m_instance.edgeGetFromId(i);
    MappingIndex to = m_instance.edgeGetToId(i);
    edges.emplace_back(from, to, i + 1);
    edges.emplace_back(to, from, i + 1);
  }
  std::sort(edges.begin(), edges.end(),
            [](const Eigen::Triplet<unsigned>& t1, const Eigen::Triplet<unsigned>& t2) {
              if(t1.col() != t2.col()) {
                return t1.col() < t2.col();
              }
              return t1.row() < t2.row();
            });
  mat.setFromTriplets(edges.begin(), edges.end());
  mat.makeCompressed();

  if(dynamic_cast<const BasicLineGraph*>(&m_instance)) {
    m_nodeOrder = this->getNodeOrderSuggestionVector(EDGE_ORDER_RELATIVE_POSITION);
  } else {
    m_nodeOrder = this->getNodeOrderSuggestionVector(NODE_ORDER_RELATIVE_POSITION);
  }
}


bool CachingAdjacencyFunctor::adjacent(
    MappingIndex node1,
    MappingIndex node2) const
{
  const auto& mat = m_cache->m_matrix;
  // perform a linear search on the neighbour indices.
  // On such low degree graphs (usually <= 4), linear searches outperform binary searches
  auto begin = mat.innerIndexPtr() + mat.outerIndexPtr()[node1];
  auto end = mat.innerIndexPtr() + mat.outerIndexPtr()[node1 + 1];
  return std::find(begin, end, node2) != end;
}

MappingIndex CachingAdjacencyFunctor::getEdgeId(
    MappingIndex nodeFromId,
    MappingIndex nodeToId) const
{
  static_assert(NO_EDGE == std::numeric_limits<MappingIndex>::max(),
                "No edge indicator must be the -1 value");
  MappingIndex edgeId = m_cache->m_matrix.coeff(nodeFromId, nodeToId);
  return edgeId - 1;
}


MappingIndex CachingAdjacencyFunctor::edgeGetFromId(
    MappingIndex edgeId) const
{
  return m_instance.edgeGetFromId(edgeId);
}


MappingIndex CachingAdjacencyFunctor::edgeGetToId(
    MappingIndex edgeId) const
{
  return m_instance.edgeGetToId(edgeId);
}


std::vector<MappingIndex> CachingAdjacencyFunctor::nodeGetEdges(
    MappingIndex nodeIdx) const
{
  const Eigen::SparseMatrix<unsigned>& mat = m_cache->m_matrix;
  Eigen::SparseMatrix<unsigned>::InnerIterator it(mat, nodeIdx);
  std::vector<MappingIndex> result;
  result.reserve(mat.outerIndexPtr()[nodeIdx + 1] - mat.outerIndexPtr()[nodeIdx]);
  for(; it; ++it) {
    result.push_back(it.value() - 1);
  }
  return result;
}


std::vector<MappingIndex> CachingAdjacencyFunctor::nodeGetNeighbours(
    MappingIndex nodeIdx) const
{
  const Eigen::SparseMatrix<unsigned>& mat = m_cache->m_matrix;
  Eigen::SparseMatrix<unsigned>::InnerIterator it(mat, nodeIdx);
  std::vector<MappingIndex> result;
  result.reserve(mat.outerIndexPtr()[nodeIdx + 1] - mat.outerIndexPtr()[nodeIdx]);
  for(; it; ++it) {
    result.push_back(it.row());
  }
  return result;
}


MappingIndex CachingAdjacencyFunctor::getNofNodes() const
{
  return m_cache->m_matrix.rows();
}


MappingIndex CachingAdjacencyFunctor::nodeGetNofEdges(
    MappingIndex nodeIdx) const
{
  return m_cache->m_matrix.outerIndexPtr()[nodeIdx + 1]
         - m_cache->m_matrix.outerIndexPtr()[nodeIdx];
}


MappingIndex CachingAdjacencyFunctor::getNofEdges() const
{
  return m_instance.getNofEdges();
}

} // namespace RIMACS
