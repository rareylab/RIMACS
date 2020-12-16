
#include <numeric>

#include <gtest/gtest.h>

#include <RIMACS/Algorithm/CompareUtils.hpp>
#include <RIMACS/MCSResult.hpp>

class MCSResultTester
    : public testing::TestWithParam<std::tuple<std::vector<RIMACS::MappingPair>,
                                               std::vector<RIMACS::MappingIndex>,
                                               std::vector<RIMACS::MappingIndex>>>
{};

TEST_P(MCSResultTester, resultIsSorted)
{
  std::vector<RIMACS::MappingPair> mapping;
  std::vector<RIMACS::MappingIndex> rawComponents;
  std::vector<RIMACS::MappingIndex> components;
  std::tie(mapping, rawComponents, components) = GetParam();

  ASSERT_FALSE(components.empty());
  EXPECT_EQ(components.size(), rawComponents.size());
  EXPECT_EQ(components.back(), mapping.size());
  for(auto it = components.begin() + 1, last = components.end(); it != last; ++it) {
    EXPECT_LE(it[-1], *it);
  }
  components.insert(components.begin(), 0);

  RIMACS::MCSResult res(mapping, rawComponents, 1.0 , 0);
  const auto& mappingResult = res.getMappings();
  for(auto it = components.begin() + 1, last = components.end(); it != last; ++it) {
    RIMACS::MappingIndex beginIdx = it[-1];
    RIMACS::MappingIndex endIdx = *it;
    for(++beginIdx; beginIdx < endIdx; ++beginIdx) {
      EXPECT_LT(mappingResult.at(beginIdx - 1), mappingResult.at(beginIdx));
    }
  }
}

using M = std::vector<RIMACS::MappingPair>;
using P = RIMACS::MappingPair;
using C = std::vector<RIMACS::MappingIndex>;


INSTANTIATE_TEST_SUITE_P(EmptyResult, MCSResultTester,
    testing::Values(std::make_tuple(M{}, C{0}, C{0})));

INSTANTIATE_TEST_SUITE_P(LinearMapping, MCSResultTester,
    testing::Values(std::make_tuple(M{P{0,0}, P{1,1}, P{2,2}, P{3,3}, P{4,4}}, C{5}, C{5}),
                    std::make_tuple(M{P{4,4}, P{3,3}, P{2,2}, P{1,1}, P{0,0}}, C{5}, C{5}),
                    std::make_tuple(M{P{4,4}, P{1,1}, P{0,0}, P{3,3}, P{2,2}}, C{5}, C{5})));

INSTANTIATE_TEST_SUITE_P(DisconnectedMapping, MCSResultTester,
    testing::Values(
                    std::make_tuple(M{P{0,0}, P{1,1}, P{2,2}, P{5,4}, P{4,5}}, C{{3, 2}}, C{{3, 5}}),
                    std::make_tuple(M{P{0,5}, P{2,4}, P{1,0}, P{3,3}, P{4,2}}, C{{2, 3}}, C{{2, 5}})));

class EquivalenceResultTester
    : public testing::TestWithParam<std::tuple<std::vector<RIMACS::MappingPair>,
                                               std::vector<RIMACS::MappingIndex>,
                                               RIMACS::MappingIndex,
                                               std::vector<RIMACS::MappingIndex>,
                                               std::vector<RIMACS::MappingIndex>,
                                               std::vector<RIMACS::MappingIndex>>> {};

namespace {

void ensure_equivalence_class_is_valid(
    const std::vector<RIMACS::MappingPair>& mapping,
    std::vector<RIMACS::MappingIndex>& equivalence_class)
{
  if(!equivalence_class.empty()) {
    return;
  }
  RIMACS::MappingIndex maxElem =
      std::max(std::max_element(mapping.begin(), mapping.end(), RIMACS::mapping_from_compare())->m_from,
               std::max_element(mapping.begin(), mapping.end(), RIMACS::mapping_to_compare())->m_to);
  equivalence_class.resize(maxElem + 1);
  std::iota(equivalence_class.begin(), equivalence_class.end(), 0);
}

} // namespace

TEST_P(EquivalenceResultTester, resultIsSorted)
{
  std::vector<RIMACS::MappingPair> mapping;
  std::vector<RIMACS::MappingIndex> rawComponents;
  RIMACS::MappingIndex initialMappingSize;
  std::vector<RIMACS::MappingIndex> queryEquivalenceClasses;
  std::vector<RIMACS::MappingIndex> targetEquivalenceClasses;
  std::vector<RIMACS::MappingIndex> initialComponentNodes;

  std::tie(mapping, rawComponents, initialMappingSize, queryEquivalenceClasses,
           targetEquivalenceClasses, initialComponentNodes) = GetParam();

  ensure_equivalence_class_is_valid(mapping, queryEquivalenceClasses);
  ensure_equivalence_class_is_valid(mapping, targetEquivalenceClasses);

  RIMACS::Compare::equivalence_class_compare<std::less<RIMACS::MappingIndex>,
                                             std::equal_to<RIMACS::MappingIndex>> compare(
      &queryEquivalenceClasses, &targetEquivalenceClasses);
  auto notCompare = [&compare](const RIMACS::MappingPair& lhs, const RIMACS::MappingPair& rhs) {
    return !compare(lhs, rhs);
  };

  RIMACS::MCSResult res(mapping, rawComponents, 1.0 , initialMappingSize,
                     queryEquivalenceClasses.empty() ? nullptr : &queryEquivalenceClasses,
                     targetEquivalenceClasses.empty() ? nullptr : &targetEquivalenceClasses);
  std::vector<RIMACS::MappingIndex> components = res.getComponentSizes();
  components.insert(components.begin(), 0);

  const auto& mappingResult = res.getMappings();
  for(auto it = components.begin() + 1, last = components.end(); it != last; ++it) {
    RIMACS::MappingIndex beginIdx = it[-1];
    beginIdx += it == components.begin() + 1 ? initialMappingSize : 0;
    RIMACS::MappingIndex endIdx = *it;
    for(++beginIdx; beginIdx < endIdx; ++beginIdx) {
      const RIMACS::MappingPair& prev = mappingResult.at(beginIdx - 1);
      const RIMACS::MappingPair& current = mappingResult.at(beginIdx);
      std::stringstream compareMessage;
      compareMessage << "Compare [" << prev.m_from << ", " << prev.m_to << "] and [" << current.m_from << ", " << current.m_to << "]";
      EXPECT_PRED2(notCompare, mappingResult.at(beginIdx), mappingResult.at(beginIdx - 1));
      if(!compare(mappingResult.at(beginIdx - 1), mappingResult.at(beginIdx))) {
        EXPECT_TRUE(!queryEquivalenceClasses.empty() || !targetEquivalenceClasses.empty());
      }
    }
  }
  for(size_t i = 0, last = components.size() - 1; i < last; ++i) {
    EXPECT_EQ(mappingResult.at(components.at(i)).m_from, initialComponentNodes.at(i));
  }
}

namespace {


using M = std::vector<RIMACS::MappingPair>;
using P = RIMACS::MappingPair;
using C = std::vector<RIMACS::MappingIndex>;
using EQ = std::vector<RIMACS::MappingIndex>;

INSTANTIATE_TEST_SUITE_P(LinearMapping, EquivalenceResultTester,
    testing::Values(std::make_tuple(M{P{0,0}, P{1,1}, P{2,2}, P{3,3}, P{4,4}},
                                    C({5}), 0u, EQ{}, EQ{}, C({0})),
                    std::make_tuple(M{P{4,4},P{1,1}, P{0,0}, P{3,3}, P{2,2}},
                                    C({5}), 1u, EQ{}, EQ{}, C({4})),
                    std::make_tuple(M{P{0,0}, P{1,1}, P{2,2}, P{3,3}, P{4,4}},
                                    C{{2, 3}}, 0u, EQ{}, EQ{}, C{{0, 3}}),
                    std::make_tuple(M{P{4,4},P{1,1}, P{0,0}, P{3,3}, P{2,2}},
                                    C{{2, 3}}, 1u, EQ{}, EQ{}, C{{4, 2}})));

INSTANTIATE_TEST_SUITE_P(DisconnectedLinear, EquivalenceResultTester,
    testing::Values(std::make_tuple(M{{4,4}, {1,9}, {0,8}, {5,2}, {7,3}, {8, 0}, {3, 7}, {6, 1}},
                                    C{{3, 2, 3}}, 1u,
                                    EQ{}, EQ{}, C{{4, 3, 5}}),
                    std::make_tuple(M{{4,4}, {1,9}, {0,8}, {5,2}, {7,3}, {8, 0}, {3, 7}, {6, 1}},
                                    C{{3, 2, 3}}, 1u,
                                    EQ(size_t{10}, RIMACS::MappingIndex{0}),
                                    EQ(size_t{10}, RIMACS::MappingIndex{0}), C{{4, 8, 5}}),
                    std::make_tuple(M{{4,4}, {1,9}, {0,8}, {5,2}, {7,3}, {8, 0}, {3, 7}, {6, 1}},
                                    C{{3, 2, 3}}, 0u,
                                    EQ(size_t{10}, RIMACS::MappingIndex{0}),
                                    EQ(), C{{8, 5, 4}}),
                    std::make_tuple(M{{4,4}, {1,9}, {0,8}, {5,2}, {7,3}, {8, 0}, {3, 7}, {6, 1}},
                                    C{{3, 2, 3}}, 0u,
                                    EQ(), EQ(size_t{10}, RIMACS::MappingIndex{0}),
                                    C{{0, 3, 5}})));
} // namespace
