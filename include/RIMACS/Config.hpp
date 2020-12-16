
#pragma once

#include <algorithm>

#include <RIMACS/ASSERT.hpp>

#include <RIMACS/Forward.hpp>

namespace RIMACS {

namespace ConfigTypes {
  enum class ResultType {
    AllMaximum,
    Maximum
  };

  enum class MCSType {
    Connected,
    Disconnected
  };

  enum Defaults : Index {
    ComponentSize = 3,
    NofComponents = 3,
    RecursionLimit = 100000000,
    MaxNofEquivalentMappings = 1
  };
} // namespace ConfigTypes

const Config& getDefaultConnectedConfig();

const Config& getDefaultDisconnectedConfig();

inline const Config& getDefaultConfig()
{
  return getDefaultConnectedConfig();
}

Config getConnectedConfig(
    ConfigTypes::ResultType,
    ComponentIndex,
    Index);

Config getDisconnectedConfig(
    ConfigTypes::ResultType,
    MappingIndex,
    MappingIndex,
    ComponentIndex,
    Index);

class Config {
public:

  using ResultType = ConfigTypes::ResultType;
  using MCSType = ConfigTypes::MCSType;
  using Defaults = ConfigTypes::Defaults;

  Config()
      : Config(getDefaultConfig())
  {}

  bool isDisconnected() const
  {
    return m_components > 1;
  }

  Config::ResultType resultType() const
  {
    return m_result;
  }

  Config::MCSType type() const
  {
    return m_type;
  }

  MappingIndex nofComponents() const
  {
    return m_components;
  }

  MappingIndex componentSize() const
  {
    return m_componentSize;
  }

  Index recursionLimit() const
  {
    return m_recursionLimit;
  }

  ComponentIndex equivalentResults() const
  {
    return m_maxNofEquivalentResults;
  }

  MappingIndex getMinimumResultSize() const
  {
    return m_minimumResultSize;
  }

  const double& getMinimumResultWeight() const
  {
    return m_minimumResultWeight;
  }

  void setMinimumSize(
      MappingIndex minimumSize)
  {
    m_minimumResultSize = minimumSize;
  }

  void setMinimumWeight(
      double minimumWeight)
  {
    m_minimumResultWeight = minimumWeight;
  }


private:

  friend const Config& getDefaultConnectedConfig();
  friend const Config& getDefaultDisconnectedConfig();
  friend Config getConnectedConfig(Config::ResultType, ComponentIndex, Index);
  friend Config getDisconnectedConfig(Config::ResultType, MappingIndex, MappingIndex, ComponentIndex, Index);

  Config(
      Config::ResultType res,
      ComponentIndex maxEquivalentRes,
      Index recursionLimit)
  : m_result(res)
  , m_type(Config::MCSType::Connected)
  , m_components(1)
  , m_componentSize(std::numeric_limits<MappingIndex>::max())
  , m_recursionLimit(recursionLimit)
  , m_maxNofEquivalentResults(maxEquivalentRes)
  {}

  Config(
      Config::ResultType res,
      MappingIndex components,
      MappingIndex componentSize,
      ComponentIndex maxEquivalentRes,
      Index recursionLimit)
  : m_result(res)
  , m_type(Config::MCSType::Disconnected)
  , m_components(components)
  , m_componentSize(componentSize)
  , m_recursionLimit(recursionLimit)
  , m_maxNofEquivalentResults(maxEquivalentRes)
  {}

  Config::ResultType m_result;
  Config::MCSType m_type;
  MappingIndex m_components;
  MappingIndex m_componentSize;
  Index m_recursionLimit;
  ComponentIndex m_maxNofEquivalentResults;
  double m_minimumResultWeight = 0.0;
  MappingIndex m_minimumResultSize = 0;
};

inline const Config& getDefaultConnectedConfig()
{
  static Config config{Config::ResultType::Maximum, Config::Defaults::MaxNofEquivalentMappings, Config::Defaults::RecursionLimit};
  return config;
}

inline const Config& getDefaultDisconnectedConfig()
{
  static Config config{Config::ResultType::Maximum, Config::Defaults::NofComponents, Config::Defaults::ComponentSize,
                       Config::Defaults::MaxNofEquivalentMappings, Config::Defaults::RecursionLimit};
  return config;
}

inline Config getConnectedConfig(
    Config::ResultType rt = getDefaultConnectedConfig().resultType(),
    ComponentIndex maxEquivalentResults = Config::Defaults::MaxNofEquivalentMappings,
    Index recursionLimit = Config::Defaults::RecursionLimit)
{
  Config config{rt, maxEquivalentResults, recursionLimit};
  return config;
}

inline Config getDisconnectedConfig(
    Config::ResultType rt = getDefaultDisconnectedConfig().resultType(),
    MappingIndex components = Config::Defaults::NofComponents,
    MappingIndex componentSize = Config::Defaults::ComponentSize,
    ComponentIndex maxEquivalentResults = Config::Defaults::MaxNofEquivalentMappings,
    Index recursionLimit = Config::Defaults::RecursionLimit)
{
  RIMACS_ASSERT(components > 1 && "Disconnected config requires at least two components");
  Config config{rt, components, std::max(MappingIndex{1}, componentSize), maxEquivalentResults,
                recursionLimit};
  return config;
}

inline Config getConfig(
    Config::ResultType rt = getDefaultDisconnectedConfig().resultType(),
    MappingIndex components = Index{1},
    MappingIndex componentSize = Config::Defaults::ComponentSize,
    ComponentIndex maxEquivalentResults = Config::Defaults::MaxNofEquivalentMappings,
    Index recursionLimit = Config::Defaults::RecursionLimit)
{
  if(components <= 1 || componentSize == std::numeric_limits<MappingIndex>::max()) {
    return getConnectedConfig(rt, maxEquivalentResults, recursionLimit);
  }
  return getDisconnectedConfig(rt, components, componentSize, maxEquivalentResults, recursionLimit);
}

} // namespace RIMACS
