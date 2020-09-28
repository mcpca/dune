//***************************************************************************
// Copyright 2020 KTH Royal Institute of Technology                         *
//***************************************************************************
// Author: Miguel Aguiar (aguiar@kth.se)                                    *
//***************************************************************************

#pragma once

#include <algorithm>
#include <memory>
#include <stdexcept>

#include <DUNE/Tasks/Task.hpp>
#include <DUNE/Parsers/HDF5Reader.hpp>
#include <DUNE/Math/Grid.hpp>

namespace KTH
{
  //! Insert short task description here.
  //!
  //! Insert explanation on task behaviour here.
  //! @author Miguel Aguiar
  namespace AlgaeBloom
  {
    namespace Simulators
    {
      namespace Density
      {
        inline double
        interpolateLinear2d(double const* values, double const* delta)
        {
          return (1 - delta[0]) * (1 - delta[1]) * values[0]
                 + delta[0] * (1 - delta[1]) * values[1]
                 + delta[0] * delta[1] * values[2]
                 + (1 - delta[0]) * delta[1] * values[3];
        }

        enum Dimensions
        {
          DIM_TIME = 0,
          DIM_LON,
          DIM_LAT
        };

        class DensityField
        {
        public:
          DensityField(std::string const& path, std::string const& dataset,
                       std::size_t time_slice, double default_val,
                       DUNE::Tasks::Task* owner)
          : m_default(default_val), m_owner(owner)
          {
            DUNE::Parsers::HDF5Reader f(path);

            auto data = f.getDataset<double>(dataset);
            m_dims = std::move(data.dimensions);
            auto field_values = std::move(data.data);

            std::size_t const space_dims = m_dims[DIM_LAT] * m_dims[DIM_LON];

            if ((time_slice + 1u) * (space_dims) > field_values.size())
              throw std::runtime_error("Bad dataset dimension");

            m_data = std::vector<double>{
              std::begin(field_values) + time_slice * space_dims,
              std::begin(field_values) + (time_slice + 1) * space_dims
            };

            m_lon = f.getDataset<double>("lon").data;
            m_lat = f.getDataset<double>("lat").data;

            std::vector<double> min = { m_lon[0], m_lat[0] };
            std::vector<double> max = { m_lon.back(), m_lat.back() };
            std::vector<std::size_t> dims = { m_dims[DIM_LON], m_dims[DIM_LAT] };

            m_grid = std::make_unique<DUNE::Math::Grid<2>>(min, max, dims);
          }

          float
          evaluate(double lat, double lon)
          try
          {
            auto neighbor = m_grid->getCorner({ lon, lat });

            auto corner = m_grid->getCoordinates(neighbor);

            double delta[2] = { (lon - corner[0]) / m_grid->getSpacing(0),
                                (lat - corner[1]) / m_grid->getSpacing(1) };

            double vals[4];

            vals[0] = m_data[m_grid->getOffset(neighbor)];

            neighbor[0] += 1;
            vals[1] = m_data[m_grid->getOffset(neighbor)];

            neighbor[1] += 1;
            vals[2] = m_data[m_grid->getOffset(neighbor)];

            neighbor[0] -= 1;
            vals[3] = m_data[m_grid->getOffset(neighbor)];

            m_owner->spew("%.2f, %.2f, %.2f, %.2f", vals[0], vals[1], vals[2], vals[3]);

            return interpolateLinear2d(vals, delta);
          }
          catch (std::exception const& e)
          {
            m_owner->spew("Error getting field value: %s", e.what());
            return m_default;
          }

        private:
          std::vector<double> m_data;
          std::vector<std::size_t> m_dims;
          std::vector<double> m_lat;
          std::vector<double> m_lon;
          std::unique_ptr<DUNE::Math::Grid<2>> m_grid;

          double m_default;
          DUNE::Tasks::Task* m_owner;
        };
      } // namespace Density
    }   // namespace Simulators
  }     // namespace AlgaeBloom
} // namespace KTH
