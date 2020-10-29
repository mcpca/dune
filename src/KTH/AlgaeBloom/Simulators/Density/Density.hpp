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

        inline void
        normalizeTimeAxis(std::vector<double>* time)
        {
          const double origin = (*time)[0];
          for (double& t : *time)
            t = (24.0 * 3600.0) * (t - origin);
        }

        class DensityField
        {
        public:
          DensityField(std::string const& path, std::string const& dataset,
                       double default_val, DUNE::Tasks::Task* owner)
          : m_default(default_val), m_owner(owner)
          {
            DUNE::Parsers::HDF5Reader f(path);

            auto data = f.getDataset<double>(dataset);
            m_dims = std::move(data.dimensions);
            m_data = std::move(data.data);

            m_time = f.getDataset<double>("time").data;
            m_lon = f.getDataset<double>("lon").data;
            m_lat = f.getDataset<double>("lat").data;

            normalizeTimeAxis(&m_time);

            std::vector<double> min = { m_time[0], m_lon[0], m_lat[0] };
            std::vector<double> max
            = { m_time.back(), m_lon.back(), m_lat.back() };
            std::vector<std::size_t> dims
            = { m_dims[DIM_TIME], m_dims[DIM_LON], m_dims[DIM_LAT] };

            m_grid = std::make_unique<DUNE::Math::Grid<3>>(min, max, dims);
          }

          double
          interpSpatial(std::array<size_t, 3> corner,
                        double const delta[2]) const
          {
            double vals[4];

            vals[0] = m_data[m_grid->getOffset(corner)];

            corner[DIM_LON] += 1;
            vals[1] = m_data[m_grid->getOffset(corner)];

            corner[DIM_LAT] += 1;
            vals[2] = m_data[m_grid->getOffset(corner)];

            corner[DIM_LON] -= 1;
            vals[3] = m_data[m_grid->getOffset(corner)];

            return interpolateLinear2d(vals, delta);
          }

          double
          evaluate(double time, double lat, double lon) const
          try
          {
            m_owner->spew("Point: %.6f %.6f %.6f", time, lon, lat);
            auto corner = m_grid->getCorner({ time, lon, lat });

            auto const corner_coords = m_grid->getCoordinates(corner);

            m_owner->spew("time=%lu [%.6f], lon=%lu [%.6f], lat=%lu [%.6f]",
                          corner[DIM_TIME], corner_coords[DIM_TIME],
                          corner[DIM_LON], corner_coords[DIM_LON],
                          corner[DIM_LAT], corner_coords[DIM_LAT]);

            double const delta[2]
            = { (lon - corner_coords[DIM_LON]) / m_grid->getSpacing(DIM_LON),
                (lat - corner_coords[DIM_LAT]) / m_grid->getSpacing(DIM_LAT) };

            double const left = interpSpatial(corner, delta);

            if (corner[DIM_TIME] == m_time.size() - 1)
              return left;

            double const weight
            = (time - m_time[corner[DIM_TIME]])
              / (m_time[corner[DIM_TIME] + 1] - m_time[corner[DIM_TIME]]);

            corner[DIM_TIME] += 1;

            double const right = interpSpatial(corner, delta);

            m_owner->spew("left=%.6f, right=%.6f, weight=%.6f", left, right,
                          weight);

            return (1.0 - weight) * left + weight * right;
          }
          catch (std::exception const& e)
          {
            m_owner->war("Error getting field value: %s", e.what());
            return m_default;
          }

        private:
          std::vector<double> m_data;
          std::vector<std::size_t> m_dims;
          std::vector<double> m_time;
          std::vector<double> m_lat;
          std::vector<double> m_lon;
          std::unique_ptr<DUNE::Math::Grid<3>> m_grid;

          double m_default;
          DUNE::Tasks::Task* m_owner;
        };
      } // namespace Density
    }   // namespace Simulators
  }     // namespace AlgaeBloom
} // namespace KTH
