//***************************************************************************
// Copyright 2020 KTH Royal Institute of Technology                         *
//***************************************************************************
// Author: Miguel Aguiar (aguiar@kth.se)                                    *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>

// Local headers.
#include "Density.hpp"

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
        using DUNE_NAMESPACES;

        struct Arguments
        {
          std::string path;
          std::string dataset_name;
          unsigned time_slice;
          float default_value;
          unsigned n_ma_pts;
        };

        struct Task : public DUNE::Tasks::Task
        {
          Arguments m_args;
          DensityField* m_field;
          std::unique_ptr<Math::MovingAverage<float>> m_data;

          //! Constructor.
          //! @param[in] name task name.
          //! @param[in] ctx context.
          Task(const std::string& name, Tasks::Context& ctx)
          : DUNE::Tasks::Task(name, ctx), m_field(NULL)
          {
            bind<IMC::EstimatedState>(this);

            param("Data File Path", m_args.path);
            param("Dataset Name", m_args.dataset_name);
            param("Time Slice", m_args.time_slice).defaultValue("0");
            param("Default Value", m_args.default_value).defaultValue("0.0");
            param("Moving Average Size", m_args.n_ma_pts).defaultValue("60");
          }

          void
          onResourceAcquisition(void)
          {
            m_data
            = std::make_unique<Math::MovingAverage<float>>(m_args.n_ma_pts);

            try
            {
              m_field
              = new DensityField(m_args.path, m_args.dataset_name,
                                 m_args.time_slice, m_args.default_value, this);
            }
            catch (std::exception const& e)
            {
              err("Could not create DensityField: %s", e.what());
              m_field = NULL;
            }
          }

          void
          consume(IMC::EstimatedState const* es)
          {
            // Don't consume other vehicles' states.
            if (es->getSource() != getSystemId())
            {
              spew("Ignoring ES msg from %u (I\'m %u)", es->getSource(),
                   getSystemId());
              return;
            }

            // Get vehicle pos
            double lat = es->lat;
            double lon = es->lon;
            WGS84::displace(es->x, es->y, &lat, &lon);

            lat = Angles::degrees(lat);
            lon = Angles::degrees(lon);

            float val
            = m_field ? m_field->evaluate(lat, lon) : m_args.default_value;

            val = m_data->update(val);

            if (m_data->sampleSize() != m_args.n_ma_pts)
              return;

            trace("Chl: %.4f", val);

            IMC::Chlorophyll msg;
            msg.value = val;
            dispatch(msg);

            m_data->clear();
          }

          //! Release resources.
          void
          onResourceRelease(void)
          {
            Memory::clear(m_field);
          }

          //! Main loop.
          void
          onMain(void)
          {
            while (!stopping())
            {
              waitForMessages(1.0);
            }
          }
        };
      } // namespace Density
    }   // namespace Simulators
  }     // namespace AlgaeBloom
} // namespace KTH

DUNE_TASK
