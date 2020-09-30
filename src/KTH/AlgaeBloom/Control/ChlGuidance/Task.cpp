//***************************************************************************
// Copyright 2020 KTH Royal Institute of Technology                         *
//***************************************************************************
// Author: Miguel Aguiar (aguiar@kth.se)                                    *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>

namespace KTH
{
  namespace AlgaeBloom
  {
    namespace Control
    {
      namespace ChlGuidance
      {
        using DUNE_NAMESPACES;

        struct Arguments
        {
          float initial_speed;
          float initial_heading;
          float target_value;
          float along_track_gain;
          float cross_track_gain;
          unsigned n_ma_pts;
        };

        struct Task : public DUNE::Tasks::Task
        {
          Arguments m_args;
          double m_last_data_timestamp;
          IMC::DesiredHeading m_heading;
          IMC::DesiredSpeed m_speed;
          std::unique_ptr<Math::MovingAverage<float>> m_data;

          float m_chl_val = 0.0;
          float m_chl_grad = 0.0;

          //! Constructor.
          //! @param[in] name task name.
          //! @param[in] ctx context.
          Task(const std::string& name, Tasks::Context& ctx)
          : DUNE::Tasks::Task(name, ctx)
          {
            bind<IMC::Chlorophyll>(this);

            param("Initial Speed", m_args.initial_speed)
            .defaultValue("1")
            .units(Units::MeterPerSecond);

            param("Initial Heading", m_args.initial_heading)
            .defaultValue("-45")
            .units(Units::Degree);

            param("Target Chl Value", m_args.target_value).defaultValue("0.5");

            param("Along Track Gain", m_args.along_track_gain)
            .defaultValue("1.0");

            param("Cross Track Gain", m_args.cross_track_gain)
            .defaultValue("1.0");

            param("Moving Average Size", m_args.n_ma_pts).defaultValue("60");

            m_last_data_timestamp = Clock::get();
          }

          void
          onResourceAcquisition(void)
          {
            m_data
            = std::make_unique<Math::MovingAverage<float>>(m_args.n_ma_pts);
          };

          void
          consume(IMC::Chlorophyll const* msg)
          {
            float new_val = m_data->update(msg->value);

            if (m_data->sampleSize() == m_args.n_ma_pts)
            {
              float old_time = m_last_data_timestamp;
              m_last_data_timestamp = Clock::get();

              m_chl_grad
              = (new_val - m_chl_val) / (m_last_data_timestamp - old_time);

              spew("%.8f %.8f %.12f", m_chl_val, new_val, m_chl_grad);

              m_chl_val = new_val;

              m_data->clear();
            }
          }

          void
          onUpdateParameters(void)
          {
            if (paramChanged(m_args.initial_heading))
            {
              m_args.initial_heading = Angles::radians(m_args.initial_heading);
              spew("initial heading: %.6f", m_args.initial_heading);
            }
          }

          void
          initializeLoops(void)
          {
            IMC::ControlLoops msg;

            msg.enable = IMC::ControlLoops::CL_ENABLE;
            msg.mask = IMC::CL_YAW | IMC::CL_SPEED;

            dispatch(msg);
          }

          void
          initializeState(void)
          {
            m_heading.value = m_args.initial_heading;

            m_speed.speed_units = IMC::SpeedUnits::SUNITS_METERS_PS;
            m_speed.value = m_args.initial_speed;
          }

          void
          update(void)
          {
          }

          //! Main loop.
          void
          onMain(void)
          {
            initializeLoops();
            initializeState();

            IMC::DesiredZ msg;
            msg.value = 0.0;
            msg.z_units = IMC::Z_DEPTH;
            dispatch(msg);

            while (!stopping())
            {
              waitForMessages(1.0);

              update();

              // Dispatch heading and speed references
              dispatch(m_heading);
              dispatch(m_speed);
            }
          }
        };
      } // namespace ChlGuidance
    }   // namespace Control
  }     // namespace AlgaeBloom
} // namespace KTH

DUNE_TASK
