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

        struct AbsolutePosition
        {
          double lat;
          double lon;
        };

        struct RelativePosition
        {
          double x;
          double y;
        };

        struct ControllerState
        {
          // Position of virtual vehicle
          AbsolutePosition vp;

          // Track bearing (rad)
          float direction;

          // Speed in mps
          float speed;

          // How many waypoints since last update
          int n_wpts;
        };

        struct ControllerParameters
        {
          // Angle with track bearing (rad)
          float angle;
          // Horizontal distance travelled by the vehicle.
          float distance;
        };

        struct Sample
        {
          double lat;
          double lon;
          float chl_value;
        };

        struct Arguments
        {
          unsigned n_ma_pts;
          float initial_speed;
          float initial_heading;
          float target_value;
          float along_track_gain;
          float cross_track_gain;
        };

        struct Task : public DUNE::Tasks::Task
        {
          Arguments m_args;
          double m_last_data_timestamp;
          std::unique_ptr<Math::MovingAverage<float>> m_data;

          ControllerState m_state;
          ControllerParameters m_params;

          RelativePosition m_last_pos;

          float m_chl_val = 0.0;
          float m_chl_grad = 0.0;

          bool m_inited = false;

          IMC::Reference m_ref;

          //! Constructor.
          //! @param[in] name task name.
          //! @param[in] ctx context.
          Task(const std::string& name, Tasks::Context& ctx)
          : DUNE::Tasks::Task(name, ctx)
          {
            bind<IMC::Chlorophyll>(this);
            bind<IMC::EstimatedState>(this);
            bind<IMC::FollowRefState>(this);

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

            param("Zigzag Angle", m_params.angle)
            .defaultValue("45")
            .units(Units::Degree);

            param("Horizontal Distance", m_params.distance).defaultValue("100");

            m_last_data_timestamp = Clock::get();

            m_ref.flags = IMC::Reference::FLAG_SPEED
                          | IMC::Reference::FLAG_LOCATION
                          | IMC::Reference::FLAG_DIRECT;
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
          consume(IMC::EstimatedState const* msg)
          {
            if (msg->getSource() != getSystemId())
              return;

            double lat;
            double lon;

            Coordinates::toWGS84(*msg, lat, lon);

            if (!m_inited)
            {
              m_state.vp.lat = lat;
              m_state.vp.lon = lon;
            }

            WGS84::displacement(m_state.vp.lat, m_state.vp.lon, msg->height,
                                lat, lon, msg->height, &m_last_pos.x,
                                &m_last_pos.y);

            if (!m_inited)
            {
              updateRef();
              m_inited = true;
            }
          }

          void
          updateVirtualPos(void)
          {
            static constexpr auto origin = RelativePosition{ 0.0, 0.0 };
            double x;

            Coordinates::getTrackPosition(origin, m_state.direction, m_last_pos,
                                          &x);

            spew("along track pos: %.4f", x);

            WGS84::displace(x * std::cos(m_state.direction),
                            x * std::sin(m_state.direction), &m_state.vp.lat,
                            &m_state.vp.lon);

            debug("Lat: %.6f Lon: %.6f", m_state.vp.lat, m_state.vp.lon);
          }

          void
          consume(IMC::FollowRefState const* msg)
          {
            if (!(msg->proximity & IMC::FollowRefState::PROX_XY_NEAR))
              return;

            inf("near waypoint");

            m_state.n_wpts++;

            updateVirtualPos();
            updateRef();
          }

          void
          updateRef(void)
          {
            int const sign = 2 * (m_state.n_wpts % 2) - 1;

            double const along_track_displacement
            = m_params.distance / std::tan(m_params.angle);
            RelativePosition next_wpt
            = { along_track_displacement, sign * m_params.distance };

            debug("next waypont is at %.6f, %.6f", next_wpt.x, next_wpt.y);

            double bearing, range;
            Coordinates::toPolar(next_wpt, &bearing, &range);

            debug("bearing: %.2f, range: %.2f", Angles::degrees(bearing), range);

            m_ref.lat = m_state.vp.lat;
            m_ref.lon = m_state.vp.lon;

            WGS84::displace(range * std::cos(bearing),
                            range * std::sin(bearing), &m_ref.lat, &m_ref.lon);
          }

          void
          onUpdateParameters(void)
          {
            if (paramChanged(m_args.initial_heading))
            {
              m_args.initial_heading = Angles::radians(m_args.initial_heading);
              spew("initial heading: %.6f", m_args.initial_heading);
            }

            if (paramChanged(m_params.angle))
            {
              m_params.angle = Angles::radians(m_params.angle);
              spew("angle: %.6f", m_params.angle);
            }
          }

          void
          initializeReferenceFollower(void)
          {
            IMC::FollowReference msg;

            msg.control_src = getSystemId();
            msg.control_ent = getEntityId();
            msg.timeout = 1;
            msg.loiter_radius = 10;
            msg.altitude_interval = 2;

            dispatch(msg);

            IMC::DesiredSpeed speed_ref;
            speed_ref.value = m_args.initial_speed;
            speed_ref.speed_units = IMC::SpeedUnits::SUNITS_METERS_PS;
            m_ref.speed.set(speed_ref);
          }

          void
          update(void)
          {
          }

          //! Main loop.
          void
          onMain(void)
          {
            m_last_pos.x = 0.0;
            m_last_pos.y = 0.0;

            // Initialize controller state
            m_state.n_wpts = 0;
            m_state.direction = m_args.initial_heading;

            // Wait for boot
            Delay::wait(10.0);

            initializeReferenceFollower();

            while (!stopping())
            {
              waitForMessages(1.0);

              if (!m_inited)
                continue;

              dispatch(m_ref);
            }
          }
        };
      } // namespace ChlGuidance
    }   // namespace Control
  }     // namespace AlgaeBloom
} // namespace KTH

DUNE_TASK
