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

          // Position of vehicle relative to virtual vehicle
          RelativePosition rp;

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
          // Front seeking gain
          float seeking_gain;
          // Front following gain
          float following_gain;
        };

        struct Sample
        {
          double lat;
          double lon;
          float chl_value;
        };

        struct Arguments
        {
          float initial_speed;
          float initial_heading;
          float target_value;
          float threshold;
        };

        struct Task : public DUNE::Tasks::Task
        {
          Arguments m_args;
          double m_last_data_timestamp;

          std::vector<Sample> m_samples;

          ControllerState m_state;
          ControllerParameters m_params;

          bool m_inited = false;
          bool m_wpt_cleared = true;
          bool m_front_crossed = false;

          IMC::Reference m_ref;
          IMC::Reference m_prev_ref;

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

            param("Front Crossing Threshold", m_args.threshold)
            .defaultValue("0.1");

            param("Following Gain", m_params.following_gain)
            .defaultValue("3.0");

            param("Seeking Gain", m_params.seeking_gain).defaultValue("1.0");

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
          consume(IMC::Chlorophyll const* msg)
          {
            double lat = m_state.vp.lat;
            double lon = m_state.vp.lon;

            WGS84::displace(m_state.rp.x, m_state.rp.y, &lat, &lon);

            m_samples.push_back({ lat, lon, msg->value });
          }

          void
          consume(IMC::EstimatedState const* msg)
          {
            if (msg->getSource() != getSystemId())
              return;

            double lat;
            double lon;

            Coordinates::toWGS84(*msg, lat, lon);

            if (lat < 0.01 || lon < 0.01)
            {
              war("Invalid initial position %.6f, %.6f", lat, lon);
              return;
            }

            if (!m_inited)
            {
              war("init pos = %.6f, %.6f", lat, lon);
              m_state.vp.lat = lat;
              m_state.vp.lon = lon;
            }

            WGS84::displacement(m_state.vp.lat, m_state.vp.lon, msg->height,
                                lat, lon, msg->height, &m_state.rp.x,
                                &m_state.rp.y);

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

            Coordinates::getTrackPosition(origin, m_state.direction, m_state.rp,
                                          &x);

            spew("along track pos: %.4f", x);

            WGS84::displace(x * std::cos(m_state.direction),
                            x * std::sin(m_state.direction), &m_state.vp.lat,
                            &m_state.vp.lon);

            debug("Lat: %.6f Lon: %.6f", m_state.vp.lat, m_state.vp.lon);
          }

          void
          resetVirtualPos(void)
          {
            WGS84::displace(m_state.rp.x, m_state.rp.y, &m_state.vp.lat,
                            &m_state.vp.lon);
          }

          void
          consume(IMC::FollowRefState const* msg)
          {
            if (msg->reference.isNull())
              return;

            // Assumes we are not working near (lat, lon) = (0, 0)
            if (msg->reference->lat != m_prev_ref.lat
                || msg->reference->lon != m_prev_ref.lon)
              m_wpt_cleared = true;

            if (!m_wpt_cleared)
              return;

            if (!(msg->proximity & IMC::FollowRefState::PROX_XY_NEAR))
              return;

            inf("near waypoint");

            m_state.n_wpts++;

            if (m_state.n_wpts >= 2)
            {
              war("Switching direction");
              m_state.n_wpts = 0;

              updateDirection();

              // resetVirtualPos();
            }

            updateVirtualPos();

            m_prev_ref = m_ref;
            updateRef();

            m_wpt_cleared = false;
          }

          void
          updateRef(void)
          {
            int const sign = 2 * (m_state.n_wpts % 2) - 1;

            double const distance = m_front_crossed ? m_params.distance : 0.0;

            double const along_track_displacement
            = m_front_crossed ? distance / std::tan(m_params.angle)
                              : m_params.distance;
            RelativePosition next_wpt
            = { along_track_displacement, sign * distance };

            debug("next waypont is at %.6f, %.6f", next_wpt.x, next_wpt.y);

            double bearing, range;
            Coordinates::toPolar(next_wpt, &bearing, &range);

            m_ref.lat = m_state.vp.lat;
            m_ref.lon = m_state.vp.lon;

            bearing += m_state.direction;

            debug("bearing: %.2f, range: %.2f", Angles::degrees(bearing),
                  range);

            WGS84::displace(range * std::cos(bearing),
                            range * std::sin(bearing), &m_ref.lat, &m_ref.lon);
          }

          static std::vector<double>
          gradientEstimator(const std::vector<Sample>& samples)
          {
            Matrix regressor(samples.size(), 3);
            Matrix outputs(samples.size(), 1);

            for (unsigned row = 0; row < samples.size(); ++row)
            {
              regressor(row, 0) = 1.0;
              regressor(row, 1) = samples[row].lat;
              regressor(row, 2) = samples[row].lon;

              outputs(row) = samples[row].chl_value;

              // std::printf("%.6f %.6f %.6f %.6f\n", regressor(row, 0),
              //             regressor(row, 1), regressor(row, 2),
                          // outputs(row, 0));
            }

            Matrix const regressor_t = transpose(regressor);
            Matrix const parameters
            = inverse(regressor_t * regressor) * regressor_t * outputs;

            std::cout << "Parameters: " << parameters;

            return { sum(outputs) / outputs.size(), parameters(1), parameters(2) };
          }

          void
          updateDirection(void)
          {
            Matrix psi(2, 1);
            float chl_value;

            if (m_front_crossed)
            {
              debug("Estimating gradient with %lu samples", m_samples.size());

              auto const parameters = gradientEstimator(m_samples);
              m_samples.clear();

              debug("Parameters: %.6f %.6f %.6f", parameters[0], parameters[1],
                    parameters[2]);

              psi(0) = parameters[1];
              psi(1) = parameters[2];

              double const norm = psi.norm_2();

              if (norm < 1e-6)
              {
                war("Small gradient signal, using default heading");
                m_state.direction = m_args.initial_heading;
                return;
              }

              psi = psi / norm;

              // Average value of measurements
              chl_value = parameters[0];
            }
            else
            {
              psi(0) = std::cos(m_args.initial_heading);
              psi(1) = std::sin(m_args.initial_heading);

              chl_value
              = std::accumulate(m_samples.begin(), m_samples.end(), 0.0,
                                [](auto const& a, auto const& b) {
                                  return a + b.chl_value;
                                })
                / m_samples.size();

              m_samples.clear();

              if (std::fabs(chl_value - m_args.target_value) < m_args.threshold)
              {
                m_front_crossed = true;
                war("Front has been found!");
              }
            }

            // orthogonal to the gradient direction
            Matrix epsi(2, 1);
            epsi(0) = -psi(1);
            epsi(1) = psi(0);

            double const error = m_args.target_value - chl_value;
            Matrix const u_seek = m_params.seeking_gain * error * psi;
            Matrix const u_follow = m_params.following_gain * epsi;
            Matrix const u = u_seek + u_follow;

            spew("seek control: %.4f, %.4f", u_seek(0), u_seek(1));
            spew("follow control: %.4f, %.4f", u_follow(0), u_follow(1));

            spew("Chl error: %.4f", error);

            double const heading = std::atan2(u(1), u(0));
            spew("Heading: %.2f", Angles::degrees(heading));
            m_state.direction = heading;
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
            m_state.rp.x = 0.0;
            m_state.rp.y = 0.0;

            // Initialize controller state
            m_state.n_wpts = 0;
            m_state.speed = m_args.initial_speed;
            m_state.direction = m_args.initial_heading;

            // Virtual position will be initialized on first ES message.
            m_state.vp.lat = 0.0;
            m_state.vp.lon = 0.0;

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
