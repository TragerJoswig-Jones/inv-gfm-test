// #![no_std]

use std::*;
use plotters::prelude::*;
use unifi_gfm::dynamics::*;
use unifi_gfm::osg::*;
use unifi_gfm::reference_frames::*;
use unifi_gfm::simulations::*;
use unifi_gfm::constants::*;

const VOLTAGE_FILE_NAME: &'static str = "images/osg_sim_voltage.png";
const THETA_FILE_NAME: &'static str = "images/osg_sim_thetas.png";
const DELTA_OUT_FILE_NAME: &'static str = "images/osg_sim_deltas.png";
const AB_OUT_FILE_NAME: &'static str = "images/osg_sim_alpha_beta_voltage.png";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env::set_var("RUST_BACKTRACE", "1");  // Enable backtrace for identifying overflow errors //TODO: Remove this after testing
    /*
    DEFINE SYSTEM PARAMETERS & CONSTRUCT OBJECTS
    */
    let v_nom: f32 = 80.;
    let f_nom: f32 = 60.;
    let w_nom: f32 = f_nom * 2.*PI;
    let fs: f32 = 10e3_f32; // Hz
    let dt: f32 = 1. / fs;  // s

    // Orthogonal system generator parameters
    let k_osg = 2.0;

    let mut osg = OrthSigGenSogi::new(w_nom, k_osg, rk2_step);

    let mut bus: AcVoltSrc<f32> = build_ac_volt_src(v_nom, w_nom);

    /*
    RUNNING DYNAMICAL SIMULATION
    */
    // Simulation settings
    let t_end = 0.5;  // Simulate time in seconds
    let t_v_step = t_end * 0.5;
    let t_phase_step = t_end * 0.75;
    let n_steps: u32 = (t_end / dt).ceil() as u32;
    let cont_n_steps = 20;  // Number of steps taken for 'continuous' dynamics for each digital step
    let cont_dt = dt / (cont_n_steps as f32);
    let steps: Vec<u32> = (0..n_steps+1).collect();
    // Initialize simulation variables
    let mut t: f32;
    let mut delta: f32; 
    let mut v_grid_sample: f32; let mut theta_grid_sample: f32;
    let grid_voltage = bus.get_x();
    let mut v_grid_sample_alpha_beta: AlphaBeta<f32> = AlphaBeta::from_polar(grid_voltage[(0)], grid_voltage[(1)]);
    // Simulation data vectors
    let mut v_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut theta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut v_alpha_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut v_beta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut vg_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut thetag_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut delta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    for step in steps {
        t = step as f32 * dt;
        if t > t_v_step {
            bus.x[(0)] = 1.1
        }
        if (t >= t_phase_step) & (t < t_phase_step + dt)  {
            bus.x[(1)] += (5. * PI / 180.) * (1. / f_nom)  // 5 deg phase jump in bus voltage
        }

        // Collect voltage values
        let osg_output = osg.get_signals();
        let osg_polar = Polar::from_ab(osg_output[0], osg_output[1], 0.);
        v_alpha_values[step as usize] = (t, osg_output[0]);
        v_beta_values[step as usize] = (t, osg_output[1]);
        v_values[step as usize] = (t, osg_polar.r);
        theta_values[step as usize] = (t, osg_polar.theta / w_nom);
        vg_values[step as usize] = (t, bus.x[(0)]);
        thetag_values[step as usize] = (t, bus.x[(1)]);
        delta = osg_polar.theta / w_nom - bus.x[(1)];
        if delta > (1. / f_nom) - 5e-4 {
            delta = -(1. / f_nom) + delta;
        } else if delta < -(1. / f_nom) + 5e-4  {
            delta = (1. / f_nom) + delta;
        }
        delta_values[step as usize] = (t, delta);

        // Sample the grid
        v_grid_sample = bus.x[(0)];
        theta_grid_sample = bus.x[(1)] * bus.w_nom;
        v_grid_sample_alpha_beta = AlphaBeta::from_polar(v_grid_sample, theta_grid_sample);

        // Step the system
        for _ in 0..cont_n_steps {
            bus.step(cont_dt, []);
        }

        // Step the controller after a z^-1 delay
        let voltage_sample = SQRT_2 * v_grid_sample * libm::cosf(theta_grid_sample);  // Add delay compensation here? theta + dt / 2. * w_nom
        osg.step(dt, [voltage_sample]);
    }
    let osg_output = osg.get_signals();
    let osg_polar = Polar::from_ab(osg_output[0], osg_output[1], 0.);
    println!("v: {}, theta: {}", osg_polar.r, osg_polar.theta / w_nom);
    println!("vg: {}, thetag: {}", bus.x[(0)], bus.x[(1)]);

    /*
    PLOTTING THE RESULTS
    */
    
    let v_alpha_values_ = v_alpha_values.to_vec();
    let v_beta_values_ = v_beta_values.to_vec();
    let v_values_ = v_values.to_vec();
    let theta_values_ = theta_values.to_vec();
    let vg_values_ = vg_values.to_vec();
    let thetag_values_ = thetag_values.to_vec();
    let delta_values_ = delta_values.to_vec();

    /* Plot voltage magnitude data */
    let (_,mut vs): (Vec<f32>, Vec<f32>) = v_values.into_iter().unzip();
    let (_,mut vgs): (Vec<_>, Vec<_>) = vg_values.into_iter().unzip();
    vs.append(&mut vgs);
    let min_v: f32 = vs.iter().fold(f32::INFINITY, |a, &b| a.min(b));
    let max_v: f32 = vs.iter().fold(0.0f32, |a, &b| a.max(b));
    let root = BitMapBackend::new(VOLTAGE_FILE_NAME, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Voltage Magnitudes", ("sans-serif", 50).into_font())
        .margin(5 as u32)
        .x_label_area_size(30 as u32)
        .y_label_area_size(30 as u32)
        .build_cartesian_2d(0.0f32..t_end, min_v..max_v)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            v_values_,
            &RED,
        ))?
        .label("OSG Voltage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    
    chart
        .draw_series(LineSeries::new(
            vg_values_,
            &GREEN,
        ))?
        .label("Bus Voltage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    /* Plot voltage angle data */
    let (_,mut ths): (Vec<f32>, Vec<f32>) = theta_values.into_iter().unzip();
    let (_,mut thgs): (Vec<_>, Vec<_>) = thetag_values.into_iter().unzip();
    ths.append(&mut thgs);
    let min_th: f32 = ths.iter().fold(f32::INFINITY, |a, &b| a.min(b));
    let max_th: f32 = ths.iter().fold(0.0f32, |a, &b| a.max(b));
    let root = BitMapBackend::new(THETA_FILE_NAME, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Voltage Angles", ("sans-serif", 50).into_font())
        .margin(5 as u32)
        .x_label_area_size(30 as u32)
        .y_label_area_size(30 as u32)
        .build_cartesian_2d(0.0f32..t_end, min_th..max_th)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            theta_values_,
            &BLUE,
        ))?
        .label("OSG Theta")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .draw_series(LineSeries::new(
            thetag_values_,
            &BLACK,
        ))?
        .label("Bus Theta")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    /* Plot voltage angle difference data */
    // Find the maximum and minimum values of delta for plotting ylims
    let (_, deltas): (Vec<f32>, Vec<f32>) = delta_values.into_iter().unzip();
    let min_d: f32 = deltas.iter().fold(0.0f32, |a, &b| a.min(b));
    let max_d: f32 = deltas.iter().fold(0.0f32, |a, &b| a.max(b));
    let root = BitMapBackend::new(DELTA_OUT_FILE_NAME, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Voltage Angle Difference", ("sans-serif", 50).into_font())
        .margin(5 as u32)
        .x_label_area_size(30 as u32)
        .y_label_area_size(30 as u32)
        .build_cartesian_2d(0.0f32..t_end, (min_d)..(max_d))?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            delta_values_,
            &RED,
        ))?
        .label("delta")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    /* Plot alpha-beta osg data */
    // Find the maximum and minimum values of ia and ib for plotting ylims
    let (_,mut vas): (Vec<f32>, Vec<f32>) = v_alpha_values.into_iter().unzip();
    let (_,mut vbs): (Vec<_>, Vec<_>) = v_beta_values.into_iter().unzip();
    vas.append(&mut vbs);
    let min_i: f32 = vas.iter().fold(0.0f32, |a, &b| a.min(b));
    let max_i: f32 = vas.iter().fold(0.0f32, |a, &b| a.max(b));
    let root = BitMapBackend::new(AB_OUT_FILE_NAME, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("OSG Signals", ("sans-serif", 50).into_font())
        .margin(5 as u32)
        .x_label_area_size(30 as u32)
        .y_label_area_size(30 as u32)
        .build_cartesian_2d(0.0f32..t_end, (min_i-0.1)..(max_i+0.1))?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            v_alpha_values_,
            &RED,
        ))?
        .label("va")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(
            v_beta_values_,
            &BLUE,
        ))?
        .label("vb")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}
