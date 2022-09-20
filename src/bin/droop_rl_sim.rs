// #![no_std]

use std::*;
use plotters::prelude::*;
use unifi_gfm::calculations::*;
use unifi_gfm::dynamics::*;
use unifi_gfm::gfm::*;
use unifi_gfm::inverter::*;
use unifi_gfm::reference_frames::*;
use unifi_gfm::simulations::*;
use unifi_gfm::constants::*;

const VOLTAGE_FILE_NAME: &'static str = "images/droop_rl_sim_voltage.png";
const THETA_FILE_NAME: &'static str = "images/droop_rl_sim_thetas.png";
const POWER_OUT_FILE_NAME: &'static str = "images/droop_rl_sim_powers.png";
const CURRENT_OUT_FILE_NAME: &'static str = "images/droop_rl_sim_currents.png";
const DELTA_OUT_FILE_NAME: &'static str = "images/droop_rl_sim_deltas.png";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    env::set_var("RUST_BACKTRACE", "1");  // Enable backtrace for identifying overflow errors //TODO: Remove this after testing
    /*
    DEFINE SYSTEM PARAMETERS & CONSTRUCT OBJECTS
    */
    let n_phases = 3.;
    let v_nom: f32 = 80.;
    let f_nom: f32 = 60.;
    let w_nom: f32 = f_nom * 2.*PI;
    let s_rated: f32 = 1000.;
    let fs: f32 = 10e3_f32; // Hz
    let dt: f32 = 1. / fs;  // s

    let rf = 0.8;  // filter-side resistance
    let lf = 1.5e-3;  // filter-side inductance

    let i_base = 3. * v_nom / s_rated;
    let z_base = 3. * v_nom * v_nom / s_rated;

    // droop controller parameters
    let mp: f32 = 0.0026; // / w_nom;  // TODO: Does this coefficient need to be per-unitized? Current values seems to make the response sluggish
    let mq: f32 = 0.005; // / v_nom;   // TODO: Does this coefficient need to be per-unitized?
    let w_c: f32 = 30.*2.*PI;  // TODO: Does this filter frequency need to be per-unitized?
    
    // presynchronization parameteres
    let gamma: f32 = 35.;  // TODO: Determine what value should be used for gamma. Too high causes instability, but too low causes sluggish presync

    // construct controllers and simulation elements
    let mut inv = build_droop_controller(v_nom, w_nom, mp, mq, w_c, n_phases);
    inv.x[(0)] = 0.005;  // Initialize inverter angle to be off from the grid to test presync
    let mut gfm = add_presynch(&mut inv, gamma);  // Place the droop controller within a GFM interface object

    let mut line: RlBranch<f32> = build_rl_branch(i_base, w_nom, rf / z_base, lf / z_base);
    line.open_switch();  // Start with the line disconnected from the ac voltage source
    let bus: AcVoltSrc<f32> = build_ac_volt_src(v_nom, w_nom);
    let mut line_to_bus: LineToBus<f32, 4, 2> = build_line_to_bus(&mut line, bus);  // TODO: Determine if these references need to be mutable and try to make it so that LTB can build its own components


    /*
    RUNNING DYNAMICAL SIMULATION
    */
    // Simulation settings
    let t_end = 0.5;  // Simulate time in seconds
    let t_step = t_end / 2.; // Active power reference step time
    let t_switch = t_end / 5.;  // Grid-side switch time
    let n_steps: u32 = (t_end / dt).ceil() as u32;
    let cont_n_steps = 20;  // Number of steps taken for 'continuous' dynamics for each digital step
    let cont_dt = dt / (cont_n_steps as f32);
    let steps: Vec<u32> = (0..n_steps+1).collect();
    // Initialize simulation variables
    let mut t: f32;
    let mut p: f32 = 0.; let mut q: f32 = 0.; let mut delta: f32; 
    let mut i_alpha_sample: f32; let mut i_beta_sample: f32;
    let mut v_grid_sample: f32; let mut theta_grid_sample: f32;
    let mut v_grid: f32; let mut theta_grid: f32; let mut v_inv: [f32; 2];
    let mut v_to_alpha_beta: AlphaBeta<f32>;
    // Simulation data vectors
    let mut v_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut theta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut vg_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut thetag_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut delta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut p_values: Vec<(f32, f32)> = vec![(0., 0.); ((n_steps+1)*cont_n_steps) as usize];
    let mut q_values: Vec<(f32, f32)> = vec![(0., 0.); ((n_steps+1)*cont_n_steps) as usize];
    let mut ia_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut ib_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    for step in steps {
        t = step as f32 * dt;
        if t > t_step {
            gfm.set_p_ref(1.0);
        }

        // Collect voltage values
        v_inv = gfm.get_pu_voltage();
        v_values[step as usize] = (t, v_inv[0]);
        theta_values[step as usize] = (t, v_inv[1]);
        vg_values[step as usize] = (t, line_to_bus.x[(0)]);
        thetag_values[step as usize] = (t, line_to_bus.x[(1)]);
        ia_values[step as usize] = (t, line_to_bus.x[(2)]);
        ib_values[step as usize] = (t, line_to_bus.x[(3)]);
        delta = v_inv[1] - line_to_bus.x[(1)];
        if delta > (1. / f_nom) {
            delta = -(1. / f_nom) + delta;
        } else if delta < -(1. / f_nom - 1.0e-4) {
            delta = (1. / f_nom) + delta;
        }
        delta_values[step as usize] = (t, delta);

        // Sample the current
        i_alpha_sample = line_to_bus.x[(2)];
        i_beta_sample = line_to_bus.x[(3)];
        // Sample the grid
        v_grid_sample = line_to_bus.x[(0)];
        theta_grid_sample = line_to_bus.x[(1)] * line_to_bus.bus.w_nom;
        let v_grid_sample_alpha_beta = AlphaBeta::from_polar(v_grid_sample, theta_grid_sample);
        // Get alpha-beta gfm voltage
        let v_inv_alpha_beta = AlphaBeta::from_polar(v_inv[0], v_inv[1] * gfm.ctrl.get_w_nom());

        // Step the system
        for n in 0..cont_n_steps {
            // Calculate power
            let i = AlphaBeta::from_ab_(line_to_bus.x[(2)], line_to_bus.x[(3)]);
            (p, q) = calc_ab_power(&v_inv_alpha_beta, &i, n_phases);
            p_values[(cont_n_steps*step + n) as usize] = (t + (n as f32)*cont_dt, p);
            q_values[(cont_n_steps*step + n) as usize] = (t + (n as f32)*cont_dt, q);
            
            // Step the system
            if (t > t_switch) & !line_to_bus.switch_is_closed()  {  // Get inductor grid-side voltage based on switch state
                gfm.disable_presync();
                line_to_bus.close_switch();
            }
            line_to_bus.step(cont_dt, [v_inv_alpha_beta.alpha, v_inv_alpha_beta.beta]);
        }

        // Step the controller after a z^-1 delay
        gfm.inv_step(dt, [i_alpha_sample, i_beta_sample], [v_grid_sample_alpha_beta.alpha, v_grid_sample_alpha_beta.beta]);
    }
    println!("v: {}, theta: {}", inv.x[(0)], inv.x[(1)]);
    println!("ia: {}, ib: {}", line_to_bus.x[(2)], line_to_bus.x[(3)]);
    println!("vg: {}, thetag: {}", line_to_bus.x[(0)], line_to_bus.x[(1)]);
    println!("p: {}, q: {}", p, q);

    /*
    PLOTTING THE RESULTS
    */

    let v_values_ = v_values.to_vec();
    let theta_values_ = theta_values.to_vec();
    let vg_values_ = vg_values.to_vec();
    let thetag_values_ = thetag_values.to_vec();
    let delta_values_ = delta_values.to_vec();
    let p_values_ = p_values.to_vec();
    let q_values_ = q_values.to_vec();
    let ia_values_ = ia_values.to_vec();
    let ib_values_ = ib_values.to_vec();

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
        .label("Inv Voltage")
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
        .label("Inv Theta")
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


    /* Plot the power data */
    // Find the maximum and minimum values of P and Q for plotting ylims
    let (_,mut ps): (Vec<f32>, Vec<f32>) = p_values.into_iter().unzip();
    let (_,mut qs): (Vec<_>, Vec<_>) = q_values.into_iter().unzip();
    ps.append(&mut qs);
    let min_pq: f32 = ps.iter().fold(0.0f32, |a, &b| a.min(b));
    let max_pq: f32 = ps.iter().fold(0.0f32, |a, &b| a.max(b));
    let root = BitMapBackend::new(POWER_OUT_FILE_NAME, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Active & Reactive Powers", ("sans-serif", 50).into_font())
        .margin(5 as u32)
        .x_label_area_size(30 as u32)
        .y_label_area_size(30 as u32)
        .build_cartesian_2d(0.0f32..t_end, (min_pq-1.)..(max_pq+1.))?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            p_values_,
            &RED,
        ))?
        .label("P")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(
            q_values_,
            &BLUE,
        ))?
        .label("Q")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    /* Plot line current data */
    // Find the maximum and minimum values of ia and ib for plotting ylims
    let (_,mut ias): (Vec<f32>, Vec<f32>) = ia_values.into_iter().unzip();
    let (_,mut ibs): (Vec<_>, Vec<_>) = ib_values.into_iter().unzip();
    ias.append(&mut ibs);
    let min_i: f32 = ias.iter().fold(0.0f32, |a, &b| a.min(b));
    let max_i: f32 = ias.iter().fold(0.0f32, |a, &b| a.max(b));
    let root = BitMapBackend::new(CURRENT_OUT_FILE_NAME, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Line Currents", ("sans-serif", 50).into_font())
        .margin(5 as u32)
        .x_label_area_size(30 as u32)
        .y_label_area_size(30 as u32)
        .build_cartesian_2d(0.0f32..t_end, (min_i-0.1)..(max_i+0.1))?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            ia_values_,
            &RED,
        ))?
        .label("ia")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(
            ib_values_,
            &BLUE,
        ))?
        .label("ib")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

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

    Ok(())
}
