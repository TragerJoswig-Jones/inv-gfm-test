// #![no_std]

use std::*;
use plotters::prelude::*;
use unifi_gfm::dvoc::*;
use unifi_gfm::sims::*;
use unifi_gfm::calc;
use unifi_gfm::refs::{alpha_beta_fr_polar, alpha_beta_fr_ab};
use unifi_gfm::constants::PI;
//use unifi_gfm::inv::test;

const OUT_FILE_NAME: &'static str = "images/voltages.png";
const POWER_OUT_FILE_NAME: &'static str = "images/powers.png";
const CURRENT_OUT_FILE_NAME: &'static str = "images/currents.png";
const DELTA_OUT_FILE_NAME: &'static str = "images/deltas.png";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    /* 
    DEFINE SYSTEM PARAMETERS & CONSTRUCT OBJECTS
    */
    let v_nom: f32 = 120.;
    let f_nom: f32 = 60.;
    let w_nom: f32 = f_nom * 2. * PI;
    let s_rated: f32 = 500.;
    let dt: f32 = 1.0e-4_f32;
    let xi: f32 = 15.;
    let c: f32 = 0.2679;
    let mut inv = build_dvoc_controller(v_nom, w_nom, s_rated, xi, c);
    inv.theta = -dt * 0.5;  // Initialize inverter angle leading the grid angle

    let rf = 0.8;
    let lf = 1.5e-3;
    let mut line: RLFilter = build_rl_line(w_nom, s_rated, rf, lf);
    let mut bus: ACVoltSrc = build_ac_volt_src(v_nom, w_nom, s_rated);

    // TODO: Debug simulation results being weird for cont_n_steps > 5. Seems to be something happening in the grid or line dynamics every some period

    /* 
    RUNNING DYNAMICAL SIMULATION
    */
    let t_end = 0.2;  // Simulate time in seconds
    let n_steps: u32 = (t_end / dt).ceil() as u32;
    let cont_n_steps = 50;  // Number of steps taken for 'continuous' dynamics for each digital step
    let cont_dt = dt / (cont_n_steps as f32);
    let steps: Vec<u32> = (0..n_steps+1).collect();
    let mut t: f32;
    let mut p: f32 = 0.; let mut q: f32 = 0.;
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
        if t > (t_end / 4.) {
            inv.set_p_ref(0.0);
        }
        
        if t > (t_end / 2.) {
            inv.set_p_ref(100.0);
        }

        // Collect voltage values
        v_values[step as usize] = (t, inv.v);
        theta_values[step as usize] = (t, inv.theta);
        vg_values[step as usize] = (t, bus.v);
        thetag_values[step as usize] = (t, bus.theta);
        ia_values[step as usize] = (t, line.i_alpha);
        ib_values[step as usize] = (t, line.i_beta);
        let delta = inv.theta - bus.theta;
        delta_values[step as usize] = (t, delta);

        // Step the controller
        inv.step(dt, [line.i_alpha, line.i_beta]);

        // Step the system
        for n in 0..cont_n_steps {
            // Calculate power
            let v = alpha_beta_fr_polar(inv.kv * inv.v, (inv.w_nom * inv.theta) % (2.*PI));
            let i = alpha_beta_fr_ab(line.i_alpha, line.i_beta);
            (p, q) = calc::calc_power(v, i);
            p_values[(cont_n_steps*step + n) as usize] = (t + (n as f32)*cont_dt, p);
            q_values[(cont_n_steps*step + n) as usize] = (t + (n as f32)*cont_dt, q);
            bus.step_(cont_dt);
            line.step(cont_dt, [inv.kv * inv.v, (inv.w_nom * inv.theta) % (2.*PI), 
                           bus.v_nom * bus.v, (bus.w_nom * bus.theta) % (2.*PI)]);
        }
    }
    println!("v: {}, theta: {}", inv.v * inv.v_nom, (inv.w_nom * inv.theta) % (2.*PI));
    println!("ia: {}, ib: {}", line.i_alpha, line.i_beta);
    println!("vg: {}, thetag: {}", bus.v * bus.v_nom, (bus.w_nom * bus.theta) % (2.*PI));
    println!("p: {}, q: {}", p, q);

    /* 
    PLOTTING THE RESULTS
    */

    let v_values_ = v_values.to_vec();
    let theta_values_ = theta_values.to_vec();
    let delta_values_ = delta_values.to_vec();
    let p_values_ = p_values.to_vec();
    let q_values_ = q_values.to_vec();
    let ia_values_ = ia_values.to_vec();
    let ib_values_ = ib_values.to_vec();

    // Plot the data
    let root = BitMapBackend::new(OUT_FILE_NAME, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Voltage & Angle", ("sans-serif", 50).into_font())
        .margin(5 as u32)
        .x_label_area_size(30 as u32)
        .y_label_area_size(30 as u32)
        .build_cartesian_2d(0.0f32..t_end, -0.1f32..1.1f32)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            v_values_,
            &RED,
        ))?
        .label("Voltage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(
            theta_values_,
            &BLUE,
        ))?
        .label("Angle")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    
    chart
        .draw_series(LineSeries::new(
            vg_values,
            &GREEN,
        ))?
        .label("bus Voltage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(
            thetag_values,
            &BLACK,
        ))?
        .label("bus Angle")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
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
