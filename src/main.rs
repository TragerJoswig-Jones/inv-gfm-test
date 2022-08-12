// #![no_std]

use std::*;
use fixed::traits::LossyInto;
use plotters::prelude::*;
use unifi_gfm::calc::*;
use unifi_gfm::dvoc::*;
use unifi_gfm::refs::*;
use unifi_gfm::sims::*;
use unifi_gfm::constants::*;
use fixed::traits::FromFixed;
type FxdSim = fixed::types::I38F26;
type FxdNum = fixed::types::I10F22; // TODO: Test with 32-bit fixed-point number and figure out what is overflowing (Seems to be related to current dynamics)
// TODO: Test how fast this runs with the package having a single fixed-point value selected (No / fewer conversions to fixed). 
// Currently running this sim with I32F32 values takes ~20s

// NOTE:    A 26-bit fractional seems to be the minimal number that results in fairly accurate power tracking / smaller oscillations at steady-state.
//          Oscillations seem to be due to lossy conversions in the simulation, and power tracking errors are due to inaccuracies in current dynamics.

const VOLTAGE_FILE_NAME: &'static str = "images/dvoc_sim_voltage.png";
const THETA_FILE_NAME: &'static str = "images/dvoc_sim_thetas.png";
const POWER_OUT_FILE_NAME: &'static str = "images/dvoc_sim_powers.png";
const CURRENT_OUT_FILE_NAME: &'static str = "images/dvoc_sim_currents.png";
const DELTA_OUT_FILE_NAME: &'static str = "images/dvoc_sim_deltas.png";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    env::set_var("RUST_BACKTRACE", "1");  // Enable backtrace for identifying overflow errors 

    /*
    DEFINE SYSTEM PARAMETERS & CONSTRUCT OBJECTS
    */
    let v_nom: f32 = 80.;
    let f_nom: f32 = 60.;
    let f_nom: f32 = f_nom;
    let s_rated: f32 = 1000.;
    let dt: f32 = 2.0e-4_f32;
    let xi: f32 = 15.;
    let c: f32 = 0.2679;
    let _pi: f32 = PI.lossy_into();

    let z_base = 3. * v_nom * v_nom / s_rated;

    let mut inv = build_dvoc_controller_from_flt::<FxdNum>(v_nom, f_nom, xi, c);
    inv.x[(1)] = FxdNum::from_num(dt * 0.53);  // Initialize inverter angle leading the grid angle by ~half a cycle to start closer to the digital equalibria
    inv.x[(0)] = FxdNum::from_num(0.999965);  // Initialize inverter voltage slightly lower than nominal to start closer to the digital equalibria

    let rf = 0.8;
    let lf = 1.5e-3;
    let mut line: RLFilter<FxdSim> = build_rl_line_from_flt(f_nom, rf / z_base, lf / z_base);
    let mut bus: ACVoltSrc<FxdSim> = build_ac_volt_src_from_flt(v_nom, f_nom);

    /*
    RUNNING DYNAMICAL SIMULATION
    */
    let t_end = 0.5;  // Simulate time in seconds
    let n_steps: u32 = (t_end / dt).ceil() as u32;
    let cont_n_steps = 20;  // Number of steps taken for 'continuous' dynamics for each digital step
    let cont_dt = dt / (cont_n_steps as f32);
    let steps: Vec<u32> = (0..n_steps+1).collect();
    let mut t: f32; let dt_: FxdNum = FxdNum::from_num(dt); let cont_dt_: FxdSim = FxdSim::from_num(cont_dt);
    let mut p: FxdSim = FxdSim::from_fixed(ZERO); let mut q: FxdSim = FxdSim::from_fixed(ZERO);
    let mut delta: f32; let mut i_alpha_sample: FxdNum; let mut i_beta_sample: FxdNum;
    let mut v_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut theta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut vg_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut thetag_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut delta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut p_values: Vec<(f32, f32)> = vec![(0., 0.); ((n_steps+1)*cont_n_steps) as usize];
    let mut q_values: Vec<(f32, f32)> = vec![(0., 0.); ((n_steps+1)*cont_n_steps) as usize];
    let mut ia_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut ib_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    for step in steps {  // TODO: Debug this to see where the overflow occurs...
        t = step as f32 * dt;
        if t > (t_end / 2.) {
            inv.set_p_ref(FxdNum::from_num(1.0));
        }

        // Collect voltage values
        v_values[step as usize] = (t, inv.x[(0)].lossy_into());
        theta_values[step as usize] = (t, inv.x[(1)].lossy_into());
        vg_values[step as usize] = (t, bus.x[(0)].lossy_into());
        thetag_values[step as usize] = (t, bus.x[(1)].lossy_into());
        ia_values[step as usize] = (t, line.x[(0)].lossy_into());
        ib_values[step as usize] = (t, line.x[(1)].lossy_into());
        delta = (inv.x[(1)] - FxdNum::from_fixed(bus.x[(1)])).lossy_into();
        if delta > (1. / f_nom) {
            delta = -(1. / f_nom) + delta;
        } else if delta < -(1. / f_nom - 5.0e-4) {
            delta = (1. / f_nom) + delta;
        }
        delta_values[step as usize] = (t, delta);

        // Sample the current
        i_alpha_sample = FxdNum::from_fixed(line.x[(0)]);
        i_beta_sample = FxdNum::from_fixed(line.x[(1)]);

        // Step the system
        for n in 0..cont_n_steps {
            // Calculate power
            let v = AlphaBeta::from_polar(FxdSim::from_fixed(inv.x[(0)]), FxdSim::from_fixed(inv.x[(1)] * inv.w_nom));
            let i = AlphaBeta::<FxdSim>::from_ab_(line.x[(0)], line.x[(1)]);
            (p, q) = calc_ab_power(v, i);
            p_values[(cont_n_steps*step + n) as usize] = (t + (n as f32)*cont_dt, p.lossy_into());
            q_values[(cont_n_steps*step + n) as usize] = (t + (n as f32)*cont_dt, q.lossy_into());
            bus.step_(cont_dt_);
            line.step(cont_dt_, [FxdSim::from_fixed(inv.x[(0)]), FxdSim::from_fixed(inv.x[(1)] * inv.w_nom), 
                                  bus.x[(0)], bus.x[(1)] * bus.w_nom]);
        }

        // Step the controller after a z^-1 delay
        inv.step(dt_, [i_alpha_sample, i_beta_sample]);
    }
    println!("v: {}, theta: {}", inv.x[(0)], inv.x[(1)]);
    println!("ia: {}, ib: {}", line.x[(0)], line.x[(1)]);
    println!("vg: {}, thetag: {}", bus.x[(0)], bus.x[(1)]);
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
