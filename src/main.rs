// #![no_std]

use std::*;
use plotters::prelude::*;
use unifi_gfm::dvoc::build_dvoc_controller as dvoc_controller;
use unifi_gfm::sims::build_rl_line as rl_line;
use unifi_gfm::sims::build_ac_volt_src as ac_vs;
use unifi_gfm::sims::Dynamics;
use unifi_gfm::sims::NoInputStep;
use unifi_gfm::sims;
use unifi_gfm::calc;
use unifi_gfm::refs::{alpha_beta_fr_polar, alpha_beta_fr_ab};
use unifi_gfm::constants::PI;
//use unifi_gfm::dvoc::test;

const OUT_FILE_NAME: &'static str = "images/sample.png";
const POWER_OUT_FILE_NAME: &'static str = "images/powers.png";
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
    let mut dvoc = dvoc_controller(v_nom, w_nom, s_rated, dt, xi, c);

    let rf = 0.8;
    let lf = 1.5e-3;
    let mut line: sims::RLFilter = rl_line(w_nom, s_rated, rf, lf);
    let mut grid: sims::ACVoltSrc = ac_vs(v_nom, w_nom, s_rated);

    /* 
    RUNNING DYNAMICAL SIMULATION
    */
    const n_steps: u32 = 40000;  // TODO: Change this to UPPERCASE or make it a variable
    let t_end = n_steps as f32 * dt;
    let steps: Vec<u32> = (0..n_steps+1).collect();
    let mut t: f32;
    let mut p: f32; let mut q: f32;
    let mut v_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut theta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut vg_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut thetag_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut p_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut q_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    for step in steps {  // TODO: Debug this simulation. Power / Current values seem wonky
        t = step as f32 * dt;
        if t > (t_end / 2.) {
            dvoc.set_p_ref(1.0);
        }

        // Collect voltage values
        v_values[step as usize] = (t, dvoc.v);
        theta_values[step as usize] = (t, dvoc.theta);
        vg_values[step as usize] = (t, grid.v);
        thetag_values[step as usize] = (t, grid.theta);

        dvoc.step((line.i_alpha, line.i_beta));
        
        let v = alpha_beta_fr_polar(dvoc.kv * dvoc.v, (dvoc.w_nom * dvoc.theta) % (2.*PI));
        let i = alpha_beta_fr_ab(line.i_alpha, line.i_beta);
        (p, q) = calc::calc_power(v, i);
        p_values[step as usize] = (t, p);
        q_values[step as usize] = (t, q);

        grid.step_(dt);
        line.step(dt, [dvoc.kv * dvoc.v, (dvoc.w_nom * dvoc.theta) % (2.*PI), grid.v, grid.theta]);
    }
    println!("v: {}, theta: {}", dvoc.v * dvoc.v_nom, (dvoc.w_nom * dvoc.theta) % (2.*PI));
    println!("ia: {}, ib: {}", line.i_alpha, line.i_beta);
    println!("vg: {}, thetag: {}", grid.v, grid.theta);

    /* 
    PLOTTING THE RESULTS
    */

    let v_values_ = v_values.to_vec();
    let theta_values_ = theta_values.to_vec();
    let p_values_ = p_values.to_vec();
    let q_values_ = q_values.to_vec();

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
            &RED,
        ))?
        .label("Grid Voltage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(
            thetag_values,
            &BLUE,
        ))?
        .label("Grid Angle")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;


    // Plot the power data
    let root = BitMapBackend::new(POWER_OUT_FILE_NAME, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Active & Reactive Powers", ("sans-serif", 50).into_font())
        .margin(5 as u32)
        .x_label_area_size(30 as u32)
        .y_label_area_size(30 as u32)
        .build_cartesian_2d(0.0f32..t_end, -10.0f32..110.0f32)?;

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

    Ok(())
}
