// #![no_std]

use std::*;
use plotters::prelude::*;
use unifi_gfm::dvoc::build_dvoc_controller as dvoc_controller;
use unifi_gfm::refs::alpha_beta_fr_polar;
use unifi_gfm::constants::PI;
//use unifi_gfm::dvoc::test;

const OUT_FILE_NAME: &'static str = "images/sample.png";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let v_nom: f32 = 120.;
    let f_nom: f32 = 60.;
    let w_nom: f32 = f_nom * 2. * PI;
    let s_rated: f32 = 500.;
    let dt: f32 = 1.0e-4_f32;
    let xi: f32 = 15.;
    let c: f32 = 0.2679;
    let mut dvoc = dvoc_controller(v_nom, w_nom, s_rated, dt, xi, c);
    let u = alpha_beta_fr_polar(0., 0.);
    // dvoc.step((u.alpha, u.beta));
    // println!("v: {}, theta: {}", dvoc.v, dvoc.theta);

    const n_steps: u32 = 20000;
    let t_end = n_steps as f32 * dt;
    let steps = (0..n_steps+1);
    let mut v_values: [(f32, f32); (n_steps + 1) as usize] = [(0., 0.); (n_steps + 1) as usize];
    let mut theta_values: [(f32, f32); (n_steps + 1) as usize] = [(0., 0.); (n_steps + 1) as usize];
    for step in steps {
        v_values[step as usize] = (step as f32, dvoc.v);
        theta_values[step as usize] = (step as f32, dvoc.theta);
        dvoc.step((u.alpha, u.beta))
    }
    println!("{}, {}", dvoc.v * dvoc.v_nom, dvoc.theta * dvoc.w_nom);

    let v_values_ = v_values.to_vec();
    let theta_values_ = theta_values.to_vec();

    // Plot the data
    let root = BitMapBackend::new(OUT_FILE_NAME, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Voltage & Angle", ("sans-serif", 50).into_font())
        .margin(5 as u32)
        .x_label_area_size(30 as u32)
        .y_label_area_size(30 as u32)
        .build_cartesian_2d(0.0f32..n_steps as f32, -0.1f32..1.1f32)?;

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
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}
