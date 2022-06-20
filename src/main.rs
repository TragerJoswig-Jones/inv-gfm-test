use unifi_gfm::dvoc::build_dvoc_controller as dvoc_controller;
use unifi_gfm::refs::alpha_beta_fr_polar;
use unifi_gfm::constants::PI;
//use unifi_gfm::dvoc::test;

fn main() {
    let v_nom: f32 = 120.;
    let f_nom: f32 = 60.;
    let w_nom: f32 = f_nom * 2. * PI;
    let s_rated: f32 = 500.;
    let dt: f32 = 1.0e-4_f32;
    let xi: f32 = 15.;
    let c: f32 = 0.2679;
    let mut dvoc = dvoc_controller(v_nom, w_nom, s_rated, dt, xi, c);
    let u = alpha_beta_fr_polar(0., 0.);
    dvoc.step((u.alpha, u.beta));
    println!("v: {}, theta: {}", dvoc.v, dvoc.theta)
}
