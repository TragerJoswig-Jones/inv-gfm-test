// #![no_std]

use std::*;
use plotters::prelude::*;
use ode_solvers::*;
use ode_solvers::dop853::*;
use unifi_gfm::dvoc::DvocController;
use unifi_gfm::dvoc::build_dvoc_controller as dvoc_controller;
use unifi_gfm::sims::*;
use unifi_gfm::calc;
use unifi_gfm::refs::{alpha_beta_fr_polar, alpha_beta_fr_ab};
use unifi_gfm::constants::PI;
//use unifi_gfm::dvoc::test;

const OUT_FILE_NAME: &'static str = "images/sample.png";
const POWER_OUT_FILE_NAME: &'static str = "images/powers.png";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    /* 
    DEFINE sys PARAMETERS & CONSTRUCT OBJECTS
    */
    let v_nom: f32 = 120.;
    let f_nom: f32 = 60.;
    let w_nom: f32 = f_nom * 2. * PI;
    let s_rated: f32 = 500.;
    let dt: f32 = 1.0e-4_f32;
    let xi: f32 = 15.;
    let c: f32 = 0.2679;
    let mut inv = dvoc_controller(v_nom, w_nom, s_rated, dt, xi, c);
    inv.theta = dt * 0.5;

    let rf = 0.8;
    let lf = 1.5e-3;
    let mut line: RLFilter = build_rl_line(w_nom, s_rated, rf, lf);
    let mut bus: ACVoltSrc = build_ac_volt_src(v_nom, w_nom, s_rated);

    // ODE solver settings
    let rtol = 1.0e-5; let atol = 1.0e-6;
    let ode_steps = 20.;
    
    /* 
    RUNNING DYNAMICAL SIMULATION
    */
    //const n_steps: u32 = 10000;  // TODO: Change this to UPPERCASE or make it a variable
    let t_end = 0.4;  // Simulate time in seconds
    let n_steps: u32 = (t_end / dt).ceil() as u32;
    let steps: Vec<u32> = (0..n_steps+1).collect();
    let mut t: f32;
    let mut p: f32; let mut q: f32;
    let mut v_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut theta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut vg_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut thetag_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut p_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut q_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    'sim: for step in steps {  // TODO: Debug this simulation. Power / Current values seem wonky. Something to do with the per unitization of theta
        t = step as f32 * dt;
        if t > (t_end / 2.) {
            inv.set_p_ref(100.0);
        }

        // Collect voltage values
        v_values[step as usize] = (t, inv.v);
        theta_values[step as usize] = (t, inv.theta);
        vg_values[step as usize] = (t, bus.v);
        thetag_values[step as usize] = (t, bus.theta);

        inv.step((line.i_alpha, line.i_beta));
        
        let v = alpha_beta_fr_polar(inv.kv * inv.v, (inv.w_nom * inv.theta) % (2.*PI));
        let i = alpha_beta_fr_ab(line.i_alpha, line.i_beta);
        (p, q) = calc::calc_power(v, i);
        p_values[step as usize] = (t, p);
        q_values[step as usize] = (t, q);

        /* ODE-solver step */
        // Create a new system object and copy over states from existing objects
        let mut ode_line: RLFilter = build_rl_line(w_nom, s_rated, rf, lf);
        let mut ode_bus: ACVoltSrc = build_ac_volt_src(v_nom, w_nom, s_rated);
        let mut ode_inv = dvoc_controller(v_nom, w_nom, s_rated, dt, xi, c);
        ode_line.i_alpha = line.i_alpha; ode_line.i_beta = line.i_beta;
        ode_bus.v = bus.v; ode_bus.theta = bus.theta;
        ode_inv.v = inv.v; ode_inv.theta = inv.theta; ode_inv.p_ref = inv.p_ref;
        let sys = DigitalInverterSMIB { 
            line: ode_line,
            bus: ode_bus,
            inv: ode_inv
        };
        let y0 = State::new(sys.bus.v, sys.bus.theta, sys.line.i_alpha, sys.line.i_beta);
        // Execute ODE solver
        let mut stepper = Dop853::new(sys, t as f64, (t+dt) as f64, (dt / ode_steps) as f64, y0, rtol, atol);
        let res = stepper.integrate();
        // Handle result
        match res {
            Ok(_stats) => {
                // let ts = stepper.x_out();
                let ys = stepper.y_out();
                let k = ys.len() - 1;
                bus.v = ys[k][(0, 0)] as f32; bus.theta = ys[k][(1, 0)] as f32;
                line.i_alpha = ys[k][(2, 0)] as f32; line.i_beta = ys[k][(3, 0)] as f32
            },
            Err(e) => {
                println!("An error occured: {}", e);
                break 'sim;
            },
        }

        
        /* Non-ODE-solver stepper
        // Uncomment to use RK2 with the same step size as the controller for the system
        grid.step_(dt);
        line.step(dt, [dvoc.kv * dvoc.v, (dvoc.w_nom * dvoc.theta) % (2.*PI), 
                       grid.v_nom * grid.v, (grid.w_nom * grid.theta) % (2.*PI)]);
        */

    }
    println!("v: {}, theta: {}", inv.v * inv.v_nom, (inv.w_nom * inv.theta) % (2.*PI));
    println!("ia: {}, ib: {}", line.i_alpha, line.i_beta);
    println!("vg: {}, thetag: {}", bus.v * bus.v_nom, (bus.w_nom * bus.theta) % (2.*PI));

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
        .build_cartesian_2d(0.0f32..t_end, -0.1f32..1.5f32)?;

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
        .label("Grid Voltage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .draw_series(LineSeries::new(
            thetag_values,
            &BLACK,
        ))?
        .label("Grid Angle")
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

    Ok(())
}

/* ODE Solver sys */
//#[derive(Copy, Clone)] // Add this to implement Copy for structs
struct DigitalInverterSMIB {  // TODO: Implement ODE Solver for single digital inverter to bus
    line: RLFilter,  // Line object 
    bus: ACVoltSrc,  // Grid object
    inv: DvocController, // Inverter object
}

type State = Vector4<f32>;

impl ode_solvers::System<State> for DigitalInverterSMIB {
    fn system(&self, _t: f64, y: &State, dy: &mut State) {
        // y: [vg, thetag, ia, ib]
        let (dvg_dt, dthetag_dt) = self.bus.dynamics([y[0], y[1]], []);
        let (di_alpha_dt, di_beta_dt) = self.line.dynamics([y[2], y[3]], [self.inv.kv * self.inv.v, (self.inv.w_nom * self.inv.theta) % (2.*PI), 
                                                                                         self.bus.v_nom * y[0], (self.bus.w_nom * y[1]) % (2.*PI)]);
        dy[0] = 0.; // dvg_dt;
        dy[1] = 1.; // dthetag_dt;
        dy[2] = di_alpha_dt;
        dy[3] = di_beta_dt;
    }
}
