// #![no_std]

use std::*;
use plotters::prelude::*;
use unifi_gfm::calc::*;
use unifi_gfm::gfm::*;
use unifi_gfm::refs::*;
use unifi_gfm::sims::*;
use unifi_gfm::constants::*;

const VOLTAGE_FILE_NAME: &'static str = "images/dvoc_lcl_dlcv_sim_voltage.png";
const THETA_FILE_NAME: &'static str = "images/dvoc_lcl_dlcv_sim_thetas.png";
const POWER_OUT_FILE_NAME: &'static str = "images/dvoc_lcl_dlcv_sim_powers.png";
const CURRENT_OUT_FILE_NAME: &'static str = "images/dvoc_lcl_dlcv_sim_currents.png";
const DELTA_OUT_FILE_NAME: &'static str = "images/dvoc_lcl_dlcv_sim_deltas.png";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    /*
    DEFINE SYSTEM PARAMETERS & CONSTRUCT OBJECTS
    */
    let v_nom: f32 = 80.;
    let f_nom: f32 = 60.;
    let w_nom: f32 = f_nom * 2.*PI;
    let s_rated: f32 = 1000.;
    let fs: f32 = 100e3_f32; // Hz
    let dt: f32 = 1. / fs;  // s
    let xi: f32 = 15.;
    let c: f32 = 0.2679;
    let gamma: f32 = 1.; 

    let i_base = 3. * v_nom / s_rated;
    let z_base = 3. * v_nom * v_nom / s_rated;

    let rf = 0.4;  // filter-side resistance
    let lf = 1.5e-3;  // filter-side inductance
    let cf = 10e-6;  // filter capacitance
    let rc = 0.05;  // filter capacitor parasitic resistance
    let rg = 0.4;  // grid-side resistance
    let lg = 1.5e-3;  // grid-side inductance
    //let rv = 0;  // virtual impedance
    //let zf = libm::sqrtf(rf*rf+(lf*w_nom)*(lf*w_nom));  // filter nominal inductance

    let w_cur = 2.*PI*2000.;  // TODO: Determine why these frequencies needed to be cranked up this high for tracking. Guessing that per-unitization is the underlying factor
    let w_vol = 2.*PI*500.;   // Originally was using 2.*PI*5000 and 2.*PI*800, but found 1e6 and 3e5 work well (Seperation is a bit low though <10x). Possibly multiply by w_nom, so remove 1/w_nom below?
    
    let kp_v = 1.*w_vol*cf * (z_base);
    let ki_v = 1.*kp_v*w_vol*w_vol/w_cur;
    let kp_i = 1.*lf*w_cur * (1. / z_base);
    let ki_i = 1.*rf*w_cur * (1. / z_base);  // TODO: Is this per-unitization of scalars here correct? Most concerned about w_nom scaling

    let mut inv = build_dvoc_controller(v_nom, w_nom, xi, c);
    let theta0 = 0.0;  // 0.003 for initializing off from grid
    inv.x[(1)] = theta0;  // Initialize inverter angle to be off from the grid to test presync
    inv.x[(0)] = 1.0;  // Initialize inverter voltage to be off from v_nom to test presync
    let inv_ab = AlphaBeta::from_polar(inv.x[(0)], inv.x[(1)] * w_nom);  // Grab alpha-beta inv voltage for initializing the LCL filter
    let mut gfm = build_gfm(&mut inv, gamma);  // Place the dVOC controller within a GFM interface object
    let mut voltage_loop = build_double_loop_voltage_controller(v_nom, kp_v, ki_v, kp_i, ki_i, lf / z_base, cf * z_base, i_base, -i_base);

    let mut line: LclFilter<f32> = build_lcl_filter(w_nom, v_nom, rf / z_base, lf / z_base, rc / z_base, cf * z_base, rg / z_base, lg / z_base);
    line.x[(2)] = inv_ab.alpha;  // Initialize capacitor voltage to align with the inverter voltage
    line.x[(3)] = inv_ab.beta;
    let mut bus: AcVoltSrc<f32> = build_ac_volt_src(v_nom, w_nom);

    /*
    RUNNING DYNAMICAL SIMULATION
    */
    // Simulation settings
    let t_end = 5.0;  // Simulate time in seconds
    let t_step = t_end; //t_end / 2.; // Active power reference step time
    let t_switch = 0.2;  // Grid-side switch time
    let n_steps: u32 = (t_end / dt).ceil() as u32;
    let cont_n_steps = 20;  // Number of steps taken for 'continuous' dynamics for each digital step
    let cont_dt = dt / (cont_n_steps as f32);
    let steps: Vec<u32> = (0..n_steps+1).collect();
    // Initialize simulation variables
    let mut t: f32;
    let mut p: f32 = 0.; let mut q: f32 = 0.; let mut delta: f32; 
    let mut i_alpha_sample: f32; let mut i_beta_sample: f32;
    let mut v_grid_sample: f32; let mut theta_grid_sample: f32;
    let mut v_grid: f32; let mut theta_grid: f32; 
    let mut v_inv: [f32; 2]; let mut v_cap: [f32; 2] = [inv_ab.alpha, inv_ab.beta]; let mut v_gfm: [f32; 2];
    let mut dv_dt_gfm:  nalgebra::SVector<f32, 2> = nalgebra::zero();
    let mut v_inv_alpha_beta: AlphaBeta<f32>; let mut v_cap_alpha_beta: AlphaBeta<f32>; 
    let mut v_grid_sample_alpha_beta: AlphaBeta<f32>; let mut v_to_alpha_beta: AlphaBeta<f32>;
    let mut i_fr: [f32; 2]; let mut i_to: [f32; 2]; let mut v_inv_dq: [f32; 2]; let mut v_inv_polar: Polar<f32>; 
    let mut sin_cos = SinCos::from_theta(theta0 * w_nom); let mut vc_dq = inv_ab.to_dqz(&sin_cos); 
    let mut if_dq = DQZ{d:0., q:0., z:0.}; let mut ig_dq: DQZ<f32> = DQZ{d:0., q:0., z:0.};
    // Simulation data vectors
    let mut v_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut theta_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut vc_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
    let mut thetac_values: Vec<(f32, f32)> = vec![(0., 0.); (n_steps+1) as usize];
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
        // Get inverter output voltage from voltage-loop and convert to polar
        v_gfm = gfm.get_pu_voltage();
        v_inv_dq = voltage_loop.output([v_gfm[0] * SQRT_2, dv_dt_gfm[(1)] * w_nom, vc_dq.d, vc_dq.q, if_dq.d, if_dq.q, ig_dq.d, ig_dq.q]);
        v_inv_polar = Polar::from_dqz(v_inv_dq[0], v_inv_dq[1], 0., &sin_cos);
        v_inv = [v_inv_polar.r, v_inv_polar.theta / w_nom];
        
        // Collect simulation values
        let vc = Polar::from_ab(v_cap[0], v_cap[1], 0.);
        let mut vc_theta = vc.theta;
        if vc_theta > (2.*PI) {
            vc_theta = -2.*PI + vc_theta;
        } else if vc_theta < 0. {
            vc_theta = 2.*PI + vc_theta;
        }
        v_values[step as usize] = (t, v_gfm[0]);
        theta_values[step as usize] = (t, v_gfm[1]);
        vc_values[step as usize] = (t, vc.r);
        thetac_values[step as usize] = (t, vc_theta / w_nom);
        vg_values[step as usize] = (t, bus.x[(0)]);
        thetag_values[step as usize] = (t, bus.x[(1)]);
        ia_values[step as usize] = (t, line.x[(0)]);
        ib_values[step as usize] = (t, line.x[(1)]);
        delta = vc_theta / w_nom - bus.x[(1)];
        if delta > (1. / f_nom - 2.0e-4) {
            delta = -(1. / f_nom) + delta;
        } else if delta < -(1. / f_nom - 1.0e-4) {
            delta = (1. / f_nom) + delta;
        }
        delta_values[step as usize] = (t, delta);

        // Sample the current on the grid side of the LCL filter
        i_alpha_sample = line.x[(4)];
        i_beta_sample = line.x[(5)];
        // Sample the grid
        v_grid_sample = bus.x[(0)];
        theta_grid_sample = bus.x[(1)] * bus.w_nom;
        v_grid_sample_alpha_beta = AlphaBeta::from_polar(v_grid_sample, theta_grid_sample);
        // Get alpha-beta inverter voltage
        v_inv_alpha_beta = AlphaBeta::from_polar(v_inv[0], v_inv[1] * gfm.ctrl.get_w_nom());
        // Get alpha-beta LCL capacitor voltage
        v_cap_alpha_beta = AlphaBeta::from_ab_(v_cap[0], v_cap[1]);

        // Step the system
        for n in 0..cont_n_steps {
            // Calculate power out of the capacitor node 
            let i = AlphaBeta::from_ab_(line.x[(4)], line.x[(5)]);
            (p, q) = calc_ab_power(&v_cap_alpha_beta, &i);  // TODO: Change this to the grid-side power once virtual impedance is implemented
            p_values[(cont_n_steps*step + n) as usize] = (t + (n as f32)*cont_dt, p);
            q_values[(cont_n_steps*step + n) as usize] = (t + (n as f32)*cont_dt, q);
            
            // Step the system
            if t > t_switch {  // Get inductor grid-side voltage based on switch state
                gfm.disable_presync();
                v_grid = bus.x[(0)];  // TODO: rename this for clarity
                theta_grid = bus.x[(1)] * bus.w_nom;
                v_to_alpha_beta = AlphaBeta::from_polar(v_grid, theta_grid);
            } else {
                v_cap = line.get_voltage();
                v_to_alpha_beta = AlphaBeta::from_ab_(v_cap[0], v_cap[1]);
            }
            line.step(cont_dt, [v_inv_alpha_beta.alpha, v_inv_alpha_beta.beta, 
                                   v_to_alpha_beta.alpha, v_to_alpha_beta.beta]);
            bus.step_(cont_dt);
        }

        // Step the controller after a z^-1 delay
        dv_dt_gfm = gfm.gfm_step(dt, [i_alpha_sample, i_beta_sample], [v_grid_sample_alpha_beta.alpha, v_grid_sample_alpha_beta.beta]);
        v_cap = line.get_voltage();
        i_fr = line.get_from_current();
        i_to = line.get_to_current();
        sin_cos = SinCos::from_theta(v_gfm[(1)] * w_nom);
        vc_dq = DQZ::from_ab_(v_cap[0], v_cap[1], &sin_cos);
        if_dq = DQZ::from_ab_(i_fr[0], i_fr[1], &sin_cos);
        ig_dq = DQZ::from_ab_(i_to[0], i_to[1], &sin_cos);
        // u=[E, w_inv, Vd_c, Vq_c, id_f, iq_f, id_g, iq_g]
        voltage_loop.step(dt, [v_gfm[0] * SQRT_2, dv_dt_gfm[(1)] * w_nom, vc_dq.d, vc_dq.q, if_dq.d, if_dq.q, ig_dq.d, ig_dq.q]);  // TODO: Should w_gfm by in per-unit or rad/s?
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
    let vc_values_ = vc_values.to_vec();
    let thetac_values_ = thetac_values.to_vec();
    let delta_values_ = delta_values.to_vec();
    let p_values_ = p_values.to_vec();
    let q_values_ = q_values.to_vec();
    let ia_values_ = ia_values.to_vec();
    let ib_values_ = ib_values.to_vec();

    /* Plot voltage magnitude data */
    let (_,mut vs): (Vec<f32>, Vec<f32>) = v_values.into_iter().unzip();
    let (_,mut vgs): (Vec<_>, Vec<_>) = vg_values.into_iter().unzip();
    let (_,mut vcs): (Vec<_>, Vec<_>) = vc_values.into_iter().unzip();
    vs.append(&mut vgs);
    vs.append(&mut vcs);
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
        .label("GFM Voltage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    
    chart
        .draw_series(LineSeries::new(
            vg_values_,
            &GREEN,
        ))?
        .label("Bus Voltage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));
    chart
        .draw_series(LineSeries::new(
            vc_values_,
            &BLUE,
        ))?
        .label("Cap Voltage")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    /* Plot voltage angle data */
    let (_,mut ths): (Vec<f32>, Vec<f32>) = theta_values.into_iter().unzip();
    let (_,mut thgs): (Vec<_>, Vec<_>) = thetag_values.into_iter().unzip();
    let (_,mut thcs): (Vec<_>, Vec<_>) = thetac_values.into_iter().unzip();
    ths.append(&mut thgs);
    ths.append(&mut thcs);
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
        .label("GFM Theta")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .draw_series(LineSeries::new(
            thetag_values_,
            &BLACK,
        ))?
        .label("Bus Theta")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    chart
        .draw_series(LineSeries::new(
            thetac_values_,
            &RED,
        ))?
        .label("Cap Theta")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

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
        .caption("Angle Difference, Vcap-Vgrid", ("sans-serif", 50).into_font())
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
