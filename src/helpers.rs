use nalgebra::Matrix1xX;
use textplots::{Chart, Plot, Shape};

pub struct ChartDetails {
    pub width: u32,
    pub height: u32,
    pub x_min: f32,
    pub x_max: f32,
    pub y_min: f32,
    pub y_max: f32,
    ///the x values paired with whatever `plot` was last given, filled once and then reused
    pub points: Vec<(f32, f32)>,
}

/// Draws one pipe's field over the chart's x values.
/// The point buffer is refilled rather than rebuilt so the time loop stays allocation free.
pub fn plot(y: &[f64], iteration: u32, pipe_id: usize, chart: &mut ChartDetails) {
    println!("Pipe #{}, Iteration {}", pipe_id, iteration);
    for (point, &value) in chart.points.iter_mut().zip(y.iter()) {
        point.1 = value as f32;
    }
    Chart::new_with_y_range(
        chart.width,
        chart.height,
        chart.x_min,
        chart.x_max,
        chart.y_min,
        chart.y_max,
    )
    .lineplot(&Shape::Points(chart.points.as_slice()))
    .display();
}

/// True when a density or pressure is not a usable physical value.
pub fn unphysical(x: f64) -> bool {
    !x.is_finite() || x <= 0.0
}

///Density, velocity and pressure across a pipe, which is what a pipe is built from.
pub type Initial = (Matrix1xX<f64>, Matrix1xX<f64>, Matrix1xX<f64>);

///Fills the left half of a pipe with one primitive state and the right half with another.
/// Each state is (density, velocity, pressure).
#[allow(dead_code)] //only reached through the cases below, which main.rs selects by name
pub fn two_state(n_cells: usize, left: [f64; 3], right: [f64; 3]) -> Initial {
    let half = n_cells / 2;
    let fill = |l: f64, r: f64| {
        let mut row: Matrix1xX<f64> = Matrix1xX::zeros(n_cells);
        row.columns_mut(0, half).fill(l);
        row.columns_mut(half, n_cells - half).fill(r);
        row
    };

    (
        fill(left[0], right[0]),
        fill(left[1], right[1]),
        fill(left[2], right[2]),
    )
}

// The four two-state cases below keep their leading underscore because only the one named
// in main.rs is ever called, and section 5.1 of the RoeM paper is where the last three
// come from. https://doi.org/10.1016/S0021-9991(02)00037-2

///Sod's problem, the standard shock tube at rest on both sides.
pub fn _sods_problem(n_cells: usize) -> Initial {
    println!("Configuration: Sod's problem.");

    two_state(n_cells, [1.0, 0.0, 1.0], [0.125, 0.0, 0.1])
}

///A shock tube whose two sides already move at the same speed.
pub fn _shock_tube(n_cells: usize) -> Initial {
    println!("Configuration: shock tube.");

    two_state(n_cells, [3.0, 0.9, 3.0], [1.0, 0.9, 1.0])
}

///A pure contact, with pressure and velocity continuous across a density jump of 80.
pub fn _contact_discontinuity(n_cells: usize) -> Initial {
    println!("Configuration: contact discontinuity.");

    two_state(n_cells, [10.0, 0.1125, 1.0], [0.125, 0.1125, 1.0])
}

///Two supersonic streams pulling apart, which is what admits an expansion shock.
pub fn _supersonic_expansion_test(n_cells: usize) -> Initial {
    println!("Configuration: supersonic expansion.");

    two_state(n_cells, [1.0, -2.0, 3.0], [1.0, 2.0, 3.0])
}

///Custom test to confirm MUSCL + Limiter functionality
/// Pulse should stay thin and sharp and remain the same shape and size for the entire simulation.
/// Basic RoeM smears this horizontally, and the conservation of area under the curve causes
/// height to decrease as well. This is INCORRECT behavior
pub fn density_pulse_test(n_cells: usize) -> Initial {
    println!("Configuration: density pulse advection test.");

    let mut rho0: Matrix1xX<f64> = Matrix1xX::from_element(n_cells, 1.0); // rho_bg
    let u0: Matrix1xX<f64> = Matrix1xX::from_element(n_cells, 0.5);
    let p0: Matrix1xX<f64> = Matrix1xX::from_element(n_cells, 1.0);

    // top-hat: 2% of domain, starting near the inlet
    let pulse_start = (0.3 * n_cells as f64) as usize;
    let pulse_end = (0.32 * n_cells as f64) as usize;
    rho0.columns_mut(pulse_start, pulse_end - pulse_start)
        .fill(2.0);

    (rho0, u0, p0)
}

/// At rest
pub fn at_rest(n_cells: usize) -> Initial {
    println!("Pipe at rest.");

    let rho0: Matrix1xX<f64> = Matrix1xX::from_element(n_cells, 1.0);
    let u0: Matrix1xX<f64> = Matrix1xX::zeros(n_cells);
    let p0: Matrix1xX<f64> = Matrix1xX::from_element(n_cells, 1.0);

    (rho0, u0, p0)
}
