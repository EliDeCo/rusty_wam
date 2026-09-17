use nalgebra::Vector3;
use textplots::{Chart, Plot, Shape};

pub struct ChartDetails {
    pub width: u32,
    pub height: u32,
    pub x_min: f32,
    pub x_max: f32,
    pub y_min: f32,
    pub y_max: f32,
    pub x: Vec<f64>,
}

pub fn plot(y: &[f64], iteration: u32, pipe_id: usize, chart: &ChartDetails) {
    println!("Pipe #{}, Iteration {}", pipe_id, iteration);
    let points: Vec<(f32, f32)> = chart
        .x
        .iter()
        .copied()
        .map(|y| y as f32)
        .zip(y.iter().copied().map(|y| y as f32))
        .collect();
    Chart::new_with_y_range(
        chart.width,
        chart.height,
        chart.x_min,
        chart.x_max,
        chart.y_min,
        chart.y_max,
    )
    .lineplot(&Shape::Points(points.as_slice()))
    .display();
}

/// Generates `n` unit vectors spread roughly evenly over the sphere.
pub fn spread_directions(n: usize) -> Vec<Vector3<f64>> {
    // ~2.399963 rad — the golden angle
    let golden_angle = std::f64::consts::PI * (3.0 - 5.0f64.sqrt());

    (0..n)
        .map(|i| {
            let i = i as f64;
            // Walk y evenly from +1 to -1 so bands have equal area
            let y = 1.0 - (i + 0.5) / n as f64 * 2.0;
            let radius = (1.0 - y * y).sqrt();
            let theta = golden_angle * i;

            Vector3::new(theta.cos() * radius, y, theta.sin() * radius)
        })
        .collect()
}

/// True when a density or pressure is not a usable physical value.
/// Written as `!is_finite() || <= 0.0` rather than `< 0.0` so that NaN and infinity
/// are caught too - NaN fails every ordinary comparison, so `x < 0.0` silently
/// passes it through.
pub fn unphysical(x: f64) -> bool {
    !x.is_finite() || x <= 0.0
}
