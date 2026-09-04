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
