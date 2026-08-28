use glam::DVec3;

/// Generates `n` unit vectors spread roughly evenly over the sphere.
pub fn spread_directions(n: usize) -> Vec<DVec3> {
    // ~2.399963 rad — the golden angle
    let golden_angle = std::f64::consts::PI * (3.0 - 5.0f64.sqrt());

    (0..n)
        .map(|i| {
            let i = i as f64;
            // Walk y evenly from +1 to -1 so bands have equal area
            let y = 1.0 - (i + 0.5) / n as f64 * 2.0;
            let radius = (1.0 - y * y).sqrt();
            let theta = golden_angle * i;

            DVec3::new(theta.cos() * radius, y, theta.sin() * radius)
        })
        .collect()
}