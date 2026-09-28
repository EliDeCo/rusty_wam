// The popup the benchmark ends in: one panel per graphed field, each carrying the computed
// cells as markers and the analytic solution as a line, which is how the Sod survey draws
// its own figures.

use super::{FIELDS, Problem, Solved, UNITS};
use eframe::egui;
use egui_plot::{Legend, Line, Plot, Points};

///Samples taken across a rarefaction fan, which is the only part of the analytic curve whose
/// state varies. Everything else is drawn from the two ends of each constant region.
const FAN: usize = 160;

///Size of one panel, fixed rather than grown to fill the window so that three of them
/// always leave room to spare on screen. Width is twice height.
const PLOT_HEIGHT: f32 = 200.0;
const PLOT_WIDTH: f32 = 2.0 * PLOT_HEIGHT;

///Cell count past which the computed solution is drawn as a line rather than as markers.
/// Beyond about one cell per pixel the markers merge into a band anyway, and their geometry
/// is what the renderer has to carry: a marker costs far more of it than a line segment.
const MARKER_LIMIT: usize = PLOT_WIDTH as usize;

///Opens the window and blocks until it is closed.
pub fn show(problem: Problem, t_end: f64, n_cells: usize, gamma: f64, solved: Solved) {
    let length = solved.x.last().copied().unwrap_or(1.0) + solved.x[0];
    let curves = problem.riemann(gamma).curve(length, t_end, FAN);

    let app = Graphs {
        title: format!("{} at t = {t_end} s, {n_cells} cells", problem.name()),
        as_markers: n_cells <= MARKER_LIMIT,
        solved,
        curves,
    };

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([PLOT_WIDTH + 60.0, 3.0 * PLOT_HEIGHT + 80.0]),
        ..Default::default()
    };

    if let Err(e) = eframe::run_native(
        "rusty_wam validation",
        options,
        Box::new(|_cc| Ok(Box::new(app))),
    ) {
        println!("could not open the graph window: {e}");
    }
}

struct Graphs {
    title: String,
    solved: Solved,
    ///the analytic curve of each field, sampled by region rather than uniformly
    curves: [Vec<[f64; 2]>; 3],
    ///whether the computed cells are drawn as individual markers
    as_markers: bool,
}

impl eframe::App for Graphs {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        ui.heading(&self.title);

        for k in 0..3 {
            let cells: Vec<[f64; 2]> = self
                .solved
                .x
                .iter()
                .zip(self.solved.numerical[k].iter())
                .map(|(&x, &y)| [x, y])
                .collect();

            Plot::new(FIELDS[k])
                .height(PLOT_HEIGHT)
                .width(PLOT_WIDTH)
                .legend(Legend::default())
                .x_axis_label("x (m)")
                .y_axis_label(format!("{} ({})", FIELDS[k], UNITS[k]))
                .link_axis(egui::Id::new("distance"), [true, false])
                .show(ui, |plot| {
                    plot.line(Line::new("Exact", self.curves[k].clone()).width(2.0));

                    let name = format!("Computed, L1 {:.3}%", 100.0 * self.solved.errors[k]);
                    match self.as_markers {
                        true => plot.points(Points::new(name, cells).radius(1.0)),
                        false => plot.line(Line::new(name, cells).width(1.0)),
                    }
                });
        }
    }
}
