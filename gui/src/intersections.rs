use std::collections::HashMap;

use core::points_uniform;

use super::{plot_line, pt_egui, Demo};
use eframe::egui::{self, epaint::Color32, Key, Ui};
use egui_plot::{self, *};
use rgeometry::data::Point;

fn gen_intersections(
    view: f64,
    count: usize,
    graph: bool,
) -> (
    Vec<(Point<f64>, Point<f64>)>,
    HashMap<(usize, usize), Point<f64>>,
) {
    let mut rng = rand::thread_rng();

    let lines = if graph {
        let points = points_uniform(&mut rng, view, count);
        let mut lines = Vec::new();
        for idx0 in 0..points.len() {
            for idx1 in (idx0 + 1)..points.len() {
                let p0 = points[idx0];
                let p1 = points[idx1];
                if p0.array[0] < p1.array[0] {
                    lines.push((p0, p1));
                } else {
                    lines.push((p1, p0));
                }
            }
        }
        lines
    } else {
        let p0 = points_uniform(&mut rng, view, count);
        let p1 = points_uniform(&mut rng, view, count);
        (0..p0.len())
            .into_iter()
            .map(|i| match p0[i].cmp(&p1[i]) {
                std::cmp::Ordering::Less => (p0[i], p1[i]),
                _ => (p1[i], p0[i]),
            })
            .collect::<Vec<_>>()
    };

    let intersections = core::intersections::Intersections::new(&lines).sweep();

    (lines, intersections)
}

pub struct DemoIntersections {
    view: f64,
    count: usize,
    graph: bool,

    lines: Vec<(Point<f64>, Point<f64>)>,
    intersections: HashMap<(usize, usize), Point<f64>>,
}

impl DemoIntersections {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let count = 10;
        let graph = false;

        let (lines, intersections) = gen_intersections(view, count, graph);

        Self {
            view,
            count,
            graph,

            lines,
            intersections,
        }
    }
}

impl Demo for DemoIntersections {
    fn name(&self) -> &'static str {
        "Intersections"
    }

    fn ui(&mut self, _t: f64, ctx: &egui::Context, ui: &mut Ui) {
        let mut regen = false;
        if ctx.input(|i| i.key_pressed(Key::R)) {
            regen = true;
        }

        ui.horizontal(|ui| {
            if ui.button("regenerate").clicked() {
                regen = true;
            }
            if ui.checkbox(&mut self.graph, "graph").clicked() {
                regen = true;
            }
            ui.separator();
            for size in [10, 100, 1000, 10000] {
                if ui
                    .radio_value(&mut self.count, size, size.to_string())
                    .clicked()
                {
                    regen = true;
                }
            }
        });

        if regen {
            let (lines, intersections) = gen_intersections(self.view, self.count, self.graph);
            self.lines = lines;
            self.intersections = intersections;
        }
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        for (p0, p1) in &self.lines {
            plot_line(plot_ui, &[p0, p1], Color32::RED);
        }

        let mut points = Vec::with_capacity(self.intersections.len());
        for ((_, _), p) in &self.intersections {
            points.push(pt_egui(p));
        }

        plot_ui.points(
            egui_plot::Points::new(PlotPoints::Owned(points))
                .color(Color32::BLUE)
                .radius(2.0),
        );
    }
}
