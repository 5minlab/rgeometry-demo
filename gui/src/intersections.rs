use core::raster::raster_bounds;
use std::collections::HashMap;

use core::{points_uniform, raster::raster, Rect};

use super::{plot_line, plot_polygon, pt_egui, Demo};
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, PlotPoint, PlotPoints, Polygon},
    Key, Ui,
};
use rgeometry::data::Point;

fn gen_intersections(
    view: f64,
    count: usize,
) -> (
    Vec<(Point<f64>, Point<f64>)>,
    HashMap<(usize, usize), Point<f64>>,
) {
    let mut rng = rand::thread_rng();

    let p0 = points_uniform(&mut rng, view, count);
    let p1 = points_uniform(&mut rng, view, count);
    let lines = (0..p0.len())
        .into_iter()
        .map(|i| match p0[i].cmp(&p1[i]) {
            std::cmp::Ordering::Less => (p0[i], p1[i]),
            _ => (p1[i], p0[i]),
        })
        .collect::<Vec<_>>();

    let intersections = core::intersections::Intersections::new(&lines).sweep();

    (lines, intersections)
}

pub struct DemoIntersections {
    view: f64,
    count: usize,

    lines: Vec<(Point<f64>, Point<f64>)>,
    intersections: HashMap<(usize, usize), Point<f64>>,
}

impl DemoIntersections {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let count = 10;

        let (lines, intersections) = gen_intersections(view, count);

        Self {
            view,
            count,

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
        if ctx.input().key_pressed(Key::R) {
            regen = true;
        }

        ui.horizontal(|ui| {
            if ui.button("regenerate").clicked() {
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
            let (lines, intersections) = gen_intersections(self.view, self.count);
            self.lines = lines;
            self.intersections = intersections;
        }
    }

    fn plot_ui(&self, plot_ui: &mut plot::PlotUi) {
        for (p0, p1) in &self.lines {
            plot_line(plot_ui, &[p0, p1], Color32::RED);
        }

        let mut points = Vec::with_capacity(self.intersections.len());
        for ((_, _), p) in &self.intersections {
            points.push(pt_egui(p));
        }

        plot_ui.points(
            plot::Points::new(PlotPoints::Owned(points))
                .color(Color32::BLUE)
                .radius(2.0),
        );
    }
}
