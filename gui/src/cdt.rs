use super::{plot_line, plot_net, pt_egui, Demo};
use core::delaunay::*;
use core::{points_circular, points_grid, points_uniform};
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, *},
    Key, Ui,
};
use rgeometry::data::Point;

fn gen_delaunay_points(view: f64, len: usize, square: bool) -> Vec<Point<f64>> {
    if square {
        points_grid(view * 0.9, 5)
    } else {
        let mut rng = rand::thread_rng();
        points_uniform(&mut rng, view, len)
    }
}

fn gen_delaunay(
    view: f64,
    points_constrained: &[Point<f64>],
    points: &[Point<f64>],
) -> TriangularNetwork<f64> {
    let v = view * 4.0;
    let mut t = TriangularNetwork::new(
        Point::new([-v, -v]),
        Point::new([v, -v]),
        Point::new([0.0, v]),
    );

    let mut r = std::usize::MAX;
    for p in points_constrained {
        if let Err(e) = t.insert(&p, &mut r) {
            eprintln!("TriangularNetwork::insert: {:?}", e);
            return t;
        }
    }
    for p in points {
        if let Err(e) = t.insert(&p, &mut r) {
            eprintln!("TriangularNetwork::insert: {:?}", e);
            return t;
        }
    }

    t
}

pub struct DemoCDT {
    view: f64,
    opt_render_supertri: bool,
    opt_square: bool,

    points_constrained: Vec<Point<f64>>,
    points: Vec<Point<f64>>,

    net: TriangularNetwork<f64>,
}

impl DemoCDT {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let points_constrained = points_circular(view, 8);
        let opt_square = true;
        let points_count = 50;

        let points = gen_delaunay_points(view, points_count, opt_square);
        let net = gen_delaunay(view, &points_constrained, &points);

        Self {
            view,

            opt_render_supertri: false,
            opt_square,
            points_constrained,
            points,

            net,
        }
    }

    fn constrain(&mut self) {
        let net = &mut self.net;
        let len = self.points_constrained.len();
        for i in 0..len {
            let p0 = VertIdx(3 + i);
            let p1 = VertIdx(3 + (i + 1) % len);

            if let Err(e) = net.constrain_edge(p0, p1) {
                eprintln!("constrain_edge: {:?}", e);
                continue;
            }
        }
    }
}

impl Demo for DemoCDT {
    fn name(&self) -> &'static str {
        "CDT"
    }

    fn ui(&mut self, _t: f64, ctx: &egui::Context, ui: &mut Ui) {
        let mut regen = false;
        if ctx.input().key_pressed(Key::D) {
            regen = true;
        }
        if ctx.input().key_pressed(Key::F) {
            self.opt_render_supertri = !self.opt_render_supertri;
        }
        if ctx.input().key_pressed(Key::G) {
            self.constrain();
        }

        ui.horizontal(|ui| {
            if ui.checkbox(&mut self.opt_square, "square").clicked() {
                regen = true;
            }
            ui.separator();
            ui.checkbox(&mut self.opt_render_supertri, "render super");
            ui.separator();
            if ui.button("regenerate").clicked() {
                regen = true;
            }
            ui.separator();
            if ui.button("force constraint").clicked() {
                self.constrain();
            }
        });
        ui.label("shortcuts: (D) Regenerate | (F) Toggle supertriangles | (G) Constraint");

        if regen {
            self.points = gen_delaunay_points(self.view, self.points.len(), self.opt_square);
            self.net = gen_delaunay(self.view, &self.points_constrained, &self.points);
        }
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        let net = &self.net;

        plot_net(net, plot_ui, self.opt_render_supertri);

        let l = self.points_constrained.len();
        for i in 0..l {
            let idx0 = VertIdx(3 + i);
            let idx1 = VertIdx(3 + (i + 1) % l);
            let p0 = net.vert(idx0);
            let p1 = net.vert(idx1);
            plot_line(plot_ui, &[p0, p1], Color32::RED);
        }

        plot_ui.points(plot::Points::new(PlotPoints::Owned(
            net.vertices.iter().skip(3).map(|v| pt_egui(v)).collect(),
        )));
    }
}
