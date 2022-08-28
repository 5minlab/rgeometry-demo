use super::{plot_line, plot_net, pt_egui, rotate, Demo};
use crate::delaunay::*;
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, *},
    Key, Ui,
};
use rand::{thread_rng, Rng};
use rgeometry::data::{Point, Vector};

fn gen_delaunay_points(view: f64, len: usize, square: bool) -> Vec<Point<f64>> {
    let mut v = Vec::new();
    if square {
        let grid = 5;
        for i in 0..grid {
            for j in 0..grid {
                // TODO: 2.0 crashes
                let inner = view * 1.9;
                let x = i as f64 / (grid - 1) as f64 * inner - inner / 2.0;
                let y = j as f64 / (grid - 1) as f64 * inner - inner / 2.0;
                v.push(Point::new([x, y]));
            }
        }
    } else {
        let mut rng = thread_rng();
        for _i in 0..len {
            let inner = view;
            let x = rng.gen_range(-inner..inner);
            let y = rng.gen_range(-inner..inner);

            v.push(Point::new([x, y]));
        }
    }
    v
}

fn gen_delaunay(
    view: f64,
    points_constrained: &[Point<f64>],
    points: &[Point<f64>],
) -> TriangularNetwork {
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

    net: TriangularNetwork,
}

impl DemoCDT {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let points_constrained = {
            let mut v = Vec::new();
            let radius = view;
            let p = Vector([radius, 0.0]);
            let len = 8;
            for i in 0..len {
                let theta = std::f64::consts::PI * 2.0 * i as f64 / len as f64;
                v.push(Point::from(rotate(p.clone(), theta)));
            }
            v
        };
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

    fn constraint(&mut self) {
        let net = &mut self.net;
        let len = self.points_constrained.len();
        for i in 0..len {
            let p0 = VertIdx(3 + i);
            let p1 = VertIdx(3 + (i + 1) % len);

            let cut = net.cut(p0, p1);
            if cut.cut_triangles.is_empty() {
                continue;
            }
            net.cut_apply(&cut).ok();
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
            self.constraint();
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
                self.constraint();
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
