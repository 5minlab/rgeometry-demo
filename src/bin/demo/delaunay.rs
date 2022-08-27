use crate::pt_egui;
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, *},
    Key, Ui,
};
use rand::{thread_rng, Rng};
use rgeometry::data::Point;
use rgeometry_playground::delaunay::TriangularNetwork;

fn gen_delaunay_points(view: f64) -> Vec<Point<f64>> {
    let mut rng = thread_rng();
    let mut v = Vec::new();
    for _i in 0..20 {
        let inner = view;
        let x = rng.gen_range(-inner..inner);
        let y = rng.gen_range(-inner..inner);

        v.push(Point::new([x, y]));
    }
    v
}

fn gen_delaunay(view: f64, points: &[Point<f64>], reductions: usize) -> (TriangularNetwork, usize) {
    let v = view * 3.0;
    let mut t = TriangularNetwork::new(
        Point::new([-v, -v]),
        Point::new([v, -v]),
        Point::new([0.0, v]),
    );

    let mut r = reductions;
    for p in points {
        t.insert(&p, &mut r);
    }
    let used = reductions - r;

    (t, used)
}

pub struct DemoDelaunay {
    view: f64,
    opt_render_outer_tri: bool,

    reductions: usize,
    reductions_max: usize,
    points: Vec<Point<f64>>,
    net: TriangularNetwork,
}

impl DemoDelaunay {
    pub fn new(view: f64) -> Self {
        let points = gen_delaunay_points(view);
        let (net, reductions) = gen_delaunay(view, &points, std::usize::MAX);

        Self {
            view,
            opt_render_outer_tri: true,
            reductions,
            reductions_max: reductions,
            points,
            net,
        }
    }

    pub fn ui(&mut self, ctx: &egui::Context, ui: &mut Ui) {
        if ctx.input().key_pressed(Key::D) {
            self.points = gen_delaunay_points(self.view);
            let (net, reductions) = gen_delaunay(self.view, &self.points, std::usize::MAX);
            self.net = net;
            self.reductions = reductions;
            self.reductions_max = reductions;
        }
        if ctx.input().key_pressed(Key::F) {
            self.opt_render_outer_tri = !self.opt_render_outer_tri;
        }

        let r = self.reductions;
        if ctx.input().key_pressed(Key::C) && self.reductions < self.reductions_max {
            self.reductions += 1;
        }
        if ctx.input().key_pressed(Key::X) && self.reductions > 0 {
            self.reductions -= 1;
        }

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.opt_render_outer_tri, "render super");
            ui.add(
                egui::Slider::new(&mut self.reductions, 0..=self.reductions_max).text("reductions"),
            );
        });
        ui.label("shortcuts: (D) Regenerate | (F) Toggle supertriangles | (C) Step forward | (X) Step backword");

        if r != self.reductions {
            let (net, _) = gen_delaunay(self.view, &self.points, self.reductions);
            self.net = net;
        }
    }

    pub fn plot_ui(&self, plot_ui: &mut PlotUi) {
        let net = &self.net;
        for t in &net.triangles {
            let [v0, v1, v2] = t.vertices;
            let p0 = net.vertices[v0];
            let p1 = net.vertices[v1];
            let p2 = net.vertices[v2];
            if !self.opt_render_outer_tri && (v0 < 3 || v1 < 3 || v2 < 3) {
                continue;
            }

            plot_ui.line(
                plot::Line::new(PlotPoints::Owned(vec![
                    pt_egui(&p0),
                    pt_egui(&p1),
                    pt_egui(&p2),
                    pt_egui(&p0),
                ]))
                .color(Color32::GREEN),
            );
        }

        plot_ui.points(plot::Points::new(PlotPoints::Owned(
            net.vertices.iter().skip(3).map(|v| pt_egui(v)).collect(),
        )));
    }
}
