use super::{pt_egui, Demo};
use crate::delaunay::*;
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, *},
    Key, Ui,
};
use rand::{thread_rng, Rng};
use rgeometry::data::Point;

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
            opt_render_outer_tri: false,
            reductions,
            reductions_max: reductions,
            points,
            net,
        }
    }
}

impl Demo for DemoDelaunay {
    fn name(&self) -> &'static str {
        "delaunay"
    }

    fn ui(&mut self, _t: f64, ctx: &egui::Context, ui: &mut Ui) {
        let mut regen = false;
        if ctx.input().key_pressed(Key::D) {
            regen = true;
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
            ui.separator();
            if ui.button("regenerate").clicked() {
                regen = true;
            }
            ui.separator();
            ui.add(
                egui::Slider::new(&mut self.reductions, 0..=self.reductions_max).text("reductions"),
            );
            ui.separator();
        });
        ui.label("shortcuts: (D) Regenerate | (F) Toggle supertriangles | (C) Step forward | (X) Step backword");

        if r != self.reductions {
            let (net, _) = gen_delaunay(self.view, &self.points, self.reductions);
            self.net = net;
        }

        if regen {
            self.points = gen_delaunay_points(self.view);
            let (net, reductions) = gen_delaunay(self.view, &self.points, std::usize::MAX);
            self.net = net;
            self.reductions = reductions;
            self.reductions_max = reductions;
        }
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        let net = &self.net;

        let v = net.cut(VertIdx(3), VertIdx(self.net.vertices.len() - 1));

        for (t_idx, t) in net.triangles.iter().enumerate() {
            let t_idx = TriIdx(t_idx);
            if let Some(_) = v.cut_triangles.iter().find(|t0| **t0 == t_idx) {
                continue;
            }

            let [v0, v1, v2] = t.vertices;
            let p0 = net.vert(v0);
            let p1 = net.vert(v1);
            let p2 = net.vert(v2);
            if !self.opt_render_outer_tri && (is_super(v0) || is_super(v1) || is_super(v2)) {
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

        for (from, to) in v.cuts {
            let v_from = net.vert(from);
            let v_to = net.vert(to);

            plot_ui.line(
                plot::Line::new(PlotPoints::Owned(vec![pt_egui(&v_from), pt_egui(&v_to)]))
                    .color(Color32::RED),
            );
        }

        for (t_idx, idx) in v.contour_ccw {
            let p_from = net.tri_vert(t_idx, idx.cw());
            let p_to = net.tri_vert(t_idx, idx);

            plot_ui.line(
                plot::Line::new(PlotPoints::Owned(vec![pt_egui(&p_from), pt_egui(&p_to)]))
                    .color(Color32::BLUE),
            );
        }

        for (t_idx, idx) in v.contour_cw {
            let p_from = net.tri_vert(t_idx, idx.cw());
            let p_to = net.tri_vert(t_idx, idx);

            plot_ui.line(
                plot::Line::new(PlotPoints::Owned(vec![pt_egui(&p_from), pt_egui(&p_to)]))
                    .color(Color32::YELLOW),
            );
        }

        plot_ui.points(plot::Points::new(PlotPoints::Owned(
            net.vertices.iter().skip(3).map(|v| pt_egui(v)).collect(),
        )));
    }
}
