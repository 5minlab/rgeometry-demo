use super::{plot_line, plot_net, pt_egui, Demo};
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
    let mut v = Vec::new();
    if false {
        let grid = 5;
        for i in 0..grid {
            for j in 0..grid {
                let inner = view * 2.0;
                let x = i as f64 / (grid - 1) as f64 * inner - inner / 2.0;
                let y = j as f64 / (grid - 1) as f64 * inner - inner / 2.0;
                v.push(Point::new([x, y]));
            }
        }
    } else {
        let mut rng = thread_rng();
        for _i in 0..20 {
            let inner = view;
            let x = rng.gen_range(-inner..inner);
            let y = rng.gen_range(-inner..inner);

            v.push(Point::new([x, y]));
        }
    }
    v
}

fn gen_delaunay(view: f64, points: &[Point<f64>], reductions: usize) -> (TriangularNetwork, usize) {
    let v = view * 4.0;
    let mut t = TriangularNetwork::new(
        Point::new([-v, -v]),
        Point::new([v, -v]),
        Point::new([0.0, v]),
    );

    let mut r = reductions;
    for p in points {
        if let Err(e) = t.insert(&p, &mut r) {
            eprintln!("{:?}", e);
            break;
        }
    }
    let used = reductions - r;

    (t, used)
}

pub struct DemoDelaunay {
    view: f64,
    opt_render_supertri: bool,

    reductions: usize,
    reductions_max: usize,
    points: Vec<Point<f64>>,

    net: TriangularNetwork,
    cut: Option<CutResult>,
}

impl DemoDelaunay {
    pub fn new(view: f64) -> Self {
        let points = gen_delaunay_points(view);
        let (net, reductions) = gen_delaunay(view, &points, std::usize::MAX);

        Self {
            view,
            opt_render_supertri: false,
            reductions,
            reductions_max: reductions,
            points,

            cut: None,
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
            self.opt_render_supertri = !self.opt_render_supertri;
        }

        let r = self.reductions;
        if ctx.input().key_pressed(Key::C) && self.reductions < self.reductions_max {
            self.reductions += 1;
        }
        if ctx.input().key_pressed(Key::X) && self.reductions > 0 {
            self.reductions -= 1;
        }

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.opt_render_supertri, "render super");
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
            self.cut = None;
            self.reductions = reductions;
            self.reductions_max = reductions;
        }

        let mut rng = thread_rng();
        if self.reductions != self.reductions_max {
            self.cut = None;
        }

        if self.reductions == self.reductions_max && self.cut.is_none() {
            let cut = loop {
                let idx0 = VertIdx(rng.gen_range(3..self.net.vertices.len()));
                let idx1 = VertIdx(rng.gen_range(3..self.net.vertices.len()));
                if idx0 == idx1 {
                    continue;
                }
                let v = self.net.cut(idx0, idx1);
                if v.cuts.len() < 2 {
                    continue;
                }
                break v;
            };

            self.net.cut_apply(&cut).ok();
            self.cut = Some(cut);
        }
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        let net = &self.net;

        plot_net(net, plot_ui, self.opt_render_supertri);

        if let Some(v) = &self.cut {
            for (from, to) in &v.cuts {
                let p_from = net.vert(*from);
                let p_to = net.vert(*to);
                plot_line(plot_ui, &[p_from, p_to], Color32::RED);
            }

            /*
            for (t_idx, idx) in &v.contour_ccw {
                let p_from = net.tri_vert(*t_idx, idx.cw());
                let p_to = net.tri_vert(*t_idx, *idx);
                plot_line(plot_ui, &[p_from, p_to], Color32::BLUE);
            }

            for (t_idx, idx) in &v.contour_cw {
                let p_from = net.tri_vert(*t_idx, idx.cw());
                let p_to = net.tri_vert(*t_idx, *idx);
                plot_line(plot_ui, &[p_from, p_to], Color32::YELLOW);
            }

            let restore = net.cut_restore(&v);
            for (v0, v1) in restore {
                let p0 = net.vert(v0);
                let p1 = net.vert(v1);
                plot_line(plot_ui, &[p0, p1], Color32::WHITE);
            }
            */
        }

        plot_ui.points(plot::Points::new(PlotPoints::Owned(
            net.vertices.iter().skip(3).map(|v| pt_egui(v)).collect(),
        )));
    }
}
