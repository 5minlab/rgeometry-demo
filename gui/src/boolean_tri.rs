use super::{p_rg_to_egui, plot_line, Demo};
use core::boolean::*;
use core::delaunay::{TriIdx, TriangularNetwork, VertIdx};
use core::{build_net, gen_rects, Rect};
use eframe::egui::{self, epaint::Color32, plot::*, Ui};
use rgeometry::data::Point;

fn rect_union(rects: &[Rect]) -> SimplicalChain<f64> {
    let mut sx = SimplicalChain::default();
    for r in rects {
        let p = r.polygon(1);
        let sx_r = SimplicalChain::from_polygon(&p);
        sx = sx.union(&sx_r);
    }
    sx
}

fn plot_net_inner(
    sx: &SimplicalChain<f64>,
    net: &TriangularNetwork<f64>,
    plot_ui: &mut PlotUi,
    prune: bool,
) {
    for (_t_idx, t) in net.triangles.iter().enumerate() {
        let [v0, v1, v2] = t.vertices;
        let p0 = net.vert(v0);
        let p1 = net.vert(v1);
        let p2 = net.vert(v2);
        let center = net.centroid(TriIdx(_t_idx));

        if prune && sx.characteristic(&center) != 1.0 {
            continue;
        }

        plot_line(plot_ui, &[p0, p1, p2, p0], Color32::GREEN);
    }
}

pub struct DemoBooleanTri {
    opt_render_rect: bool,
    opt_render_union: bool,
    opt_render_tri: bool,
    opt_render_vis: bool,
    opt_cut: bool,
    opt_prune: bool,
    opt_limit: bool,

    view: f64,

    rects: Vec<Rect>,
    sx: SimplicalChain<f64>,
    net: TriangularNetwork<f64>,
    constraints: Vec<(VertIdx, VertIdx)>,
    vis: Vec<(Point<f64>, Point<f64>)>,
}

impl DemoBooleanTri {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let mut rng = rand::thread_rng();

        let rects = gen_rects(&mut rng, view, 100);
        let sx = rect_union(&rects);
        let opt_cut = true;
        let (net, constraints) = build_net(view, &sx, opt_cut);

        Self {
            opt_render_rect: false,
            opt_render_union: true,
            opt_render_tri: false,
            opt_render_vis: true,
            opt_cut,
            opt_prune: true,
            opt_limit: true,

            view,

            rects,
            sx,
            net,
            constraints,

            vis: Vec::new(),
        }
    }
}

impl Demo for DemoBooleanTri {
    fn name(&self) -> &'static str {
        "boolean_tri"
    }

    fn ui(&mut self, t: f64, ctx: &egui::Context, ui: &mut Ui) {
        use rand::prelude::*;
        use rand_chacha::*;

        let mut seed: <ChaCha20Rng as SeedableRng>::Seed = Default::default();
        rand::thread_rng().fill(&mut seed);
        let mut rng = ChaCha20Rng::from_seed(seed);

        if ctx.input().key_pressed(egui::Key::G) {
            self.rects = gen_rects(&mut rng, self.view, self.rects.len());
        }

        ui.horizontal(|ui| {
            if ui.button("(R) regenerate").clicked() {
                self.rects = gen_rects(&mut rng, self.view, self.rects.len());
            }
            ui.separator();
            ui.label("render?");
            ui.checkbox(&mut self.opt_render_rect, "rect");
            ui.checkbox(&mut self.opt_render_union, "union");
            ui.checkbox(&mut self.opt_render_tri, "tri");
            ui.checkbox(&mut self.opt_render_vis, "vis");
            ui.checkbox(&mut self.opt_cut, "cut");
            ui.separator();
            ui.checkbox(&mut self.opt_prune, "prune");
            ui.checkbox(&mut self.opt_limit, "limit");
        });

        for r in &mut self.rects {
            r.rot = t;
        }

        self.sx = rect_union(&self.rects);
        let (net, c) = build_net(self.view, &self.sx, self.opt_cut);
        self.net = net;
        self.constraints = c;

        let vis = {
            let pos = Point::new([0.0, 0.0]);
            let mut vis = self.net.visibility(&self.constraints, &pos);

            if self.opt_limit {
                let limit = 15.0f64;
                let limit_sq = limit * limit;

                for i in 0..vis.len() {
                    let (ref mut p0, ref mut p1) = vis[i];

                    let dist: f64 = p0.squared_euclidean_distance(&pos);
                    if dist > limit_sq {
                        let ratio = (limit_sq / dist).sqrt();
                        let v = *p0 - pos;
                        *p0 = pos + v * ratio;
                    }

                    let dist: f64 = p1.squared_euclidean_distance(&pos);
                    if dist > limit_sq {
                        let ratio = (limit_sq / dist).sqrt();
                        let v = *p1 - pos;
                        *p1 = pos + v * ratio;
                    }
                }
            }
            vis
        };

        self.vis = vis;
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        if self.opt_render_rect {
            for r in &self.rects {
                let p = r.polygon(1);
                let pe = p_rg_to_egui(&p);

                plot_ui.polygon(pe.color(Color32::RED));
            }
        }

        if self.opt_render_union {
            for s in &self.sx.simplices {
                plot_line(plot_ui, &[&s.src, &s.dst], Color32::BLUE);
            }
        }

        if self.opt_render_tri {
            plot_net_inner(&self.sx, &self.net, plot_ui, self.opt_prune);
        }

        if self.opt_render_vis {
            let mut last = &self.vis[self.vis.len() - 1].1;
            for (p0, p1) in &self.vis {
                plot_line(plot_ui, &[last, p0], Color32::YELLOW);
                plot_line(plot_ui, &[p0, p1], Color32::YELLOW);
                last = p1;
            }
        }
    }
}
