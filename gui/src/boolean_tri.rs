use super::{p_rg_to_egui, plot_line, pt_egui, Demo};
use core::{
    boolean::*,
    build_net,
    delaunay::{TriIdx, TriangularNetwork, VertIdx},
    gen_rects, points_uniform,
    visibility::VisibilityResult,
    visibility_limit, Rect,
};
use eframe::egui::{self, epaint::Color32, Ui};
use egui_plot::{self, *};
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

fn plot_vis(plot_ui: &mut PlotUi, vis: &VisibilityResult<f64>, color: Color32) {
    if vis.pairs.is_empty() {
        return;
    }

    let mut last = if !vis.arc {
        &vis.pairs[vis.pairs.len() - 1].1
    } else {
        &vis.origin
    };

    for (p0, p1) in &vis.pairs {
        plot_line(plot_ui, &[last, p0], color);
        plot_line(plot_ui, &[p0, p1], color);
        last = p1;
    }

    if vis.arc {
        plot_line(plot_ui, &[last, &vis.origin], color);
    }
}

pub struct DemoBooleanTri {
    opt_render_rect: bool,
    opt_render_union: bool,
    opt_render_tri: bool,
    opt_render_vis: bool,
    opt_render_vis_dir: bool,
    opt_cut: bool,
    opt_prune: bool,
    opt_limit: bool,
    opt_clip: bool,
    opt_bound: bool,

    view: f64,

    rects: Vec<Rect>,
    sx: SimplicalChain<f64>,
    net: TriangularNetwork<f64>,
    constraints: Vec<(VertIdx, VertIdx)>,
    vis: VisibilityResult<f64>,
    vis_dir: VisibilityResult<f64>,

    points: Vec<Point<f64>>,
}

impl DemoBooleanTri {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let mut rng = rand::thread_rng();

        let opt_bound = true;

        eprintln!("view: {:?}", view);
        let mut rects = gen_rects(&mut rng, view, 100);

        let sx = rect_union(&rects);
        let opt_cut = true;
        let (net, constraints) = build_net(view, &sx, opt_cut);

        Self {
            opt_render_rect: false,
            opt_render_union: true,
            opt_render_tri: false,
            opt_render_vis: true,
            opt_render_vis_dir: true,
            opt_cut,
            opt_prune: true,
            opt_limit: false,
            opt_clip: false,
            opt_bound,

            view,

            rects,
            sx,
            net,
            constraints,

            vis: VisibilityResult::empty(Point::new([0.0, 0.0])),
            vis_dir: VisibilityResult::empty(Point::new([0.0, 0.0])),

            points: points_uniform(&mut rng, view, 100),
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

        if ctx.input(|i| i.key_pressed(egui::Key::G)) {
            self.rects = gen_rects(&mut rng, self.view, self.rects.len());
        }

        ui.horizontal(|ui| {
            if ui.button("(G) regenerate").clicked() {
                self.rects = gen_rects(&mut rng, self.view, self.rects.len());
            }
            ui.separator();
            ui.label("render?");
            ui.checkbox(&mut self.opt_render_rect, "rect");
            ui.checkbox(&mut self.opt_render_union, "union");
            ui.checkbox(&mut self.opt_render_tri, "tri");
            ui.checkbox(&mut self.opt_render_vis, "vis");
            ui.checkbox(&mut self.opt_render_vis_dir, "visdir");
            ui.checkbox(&mut self.opt_cut, "cut");
            ui.separator();
            ui.checkbox(&mut self.opt_prune, "prune");
            ui.checkbox(&mut self.opt_limit, "limit");
            ui.checkbox(&mut self.opt_clip, "clip");
            ui.checkbox(&mut self.opt_bound, "bound");
        });

        ui.horizontal(|ui| {
            let mut len = self.points.len();
            for count in [0, 10, 100, 1000, 10000] {
                ui.radio_value(&mut len, count, count.to_string());
            }
            if len != self.points.len() {
                let mut rng = rand::thread_rng();
                self.points = points_uniform(&mut rng, self.view, len);
            }
        });

        let mut rects = self.rects.clone();
        for r in &mut rects {
            r.rot = t;
        }

        if self.opt_bound {
            let view = self.view;
            let width = 1.0;
            rects.push(Rect::new(view + width, width / 2.0).pos(0.0, view + width / 2.0));
            rects.push(Rect::new(view + width, width / 2.0).pos(0.0, -view - width / 2.0));

            rects.push(Rect::new(width / 2.0, view + width).pos(view + width / 2.0, 0.0));
            rects.push(Rect::new(width / 2.0, view + width).pos(-view - width / 2.0, 0.0));
        }

        self.sx = rect_union(&rects);
        let (net, c) = build_net(self.view, &self.sx, self.opt_cut);
        self.net = net;
        self.constraints = c;

        self.vis = {
            let pos = Point::new([0.0, 0.0]);
            let mut vis = self.net.visibility(&self.constraints, &pos).unwrap();

            if self.opt_clip {
                let d0p = Point::new([1.0, 0.0]);
                let d1p = Point::new([0.0, 1.0]);
                let d0 = rgeometry::data::Direction::Through(&d0p);
                let d1 = rgeometry::data::Direction::Through(&d1p);

                vis = vis.clip(d0, d1);
            }

            if self.opt_limit {
                visibility_limit(&mut vis, 15.0f64);
            }
            vis
        };

        self.vis_dir = {
            let pos = Point::new([0.0, 0.0]);
            let mut vis = self
                .net
                .visibility_dir(&self.constraints, &pos, true)
                .unwrap();

            if self.opt_limit {
                visibility_limit(&mut vis, 15.0f64);
            }
            vis
        };
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        if self.opt_render_rect {
            for r in &self.rects {
                let p = r.polygon(1);
                let pe = p_rg_to_egui(&p);

                plot_ui.polygon(pe.fill_color(Color32::RED));
            }
        }

        if self.opt_render_tri {
            plot_net_inner(&self.sx, &self.net, plot_ui, self.opt_prune);
        }

        if self.opt_render_union {
            for s in &self.sx.simplices {
                plot_line(plot_ui, &[&s.src, &s.dst], Color32::BLUE);
            }
        }

        if self.opt_render_vis {
            plot_vis(plot_ui, &self.vis, Color32::YELLOW);
        }

        if self.opt_render_vis_dir {
            plot_vis(plot_ui, &self.vis_dir, Color32::LIGHT_YELLOW);
        }

        let mut points = Vec::new();
        let mut points_visibles = Vec::new();

        for p in &self.points {
            if self.vis.point_visible(p) {
                points_visibles.push(pt_egui(p));
            } else {
                points.push(pt_egui(p));
            }
        }

        plot_ui.points(
            egui_plot::Points::new(PlotPoints::Owned(points))
                .color(Color32::RED)
                .radius(2.0),
        );

        plot_ui.points(
            egui_plot::Points::new(PlotPoints::Owned(points_visibles))
                .color(Color32::LIGHT_BLUE)
                .radius(2.0),
        );
    }
}
