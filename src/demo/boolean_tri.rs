use super::{gen_rects, p_rg_to_egui, plot_line, plot_net, pt_mean, Demo, Rect};
use crate::boolean::*;
use crate::delaunay::VertIdx;
use crate::demo::{is_super, TriangularNetwork};
use eframe::egui::{self, epaint::Color32, plot::*, Ui};
use rgeometry::data::Point;

fn rect_union(rects: &[Rect]) -> SimplicalChain {
    let mut sx = SimplicalChain::default();
    for r in rects {
        let p = r.polygon();
        let sx_r = SimplicalChain::from_polygon(&p);
        sx = sx.bool_union(&sx_r);
    }
    sx
}

fn build_net(view: f64, sx: &SimplicalChain, cut: bool) -> TriangularNetwork {
    let v = view * 4.0;
    let mut net = TriangularNetwork::new(
        Point::new([-v, -v]),
        Point::new([v, -v]),
        Point::new([0.0, v]),
    );

    let mut r = std::usize::MAX;
    for s in &sx.simplices {
        if let Err(e) = net.insert(&s.dst, &mut r) {
            eprintln!("TriangularNetwork::insert: {:?}", e);
            return net;
        }
    }

    if cut {
        for s in &sx.simplices {
            let idx0 = net.vertices.iter().position(|p| *p == s.src).unwrap();
            let idx1 = net.vertices.iter().position(|p| *p == s.dst).unwrap();

            let cut = net.cut(VertIdx(idx0), VertIdx(idx1));
            if let Err(e) = net.cut_restore(&cut) {
                eprintln!("failed to cut: {:?}", e);
            }
        }
    }

    net
}

fn plot_net_inner(
    sx: &SimplicalChain,
    net: &TriangularNetwork,
    plot_ui: &mut PlotUi,
    render_supertri: bool,
) {
    for (_t_idx, t) in net.triangles.iter().enumerate() {
        let [v0, v1, v2] = t.vertices;
        let p0 = net.vert(v0);
        let p1 = net.vert(v1);
        let p2 = net.vert(v2);
        if !render_supertri && (is_super(v0) || is_super(v1) || is_super(v2)) {
            continue;
        }
        let center = pt_mean(&[p0, p1, p2]);

        if sx.characteristic(&center) != 1.0 {
            continue;
        }

        plot_line(plot_ui, &[p0, p1, p2, p0], Color32::GREEN);
    }
}

pub struct DemoBooleanTri {
    opt_render_rect: bool,
    opt_render_union: bool,
    opt_cut: bool,

    view: f64,

    rects: Vec<Rect>,
    sx: SimplicalChain,
    net: TriangularNetwork,
}

impl DemoBooleanTri {
    pub fn new(view: f64) -> Self {
        let rects = gen_rects(view, 100);
        let sx = rect_union(&rects);
        let opt_cut = false;
        let net = build_net(view, &sx, opt_cut);

        Self {
            opt_render_rect: false,
            opt_render_union: true,
            opt_cut,

            view,

            rects,
            sx,
            net,
        }
    }
}

impl Demo for DemoBooleanTri {
    fn name(&self) -> &'static str {
        "boolean_tri"
    }

    fn ui(&mut self, t: f64, _ctx: &egui::Context, ui: &mut Ui) {
        ui.horizontal(|ui| {
            if ui.button("regenerate").clicked() {
                self.rects = gen_rects(self.view, 100);
            }
            ui.separator();
            ui.checkbox(&mut self.opt_render_rect, "render rect");
            ui.separator();
            ui.checkbox(&mut self.opt_render_union, "render union");
            ui.separator();
            ui.checkbox(&mut self.opt_cut, "cut");
        });

        for r in &mut self.rects {
            r.rot = t;
        }

        self.sx = rect_union(&self.rects);
        self.net = build_net(self.view, &self.sx, self.opt_cut);
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        plot_net_inner(&self.sx, &self.net, plot_ui, false);

        if self.opt_render_rect {
            for r in &self.rects {
                let p = r.polygon();
                let pe = p_rg_to_egui(&p);

                plot_ui.polygon(pe.color(Color32::RED));
            }
        }

        if self.opt_render_union {
            for s in &self.sx.simplices {
                plot_line(plot_ui, &[&s.src, &s.dst], Color32::BLUE);
            }
        }
    }
}