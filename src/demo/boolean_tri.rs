use super::{gen_rects, p_rg_to_egui, plot_line, plot_net, Demo, Rect};
use crate::boolean::*;
use crate::demo::TriangularNetwork;
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

fn build_net(view: f64, sx: &SimplicalChain) -> TriangularNetwork {
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
    net
}

pub struct DemoBooleanTri {
    opt_render_rect: bool,
    opt_render_union: bool,

    view: f64,

    rects: Vec<Rect>,
    sx: SimplicalChain,
    net: TriangularNetwork,
}

impl DemoBooleanTri {
    pub fn new(view: f64) -> Self {
        let rects = gen_rects(view, 100);
        let sx = rect_union(&rects);
        let net = build_net(view, &sx);

        Self {
            opt_render_rect: false,
            opt_render_union: true,

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
        });

        for r in &mut self.rects {
            r.rot = t;
        }

        self.sx = rect_union(&self.rects);
        self.net = build_net(self.view, &self.sx);
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        plot_net(&self.net, plot_ui, false);

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
