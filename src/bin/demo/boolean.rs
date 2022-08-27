use crate::{gen_rects, p_rg_to_egui, pt_egui, Demo, Rect};
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, *},
    Ui,
};
use rgeometry_playground::boolean::*;

pub struct DemoBoolean {
    opt_render_rect: bool,
    opt_render_union: bool,

    rects: Vec<Rect>,
    sx: SimplicalChain,
}

fn rect_union(rects: &[Rect]) -> SimplicalChain {
    let mut sx = SimplicalChain::default();
    for r in rects {
        let p = r.polygon();
        let sx_r = SimplicalChain::from_polygon(&p);
        sx = sx.bool_union(&sx_r);
    }
    sx
}

impl DemoBoolean {
    pub fn new(view: f64) -> Self {
        let rects = gen_rects(view, 100);
        let sx = rect_union(&rects);

        Self {
            opt_render_rect: true,
            opt_render_union: true,

            rects,
            sx,
        }
    }
}

impl Demo for DemoBoolean {
    fn name(&self) -> &'static str {
        "boolean"
    }

    fn ui(&mut self, t: f64, _ctx: &egui::Context, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.opt_render_rect, "render rect");
            ui.separator();
            ui.checkbox(&mut self.opt_render_union, "render union");
        });

        for r in &mut self.rects {
            r.rot = t;
        }
        self.sx = rect_union(&self.rects);
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        if self.opt_render_rect {
            for r in &self.rects {
                let p = r.polygon();
                let pe = p_rg_to_egui(&p);

                plot_ui.polygon(pe.color(Color32::RED));
            }
        }

        if self.opt_render_union {
            for s in &self.sx.simplices {
                plot_ui.line(
                    plot::Line::new(PlotPoints::Owned(vec![pt_egui(&s.src), pt_egui(&s.dst)]))
                        .color(Color32::GREEN),
                );
            }
        }
    }
}
