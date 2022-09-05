use super::{p_rg_to_egui, plot_line, Demo};
use core::boolean::*;
use core::{gen_rects, points_circular, Rect};
use eframe::egui::{self, epaint::Color32, plot::*, Ui};
use rgeometry::data::{Point, Polygon};

#[derive(Clone, Copy, PartialEq, Eq)]
enum CircleMode {
    Disabled,
    Union,
    Intersect,
    Subtract0,
    Subtract1,
}

pub struct DemoBoolean {
    opt_render_rect: bool,
    opt_render_union: bool,

    opt_circle_mode: CircleMode,

    circle: Vec<Point<f64>>,

    view: f64,

    rects: Vec<Rect>,
    sx: SimplicalChain<f64>,
}

fn rect_union(rects: &[Rect]) -> SimplicalChain<f64> {
    let mut sx = SimplicalChain::default();
    for r in rects {
        let p = r.polygon(1);
        let sx_r = SimplicalChain::from_polygon(&p);
        sx = sx.union(&sx_r);
    }
    sx
}

impl DemoBoolean {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let mut rng = rand::thread_rng();
        let rects = gen_rects(&mut rng, view, 100);
        let sx = rect_union(&rects);

        Self {
            opt_render_rect: true,
            opt_render_union: true,

            opt_circle_mode: CircleMode::Disabled,

            circle: points_circular(view / 2.0, 32),

            view,

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
        let mut rng = rand::thread_rng();
        ui.horizontal(|ui| {
            if ui.button("regenerate").clicked() {
                self.rects = gen_rects(&mut rng, self.view, 100);
            }
            ui.separator();

            for (mode, label) in &[
                (CircleMode::Disabled, "none"),
                (CircleMode::Union, "union"),
                (CircleMode::Intersect, "intersect"),
                (CircleMode::Subtract0, "subtract0"),
                (CircleMode::Subtract1, "subtract1"),
            ] {
                ui.radio_value(&mut self.opt_circle_mode, *mode, *label);
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

        if self.opt_circle_mode == CircleMode::Disabled {
            return;
        }

        let p = Polygon::new(self.circle.clone()).unwrap();
        let sx_circle = SimplicalChain::from_polygon(&p);

        match self.opt_circle_mode {
            CircleMode::Union => {
                self.sx = self.sx.union(&sx_circle);
            }
            CircleMode::Intersect => {
                self.sx = self.sx.intersect(&sx_circle);
            }
            CircleMode::Subtract0 => {
                self.sx = self.sx.subtract(&sx_circle);
            }
            CircleMode::Subtract1 => {
                self.sx = sx_circle.subtract(&self.sx);
            }
            _ => (),
        }
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
                plot_line(plot_ui, &[&s.src, &s.dst], Color32::GREEN);
            }
        }
    }
}
