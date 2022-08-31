use crate::rester::rester;

use super::{plot_line, points_tri, Demo};
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, PlotPoint, PlotPoints, Polygon},
    Key, Ui,
};
use rgeometry::data::Point;

pub struct DemoRester {
    view: f64,

    verts: [Point<f64>; 3],
}

impl DemoRester {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let mut rng = rand::thread_rng();
        Self {
            view,
            verts: points_tri(&mut rng, view),
        }
    }
}

impl Demo for DemoRester {
    fn name(&self) -> &'static str {
        "Rester"
    }

    fn ui(&mut self, _t: f64, ctx: &egui::Context, ui: &mut Ui) {
        let mut regen = false;
        if ctx.input().key_pressed(Key::R) {
            regen = true;
        }

        ui.horizontal(|ui| {
            if ui.button("regenerate").clicked() {
                regen = true;
            }
        });

        if regen {
            let mut rng = rand::thread_rng();
            self.verts = points_tri(&mut rng, self.view);
        }
    }

    fn plot_ui(&self, plot_ui: &mut plot::PlotUi) {
        let [ref p0, ref p1, ref p2] = self.verts;
        plot_line(plot_ui, &[p0, p1, p2, p0], Color32::RED);
        rester(&self.verts, |x, y| {
            let p = Polygon::new(PlotPoints::Owned(vec![
                PlotPoint::new(x - 0.5, y - 0.5),
                PlotPoint::new(x - 0.5, y + 0.5),
                PlotPoint::new(x + 0.5, y + 0.5),
                PlotPoint::new(x + 0.5, y - 0.5),
            ]))
            .color(Color32::WHITE);
            plot_ui.polygon(p);
        });
    }
}
