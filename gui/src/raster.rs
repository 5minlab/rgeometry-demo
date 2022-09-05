use core::{points_tri, raster::raster};

use super::{plot_line, Demo};
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, PlotPoint, PlotPoints, Polygon},
    Key, Ui,
};
use rgeometry::data::Point;

pub struct DemoRaster {
    view: f64,
    grid_size: usize,

    verts: [Point<f64>; 3],
}

impl DemoRaster {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let mut rng = rand::thread_rng();
        Self {
            view,
            grid_size: 5,
            verts: points_tri(&mut rng, view),
        }
    }
}

impl Demo for DemoRaster {
    fn name(&self) -> &'static str {
        "Raster"
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
            ui.separator();
            for size in [1, 5, 10, 20] {
                ui.radio_value(&mut self.grid_size, size, size.to_string());
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
        raster(self.grid_size, &self.verts, |x, y| {
            let m = self.grid_size as f64;
            let p = Polygon::new(PlotPoints::Owned(vec![
                PlotPoint::new(m * x, m * y),
                PlotPoint::new(m * x, m * y + m),
                PlotPoint::new(m * x + m, m * y + m),
                PlotPoint::new(m * x + m, m * y),
            ]))
            .color(Color32::WHITE);
            plot_ui.polygon(p);
        });
    }
}
