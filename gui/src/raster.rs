use core::raster::raster_bounds;

use core::{points_tri, raster::raster, Rect};

use super::{plot_line, plot_polygon, Demo};
use eframe::egui::{self, epaint::Color32, Key, Ui};
use egui_plot::{self, *};
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
        if ctx.input(|i| i.key_pressed(Key::R)) {
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

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        let m = self.grid_size as f64;

        let [ref p0, ref p1, ref p2] = self.verts;
        plot_line(plot_ui, &[p0, p1, p2, p0], Color32::RED);

        let aabb = raster_bounds(self.grid_size, &self.verts);
        let bb_r = Rect::from_aabb(aabb);
        plot_polygon(
            plot_ui,
            &bb_r
                .polygon(1)
                .map_points(|p| Point::new([p.array[0] * m, p.array[1] * m])),
            Color32::BLUE,
        );

        raster(self.grid_size, &self.verts, |x, y| {
            let p = Polygon::new(PlotPoints::Owned(vec![
                PlotPoint::new(m * x, m * y),
                PlotPoint::new(m * x, m * y + m),
                PlotPoint::new(m * x + m, m * y + m),
                PlotPoint::new(m * x + m, m * y),
            ]))
            .fill_color(Color32::WHITE);
            plot_ui.polygon(p);
        });
    }
}
