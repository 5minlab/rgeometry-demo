use crate::delaunay::SubIdx;

use super::{plot_line, points_uniform, Demo};
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, PlotPoint, PlotPoints, Polygon},
    Key, Ui,
};
use rgeometry::{
    data::{Direction, Point},
    Orientation,
};

pub struct DemoRester {
    view: f64,

    verts: [Point<f64>; 3],
}

fn gen_tri(view: f64) -> [Point<f64>; 3] {
    let verts = points_uniform(view, 3);

    let p0 = verts[0];
    let p1 = verts[1];
    let p2 = verts[2];

    match Point::orient_along_direction(&p0, Direction::Through(&p1), &p2) {
        Orientation::CounterClockWise => [p0, p2, p1],
        _ => [p0, p1, p2],
    }
}

impl DemoRester {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        Self {
            view,
            verts: gen_tri(view),
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
            self.verts = gen_tri(self.view);
        }
    }

    fn plot_ui(&self, plot_ui: &mut plot::PlotUi) {
        let [ref p0, ref p1, ref p2] = self.verts;
        plot_line(plot_ui, &[p0, p1, p2, p0], Color32::RED);

        let idx_topmost = if p0.array[1] < p1.array[1] {
            if p1.array[1] < p2.array[1] {
                2
            } else {
                1
            }
        } else {
            if p0.array[1] < p2.array[1] {
                2
            } else {
                0
            }
        };

        // ccw
        let mut left = SubIdx(idx_topmost).ccw();
        let mut right = SubIdx(idx_topmost).cw();

        fn midcoord(v: f64) -> f64 {
            (v + 0.5).floor() - 0.5
        }

        fn roundx(v: f64) -> f64 {
            (v + 0.5).ceil() - 0.5
        }

        fn interp(p0: &Point<f64>, p1: &Point<f64>, y: f64) -> f64 {
            let [p0x, p0y] = p0.array;
            let [p1x, p1y] = p1.array;
            if p0y == p1y {
                return p0x;
            }
            let t = (y - p0y) / (p1y - p0y);
            return (p1x - p0x) * t + p0x;
        }

        let mut y = midcoord(self.verts[idx_topmost].array[1]);
        loop {
            if y <= self.verts[left.0].array[1] {
                if self.verts[left.ccw().0].array[1] < y {
                    left = left.ccw();
                } else {
                    break;
                }
            }
            if y <= self.verts[right.0].array[1] {
                if self.verts[right.cw().0].array[1] < y {
                    right = right.cw();
                } else {
                    break;
                }
            }

            let b = self.verts[left.0].array[1].max(self.verts[right.0].array[1]);

            while y > b {
                let xl = interp(&self.verts[right.ccw().0], &self.verts[right.0], y);
                let xr = interp(&self.verts[left.cw().0], &self.verts[left.0], y);

                // round to next 0.5
                let mut x0 = roundx(xl);
                while x0 < xr {
                    let p = Polygon::new(PlotPoints::Owned(vec![
                        PlotPoint::new(x0 - 0.5, y - 0.5),
                        PlotPoint::new(x0 - 0.5, y + 0.5),
                        PlotPoint::new(x0 + 0.5, y + 0.5),
                        PlotPoint::new(x0 + 0.5, y - 0.5),
                    ]))
                    .color(Color32::WHITE);
                    plot_ui.polygon(p);
                    x0 += 1.0;
                }

                y -= 1.0;
            }
        }
    }
}
