use std::f32::consts::SQRT_2;

use eframe::{
    egui::{
        self,
        plot::{self, Legend, Plot, PlotPoint, PlotPoints},
        Key,
    },
    epaint::Color32,
};
use rgeometry::data::{Point, PointLocation, Polygon, Vector};

type P = Vector<f32, 2>;

fn rotate(p: P, rot: f32) -> P {
    Vector([
        rot.cos() * p.0[0] + rot.sin() * p.0[1],
        rot.sin() * p.0[0] - rot.cos() * p.0[1],
    ])
}

fn p_rg_to_egui(p: &Polygon<f32>) -> plot::Polygon {
    let points: Vec<[f64; 2]> = p
        .iter()
        .map(|p| [p.array[0] as f64, p.array[1] as f64])
        .collect::<Vec<_>>();
    plot::Polygon::new(points)
}

#[derive(Clone, Copy)]
struct Rect {
    pub pos: [f32; 2],
    pub extent: [f32; 2],
    pub rot: f32,
}

impl Rect {
    fn new(extent_x: f32, extent_y: f32) -> Self {
        Self {
            pos: [0f32; 2],
            extent: [extent_x, extent_y],
            rot: 0f32,
        }
    }

    fn from_bb(t: &(Point<f32, 2>, Point<f32, 2>)) -> Self {
        let (p0, p1) = t;
        Self {
            pos: [
                (p0.array[0] + p1.array[0]) / 2.0,
                (p0.array[1] + p1.array[1]) / 2.0,
            ],
            extent: [
                (p1.array[0] - p0.array[0]) / 2.0,
                (p1.array[1] - p0.array[1]) / 2.0,
            ],
            rot: 0.0,
        }
    }

    fn pos(self, x: f32, y: f32) -> Self {
        Self {
            pos: [x, y],
            ..self
        }
    }

    fn add_extents(self, x: f32, y: f32) -> Self {
        let [x0, y0] = self.extent;
        Self {
            extent: [x0 + x, y0 + y],
            ..self
        }
    }

    fn rot(self, rot: f32) -> Self {
        Self { rot, ..self }
    }

    fn polygon(&self) -> Polygon<f32> {
        let center = Point::new(self.pos);
        let p0 = center + rotate(Vector([-self.extent[0], -self.extent[1]]), self.rot);
        let p1 = center + rotate(Vector([self.extent[0], -self.extent[1]]), self.rot);
        let p2 = center + rotate(Vector([self.extent[0], self.extent[1]]), self.rot);
        let p3 = center + rotate(Vector([-self.extent[0], self.extent[1]]), self.rot);

        Polygon::new_unchecked(vec![p0, p1, p2, p3])
    }
}

fn main() {
    let options = eframe::NativeOptions::default();
    eframe::run_native(
        "rgeometry playground",
        options,
        Box::new(|_cc| Box::new(MyApp::default())),
    );
}

struct MyApp {
    pause: bool,
    t: f32,
}

impl Default for MyApp {
    fn default() -> Self {
        Self {
            pause: false,
            t: 0.0,
        }
    }
}

struct GridIterator {
    x_min: i32,
    x_max: i32,
    #[allow(unused)]
    y_min: i32,
    y_max: i32,

    x: i32,
    y: i32,
}

impl GridIterator {
    fn new(t: &(Point<f32, 2>, Point<f32, 2>)) -> Self {
        let x_min = t.0[0].floor() as i32;
        let x_max = t.1[0].ceil() as i32;
        let y_min = t.0[1].floor() as i32;
        let y_max = t.1[1].ceil() as i32;
        Self {
            x_min,
            x_max,
            y_min,
            y_max,

            x: x_min - 1,
            y: y_min,
        }
    }
}

impl Iterator for GridIterator {
    type Item = (i32, i32);

    fn next(&mut self) -> Option<Self::Item> {
        self.x += 1;
        if self.x == self.x_max {
            self.x = self.x_min;
            self.y += 1;
        }

        if self.y == self.y_max {
            None
        } else {
            Some((self.x, self.y))
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if ctx.input().key_pressed(Key::Space) {
            self.pause = !self.pause;
        }
        if ctx.input().key_pressed(Key::R) {
            self.t = 0.0;
        }
        if ctx.input().key_pressed(Key::S) {
            self.t += 1.0 / 60.0;
        }

        if !self.pause {
            self.t += ctx.input().predicted_dt;
            ctx.request_repaint();
        }

        let view = 5f32;

        let w = 3f32;
        let h = 0.1f32;

        let r = Rect::new(w, h).rot(self.t);
        let p = r.polygon();
        let pe = p_rg_to_egui(&p);

        let r_extend = r.add_extents(SQRT_2 / 2.0, SQRT_2 / 2.0);
        let r_extend_p = r_extend.polygon();

        let bb = p.bounding_box();
        let bb_r = Rect::from_bb(&bb);
        let bb_p = bb_r.polygon();
        let bb_pe = p_rg_to_egui(&bb_p);

        egui::CentralPanel::default().show(ctx, |ui| {
            let plot = Plot::new(0)
                .legend(Legend::default())
                .data_aspect(1.0)
                .center_x_axis(true)
                .center_y_axis(true)
                .include_x(-view)
                .include_x(view)
                .include_y(-view)
                .include_y(view);

            plot.show(ui, |plot_ui| {
                plot_ui.polygon(pe.color(Color32::RED));

                plot_ui.polygon(bb_pe.color(Color32::BLUE));

                let mut points = Vec::new();
                let mut p_points = Vec::new();
                let mut pe_points = Vec::new();

                let iter = GridIterator::new(&bb);

                for (x, y) in iter {
                    // cell center
                    let x = x as f64 + 0.5;
                    let y = y as f64 + 0.5;
                    let center = Point::new([x as f32, y as f32]);
                    let point = PlotPoint::new(x, y);

                    let p = match r_extend_p.locate(&center) {
                        PointLocation::Outside => &mut points,
                        PointLocation::Inside | PointLocation::OnBoundary => {
                            match p.locate(&center) {
                                PointLocation::Inside | PointLocation::OnBoundary => &mut p_points,
                                PointLocation::Outside => &mut pe_points,
                            }
                        }
                    };
                    p.push(point);
                }

                plot_ui.points(
                    plot::Points::new(PlotPoints::Owned(points)).color(Color32::LIGHT_BLUE),
                );
                plot_ui.points(
                    plot::Points::new(PlotPoints::Owned(p_points)).color(Color32::LIGHT_RED),
                );
                plot_ui
                    .points(plot::Points::new(PlotPoints::Owned(pe_points)).color(Color32::GREEN));
            });
            //
        });
    }
}
