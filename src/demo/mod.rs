use crate::delaunay::TriangularNetwork;
use rand::{thread_rng, Rng};

use eframe::{
    egui::{
        self,
        plot::{self, Legend, Plot, PlotPoint, PlotPoints, PlotUi},
        Key,
    },
    epaint::Color32,
};
use rgeometry::data::{Point, Polygon, Vector};

mod boolean;
mod boolean_tri;
mod cdt;
mod delaunay;
mod gridnet;

pub fn plot_line(plot_ui: &mut PlotUi, points: &[&Point<f64>], color: Color32) {
    let e_points = points.iter().map(|p| pt_egui(p)).collect();
    plot_ui.line(plot::Line::new(PlotPoints::Owned(e_points)).color(color));
}

pub fn pt_mean(points: &[&Point<f64>]) -> Point<f64> {
    let mut x = 0.0;
    let mut y = 0.0;
    for p in points {
        x += p.array[0];
        y += p.array[1];
    }
    let l = points.len() as f64;
    Point::new([x / l, y / l])
}

fn plot_net(net: &TriangularNetwork, plot_ui: &mut PlotUi, render_supertri: bool) {
    for (_t_idx, t) in net.triangles.iter().enumerate() {
        let [v0, v1, v2] = t.vertices;
        let p0 = net.vert(v0);
        let p1 = net.vert(v1);
        let p2 = net.vert(v2);
        if !render_supertri && t.is_super() {
            continue;
        }

        plot_line(plot_ui, &[p0, p1, p2, p0], Color32::GREEN);

        if false {
            let center = pt_mean(&[p0, p1, p2]);
            let mut label = format!("t={}\n(", _t_idx);
            for i in 0..3 {
                let n = t.neighbors[i];
                let prefix = if i == 0 { "" } else { ", " };
                match n {
                    Some(n) => label += &format!("{}{}", prefix, n.0),
                    None => label += &format!("{}_", prefix),
                }
            }
            label += ")";

            plot_ui.text(plot::Text::new(pt_egui(&center), label));
        }
    }
}

type P = Vector<f64, 2>;

fn rotate(p: P, rot: f64) -> P {
    Vector([
        rot.cos() * p.0[0] + rot.sin() * p.0[1],
        rot.sin() * p.0[0] - rot.cos() * p.0[1],
    ])
}

fn p_rg_to_egui(p: &Polygon<f64>) -> plot::Polygon {
    let points: Vec<[f64; 2]> = p
        .iter()
        .map(|p| [p.array[0] as f64, p.array[1] as f64])
        .collect::<Vec<_>>();
    plot::Polygon::new(points)
}

pub fn pt_egui(p: &Point<f64>) -> PlotPoint {
    let [x, y] = p.array;
    PlotPoint::new(x, y)
}

#[derive(Clone, Copy, Debug)]
pub struct Rect {
    pub pos: [f64; 2],
    pub extent: [f64; 2],
    pub rot: f64,
}

impl Rect {
    fn new(extent_x: f64, extent_y: f64) -> Self {
        Self {
            pos: [0f64; 2],
            extent: [extent_x, extent_y],
            rot: 0f64,
        }
    }

    fn from_bb(t: &(Point<f64>, Point<f64>)) -> Self {
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

    fn pos(self, x: f64, y: f64) -> Self {
        Self {
            pos: [x, y],
            ..self
        }
    }

    fn add_extents(self, x: f64, y: f64) -> Self {
        let [x0, y0] = self.extent;
        Self {
            extent: [x0 + x, y0 + y],
            ..self
        }
    }

    fn rot(self, rot: f64) -> Self {
        Self { rot, ..self }
    }

    fn polygon(&self) -> Polygon<f64> {
        let center = Point::new(self.pos);
        let p0 = center + rotate(Vector([-self.extent[0], -self.extent[1]]), self.rot);
        let p1 = center + rotate(Vector([self.extent[0], -self.extent[1]]), self.rot);
        let p2 = center + rotate(Vector([self.extent[0], self.extent[1]]), self.rot);
        let p3 = center + rotate(Vector([-self.extent[0], self.extent[1]]), self.rot);

        Polygon::new(vec![p0, p1, p2, p3]).unwrap()
    }
}

pub fn gen_rects(view: f64, count: usize) -> Vec<Rect> {
    let mut rng = thread_rng();
    let mut rects = Vec::new();
    for _ in 0..count {
        let w = rng.gen_range(3.0..6.0);
        let h = rng.gen_range(0.1..3.0);

        let inner = view - (w + h) - 2.0;
        let x = rng.gen_range(-inner..inner);
        let y = rng.gen_range(-inner..inner);

        rects.push(Rect::new(w, h).pos(x, y));
    }
    rects
}

trait Demo {
    fn name(&self) -> &'static str;
    fn ui(&mut self, t: f64, ctx: &egui::Context, ui: &mut egui::Ui);
    fn plot_ui(&self, plot_ui: &mut PlotUi);
}

pub struct MyApp {
    pause: bool,
    t: f64,

    view: f64,

    selected: &'static str,
    demos: Vec<Box<dyn Demo>>,
}

impl Default for MyApp {
    fn default() -> Self {
        let view = 30f64;

        let demos: Vec<Box<dyn Demo>> = vec![
            Box::new(boolean_tri::DemoBooleanTri::new(view)),
            Box::new(delaunay::DemoDelaunay::new(view)),
            Box::new(cdt::DemoCDT::new(view)),
            Box::new(boolean::DemoBoolean::new(view)),
            Box::new(gridnet::DemoGridNet::new(view)),
        ];
        let selected = demos[0].name();

        Self {
            pause: false,
            t: 0.0,

            view,

            selected,
            demos,
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
        if ctx.input().key_pressed(Key::A) {
            self.t -= 1.0 / 60.0;
        }
        if ctx.input().key_pressed(Key::S) {
            self.t += 1.0 / 60.0;
        }

        if !self.pause {
            self.t += ctx.input().predicted_dt as f64;
            ctx.request_repaint();
        }

        let view = self.view;

        let t = self.t;

        egui::TopBottomPanel::top("top").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("time={:.2}", t));
                ui.separator();
                ui.checkbox(&mut self.pause, "pause");
                ui.separator();
                ui.label("shortcuts: (Space) pause/resume | (R) Reset time | (S) single-step");
            });

            ui.horizontal(|ui| {
                ui.visuals_mut().button_frame = false;
                for d in &self.demos {
                    let name = d.name();
                    if ui.selectable_label(self.selected == name, name).clicked() {
                        self.selected = name;
                    }
                    ui.separator();
                }
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            for demo in &mut self.demos {
                if demo.name() == self.selected {
                    demo.ui(t, ctx, ui);
                }
            }

            let plot = Plot::new(0)
                .legend(Legend::default())
                .data_aspect(1.0)
                /*
                .center_x_axis(true)
                .center_y_axis(true)
                .include_x(-view)
                .include_x(view)
                .include_y(-view)
                .include_y(view)
                */
                .allow_zoom(true)
                .allow_drag(true);

            plot.show(ui, |plot_ui| {
                for demo in &self.demos {
                    if demo.name() == self.selected {
                        demo.plot_ui(plot_ui);
                    }
                }

                let area = view;
                plot_ui.points(
                    plot::Points::new(PlotPoints::Owned(vec![
                        PlotPoint::new(-area, -area),
                        PlotPoint::new(area, -area),
                        PlotPoint::new(-area, area),
                        PlotPoint::new(area, area),
                    ]))
                    .color(Color32::WHITE),
                );
            });
        });
    }
}
