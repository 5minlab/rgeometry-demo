use std::f64::consts::SQRT_2;

use eframe::{
    egui::{
        self,
        plot::{self, Legend, Plot, PlotPoint, PlotPoints},
        Key, Label,
    },
    epaint::Color32,
};
use rgeometry::data::{DirectedEdge, Point, PointLocation, Polygon, Vector};
use rgeometry::Intersects;
use stopwatch::Stopwatch;

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

fn pt_egui(p: &Point<f64, 2>) -> PlotPoint {
    let [x, y] = p.array;
    PlotPoint::new(x, y)
}

pub fn contains(p: &Polygon<f64>, p0: &Point<f64, 2>) -> bool {
    if true {
        let p1 = p0 + &Vector([1000., 0.]);
        intersects(p, p0, &p1)
    } else {
        match p.locate(p0) {
            PointLocation::Inside => true,
            _ => false,
        }
    }
}

pub fn intersects(p: &Polygon<f64>, p0: &Point<f64, 2>, p1: &Point<f64, 2>) -> bool {
    let ray = DirectedEdge { src: p0, dst: p1 };
    let mut intersections = 0;
    for edge in p.iter_boundary_edges() {
        if let Some(rgeometry::data::ILineSegment::Crossing) = ray.intersect(edge) {
            // Only count crossing that aren't leaning to the right.
            intersections += 1;
        }
    }
    // handle polygon without holes
    intersections % 2 == 1
}

#[derive(Clone, Copy)]
struct Rect {
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

    fn from_bb(t: &(Point<f64, 2>, Point<f64, 2>)) -> Self {
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
    t: f64,
}

impl Default for MyApp {
    fn default() -> Self {
        Self {
            pause: false,
            t: 0.0,
        }
    }
}

/// for GridPos (x, y),
///  - covering area: ([x, x+1), [y, y+1)]
///  - center: (x + 0.5, y + 0.5)
///  - world to grid idx, with unit grid size: world_x.floor()
#[derive(Clone, Copy)]
struct GridPos {
    x: i32,
    y: i32,
}

impl GridPos {
    const fn new(x: i32, y: i32) -> Self {
        Self { x, y }
    }
    fn center(&self) -> Point<f64, 2> {
        Point::new([self.x as f64 + 0.5, self.y as f64 + 0.5])
    }
}

impl std::ops::Add<GridPos> for GridPos {
    type Output = GridPos;
    fn add(self, other: GridPos) -> Self::Output {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

#[allow(unused)]
const VECS_HALF: [GridPos; 4] = [
    GridPos::new(1, 0),
    GridPos::new(1, 1),
    GridPos::new(0, 1),
    GridPos::new(-1, 1),
];

#[allow(unused)]
const VECS: [GridPos; 8] = [
    GridPos::new(1, 0),
    GridPos::new(1, 1),
    GridPos::new(0, 1),
    GridPos::new(-1, 1),
    GridPos::new(-1, 0),
    GridPos::new(-1, -1),
    GridPos::new(0, -1),
    GridPos::new(1, -1),
];

struct GridNet {
    area: GridArea,

    masks: Box<[u8]>,
}

impl GridNet {
    fn new(area: GridArea) -> Self {
        let w = area.x_max - area.x_min;
        let h = area.y_max - area.y_min;
        let l = (w * h) as usize;
        let mut v = Vec::with_capacity(l);
        v.resize(l, 0u8);
        Self {
            area,
            masks: v.into_boxed_slice(),
        }
    }

    fn update(&mut self, p: &Polygon<f64>) {
        for g in self.area.iter() {
            // cell center
            let center = g.center();

            for (idx, dir) in VECS_HALF.iter().enumerate() {
                let v = Vector([dir.x as f64, dir.y as f64]);

                let t = center + v;
                if !intersects(p, &center, &t) {
                    continue;
                }

                let neighbor = g + *dir;
                if let Some(mask) = self.at(&g) {
                    *mask |= 1u8 << idx;
                }
                if let Some(mask) = self.at(&neighbor) {
                    *mask |= 1u8 << (idx + 4);
                }
            }
        }
    }

    fn at(&mut self, p: &GridPos) -> Option<&mut u8> {
        let a = &self.area;
        if p.x < a.x_min || p.x >= a.x_max || p.y <= a.y_min || p.y >= a.y_max {
            return None;
        }
        let idx = p.x - a.x_min + (p.y - a.y_min) * (a.x_max - a.x_min);
        Some(&mut self.masks[idx as usize])
    }
}

#[derive(Clone, Copy)]
struct GridArea {
    pub x_min: i32,
    pub x_max: i32,
    pub y_min: i32,
    pub y_max: i32,
}

impl GridArea {
    fn new(t: &(Point<f64, 2>, Point<f64, 2>)) -> Self {
        let x_min = t.0[0].floor() as i32;
        let x_max = t.1[0].ceil() as i32;
        let y_min = t.0[1].floor() as i32;
        let y_max = t.1[1].ceil() as i32;
        Self {
            x_min,
            x_max,
            y_min,
            y_max,
        }
    }

    fn iter(&self) -> GridIterator {
        GridIterator {
            area: *self,
            x: self.x_min - 1,
            y: self.y_min,
        }
    }
}

struct GridIterator {
    area: GridArea,

    x: i32,
    y: i32,
}

impl Iterator for GridIterator {
    type Item = GridPos;

    fn next(&mut self) -> Option<Self::Item> {
        self.x += 1;
        if self.x == self.area.x_max {
            self.x = self.area.x_min;
            self.y += 1;
        }

        if self.y == self.area.y_max {
            None
        } else {
            Some(GridPos {
                x: self.x,
                y: self.y,
            })
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
            self.t += ctx.input().predicted_dt as f64;
            ctx.request_repaint();
        }

        let view = 15f64;

        let w = 10f64;
        let h = 2f64;

        let sw = Stopwatch::start_new();

        let r = Rect::new(w, h).rot(self.t);
        let p = r.polygon();
        let pe = p_rg_to_egui(&p);

        let r_extend = r.add_extents(SQRT_2 / 2.0, SQRT_2 / 2.0);
        let r_extend_p = r_extend.polygon();

        let bb = p.bounding_box();
        let bb_r = Rect::from_bb(&bb);
        let bb_p = bb_r.polygon();
        let bb_pe = p_rg_to_egui(&bb_p);

        let elapsed_poly = sw.elapsed_ms();
        let sw = Stopwatch::start_new();

        let mut points = Vec::new();
        let mut p_points = Vec::new();
        let mut pe_points = Vec::new();

        let area = GridArea::new(&bb);

        let mut count = 0;
        for g in area.iter() {
            count += 1;
            // cell center
            let center = g.center();

            let bucket = match contains(&r_extend_p, &center) {
                false => &mut points,
                true => match contains(&p, &center) {
                    true => &mut p_points,
                    false => &mut pe_points,
                },
            };
            bucket.push(pt_egui(&center));
        }
        eprintln!("count={count}");

        let elapsed_points = sw.elapsed_ms();
        let sw = Stopwatch::start_new();

        // extended bounding box
        let area2 = {
            let mut a = area;
            a.x_min -= 1;
            a.y_min -= 1;
            a.x_max += 1;
            a.y_max += 1;
            a
        };

        let mut net = GridNet::new(area2);
        net.update(&p);

        let elapsed_net = sw.elapsed_ms();

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

            ui.label(format!(
                "timings: poly={}ms, points={}ms, net={}ms",
                elapsed_poly, elapsed_points, elapsed_net
            ));
            plot.show(ui, |plot_ui| {
                plot_ui.polygon(pe.color(Color32::RED));

                plot_ui.polygon(bb_pe.color(Color32::BLUE));

                const HEX_OFFSET: f64 = 0.2886751345948128 * 0.5;
                const C_H: f64 = 0.5;
                const HEX_VECS: [Vector<f64, 2>; 9] = [
                    Vector([C_H, -HEX_OFFSET]),
                    Vector([C_H, HEX_OFFSET]),
                    Vector([HEX_OFFSET, C_H]),
                    Vector([-HEX_OFFSET, C_H]),
                    Vector([-C_H, HEX_OFFSET]),
                    Vector([-C_H, -HEX_OFFSET]),
                    Vector([-HEX_OFFSET, -C_H]),
                    Vector([HEX_OFFSET, -C_H]),
                    Vector([C_H, -HEX_OFFSET]),
                ];

                for g in area2.iter() {
                    if let Some(mask) = net.at(&g) {
                        for i in 0..8u8 {
                            let blocked = *mask & (1 << i) != 0;
                            if !blocked {
                                continue;
                            }
                            let center = g.center();
                            let p0 = center + &HEX_VECS[i as usize];
                            let p1 = center + &HEX_VECS[(i + 1) as usize];

                            plot_ui.line(
                                plot::Line::new(PlotPoints::Owned(vec![
                                    pt_egui(&p0),
                                    pt_egui(&p1),
                                ]))
                                .color(Color32::BROWN),
                            );

                            /*
                            let dir = &VECS[i as usize];
                            let v = Vector([dir.x as f64, dir.y as f64]);
                            let t = center + v;
                            line_lists.push(vec![pt_egui(&center), pt_egui(&t)]);
                            */
                        }
                    }
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
