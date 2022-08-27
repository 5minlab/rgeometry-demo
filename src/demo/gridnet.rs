use super::{gen_rects, p_rg_to_egui, pt_egui, Demo, Rect};
use eframe::egui::{
    self,
    epaint::Color32,
    plot::{self, *},
    Key,
};
use rgeometry::data::{DirectedEdge, Point, PointLocation, Polygon, Vector};
use rgeometry::Intersects;
use std::f64::consts::SQRT_2;

pub fn contains(p: &Polygon<f64>, p0: &Point<f64>) -> bool {
    if true {
        let p1 = p0 + &Vector([1000., 0.]);
        intersects(p, p0, &p1) % 2 == 0
    } else {
        match p.locate(p0) {
            PointLocation::Inside => true,
            _ => false,
        }
    }
}

pub fn intersects(p: &Polygon<f64>, p0: &Point<f64>, p1: &Point<f64>) -> usize {
    let ray = DirectedEdge { src: p0, dst: p1 };
    let mut intersections = 0;
    for edge in p.iter_boundary_edges() {
        if let Some(rgeometry::data::ILineSegment::Crossing) = ray.intersect(edge) {
            intersections += 1;
        }
    }
    intersections
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

/// for GridPos (x, y),
///  - covering area: ([x, x+1), [y, y+1)]
///  - center: (x + 0.5, y + 0.5)
///  - world to grid idx, with unit grid size: world_x.floor()
#[derive(Clone, Copy, Debug)]
struct GridPos {
    x: i32,
    y: i32,
}

impl GridPos {
    const fn new(x: i32, y: i32) -> Self {
        Self { x, y }
    }
    fn center(&self) -> Point<f64> {
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
                if intersects(p, &center, &t) == 0 {
                    continue;
                }

                let neighbor = g + *dir;
                if let Some(mask) = self.at_mut(&g) {
                    *mask |= 1u8 << idx;
                }
                if let Some(mask) = self.at_mut(&neighbor) {
                    *mask |= 1u8 << (idx + 4);
                }
            }
        }
    }

    fn apply(&mut self, other: &Self) {
        if let Some(area) = self.area.intersects(&other.area) {
            let stripe = (area.x_max - area.x_min) as usize;

            for y in area.y_min..area.y_max {
                let g = GridPos { x: area.x_min, y };
                let idx_dst = self.idx(&g).unwrap();
                let idx_src = other.idx(&g).unwrap();

                let slice_dst = &mut self.masks[idx_dst..(idx_dst + stripe)];
                let slice_src = &other.masks[idx_src..(idx_src + stripe)];

                for i in 0..stripe {
                    slice_dst[i] |= slice_src[i];
                }
            }
        }
    }

    fn idx(&self, p: &GridPos) -> Option<usize> {
        let a = &self.area;
        if p.x < a.x_min || p.x >= a.x_max || p.y < a.y_min || p.y >= a.y_max {
            None
        } else {
            Some((p.x - a.x_min + (p.y - a.y_min) * (a.x_max - a.x_min)) as usize)
        }
    }

    fn at_mut(&mut self, p: &GridPos) -> Option<&mut u8> {
        self.idx(p).map(|idx| &mut self.masks[idx])
    }

    fn at(&self, p: &GridPos) -> Option<u8> {
        self.idx(p).map(|idx| self.masks[idx])
    }

    fn ui(&self, plot_ui: &mut PlotUi) {
        for g in self.area.iter() {
            if let Some(mask) = self.at(&g) {
                for i in 0..8u8 {
                    let blocked = mask & (1 << i) != 0;
                    if !blocked {
                        continue;
                    }
                    let center = g.center();
                    let p0 = center + &HEX_VECS[i as usize];
                    let p1 = center + &HEX_VECS[(i + 1) as usize];

                    plot_ui.line(
                        plot::Line::new(PlotPoints::Owned(vec![pt_egui(&p0), pt_egui(&p1)]))
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
    }
}

#[derive(Clone, Copy, Debug)]
struct GridArea {
    pub x_min: i32,
    pub x_max: i32,
    pub y_min: i32,
    pub y_max: i32,
}

impl GridArea {
    fn extent(x: f64, y: f64) -> Self {
        Self {
            x_min: (-x).floor() as i32,
            x_max: x.ceil() as i32,
            y_min: (-y).floor() as i32,
            y_max: y.ceil() as i32,
        }
    }

    fn new(t: &(Point<f64>, Point<f64>)) -> Self {
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

    fn intersects(&self, other: &Self) -> Option<Self> {
        let x_min = self.x_min.max(other.x_min);
        let x_max = self.x_max.min(other.x_max);
        let y_min = self.y_min.max(other.y_min);
        let y_max = self.y_max.min(other.y_max);
        if x_min >= x_max || y_min >= y_max {
            None
        } else {
            Some(Self {
                x_min,
                x_max,
                y_min,
                y_max,
            })
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

#[derive(Default)]
struct EguiRectOpts {
    render_points: bool,
    render_bb: bool,
    render_polygon: bool,
    render_rectnet: bool,
    render_net: bool,
}

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

struct EguiRect {
    rect: Rect,

    points: Vec<PlotPoint>,
    p_points: Vec<PlotPoint>,
    pe_points: Vec<PlotPoint>,

    net: GridNet,
}

impl EguiRect {
    fn new(rect: Rect) -> Self {
        let mut points = Vec::new();
        let mut p_points = Vec::new();
        let mut pe_points = Vec::new();

        let p = rect.polygon();

        let r_extend = rect.add_extents(SQRT_2 / 2.0, SQRT_2 / 2.0);
        let r_extend_p = r_extend.polygon();

        let bb = p.bounding_box();
        let area = GridArea::new(&bb);

        for g in area.iter() {
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
        let area = GridArea::new(&bb);

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

        Self {
            rect,

            points,
            p_points,
            pe_points,

            net,
        }
    }

    fn ui(&self, plot_ui: &mut PlotUi, opts: &EguiRectOpts) {
        let p = self.rect.polygon();
        let pe = p_rg_to_egui(&p);

        let bb = p.bounding_box();
        let bb_r = Rect::from_bb(&bb);
        let bb_p = bb_r.polygon();
        let bb_pe = p_rg_to_egui(&bb_p);

        if opts.render_polygon {
            plot_ui.polygon(pe.color(Color32::RED));
        }
        if opts.render_bb {
            plot_ui.polygon(bb_pe.color(Color32::BLUE));
        }

        if opts.render_points {
            plot_ui.points(
                plot::Points::new(PlotPoints::Owned(self.points.clone()))
                    .color(Color32::LIGHT_BLUE),
            );
            plot_ui.points(
                plot::Points::new(PlotPoints::Owned(self.p_points.clone()))
                    .color(Color32::LIGHT_RED),
            );
            plot_ui.points(
                plot::Points::new(PlotPoints::Owned(self.pe_points.clone())).color(Color32::GREEN),
            );
        }

        if opts.render_rectnet {
            self.net.ui(plot_ui);
        }
    }
}

pub struct DemoGridNet {
    view: f64,
    count: usize,

    opts: EguiRectOpts,
    rects: Vec<Rect>,
    egui_rects: Vec<EguiRect>,
}

impl DemoGridNet {
    pub fn new(view: f64) -> Self {
        let count = 100;

        let rects = gen_rects(view, count);
        let opts = EguiRectOpts {
            ..EguiRectOpts::default()
        };

        let t = 0.0f64;
        let egui_rects = rects
            .iter()
            .map(|r| EguiRect::new(r.rot(t)))
            .collect::<Vec<_>>();

        Self {
            view,
            count,
            opts,
            rects,
            egui_rects,
        }
    }
}

impl Demo for DemoGridNet {
    fn name(&self) -> &'static str {
        "gridnet"
    }
    fn ui(&mut self, t: f64, ctx: &egui::Context, ui: &mut egui::Ui) {
        if ctx.input().key_released(Key::G) {
            self.rects = gen_rects(self.view, self.count);
        }

        self.egui_rects = self
            .rects
            .iter()
            .map(|r| EguiRect::new(r.rot(t)))
            .collect::<Vec<_>>();

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.opts.render_polygon, "polygon");
            ui.checkbox(&mut self.opts.render_bb, "bb");
            ui.checkbox(&mut self.opts.render_points, "points");
            ui.checkbox(&mut self.opts.render_rectnet, "rectnet");
            ui.checkbox(&mut self.opts.render_net, "net");
            /*
            ui.label(format!(
                "elapsed build={}ms, net={}ms",
                elapsed_build, elapsed_net
            ));
            */
        });
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        for r in &self.egui_rects {
            r.ui(plot_ui, &self.opts);
        }

        if self.opts.render_net {
            let area = GridArea::extent(self.view, self.view);
            let mut net = GridNet::new(area);

            for rect in &self.egui_rects {
                net.apply(&rect.net);
            }

            net.ui(plot_ui);
        }
    }
}
