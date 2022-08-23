use rand::{thread_rng, Rng};
use std::f64::consts::SQRT_2;

// https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.83.6811&rep=rep1&type=pdf
pub mod boolean {
    use rgeometry::{
        data::{Direction, EndPoint, ILineSegment, Line, LineSegment},
        Orientation,
    };
    use std::cmp::Ordering;

    use super::*;

    #[derive(Debug, PartialEq, Eq)]
    pub struct SimplicalChain {
        simplices: Vec<Simplex>,
    }

    impl SimplicalChain {
        pub fn from_polygon(p: &Polygon<f64>) -> Self {
            let simplices = p
                .iter_boundary_edges()
                .map(|e| Simplex {
                    p0: e.src.clone(),
                    p1: e.dst.clone(),
                })
                .collect::<Vec<_>>();
            Self { simplices }
        }

        pub fn characteristic(&self, q: &Point<f64>) -> f64 {
            for simplex in &self.simplices {
                if simplex.on_non_original_edge(q) {
                    eprintln!("on-non-original-edge");
                    return 1.0;
                }
            }

            let mut sum = 0f64;
            for simplex in &self.simplices {
                let c = simplex.beta(q);
                eprintln!("c={:?}", c);
                sum += c;
            }
            return sum;
        }

        pub fn subdivide(&self, other: &SimplicalChain) -> SimplicalChain {
            let mut simplices = Vec::new();

            for s0 in &self.simplices {
                let l0 = LineSegment::new(EndPoint::Inclusive(s0.p0), EndPoint::Inclusive(s0.p1));
                let mut intersection_points = vec![s0.p0, s0.p1];

                for s1 in &other.simplices {
                    let l1 =
                        LineSegment::new(EndPoint::Inclusive(s1.p0), EndPoint::Inclusive(s1.p1));

                    match l0.intersect(&l1) {
                        Some(ILineSegment::Crossing) => {
                            eprintln!("crossing, l0={:?}, l1={:?}", l0, l1);
                            let line0 = Line::new_through(&s0.p0, &s0.p1);
                            let line1 = Line::new_through(&s1.p0, &s1.p1);

                            if let Some(p) = line0.intersection_point(&line1) {
                                if !intersection_points.contains(&p)
                                    && Ordering::Less == s0.p0.cmp_distance_to(&p, &s0.p1)
                                {
                                    intersection_points.push(p);
                                }
                            }
                        }
                        Some(ILineSegment::Overlap(view)) => {
                            eprintln!("{:?} / {:?} / {:?}", l0, l1, view);
                            let p = **view.min.inner();
                            if !intersection_points.contains(&p) {
                                intersection_points.push(p);
                            }
                            let p = **view.max.inner();
                            if !intersection_points.contains(&p) {
                                intersection_points.push(p);
                            }
                        }
                        _ => (),
                    }
                }

                intersection_points.sort_by(|a, b| s0.p0.cmp_distance_to(a, b));
                for i in 0..(intersection_points.len() - 1) {
                    let p0 = intersection_points[i];
                    let p1 = intersection_points[i + 1];
                    eprintln!("{:?} {:?}", p0, p1);
                    simplices.push(Simplex { p0, p1 });
                }
            }

            SimplicalChain { simplices }
        }
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    struct Simplex {
        // omit original point
        p0: Point<f64>,
        p1: Point<f64>,
    }

    enum SimplexLocation {
        NonOriginalEdge,
        OriginalEdge,
        Inside,
        Outside,
    }

    const ORIGIN: Point<f64> = Point::new([0f64, 0f64]);

    impl Simplex {
        // clockwise = 1
        // counterclockwise = -1
        // colinear = 0
        fn sign(&self) -> Orientation {
            Point::orient_along_direction(&ORIGIN, Direction::Through(&self.p0), &self.p1)
        }

        fn on_non_original_edge(&self, q: &Point<f64>) -> bool {
            use std::cmp::Ordering::*;
            use Orientation::*;

            let dir2 = Point::orient_along_direction(&self.p0, Direction::Through(&self.p1), q);
            match dir2 {
                CoLinear => {
                    // given point might be on non-original edge
                    let cmp_0 = self.p0.cmp_distance_to(q, &self.p1);
                    let cmp_1 = self.p1.cmp_distance_to(q, &self.p0);

                    match (cmp_0, cmp_1) {
                        (Greater, _) => false,
                        (_, Greater) => false,
                        // p is on non-original edge
                        _ => true,
                    }
                }
                _ => false,
            }
        }

        // assume that p is not on non-original edge of the simplex
        fn beta(&self, q: &Point<f64>) -> f64 {
            use std::cmp::Ordering::*;
            use Orientation::*;

            let sign = self.sign();
            eprintln!(
                "characteristic: sign={:?}, q={:?}, self={:?}",
                q, sign, self
            );
            let signval = match sign {
                ClockWise => -1.0,
                CounterClockWise => 1.0,
                CoLinear => 0.0,
            };

            let dir0 = Point::orient_along_direction(&ORIGIN, Direction::Through(&self.p0), q);
            let dir1 = Point::orient_along_direction(&ORIGIN, Direction::Through(&self.p1), q);

            match (dir0, dir1) {
                // Q is on some original edge
                (CoLinear, _) => {
                    if dir0 == CoLinear && ORIGIN.cmp_distance_to(q, &self.p0) != Greater {
                        signval / 2.0
                    } else {
                        0.0
                    }
                }
                (_, CoLinear) => {
                    if dir1 == CoLinear && ORIGIN.cmp_distance_to(q, &self.p1) != Greater {
                        signval / 2.0
                    } else {
                        0.0
                    }
                }
                (CounterClockWise, ClockWise) => {
                    let dir2 =
                        Point::orient_along_direction(&self.p0, Direction::Through(&self.p1), q);
                    if dir2 == CounterClockWise {
                        signval
                    } else {
                        0.0
                    }
                }
                _ => 0.0,
            }
        }
    }

    struct Polygon2 {
        rings: Vec<Point<f64>>,
    }

    #[cfg(test)]
    mod test {
        use super::*;

        #[test]
        fn simplex_sign() {
            let ccw = Simplex {
                p0: Point::new([1.0, 0.0]),
                p1: Point::new([1.0, 1.0]),
            };
            assert_eq!(ccw.sign(), Orientation::CounterClockWise);

            let cw = Simplex {
                p0: Point::new([1.0, 1.0]),
                p1: Point::new([1.0, 0.0]),
            };
            assert_eq!(cw.sign(), Orientation::ClockWise);

            let cl = Simplex {
                p0: Point::new([1.0, 0.0]),
                p1: Point::new([2.0, 0.0]),
            };
            assert_eq!(cl.sign(), Orientation::CoLinear);
        }

        #[test]
        fn characteristic0() {
            let p = Polygon::new(vec![
                Point::new([1.0, 1.0]),
                Point::new([4.0, 1.0]),
                Point::new([1.0, 4.0]),
            ])
            .unwrap();

            let s = SimplicalChain::from_polygon(&p);
            eprintln!("{:?}", s);

            // on non-original edge, middle
            assert_eq!(s.characteristic(&Point::new([2.5, 2.5])), 1.0);

            // on non-original edge, src
            assert_eq!(s.characteristic(&Point::new([4.0, 1.0])), 1.0);

            // on non-original edge, dst
            assert_eq!(s.characteristic(&Point::new([1.0, 4.0])), 1.0);

            // inside of open region
            assert_eq!(s.characteristic(&Point::new([2.0, 2.0])), 1.0);

            // outside
            assert_eq!(s.characteristic(&Point::new([0.5, 0.5])), 0.0);
        }

        fn lerp_f64(v0: f64, v1: f64, t: f64) -> f64 {
            v0 + (v1 - v0) * t
        }

        fn lerp(p0: &Point<f64>, p1: &Point<f64>, t: f64) -> Point<f64> {
            Point::new([
                lerp_f64(p0.array[0], p1.array[0], t),
                lerp_f64(p0.array[1], p1.array[1], t),
            ])
        }

        fn points_along(
            src: &Point<f64>,
            dst: &Point<f64>,
            subdivide: usize,
            out: &mut Vec<Point<f64>>,
        ) {
            for i in 0..subdivide {
                let t = ((i + 1) as f64) / (subdivide as f64);
                out.push(lerp(src, dst, t));
            }
        }

        fn polygon_cube_subdivide(pos: Point<f64>, extent: f64, subdivide: usize) -> Polygon<f64> {
            let mut points = Vec::new();

            let p0 = pos + Vector([-extent, -extent]);
            let p1 = pos + Vector([extent, -extent]);
            let p2 = pos + Vector([extent, extent]);
            let p3 = pos + Vector([-extent, extent]);

            points_along(&p3, &p0, subdivide, &mut points);
            points_along(&p0, &p1, subdivide, &mut points);
            points_along(&p1, &p2, subdivide, &mut points);
            points_along(&p2, &p3, subdivide, &mut points);

            Polygon::new(points).expect("polygon_cube_subdivide")
        }

        fn polygon_cube(pos: Point<f64>, extent: f64) -> Polygon<f64> {
            let p0 = pos + Vector([-extent, -extent]);
            let p1 = pos + Vector([extent, -extent]);
            let p2 = pos + Vector([extent, extent]);
            let p3 = pos + Vector([-extent, extent]);

            Polygon::new(vec![p0, p1, p2, p3]).expect("polygon_cube")
        }

        #[test]
        fn test_polygon_cube() {
            let p = Point::new([0.0, 0.0]);
            let extent = 1.0;

            let c0 = polygon_cube(p, extent);
            let c1 = polygon_cube_subdivide(p, extent, 1);

            let p0 = c0.iter().collect::<Vec<_>>();
            let p1 = c1.iter().collect::<Vec<_>>();

            assert_eq!(p0, p1);
        }

        #[test]
        fn subdivide_nonoverlapping() {
            let p0 = polygon_cube(Point::new([0.0, 0.0]), 1.0);
            let s0 = SimplicalChain::from_polygon(&p0);

            let p1 = polygon_cube(Point::new([0.0, 0.0]), 2.0);
            let s1 = SimplicalChain::from_polygon(&p1);

            let subdivide = s0.subdivide(&s1);
            assert_eq!(subdivide, s0);
        }

        #[test]
        fn subdivide_exact_edge() {
            let p0 = polygon_cube(Point::new([0.0, 0.0]), 1.0);
            let s0 = SimplicalChain::from_polygon(&p0);

            let p1 = polygon_cube(Point::new([2.0, 0.0]), 1.0);
            let s1 = SimplicalChain::from_polygon(&p1);

            let subdivide = s0.subdivide(&s1);
            assert_eq!(subdivide, s0);
        }

        // TODO: fix test with signed zero
        #[test]
        #[ignore]
        fn subdivide_exact_overlap_zerosign() {
            let p0 = polygon_cube(Point::new([1.0, 0.0]), 1.0);
            let s0 = SimplicalChain::from_polygon(&p0);

            let p1 = polygon_cube(Point::new([2.0, 0.0]), 1.0);
            let s1 = SimplicalChain::from_polygon(&p1);

            let subdivide = s0.subdivide(&s1);
            assert_eq!(subdivide.simplices.len(), 6);
        }

        #[test]
        fn subdivide_exact_overlap() {
            let p0 = polygon_cube(Point::new([1.0, 0.0]), 1.0);
            let s0 = SimplicalChain::from_polygon(&p0);

            let p1 = polygon_cube(Point::new([2.0, 0.0]), 1.0);
            let s1 = SimplicalChain::from_polygon(&p1);

            let subdivide = s0.subdivide(&s1);
            assert_eq!(subdivide.simplices.len(), 6);
        }

        #[test]
        fn subdivide_exact_overlap_partialedge() {
            let p0 = polygon_cube(Point::new([0.0, 0.0]), 1.0);
            let s0 = SimplicalChain::from_polygon(&p0);

            let p1 = polygon_cube(Point::new([2.0, 0.5]), 1.0);
            let s1 = SimplicalChain::from_polygon(&p1);

            let subdivide = s0.subdivide(&s1);
            assert_eq!(subdivide.simplices.len(), 5);
        }

        #[test]
        fn subdivide_exact_overlap_multisplit_partial_small() {
            let p0 = polygon_cube(Point::new([0.0, 0.0]), 3.0);
            let s0 = SimplicalChain::from_polygon(&p0);

            let p1 = polygon_cube_subdivide(Point::new([4.0, 0.0]), 1.0, 3);
            let s1 = SimplicalChain::from_polygon(&p1);

            let subdivide = s0.subdivide(&s1);
            assert_eq!(subdivide.simplices.len(), 9);
        }

        #[test]
        fn subdivide_exact_overlap_multisplit_partial_large() {
            let p0 = polygon_cube(Point::new([0.0, 0.0]), 3.0);
            let s0 = SimplicalChain::from_polygon(&p0);

            let p1 = polygon_cube_subdivide(Point::new([7.0, 0.0]), 4.0, 3);
            let s1 = SimplicalChain::from_polygon(&p1);

            let subdivide = s0.subdivide(&s1);
            eprintln!("{:?}", subdivide);
            assert_eq!(subdivide.simplices.len(), 6);
        }

        #[test]
        fn subdivide_exact_overlap_multisplit_full() {
            let p0 = polygon_cube(Point::new([0.0, 0.0]), 3.0);
            let s0 = SimplicalChain::from_polygon(&p0);

            let p1 = polygon_cube_subdivide(Point::new([6.0, 0.0]), 3.0, 3);
            let s1 = SimplicalChain::from_polygon(&p1);

            let subdivide = s0.subdivide(&s1);
            assert_eq!(subdivide.simplices.len(), 6);
        }

        #[ignore]
        #[test]
        fn test_crossing() {
            let p0 = Point::new([-3.0, -3.0]);
            let p1 = Point::new([3.0, -3.0]);

            let p2 = Point::new([3.0, -4.0]);
            let p3 = Point::new([3.0, -1.333333333333333]);

            let l0 = LineSegment::new(EndPoint::Inclusive(p0), EndPoint::Inclusive(p1));
            let l1 = LineSegment::new(EndPoint::Inclusive(p2), EndPoint::Inclusive(p3));

            assert_eq!(Some(ILineSegment::Crossing), l0.intersect(&l1));

            let l0 = Line::new_through(&p0, &p1);
            let l1 = Line::new_through(&p2, &p3);

            let res = l0.intersection_point(&l1);
            if let Some(p) = res {
                // invariant: distance?
                assert_ne!(p0.cmp_distance_to(&p, &p1), std::cmp::Ordering::Greater);

                // invariant: colinear?
                assert_eq!(
                    Orientation::CoLinear,
                    Point::orient_along_direction(&p0, Direction::Through(&p1), &p)
                );
            } else {
                todo!();
            }
        }
    }
}

use eframe::{
    egui::{
        self,
        plot::{self, Legend, Plot, PlotPoint, PlotPoints, PlotUi},
        Key,
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

fn pt_egui(p: &Point<f64>) -> PlotPoint {
    let [x, y] = p.array;
    PlotPoint::new(x, y)
}

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

#[derive(Clone, Copy, Debug)]
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

    fn ui(self, plot_ui: &mut PlotUi, opts: &EguiRectOpts) {
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
                plot::Points::new(PlotPoints::Owned(self.points)).color(Color32::LIGHT_BLUE),
            );
            plot_ui.points(
                plot::Points::new(PlotPoints::Owned(self.p_points)).color(Color32::LIGHT_RED),
            );
            plot_ui
                .points(plot::Points::new(PlotPoints::Owned(self.pe_points)).color(Color32::GREEN));
        }

        if opts.render_rectnet {
            self.net.ui(plot_ui);
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

struct MyApp {
    pause: bool,
    t: f64,

    count: usize,
    view: f64,

    opts: EguiRectOpts,

    rects: Vec<Rect>,
}

fn gen_rects(view: f64, count: usize) -> Vec<Rect> {
    if false {
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
    } else {
        vec![Rect::new(10.0, 3.0)]
    }
}

impl Default for MyApp {
    fn default() -> Self {
        let view = 30f64;
        let count = 100;

        let rects = gen_rects(view, count);
        let opts = EguiRectOpts {
            render_rectnet: true,
            ..EguiRectOpts::default()
        };

        Self {
            pause: false,
            t: 0.0,

            count,
            view,

            opts,

            rects,
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
        if ctx.input().key_released(Key::G) {
            self.rects = gen_rects(self.view, self.count);
        }
        if ctx.input().key_pressed(Key::S) {
            self.t += 1.0 / 60.0;
        }

        if !self.pause {
            // self.t += ctx.input().predicted_dt as f64;
            ctx.request_repaint();
        }

        let view = self.view;

        let t = self.t;

        let sw = Stopwatch::start_new();
        let rects = self
            .rects
            .iter()
            .map(|r| EguiRect::new(r.rot(t)))
            .collect::<Vec<_>>();

        let elapsed_build = sw.elapsed_ms();
        let sw = Stopwatch::start_new();

        let area = GridArea::extent(view, view);
        let mut net = GridNet::new(area);

        for rect in &rects {
            net.apply(&rect.net);
        }
        let elapsed_net = sw.elapsed_ms();

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.opts.render_polygon, "polygon");
                ui.checkbox(&mut self.opts.render_bb, "bb");
                ui.checkbox(&mut self.opts.render_points, "points");
                ui.checkbox(&mut self.opts.render_rectnet, "rectnet");
                ui.checkbox(&mut self.opts.render_net, "net");
                ui.label(format!(
                    "elapsed build={}ms, net={}ms",
                    elapsed_build, elapsed_net
                ));
            });

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
                for r in rects {
                    r.ui(plot_ui, &self.opts);
                }

                if self.opts.render_net {
                    net.ui(plot_ui);
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
