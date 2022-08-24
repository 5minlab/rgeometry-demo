// https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.83.6811&rep=rep1&type=pdf

use rgeometry::{data::*, Intersects, Orientation};
use std::cmp::Ordering;

fn lerp_f64(v0: f64, v1: f64, t: f64) -> f64 {
    v0 + (v1 - v0) * t
}

fn lerp_pf64(p0: &Point<f64>, p1: &Point<f64>, t: f64) -> Point<f64> {
    Point::new([
        lerp_f64(p0.array[0], p1.array[0], t),
        lerp_f64(p0.array[1], p1.array[1], t),
    ])
}

#[derive(Default, Debug, PartialEq, Eq)]
pub struct SimplicalChain {
    pub simplices: Vec<Simplex>,
}

impl SimplicalChain {
    pub fn from_polygon(p: &Polygon<f64>) -> Self {
        let simplices = p
            .iter_boundary_edges()
            .map(|e| Simplex {
                src: *e.src,
                dst: *e.dst,
            })
            .collect::<Vec<_>>();
        Self { simplices }
    }

    pub fn characteristic(&self, q: &Point<f64>) -> f64 {
        for simplex in &self.simplices {
            if simplex.on_non_original_edge(q) {
                // eprintln!("on-non-original-edge");
                return 1.0;
            }
        }

        let mut sum = 0f64;
        for simplex in &self.simplices {
            let c = simplex.beta(q);
            // eprintln!("c={:?}", c);
            sum += c;
        }
        return sum;
    }

    fn subdivide(&self, other: &SimplicalChain) -> SimplicalChain {
        let mut simplices = Vec::new();

        for s0 in &self.simplices {
            let l0 = LineSegment::new(EndPoint::Inclusive(s0.src), EndPoint::Inclusive(s0.dst));
            let mut intersection_points = vec![s0.src, s0.dst];

            for s1 in &other.simplices {
                let l1 = LineSegment::new(EndPoint::Inclusive(s1.src), EndPoint::Inclusive(s1.dst));

                match l0.intersect(&l1) {
                    Some(ILineSegment::Crossing) => {
                        // eprintln!("crossing, l0={:?}, l1={:?}", l0, l1);
                        let line0 = Line::new_through(&s0.src, &s0.dst);
                        let line1 = Line::new_through(&s1.src, &s1.dst);

                        if let Some(p) = line0.intersection_point(&line1) {
                            if !intersection_points.contains(&p)
                                && Ordering::Less == s0.src.cmp_distance_to(&p, &s0.dst)
                                && Ordering::Less == s0.dst.cmp_distance_to(&p, &s0.src)
                            {
                                intersection_points.push(p);
                            }
                        }
                    }
                    Some(ILineSegment::Overlap(view)) => {
                        // eprintln!("{:?} / {:?} / {:?}", l0, l1, view);
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

            if intersection_points.len() > 2 {
                intersection_points.sort_by(|a, b| s0.src.cmp_distance_to(a, b));
            }
            for i in 0..(intersection_points.len() - 1) {
                let p0 = intersection_points[i];
                let p1 = intersection_points[i + 1];
                // eprintln!("{:?} {:?}", p0, p1);
                simplices.push(Simplex { src: p0, dst: p1 });
            }
        }

        SimplicalChain { simplices }
    }

    fn subdivide_chracteristics(
        &self,
        other: &SimplicalChain,
    ) -> (SimplicalChain, Vec<f64>, SimplicalChain, Vec<f64>) {
        let sx0 = self;
        let sx1 = other;

        // subdivide, simplical chain
        let sx0_subdivide = sx0.subdivide(sx1);
        let sx1_subdivide = sx1.subdivide(sx0);

        // calculate edge characteristics

        // f(s_i*, P2)
        let mut ec_sx0 = Vec::with_capacity(sx0_subdivide.simplices.len());
        ec_sx0.resize(sx0_subdivide.simplices.len(), 0f64);

        // f(u_i*, P1)
        let mut ec_sx1 = Vec::with_capacity(sx1_subdivide.simplices.len());
        ec_sx1.resize(sx1_subdivide.simplices.len(), 0f64);

        for (idx0, s) in sx0_subdivide.simplices.iter().enumerate() {
            if sx1_subdivide.simplices.contains(&s) {
                ec_sx0[idx0] = 1.0;
                continue;
            }
            if sx1_subdivide.simplices.contains(&s.reverse()) {
                continue;
            }

            let q = s.midpoint();
            for u in &sx1.simplices {
                ec_sx0[idx0] += u.beta(&q);
            }
            // eprintln!("idx={}, c={}", idx0, ec_sx0[idx0]);
        }

        for (idx1, u) in sx1_subdivide.simplices.iter().enumerate() {
            let u_in_f1 = sx0_subdivide.simplices.contains(&u);
            let u_in_f2 = sx0_subdivide.simplices.contains(&u.reverse());

            if u_in_f1 || u_in_f2 {
                continue;
            }

            let q = u.midpoint();
            for s in &sx0.simplices {
                ec_sx1[idx1] += s.beta(&q);
            }
            // eprintln!("idx={}, c={}", idx1, ec_sx1[idx1]);
        }
        // eprintln!("0 {:?} {:?}", sx0_subdivide, ec_sx0);
        // eprintln!("1 {:?} {:?}", sx1_subdivide, ec_sx1);

        return (sx0_subdivide, ec_sx0, sx1_subdivide, ec_sx1);
    }

    pub fn bool_intersect(&self, other: &SimplicalChain) -> SimplicalChain {
        let sx0 = self;
        let sx1 = other;

        // subdivide, simplical chain
        let (sx0_subdivide, ec_sx0, sx1_subdivide, ec_sx1) = sx0.subdivide_chracteristics(sx1);

        let mut simplices = Vec::new();

        for (idx, s) in sx0_subdivide.simplices.into_iter().enumerate() {
            if ec_sx0[idx] == 1.0 {
                simplices.push(s);
            }
        }

        for (idx, s) in sx1_subdivide.simplices.into_iter().enumerate() {
            if ec_sx1[idx] == 1.0 {
                simplices.push(s);
            }
        }

        SimplicalChain { simplices }
    }

    pub fn bool_union(&self, other: &SimplicalChain) -> SimplicalChain {
        let sx0 = self;
        let sx1 = other;

        // subdivide, simplical chain
        let (sx0_subdivide, ec_sx0, sx1_subdivide, ec_sx1) = sx0.subdivide_chracteristics(sx1);

        let mut simplices = Vec::new();

        for (idx, s) in sx0_subdivide.simplices.into_iter().enumerate() {
            if ec_sx0[idx] != 1.0 {
                simplices.push(s);
            }
        }
        for (idx, s) in sx1_subdivide.simplices.into_iter().enumerate() {
            if ec_sx1[idx] != 1.0 {
                simplices.push(s);
            }
        }

        SimplicalChain { simplices }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Simplex {
    // omit original point
    pub src: Point<f64>,
    pub dst: Point<f64>,
}

impl Simplex {
    fn reverse(&self) -> Self {
        Simplex {
            src: self.dst,
            dst: self.src,
        }
    }
    fn midpoint(&self) -> Point<f64> {
        lerp_pf64(&self.src, &self.dst, 0.5)
    }
}

const ORIGIN: Point<f64> = Point::new([0f64, 0f64]);

impl Simplex {
    // clockwise = -1
    // counterclockwise = 1
    // colinear = 0
    fn sign(&self) -> Orientation {
        Point::orient_along_direction(&ORIGIN, Direction::Through(&self.src), &self.dst)
    }

    fn on_non_original_edge(&self, q: &Point<f64>) -> bool {
        use std::cmp::Ordering::*;
        use Orientation::*;

        let dir2 = Point::orient_along_direction(&self.src, Direction::Through(&self.dst), q);
        match dir2 {
            CoLinear => {
                // given point might be on non-original edge
                let cmp_0 = self.src.cmp_distance_to(q, &self.dst);
                let cmp_1 = self.dst.cmp_distance_to(q, &self.src);

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
        /*
        eprintln!(
            "characteristic: sign={:?}, q={:?}, self={:?}",
            q, sign, self
        );
        */
        let signval = match sign {
            ClockWise => -1.0,
            CounterClockWise => 1.0,
            CoLinear => 0.0,
        };

        let dir0 = Point::orient_along_direction(&ORIGIN, Direction::Through(&self.src), q);
        let dir1 = Point::orient_along_direction(&ORIGIN, Direction::Through(&self.dst), q);

        match (dir0, dir1) {
            // Q is on some original edge
            (CoLinear, _) => {
                if dir0 == CoLinear && ORIGIN.cmp_distance_to(q, &self.src) != Greater {
                    signval / 2.0
                } else {
                    0.0
                }
            }
            (_, CoLinear) => {
                if dir1 == CoLinear && ORIGIN.cmp_distance_to(q, &self.dst) != Greater {
                    signval / 2.0
                } else {
                    0.0
                }
            }
            (CounterClockWise, ClockWise) => {
                let dir2 =
                    Point::orient_along_direction(&self.src, Direction::Through(&self.dst), q);
                match (sign, dir2) {
                    (CounterClockWise, CounterClockWise) => signval,
                    _ => 0.0,
                }
            }
            (ClockWise, CounterClockWise) => {
                let dir2 =
                    Point::orient_along_direction(&self.src, Direction::Through(&self.dst), q);
                match (sign, dir2) {
                    (ClockWise, ClockWise) => signval,
                    _ => 0.0,
                }
            }
            _ => 0.0,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn simplex_sign() {
        let ccw = Simplex {
            src: Point::new([1.0, 0.0]),
            dst: Point::new([1.0, 1.0]),
        };
        assert_eq!(ccw.sign(), Orientation::CounterClockWise);

        let cw = Simplex {
            src: Point::new([1.0, 1.0]),
            dst: Point::new([1.0, 0.0]),
        };
        assert_eq!(cw.sign(), Orientation::ClockWise);

        let cl = Simplex {
            src: Point::new([1.0, 0.0]),
            dst: Point::new([2.0, 0.0]),
        };
        assert_eq!(cl.sign(), Orientation::CoLinear);
    }

    #[test]
    fn simplex_beta() {
        let s = Simplex {
            src: Point::new([1.0, 0.0]),
            dst: Point::new([1.0, 1.0]),
        };

        // original edge
        assert_eq!(s.beta(&Point::new([0.5, 0.0])), 0.5);

        // inside
        assert_eq!(s.beta(&Point::new([0.9, 0.5])), 1.0);

        // outside
        assert_eq!(s.beta(&Point::new([1.1, 0.5])), 0.0);
    }

    #[test]
    fn simplex_beta_neg() {
        let s = Simplex {
            src: Point::new([1.0, 1.0]),
            dst: Point::new([1.0, 0.0]),
        };

        // original edge
        assert_eq!(s.beta(&Point::new([0.5, 0.0])), -0.5);

        // inside
        assert_eq!(s.beta(&Point::new([0.9, 0.5])), -1.0);

        // outside
        assert_eq!(s.beta(&Point::new([1.1, 0.5])), 0.0);
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

    fn points_along(
        src: &Point<f64>,
        dst: &Point<f64>,
        subdivide: usize,
        out: &mut Vec<Point<f64>>,
    ) {
        for i in 0..subdivide {
            let t = ((i + 1) as f64) / (subdivide as f64);
            out.push(lerp_pf64(src, dst, t));
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
    fn intersects_contains() {
        let p0 = polygon_cube(Point::new([0.0, 0.0]), 1.0);
        let s0 = SimplicalChain::from_polygon(&p0);

        let p1 = polygon_cube(Point::new([0.0, 0.0]), 2.0);
        let s1 = SimplicalChain::from_polygon(&p1);

        let sx = s0.bool_intersect(&s1);
        assert_eq!(sx.simplices.len(), 4);
        for s in &sx.simplices {
            assert!(s0.simplices.contains(s));
        }
    }

    #[test]
    fn union_contains() {
        let p0 = polygon_cube(Point::new([0.0, 0.0]), 1.0);
        let s0 = SimplicalChain::from_polygon(&p0);

        let p1 = polygon_cube(Point::new([0.0, 0.0]), 2.0);
        let s1 = SimplicalChain::from_polygon(&p1);

        let sx = s0.bool_union(&s1);
        assert_eq!(sx.simplices.len(), 4);
        for s in &sx.simplices {
            assert!(s1.simplices.contains(s));
        }
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

    #[ignore]
    #[test]
    fn intersect_exact_edge() {
        let p0 = polygon_cube(Point::new([0.0, 0.0]), 1.0);
        let s0 = SimplicalChain::from_polygon(&p0);

        let p1 = polygon_cube(Point::new([2.0, 0.0]), 1.0);
        let s1 = SimplicalChain::from_polygon(&p1);

        let s = s0.bool_intersect(&s1);
        eprintln!("out={:?}", s);
        todo!();
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
    fn intersect_basic() {
        let p0 = polygon_cube(Point::new([2.0, 2.0]), 2.0);
        let s0 = SimplicalChain::from_polygon(&p0);

        let p1 = polygon_cube(Point::new([4.0, 4.0]), 2.0);
        let s1 = SimplicalChain::from_polygon(&p1);

        let sx = s0.bool_intersect(&s1);

        let p_expected = polygon_cube(Point::new([3.0, 3.0]), 1.0);
        let sx_expected = SimplicalChain::from_polygon(&p_expected);

        assert_eq!(sx.simplices.len(), sx_expected.simplices.len());

        for s in sx_expected.simplices {
            assert!(sx.simplices.contains(&s));
        }
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

    #[test]
    fn test_lerp() {
        assert_eq!(lerp_f64(1.0, 3.0, 0.0), 1.0);
        assert_eq!(lerp_f64(1.0, 3.0, 1.0), 3.0);

        assert_eq!(
            lerp_pf64(&Point::new([0.0, 0.0]), &Point::new([2.0, 2.0]), 0.5),
            Point::new([1.0, 1.0])
        );
    }
}
