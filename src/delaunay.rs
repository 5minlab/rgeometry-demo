// https://www.personal.psu.edu/cxc11/AERSP560/DELAUNEY/13_Two_algorithms_Delauney.pdf
use rgeometry::{data::*, Orientation};

#[derive(Debug)]
pub struct TriangularNetwork {
    pub vertices: Vec<Point<f64>>,
    pub triangles: Vec<Triangle>,
}

#[derive(Debug, PartialEq, Eq)]
enum TriangularNetworkLocation {
    InTriangle(usize),
    // colinear with edge
    Colinear(usize, usize),
    Outside(usize, usize),
}

fn det2(v: [f64; 4]) -> f64 {
    let [a, b, c, d] = v;
    a * d - b * c
}

fn det3(v: [f64; 9]) -> f64 {
    let [a, b, c, d, e, f, g, h, i] = v;
    a * det2([e, f, h, i]) - b * det2([d, f, g, i]) + c * det2([d, e, g, h])
}

fn inside_circle(a: &Point<f64>, b: &Point<f64>, c: &Point<f64>, d: &Point<f64>) -> bool {
    let [ax, ay] = a.array;
    let [bx, by] = b.array;
    let [cx, cy] = c.array;
    let [dx, dy] = d.array;

    det3([
        ax - dx,
        ay - dy,
        ax * ax - dx * dx + ay * ay - dy * dy,
        bx - dx,
        by - dy,
        bx * bx - dx * dx + by * by - dy * dy,
        cx - dx,
        cy - dy,
        cx * cx - dx * dx + cy * cy - dy * dy,
    ]) > 0.0
}

fn is_super(idx: usize) -> bool {
    idx < 3
}

impl TriangularNetwork {
    pub fn new(p0: Point<f64>, p1: Point<f64>, p2: Point<f64>) -> Self {
        let (p1, p2) = match Point::orient_along_direction(&p0, Direction::Through(&p1), &p2) {
            Orientation::CounterClockWise => (p1, p2),
            Orientation::CoLinear => todo!(),
            // TODO: colinear?
            _ => (p2, p1),
        };

        Self {
            vertices: vec![p0, p1, p2],
            triangles: vec![Triangle {
                vertices: [0, 1, 2],
                neighbors: [None, None, None],
            }],
        }
    }

    fn check_invariant_tri(&self, idx: usize) {
        let t = &self.triangles[idx];
        for i in 0..3 {
            if let Some(idx_neighbor) = t.neighbors[i] {
                let n = &self.triangles[idx_neighbor];
                if let Some(j) = n.neighbor_idx(idx) {
                    // eprintln!("{}={:?}, {}={:?}", idx, t, idx_neighbor, n);
                    assert_eq!(t.vertices[i], n.vertices[(j + 2) % 3]);
                    assert_eq!(t.vertices[(i + 2) % 3], n.vertices[j]);
                } else {
                    panic!("{:#?}", self);
                }
            }
        }
    }

    fn check_invariant(&self) {
        for idx in 0..self.triangles.len() {
            self.check_invariant_tri(idx);
        }
    }

    fn maybe_swap(&mut self, idx0: usize, reductions: &mut usize) -> bool {
        let idx1 = match self.triangles[idx0].neighbors[2] {
            Some(idx) => idx,
            None => return false,
        };

        if *reductions == 0 {
            return false;
        }
        *reductions -= 1;

        let t0 = &self.triangles[idx0];
        let t1 = &self.triangles[idx1];

        let t0_t1_idx = 2;
        let t1_t0_idx = t1.neighbor_idx(idx0).unwrap();

        let v0 = t0.vertices[t0_t1_idx];
        let v1 = t0.vertices[(t0_t1_idx + 1) % 3];
        let v2 = t0.vertices[(t0_t1_idx + 2) % 3];
        let v3 = t1.vertices[(t1_t0_idx + 1) % 3];

        let p0 = &self.vertices[v0];
        let p1 = &self.vertices[v1];
        let p2 = &self.vertices[v2];
        let p3 = &self.vertices[v3];

        // check if four points are convex in order
        let d0 = Point::orient_along_direction(p0, Direction::Through(p2), p1);
        let d1 = Point::orient_along_direction(p0, Direction::Through(p2), p3);
        let d2 = Point::orient_along_direction(p1, Direction::Through(p3), p0);
        let d3 = Point::orient_along_direction(p1, Direction::Through(p3), p2);

        if d0 == d1 || d2 == d3 {
            return false;
        }

        let should_swap = if is_super(v0) || is_super(v2) {
            true
        } else {
            if is_super(v1) || is_super(v3) {
                false
            } else {
                inside_circle(p0, p1, p2, p3)
            }
        };

        if !should_swap {
            return false;
        }

        // swap
        let n0 = t1.neighbors[(t1_t0_idx + 2) % 3];
        let n1 = t0.neighbors[(t0_t1_idx + 1) % 3];
        let n2 = t0.neighbors[(t0_t1_idx + 2) % 3];
        let n3 = t1.neighbors[(t1_t0_idx + 1) % 3];

        self.triangles[idx0] = Triangle {
            vertices: [v1, v2, v3],
            neighbors: [Some(idx1), n2, n3],
        };

        self.triangles[idx1] = Triangle {
            vertices: [v1, v3, v0],
            neighbors: [n1, Some(idx0), n0],
        };

        // n0, n2 stays same, n1, n3 changes neighbor
        if let Some(idx) = n1 {
            self.triangles[idx].update_neighbor(idx0, idx1);
        }
        if let Some(idx) = n3 {
            self.triangles[idx].update_neighbor(idx1, idx0);
        }

        self.maybe_swap(idx0, reductions);
        self.maybe_swap(idx1, reductions);

        self.check_invariant();

        true
    }

    pub fn insert(&mut self, p: &Point<f64>, reductions: &mut usize) {
        use TriangularNetworkLocation::*;

        if *reductions == 0 {
            return;
        }
        *reductions -= 1;

        match self.locate_recursive(p) {
            InTriangle(idx_t) => {
                let idx_v = self.vertices.len();
                self.vertices.push(*p);

                let t = self.triangles[idx_t].clone();

                let idx_t0 = idx_t;
                let idx_t1 = self.triangles.len();
                let idx_t2 = idx_t1 + 1;

                let [v0, v1, v2] = t.vertices;
                let [n0, n1, n2] = t.neighbors;

                self.triangles[idx_t0] = Triangle {
                    vertices: [idx_v, v0, v1],
                    neighbors: [Some(idx_t1), Some(idx_t2), n1],
                };
                self.triangles.push(Triangle {
                    vertices: [idx_v, v1, v2],
                    neighbors: [Some(idx_t2), Some(idx_t0), n2],
                });
                self.triangles.push(Triangle {
                    vertices: [idx_v, v2, v0],
                    neighbors: [Some(idx_t0), Some(idx_t1), n0],
                });

                if let Some(idx_neighbor) = n2 {
                    self.triangles[idx_neighbor].update_neighbor(idx_t0, idx_t1);
                }
                if let Some(idx_neighbor) = n0 {
                    self.triangles[idx_neighbor].update_neighbor(idx_t0, idx_t2);
                }

                self.maybe_swap(idx_t0, reductions);
                self.maybe_swap(idx_t1, reductions);
                self.maybe_swap(idx_t2, reductions);
            }

            Colinear(idx_t, idx_neighbor) => {
                let idx_v = self.vertices.len();
                self.vertices.push(*p);

                let idx_t0 = idx_t;
                let t0 = self.triangles[idx_t0].clone();
                let idx_t1 = t0.neighbors[idx_neighbor];

                let idx_t2 = self.triangles.len();

                let (t3, idx_t3) = if let Some(idx_t1) = idx_t1 {
                    let t1 = self.triangles[idx_t1].clone();
                    let idx_neighbor_t1 = t1.neighbors.iter().position(|p| p == &Some(idx_t));
                    let t3 = self.triangles[idx_t1].split(
                        idx_t1,
                        idx_neighbor_t1.unwrap(),
                        idx_v,
                        idx_t2,
                        Some(idx_t2),
                    );
                    (Some(t3), Some(idx_t2 + 1))
                } else {
                    (None, None)
                };

                let t2 = self.triangles[idx_t0].split(idx_t0, idx_neighbor, idx_v, idx_t2, idx_t3);

                self.triangles.push(t2);
                if let Some(t3) = t3 {
                    self.triangles.push(t3);
                }

                // TODO: swap
            }

            Outside(_, _) => {
                todo!();
            }
        }

        self.check_invariant();
    }

    fn locate(&self, start: usize, p: &Point<f64>) -> TriangularNetworkLocation {
        use Orientation::*;
        use TriangularNetworkLocation::*;

        let t = &self.triangles[start];
        let p0 = &self.vertices[t.vertices[0]];
        let p1 = &self.vertices[t.vertices[1]];
        let p2 = &self.vertices[t.vertices[2]];

        // 0, self <- ccw | cw -> 1
        let d0 = Point::orient_along_direction(p0, Direction::Through(p1), p);
        // 1, self <- ccw | cw -> 2
        let d1 = Point::orient_along_direction(p1, Direction::Through(p2), p);
        // 2, self <- ccw | cw -> 0
        let d2 = Point::orient_along_direction(p2, Direction::Through(p0), p);

        match (d0, d1, d2) {
            (CounterClockWise, CounterClockWise, CounterClockWise) => return InTriangle(start),
            // handle cw case first
            (ClockWise, _, _) => Outside(start, 1),
            (_, ClockWise, _) => Outside(start, 2),
            (_, _, ClockWise) => Outside(start, 0),
            // TODO: (Colinear, Colinear, _) -> exact point?
            (CoLinear, _, _) => Colinear(start, 1),
            (_, CoLinear, _) => Colinear(start, 2),
            (_, _, CoLinear) => Colinear(start, 0),
        }
    }

    fn locate_recursive(&self, p: &Point<f64>) -> TriangularNetworkLocation {
        let mut start = 0;

        loop {
            use TriangularNetworkLocation::*;

            start = match self.locate(start, p) {
                InTriangle(idx) => return InTriangle(idx),
                Colinear(idx, idx_neighbor) => return Colinear(idx, idx_neighbor),
                Outside(idx, idx_neighbor) => match self.triangles[idx].neighbors[idx_neighbor] {
                    Some(idx) => idx,
                    None => {
                        eprintln!("{:#?}, {:?}", self, p);
                        todo!();
                        // return Outside(idx, idx_neighbor);
                    }
                },
            };
        }
    }
}

#[derive(Clone, Debug)]
pub struct Triangle {
    // counterclockwise
    pub vertices: [usize; 3],
    //  - self, neighbors[0], neighbors[1] meets at vertices[0], counterclockwise
    //  - vertices[0], vertices[1] meets with neighbors[1], ...
    neighbors: [Option<usize>; 3],
}

impl Triangle {
    fn update_neighbor(&mut self, idx_from: usize, idx_to: usize) -> bool {
        for i in 0..3 {
            if self.neighbors[i] == Some(idx_from) {
                self.neighbors[i] = Some(idx_to);
                return true;
            }
        }
        todo!();
    }
    fn neighbor_idx(&self, idx: usize) -> Option<usize> {
        for i in 0..3 {
            if self.neighbors[i] == Some(idx) {
                return Some(i);
            }
        }
        None
    }

    // split triangle
    fn split(
        &mut self,
        idx_self: usize,
        idx_neighbor: usize,
        idx_vert: usize,
        idx_new_tri: usize,
        idx_neighbor_new: Option<usize>,
    ) -> Triangle {
        let mut t = self.clone();
        t.neighbors[idx_neighbor] = idx_neighbor_new;
        t.neighbors[(idx_neighbor + 1) % 3] = Some(idx_self);
        t.vertices[idx_neighbor] = idx_vert;

        self.vertices[(idx_neighbor + 2) % 3] = idx_vert;
        self.neighbors[(idx_neighbor + 2) % 3] = Some(idx_new_tri);

        t
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn locate() {
        use TriangularNetworkLocation::*;
        let p0 = Point::new([0.0, 0.0]);
        let p1 = Point::new([1.0, 0.0]);
        let p2 = Point::new([1.0, 1.0]);

        let net = TriangularNetwork::new(p0, p1, p2);

        let cases = vec![
            (0.5, 0.0, Colinear(0, 1)),
            (0.5, -0.1, Outside(0, 1)),
            (1.0, 0.5, Colinear(0, 2)),
            (1.1, 0.5, Outside(0, 2)),
            (0.5, 0.5, Colinear(0, 0)),
            (0.4, 0.6, Outside(0, 0)),
            (0.5, 0.1, InTriangle(0)),
        ];

        for (x, y, expected) in cases {
            assert_eq!(net.locate_recursive(&Point::new([x, y])), expected);
        }
    }
}
