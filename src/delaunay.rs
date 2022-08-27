// https://www.personal.psu.edu/cxc11/AERSP560/DELAUNEY/13_Two_algorithms_Delauney.pdf
use rgeometry::{data::*, Orientation};

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy)]
pub struct TriIdx(pub usize);
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy)]
pub struct VertIdx(pub usize);

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy)]
pub struct SubIdx(pub usize);

type Result<T> = anyhow::Result<T>;

impl SubIdx {
    pub fn ccw(self) -> Self {
        Self((self.0 + 1) % 3)
    }
    pub fn cw(self) -> Self {
        Self((self.0 + 2) % 3)
    }
}

pub struct CutResult {
    pub cuts: Vec<(VertIdx, VertIdx)>,
    pub cut_triangles: Vec<TriIdx>,
    pub contour_cw: Vec<(TriIdx, SubIdx)>,
    pub contour_ccw: Vec<(TriIdx, SubIdx)>,
}

#[derive(Debug, PartialEq, Eq)]
enum TriangularNetworkLocation {
    InTriangle(TriIdx),
    // colinear with edge
    Colinear(TriIdx, SubIdx),
    Outside(TriIdx, SubIdx),
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

pub fn is_super(idx: VertIdx) -> bool {
    idx.0 < 3
}

#[derive(Debug)]
pub struct TriangularNetwork {
    pub vertices: Vec<Point<f64>>,
    pub triangles: Vec<Triangle>,
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
                vertices: [VertIdx(0), VertIdx(1), VertIdx(2)],
                neighbors: [None, None, None],
            }],
        }
    }

    pub fn tri(&self, idx: TriIdx) -> &Triangle {
        &self.triangles[idx.0]
    }

    pub fn tri_mut(&mut self, idx: TriIdx) -> &mut Triangle {
        &mut self.triangles[idx.0]
    }

    fn add_tri(&mut self) -> TriIdx {
        let v = VertIdx(0);
        let idx = self.triangles.len();
        self.triangles.push(Triangle {
            vertices: [v, v, v],
            neighbors: [None, None, None],
        });
        TriIdx(idx)
    }

    pub fn vert(&self, idx: VertIdx) -> &Point<f64> {
        &self.vertices[idx.0]
    }

    fn add_vert(&mut self, p: Point<f64>) -> VertIdx {
        let idx = self.vertices.len();
        self.vertices.push(p);
        VertIdx(idx)
    }

    pub fn tri_vert(&self, tri_idx: TriIdx, idx: SubIdx) -> &Point<f64> {
        self.vert(self.tri(tri_idx).vertices[idx.0])
    }

    pub fn find_vert_dest(&self, v_from: VertIdx, v_to: VertIdx) -> Option<(TriIdx, SubIdx)> {
        use Orientation::*;

        let p_end = self.vert(v_to);

        for (t_idx, t) in self.triangles.iter().enumerate() {
            let t_idx = TriIdx(t_idx);
            if let Some(idx) = t.vertex_idx(v_from) {
                let p0 = self.tri_vert(t_idx, idx);
                let p1 = self.tri_vert(t_idx, idx.ccw());
                let p2 = self.tri_vert(t_idx, idx.cw());

                let d0 = Point::orient_along_direction(p0, Direction::Through(p1), p_end);
                let d1 = Point::orient_along_direction(p0, Direction::Through(p2), p_end);

                match (d0, d1) {
                    (CounterClockWise, ClockWise) => {
                        return Some((t_idx, idx));
                    }
                    _ => (),
                }
            }
        }
        None
    }

    pub fn cut_restore_subdevide(
        &mut self,
        idx_p: Option<TriIdx>,
        slice: &[(Option<(TriIdx, SubIdx)>, VertIdx)],
        indices: &mut Vec<TriIdx>,
        out: &mut Vec<(VertIdx, VertIdx)>,
    ) -> Option<TriIdx> {
        if slice.len() < 2 {
            return None;
        }
        if slice.len() == 2 {
            if let Some((tri, sub)) = slice[1].0 {
                self.tri_mut(tri).neighbors[sub.0] = idx_p;
                return Some(tri);
            }
            return None;
        }

        let last = slice.len() - 1;

        let p_start = self.vert(slice[0].1);
        let p_end = self.vert(slice[last].1);

        let mut curidx = 1;

        for i in 2..last {
            let cur = self.vert(slice[curidx].1);
            let next = self.vert(slice[i].1);
            if inside_circle(p_start, cur, p_end, next) {
                curidx = i;
            }
        }

        let idx_self = indices.pop().unwrap();

        let idx_t0 =
            self.cut_restore_subdevide(Some(idx_self), &slice[0..curidx + 1], indices, out);
        let idx_t1 =
            self.cut_restore_subdevide(Some(idx_self), &slice[curidx..slice.len()], indices, out);

        *self.tri_mut(idx_self) = Triangle {
            vertices: [slice[0].1, slice[curidx].1, slice[last].1],
            neighbors: [idx_p, idx_t0, idx_t1],
        };

        out.push((slice[0].1, slice[curidx].1));
        out.push((slice[curidx].1, slice[last].1));

        Some(idx_self)
    }

    pub fn cut_restore(&mut self, res: &CutResult) -> Result<Vec<(VertIdx, VertIdx)>> {
        let mut out = Vec::new();

        let mut verts_ccw = Vec::new();
        for (idx, (tri, sub)) in res.contour_ccw.iter().rev().enumerate() {
            let t = self.tri(*tri);
            if idx == 0 {
                verts_ccw.push((None, t.vert(sub.cw())));
            }
            let n = t
                .neighbor(*sub)
                .map(|n| (n, self.tri(n).neighbor_idx(*tri).unwrap()));
            verts_ccw.push((n, t.vert(*sub)));
        }

        let mut verts_cw = Vec::new();
        for (idx, (tri, sub)) in res.contour_cw.iter().enumerate() {
            let t = self.tri(*tri);
            if idx == 0 {
                verts_cw.push((None, t.vert(sub.cw())));
            }
            let n = t
                .neighbor(*sub)
                .map(|n| (n, self.tri(n).neighbor_idx(*tri).unwrap()));
            verts_cw.push((n, t.vert(*sub)));
        }

        let mut triangles = res.cut_triangles.clone();

        let idx0 = self.cut_restore_subdevide(None, &verts_ccw, &mut triangles, &mut out);
        let idx1 = self.cut_restore_subdevide(None, &verts_cw, &mut triangles, &mut out);

        if let Some(idx0) = idx0 {
            self.tri_mut(idx0).neighbors[0] = idx1;
        }
        if let Some(idx1) = idx1 {
            self.tri_mut(idx1).neighbors[0] = idx0;
        }
        self.check_invariant()?;

        Ok(out)
    }

    pub fn cut(&self, v_from: VertIdx, v_to: VertIdx) -> CutResult {
        use Orientation::*;

        let mut cuts = vec![(v_from, v_to)];
        let mut cut_triangles = vec![];
        let mut contour_cw = vec![];
        let mut contour_ccw = vec![];

        let p_start = self.vert(v_from);
        let p_end = self.vert(v_to);

        #[derive(Debug)]
        enum CutIter {
            FromVertex(TriIdx, SubIdx),
            ToEdge(TriIdx, SubIdx),
        }

        let mut cur = {
            match self.find_vert_dest(v_from, v_to) {
                Some((t_idx, idx)) => Some(CutIter::FromVertex(t_idx, idx)),
                _ => None,
            }
        };

        while let Some(iter) = cur.take() {
            match iter {
                // ray from vertex idx
                CutIter::FromVertex(t_idx, idx) => {
                    cut_triangles.push(t_idx);

                    let t = self.tri(t_idx);

                    let v0 = t.vert(idx);
                    contour_ccw.push((t_idx, idx));
                    contour_cw.push((t_idx, idx.ccw()));
                    if v0 == v_to {
                        break;
                    }
                    let v1 = t.vert(idx.ccw());
                    let v2 = t.vert(idx.cw());

                    cuts.push((v1, v2));

                    match t.neighbor(idx.cw()) {
                        Some(neighbor) => {
                            let t_neighbor = self.tri(neighbor);
                            if let Some(idx_neighbor) = t_neighbor.neighbor_idx(t_idx) {
                                cur = Some(CutIter::ToEdge(neighbor, idx_neighbor));
                            } else {
                                todo!();
                            }
                        }
                        None => todo!(),
                    }
                }
                CutIter::ToEdge(t_idx, idx) => {
                    cut_triangles.push(t_idx);

                    let t = self.tri(t_idx);

                    let v1 = t.vert(idx.ccw());

                    let p0 = self.tri_vert(t_idx, idx);
                    let p1 = self.tri_vert(t_idx, idx.ccw());
                    let p2 = self.tri_vert(t_idx, idx.cw());

                    // should not colinear
                    let d0 = Point::orient_along_direction(p_start, Direction::Through(p_end), p0);
                    let d1 = Point::orient_along_direction(p_start, Direction::Through(p_end), p1);
                    // should not colineaer
                    let d2 = Point::orient_along_direction(p_start, Direction::Through(p_end), p2);

                    let idx_n = if d1 == CoLinear {
                        contour_ccw.push((t_idx, idx.cw()));
                        contour_cw.push((t_idx, idx.ccw()));
                        cur = self
                            .find_vert_dest(v1, v_to)
                            .map(|(t_idx, idx)| CutIter::FromVertex(t_idx, idx));
                        continue;
                    } else if d1.reverse() == d0 {
                        contour_ccw.push((t_idx, idx.cw()));
                        idx.ccw()
                    } else if d1.reverse() == d2 {
                        contour_cw.push((t_idx, idx.ccw()));
                        idx.cw()
                    } else {
                        todo!();
                    };

                    let v_t_from = t.vert(idx_n.cw());
                    let v_t_to = t.vert(idx_n);

                    cuts.push((v_t_from, v_t_to));

                    if let Some(idx_neighbor) = t.neighbor(idx_n) {
                        let t_neighbor = self.tri(idx_neighbor);
                        let idx = t_neighbor.neighbor_idx(t_idx).unwrap();
                        cur = Some(CutIter::ToEdge(idx_neighbor, idx));
                    } else {
                        break;
                    }
                }
            }
        }

        CutResult {
            cuts,
            cut_triangles,
            contour_cw,
            contour_ccw,
        }
    }

    fn check_invariant_tri(&self, idx: TriIdx) -> Result<()> {
        let t = self.tri(idx);
        for i in 0..3 {
            let i = SubIdx(i);
            if let Some(idx_neighbor) = t.neighbor(i) {
                let n = self.tri(idx_neighbor);
                let violated = if let Some(j) = n.neighbor_idx(idx) {
                    t.vert(i) != n.vert(j.cw()) || t.vert(i.cw()) != n.vert(j)
                } else {
                    false
                };

                if violated {
                    anyhow::bail!(
                        "invariant violated: {:?}={:?}, {:?}={:?}",
                        idx,
                        t,
                        idx_neighbor,
                        n
                    );
                }
            }
        }
        Ok(())
    }

    fn check_invariant(&self) -> Result<()> {
        for idx in 0..self.triangles.len() {
            self.check_invariant_tri(TriIdx(idx))?;
        }
        Ok(())
    }

    fn maybe_swap(&mut self, idx0: TriIdx, reductions: &mut usize) -> Result<bool> {
        let idx1 = match self.tri(idx0).neighbor(SubIdx(2)) {
            Some(idx) => idx,
            None => return Ok(false),
        };

        if *reductions == 0 {
            return Ok(false);
        }
        *reductions -= 1;

        let t0 = self.tri(idx0);
        let t1 = self.tri(idx1);

        let t0_t1_idx = SubIdx(2);
        let t1_t0_idx = t1.neighbor_idx(idx0).unwrap();

        let v0 = t0.vert(t0_t1_idx);
        let v1 = t0.vert(t0_t1_idx.ccw());
        let v2 = t0.vert(t0_t1_idx.cw());
        let v3 = t1.vert(t1_t0_idx.ccw());

        let p0 = self.vert(v0);
        let p1 = self.vert(v1);
        let p2 = self.vert(v2);
        let p3 = self.vert(v3);

        // check if four points are convex in order
        let d0 = Point::orient_along_direction(p0, Direction::Through(p2), p1);
        let d1 = Point::orient_along_direction(p0, Direction::Through(p2), p3);
        let d2 = Point::orient_along_direction(p1, Direction::Through(p3), p0);
        let d3 = Point::orient_along_direction(p1, Direction::Through(p3), p2);

        if d0 == d1 || d2 == d3 {
            return Ok(false);
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
            return Ok(false);
        }

        // swap
        let n0 = t1.neighbor(t1_t0_idx.cw());
        let n1 = t0.neighbor(t0_t1_idx.ccw());
        let n2 = t0.neighbor(t0_t1_idx.cw());
        let n3 = t1.neighbor(t1_t0_idx.ccw());

        *self.tri_mut(idx0) = Triangle {
            vertices: [v1, v2, v3],
            neighbors: [Some(idx1), n2, n3],
        };

        *self.tri_mut(idx1) = Triangle {
            vertices: [v1, v3, v0],
            neighbors: [n1, Some(idx0), n0],
        };

        // n0, n2 stays same, n1, n3 changes neighbor
        if let Some(idx) = n1 {
            self.tri_mut(idx).update_neighbor(idx0, idx1);
        }
        if let Some(idx) = n3 {
            self.tri_mut(idx).update_neighbor(idx1, idx0);
        }

        self.check_invariant()?;

        self.maybe_swap(idx0, reductions)?;
        self.maybe_swap(idx1, reductions)?;

        self.check_invariant()?;

        Ok(true)
    }

    pub fn insert(&mut self, p: &Point<f64>, reductions: &mut usize) -> Result<()> {
        use TriangularNetworkLocation::*;

        if *reductions == 0 {
            return Ok(());
        }
        *reductions -= 1;

        match self.locate_recursive(p) {
            InTriangle(idx_t) => {
                let idx_v = self.add_vert(*p);
                let t = self.tri(idx_t).clone();

                let idx_t0 = idx_t;
                let idx_t1 = self.add_tri();
                let idx_t2 = self.add_tri();

                let [v0, v1, v2] = t.vertices;
                let [n0, n1, n2] = t.neighbors;

                *self.tri_mut(idx_t0) = Triangle {
                    vertices: [idx_v, v0, v1],
                    neighbors: [Some(idx_t1), Some(idx_t2), n1],
                };
                *self.tri_mut(idx_t1) = Triangle {
                    vertices: [idx_v, v1, v2],
                    neighbors: [Some(idx_t2), Some(idx_t0), n2],
                };
                *self.tri_mut(idx_t2) = Triangle {
                    vertices: [idx_v, v2, v0],
                    neighbors: [Some(idx_t0), Some(idx_t1), n0],
                };

                if let Some(idx_neighbor) = n2 {
                    self.tri_mut(idx_neighbor).update_neighbor(idx_t0, idx_t1);
                }
                if let Some(idx_neighbor) = n0 {
                    self.tri_mut(idx_neighbor).update_neighbor(idx_t0, idx_t2);
                }

                self.check_invariant()?;

                self.maybe_swap(idx_t0, reductions)?;
                self.maybe_swap(idx_t1, reductions)?;
                self.maybe_swap(idx_t2, reductions)?;

                self.check_invariant()?;
            }

            Colinear(idx_t, idx_neighbor) => {
                if p == self.tri_vert(idx_t, SubIdx(0))
                    || p == self.tri_vert(idx_t, SubIdx(1))
                    || p == self.tri_vert(idx_t, SubIdx(2))
                {
                    // ignore duplicated points
                    return Ok(());
                }

                let idx_t0 = idx_t;
                let t0 = self.tri(idx_t0).clone();

                let idx_t1 = t0.neighbor(idx_neighbor);

                let idx_t2 = self.add_tri();
                let idx_t3 = if idx_t1.is_some() {
                    Some(self.add_tri())
                } else {
                    None
                };

                let v0 = t0.vert(idx_neighbor);
                let v1 = t0.vert(idx_neighbor.ccw());
                let v2 = t0.vert(idx_neighbor.cw());
                let idx_v = self.add_vert(*p);

                *self.tri_mut(idx_t0) = Triangle {
                    vertices: [idx_v, v0, v1],
                    neighbors: [Some(idx_t2), idx_t1, t0.neighbor(idx_neighbor.ccw())],
                };

                *self.tri_mut(idx_t2) = Triangle {
                    vertices: [idx_v, v1, v2],
                    neighbors: [idx_t3, Some(idx_t0), t0.neighbor(idx_neighbor.cw())],
                };

                if let Some(idx_t1) = idx_t1 {
                    let idx_t3 = idx_t3.unwrap();

                    let t1 = self.tri(idx_t1).clone();
                    let idx_neighbor = t1.neighbor_idx(idx_t).unwrap();

                    let v0 = t1.vert(idx_neighbor);
                    let v1 = t1.vert(idx_neighbor.ccw());
                    let v2 = t1.vert(idx_neighbor.cw());

                    *self.tri_mut(idx_t1) = Triangle {
                        vertices: [idx_v, v1, v2],
                        neighbors: [Some(idx_t0), Some(idx_t3), t1.neighbor(idx_neighbor.cw())],
                    };

                    *self.tri_mut(idx_t3) = Triangle {
                        vertices: [idx_v, v0, v1],
                        neighbors: [Some(idx_t1), Some(idx_t2), t1.neighbor(idx_neighbor.ccw())],
                    };
                }

                self.check_invariant_tri(idx_t0)?;
                if let Some(idx_t1) = idx_t1 {
                    self.check_invariant_tri(idx_t1)?;
                }
                self.check_invariant_tri(idx_t2)?;
                if let Some(idx_t3) = idx_t1 {
                    self.check_invariant_tri(idx_t3)?;
                }
            }

            Outside(_, _) => {
                todo!();
            }
        }

        self.check_invariant()
    }

    fn locate(&self, start: TriIdx, p: &Point<f64>) -> TriangularNetworkLocation {
        use Orientation::*;
        use TriangularNetworkLocation::*;

        let p0 = self.tri_vert(start, SubIdx(0));
        let p1 = self.tri_vert(start, SubIdx(1));
        let p2 = self.tri_vert(start, SubIdx(2));

        // 0, self <- ccw | cw -> 1
        let d0 = Point::orient_along_direction(p0, Direction::Through(p1), p);
        // 1, self <- ccw | cw -> 2
        let d1 = Point::orient_along_direction(p1, Direction::Through(p2), p);
        // 2, self <- ccw | cw -> 0
        let d2 = Point::orient_along_direction(p2, Direction::Through(p0), p);

        match (d0, d1, d2) {
            (CounterClockWise, CounterClockWise, CounterClockWise) => return InTriangle(start),
            // handle cw case first
            (ClockWise, _, _) => Outside(start, SubIdx(1)),
            (_, ClockWise, _) => Outside(start, SubIdx(2)),
            (_, _, ClockWise) => Outside(start, SubIdx(0)),
            // TODO: (Colinear, Colinear, _) -> exact point?
            (CoLinear, _, _) => Colinear(start, SubIdx(1)),
            (_, CoLinear, _) => Colinear(start, SubIdx(2)),
            (_, _, CoLinear) => Colinear(start, SubIdx(0)),
        }
    }

    fn locate_recursive(&self, p: &Point<f64>) -> TriangularNetworkLocation {
        let mut start = TriIdx(0);

        loop {
            use TriangularNetworkLocation::*;

            start = match self.locate(start, p) {
                InTriangle(idx) => return InTriangle(idx),
                Colinear(idx, idx_neighbor) => return Colinear(idx, idx_neighbor),
                Outside(idx, idx_neighbor) => match self.tri(idx).neighbor(idx_neighbor) {
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
    pub vertices: [VertIdx; 3],
    //  - self, neighbors[0], neighbors[1] meets at vertices[0], counterclockwise
    //  - vertices[0], vertices[1] meets with neighbors[1], ...
    neighbors: [Option<TriIdx>; 3],
}

impl Triangle {
    pub fn vert(&self, idx: SubIdx) -> VertIdx {
        self.vertices[idx.0]
    }

    fn neighbor(&self, idx: SubIdx) -> Option<TriIdx> {
        self.neighbors[idx.0]
    }

    fn update_neighbor(&mut self, idx_from: TriIdx, idx_to: TriIdx) -> bool {
        for i in 0..3 {
            if self.neighbors[i] == Some(idx_from) {
                self.neighbors[i] = Some(idx_to);
                return true;
            }
        }
        todo!();
    }

    fn vertex_idx(&self, v_idx: VertIdx) -> Option<SubIdx> {
        self.vertices.iter().position(|p| *p == v_idx).map(SubIdx)
    }

    fn neighbor_idx(&self, idx: TriIdx) -> Option<SubIdx> {
        self.neighbors
            .iter()
            .position(|p| *p == Some(idx))
            .map(SubIdx)
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
            (0.5, 0.0, Colinear(TriIdx(0), SubIdx(1))),
            (0.5, -0.1, Outside(TriIdx(0), SubIdx(1))),
            (1.0, 0.5, Colinear(TriIdx(0), SubIdx(2))),
            (1.1, 0.5, Outside(TriIdx(0), SubIdx(2))),
            (0.5, 0.5, Colinear(TriIdx(0), SubIdx(0))),
            (0.4, 0.6, Outside(TriIdx(0), SubIdx(0))),
            (0.5, 0.1, InTriangle(TriIdx(0))),
        ];

        for (x, y, expected) in cases {
            assert_eq!(net.locate(TriIdx(0), &Point::new([x, y])), expected);
        }
    }
}
