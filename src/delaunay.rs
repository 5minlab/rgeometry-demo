// https://www.personal.psu.edu/cxc11/AERSP560/DELAUNEY/13_Two_algorithms_Delauney.pdf
use rgeometry::{data::*, Orientation, PolygonScalar};

use crate::boolean::SimplicalChain;

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct TriIdx(pub usize);
impl std::fmt::Debug for TriIdx {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(fmt, "t{}", self.0)
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct VertIdx(pub usize);
impl std::fmt::Debug for VertIdx {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(fmt, "v{}", self.0)
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct SubIdx(pub usize);
impl std::fmt::Debug for SubIdx {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(fmt, "s{}", self.0)
    }
}

type Result<T> = anyhow::Result<T>;

impl SubIdx {
    pub fn ccw(self) -> Self {
        Self((self.0 + 1) % 3)
    }
    pub fn cw(self) -> Self {
        Self((self.0 + 2) % 3)
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug)]
pub struct Edge {
    tri: TriIdx,
    sub: SubIdx,
}
impl Edge {
    fn new(tri: TriIdx, sub: SubIdx) -> Self {
        Self { tri, sub }
    }
}

#[derive(Debug)]
pub struct CutResult {
    pub from: VertIdx,
    pub to: VertIdx,
    pub cuts: Vec<(VertIdx, VertIdx)>,
    pub cut_triangles: Vec<TriIdx>,
    pub contour_cw: Vec<Edge>,
    pub contour_ccw: Vec<Edge>,
}

#[derive(Debug, PartialEq, Eq)]
pub enum TriangularNetworkLocation {
    InTriangle(TriIdx),
    // colinear with edge
    OnVertex(TriIdx, SubIdx),
    // colinear with edge
    OnEdge(Edge),
    Outside(Edge),
}

fn det2<T>(v: [T; 4]) -> T
where
    T: PolygonScalar,
{
    let [a, b, c, d] = v;
    a * d - b * c
}

fn det3<T>(v: [T; 9]) -> T
where
    T: PolygonScalar + Copy,
{
    let [a, b, c, d, e, f, g, h, i] = v;
    a * det2([e, f, h, i]) - b * det2([d, f, g, i]) + c * det2([d, e, g, h])
}

fn inside_circle<T>(a: &Point<T>, b: &Point<T>, c: &Point<T>, d: &Point<T>) -> bool
where
    T: PolygonScalar + Copy,
{
    let [ax, ay] = a.array;
    let [bx, by] = b.array;
    let [cx, cy] = c.array;
    let [dx, dy] = d.array;

    let d = det3([
        ax - dx,
        ay - dy,
        ax * ax - dx * dx + ay * ay - dy * dy,
        bx - dx,
        by - dy,
        bx * bx - dx * dx + by * by - dy * dy,
        cx - dx,
        cy - dy,
        cx * cx - dx * dx + cy * cy - dy * dy,
    ]);
    d > T::from_constant(0)
}

fn is_super(idx: VertIdx) -> bool {
    idx.0 < 3
}

#[derive(Debug)]
struct CutEdge {
    inner: Edge,
    outer: Option<Edge>,
    vert: VertIdx,
}

#[derive(Debug)]
enum CutIter {
    FromVertex(TriIdx, SubIdx),
    CoLinear {
        tri: TriIdx,
        src: SubIdx,
        dst: SubIdx,
    },
    ToEdge(Edge),
}

fn pt_mean<T>(points: &[&Point<T>]) -> Point<T>
where
    T: PolygonScalar + Copy,
{
    let mut x = T::from_constant(0);
    let mut y = T::from_constant(0);
    for p in points {
        x += p.array[0];
        y += p.array[1];
    }
    // TODO
    let l = T::from_constant(points.len() as i8);
    Point::new([x / l, y / l])
}

#[derive(Debug, Clone)]
pub struct TriangularNetwork<T> {
    pub vertices: Vec<Point<T>>,
    pub triangles: Vec<Triangle>,
}

impl<T: PolygonScalar + Copy> TriangularNetwork<T> {
    pub fn new(p0: Point<T>, p1: Point<T>, p2: Point<T>) -> Self {
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

    pub fn vert(&self, idx: VertIdx) -> &Point<T> {
        &self.vertices[idx.0]
    }

    fn add_vert(&mut self, p: Point<T>) -> VertIdx {
        let idx = self.vertices.len();
        self.vertices.push(p);
        VertIdx(idx)
    }

    pub fn tri_vert(&self, tri_idx: TriIdx, idx: SubIdx) -> &Point<T> {
        self.vert(self.tri(tri_idx).vertices[idx.0])
    }

    pub fn edge_duel(&self, edge: &Edge) -> Option<Edge> {
        let t = self.tri(edge.tri);
        let idx_neighbor = t.neighbor(edge.sub)?;
        let t_neighbor = self.tri(idx_neighbor);
        let sub_neighbor = match t_neighbor.neighbor_idx(edge.tri) {
            Some(idx) => idx,
            None => {
                panic!(
                    "invariant: t1={:?}={:?}, t2={:?}={:?}",
                    edge.tri, t, idx_neighbor, t_neighbor
                );
            }
        };
        Some(Edge::new(idx_neighbor, sub_neighbor))
    }

    fn find_vert_dest(&self, v_from: VertIdx, v_to: VertIdx) -> Option<CutIter> {
        use Orientation::*;

        let p_from = self.vert(v_from);
        let l = self.locate_recursive(&p_from);

        let (tri0, sub0) = match l {
            TriangularNetworkLocation::OnVertex(tri, sub) => (tri, sub),
            _ => return None,
        };

        let p_end = self.vert(v_to);
        let mut candidates = Vec::new();

        // TODO: to iterator
        // ccw triangles
        let mut curtri = tri0;
        let mut cursub = sub0;
        loop {
            let t = self.tri(curtri);
            assert_eq!(t.vert(cursub), v_from);
            if t.vert(cursub.cw()) == v_to {
                return None;
            }
            if candidates.contains(&(curtri, cursub)) {
                break;
            }
            candidates.push((curtri, cursub));

            match self.edge_duel(&Edge::new(curtri, cursub)) {
                Some(Edge { tri, sub }) => {
                    curtri = tri;
                    cursub = sub.cw();
                }
                None => break,
            }
        }

        // TODO: to iterator
        // cw triangles
        let mut curtri = tri0;
        let mut cursub = sub0;
        loop {
            let t = self.tri(curtri);
            assert_eq!(t.vert(cursub), v_from);
            if t.vert(cursub.ccw()) == v_to {
                return None;
            }
            if candidates.contains(&(curtri, cursub)) {
                break;
            }
            candidates.push((curtri, cursub));

            match self.edge_duel(&Edge::new(curtri, cursub.ccw())) {
                Some(Edge { tri, sub }) => {
                    if tri == tri0 {
                        break;
                    }
                    curtri = tri;
                    cursub = sub;
                }
                None => break,
            }
        }

        for (t_idx, idx) in candidates {
            let p0 = self.tri_vert(t_idx, idx);
            let p1 = self.tri_vert(t_idx, idx.ccw());
            let p2 = self.tri_vert(t_idx, idx.cw());

            let d0 = Point::orient_along_direction(p0, Direction::Through(p1), p_end);
            let d1 = Point::orient_along_direction(p1, Direction::Through(p2), p_end);
            let d2 = Point::orient_along_direction(p2, Direction::Through(p0), p_end);

            match (d0, d1, d2) {
                (CounterClockWise, ClockWise, CounterClockWise) => {
                    return Some(CutIter::FromVertex(t_idx, idx));
                }
                (CoLinear, ClockWise, CounterClockWise) => {
                    return Some(CutIter::CoLinear {
                        tri: t_idx,
                        src: idx,
                        dst: idx.ccw(),
                    });
                }
                (CounterClockWise, ClockWise, CoLinear) => {
                    return Some(CutIter::CoLinear {
                        tri: t_idx,
                        src: idx,
                        dst: idx.cw(),
                    });
                }
                _ => (),
            }
        }
        None
    }

    pub fn cut(&self, v_from: VertIdx, v_to: VertIdx) -> CutResult {
        use Orientation::*;

        let mut cuts = vec![];
        let mut cut_triangles = vec![];
        let mut contour_ccw = vec![];
        let mut contour_cw = vec![];

        let p_start = self.vert(v_from);
        let p_end = self.vert(v_to);

        let mut cur = self.find_vert_dest(v_from, v_to);
        while let Some(iter) = cur.take() {
            match iter {
                // ray from vertex idx
                CutIter::FromVertex(t_idx, idx) => {
                    cut_triangles.push(t_idx);

                    let t = self.tri(t_idx);

                    let v0 = t.vert(idx);
                    assert!(v0 != v_to);

                    contour_ccw.push(Edge::new(t_idx, idx));
                    contour_cw.push(Edge::new(t_idx, idx.ccw()));
                    let v1 = t.vert(idx.ccw());
                    let v2 = t.vert(idx.cw());

                    cuts.push((v1, v2));

                    let next = self.edge_duel(&Edge::new(t_idx, idx.cw())).unwrap();
                    cur = Some(CutIter::ToEdge(next));
                }

                CutIter::ToEdge(Edge {
                    tri: t_idx,
                    sub: idx,
                }) => {
                    cut_triangles.push(t_idx);

                    let t = self.tri(t_idx);

                    let p0 = self.tri_vert(t_idx, idx);
                    let p1 = self.tri_vert(t_idx, idx.ccw());
                    let p2 = self.tri_vert(t_idx, idx.cw());

                    // should not colinear
                    let d0 = Point::orient_along_direction(p_start, Direction::Through(p_end), p0);
                    let d1 = Point::orient_along_direction(p_start, Direction::Through(p_end), p1);
                    // should not colineaer
                    let d2 = Point::orient_along_direction(p_start, Direction::Through(p_end), p2);

                    let idx_n = if d1 == CoLinear {
                        contour_ccw.push(Edge::new(t_idx, idx.cw()));
                        contour_cw.push(Edge::new(t_idx, idx.ccw()));

                        cur = self.find_vert_dest(t.vert(idx.ccw()), v_to);
                        continue;
                    } else if d1.reverse() == d0 {
                        contour_ccw.push(Edge::new(t_idx, idx.cw()));
                        idx.ccw()
                    } else if d1.reverse() == d2 {
                        contour_cw.push(Edge::new(t_idx, idx.ccw()));
                        idx.cw()
                    } else {
                        dbg!((d0, d1, d2));
                        todo!();
                    };

                    let v_t_from = t.vert(idx_n.cw());
                    let v_t_to = t.vert(idx_n);

                    cuts.push((v_t_from, v_t_to));

                    let next = self.edge_duel(&Edge::new(t_idx, idx_n)).unwrap();
                    cur = Some(CutIter::ToEdge(next));
                }

                CutIter::CoLinear { tri, src, dst } => {
                    let (edge_cw, edge_ccw) = if src.ccw() == dst {
                        let edge_cw = Edge::new(tri, dst);
                        let edge_ccw = self.edge_duel(&edge_cw).unwrap();
                        (edge_cw, edge_ccw)
                    } else {
                        let edge_ccw = Edge::new(tri, src);
                        let edge_cw = self.edge_duel(&edge_ccw).unwrap();
                        (edge_cw, edge_ccw)
                    };
                    contour_cw.push(edge_cw);
                    contour_ccw.push(edge_ccw);

                    cur = self.find_vert_dest(self.tri(tri).vert(dst), v_to);
                }
            }
        }

        CutResult {
            from: v_from,
            to: v_to,
            cut_triangles,
            cuts,
            contour_cw,
            contour_ccw,
        }
    }

    fn cut_apply_subdivide(
        &mut self,
        idx_p: Option<Edge>,
        slice: &[CutEdge],
        indices: &mut Vec<TriIdx>,
        proj: &mut Vec<(Edge, Edge)>,
        out: &mut Vec<(VertIdx, VertIdx)>,
    ) -> Option<TriIdx> {
        assert!(slice.len() > 1);
        if slice.len() == 2 {
            if let CutEdge {
                inner,
                outer: Some(mut outer),
                vert: _,
            } = slice[1]
            {
                if let Some(p) = idx_p {
                    proj.push((inner, p));
                }

                if let Some((_before, after)) = proj.iter().find(|t| t.0 == outer) {
                    outer = *after;
                }

                *self.tri_mut(outer.tri).neighbor_mut(outer.sub) = idx_p.map(|e| e.tri);
                return Some(outer.tri);
            }
            return None;
        }

        let last = slice.len() - 1;

        let v_start = slice[0].vert;
        let v_end = slice[last].vert;

        let p_start = self.vert(v_start);
        let p_end = self.vert(v_end);

        let mut i_mid = 1;
        for i in 2..last {
            let cur = self.vert(slice[i_mid].vert);
            let next = self.vert(slice[i].vert);

            if inside_circle(p_start, cur, p_end, next) {
                i_mid = i;
            }
        }
        let v_mid = slice[i_mid].vert;

        let idx_self = indices.pop().unwrap();

        let idx_t0 = self.cut_apply_subdivide(
            Some(Edge::new(idx_self, SubIdx(1))),
            &slice[..i_mid + 1],
            indices,
            proj,
            out,
        );
        let idx_t1 = self.cut_apply_subdivide(
            Some(Edge::new(idx_self, SubIdx(2))),
            &slice[i_mid..],
            indices,
            proj,
            out,
        );

        *self.tri_mut(idx_self) = Triangle {
            vertices: [v_start, v_mid, v_end],
            neighbors: [idx_p.map(|e| e.tri), idx_t0, idx_t1],
        };

        out.push((v_start, v_mid));
        out.push((v_mid, v_end));

        Some(idx_self)
    }

    fn cut_apply_prepare(&mut self, res: &CutResult) -> Vec<CutEdge> {
        let first = res.contour_cw[0];
        let mut verts = vec![CutEdge {
            inner: first,
            outer: None,
            vert: self.tri(first.tri).vert(first.sub.cw()),
        }];

        for edge in res.contour_cw.iter() {
            verts.push(CutEdge {
                inner: *edge,
                outer: self.edge_duel(edge),
                vert: self.tri(edge.tri).vert(edge.sub),
            });
        }
        for edge in res.contour_ccw.iter().rev() {
            verts.push(CutEdge {
                inner: *edge,
                outer: self.edge_duel(edge),
                vert: self.tri(edge.tri).vert(edge.sub),
            });
        }
        verts
    }

    pub fn cut_apply(&mut self, res: &CutResult) -> Result<Vec<(VertIdx, VertIdx)>> {
        if res.cut_triangles.len() == 0 {
            return Ok(vec![]);
        }

        let mut indices = res.cut_triangles.clone();
        let mut proj = Vec::new();

        let mut out = Vec::new();
        let verts = self.cut_apply_prepare(res);
        let mut slice = verts.as_slice();

        let p_start = self.vert(res.from).clone();
        let p_end = self.vert(res.to).clone();

        let mut out_triangles = Vec::new();
        while slice.len() > 1 {
            let split = slice.iter().skip(1).position(|e| {
                Point::orient_along_direction(
                    &p_start,
                    Direction::Through(&p_end),
                    self.vert(e.vert),
                ) == Orientation::CoLinear
            });
            let i = match split {
                Some(i) => i + 1,
                None => slice.len() - 1,
            };
            let idx = self
                .cut_apply_subdivide(None, &slice[..i + 1], &mut indices, &mut proj, &mut out)
                .unwrap();
            out_triangles.push(idx);
            slice = &slice[i..];
        }
        assert!(out_triangles.len() % 2 == 0);

        for i in 0..(out_triangles.len() / 2) {
            let t_ccw = out_triangles[i];
            let t_cw = out_triangles[out_triangles.len() - 1 - i];

            *self.tri_mut(t_ccw).neighbor_mut(SubIdx(0)) = Some(t_cw);
            *self.tri_mut(t_cw).neighbor_mut(SubIdx(0)) = Some(t_ccw);
        }

        assert!(indices.is_empty());
        self.check_invariant("post-cut_resolve")?;

        Ok(out)
    }

    #[allow(unused)]
    #[cfg(not(debug_assertions))]
    fn check_invariant_tri_opt(&self, _idx: Option<TriIdx>, msg: &str) -> Result<()> {
        Ok(())
    }

    #[cfg(debug_assertions)]
    fn check_invariant_tri_opt(&self, idx: Option<TriIdx>, msg: &str) -> Result<()> {
        if let Some(idx) = idx {
            self.check_invariant_tri(idx, msg)
        } else {
            Ok(())
        }
    }

    #[allow(unused)]
    #[cfg(not(debug_assertions))]
    fn check_invariant_tri(&self, _idx: TriIdx, _msg: &str) -> Result<()> {
        Ok(())
    }

    #[cfg(debug_assertions)]
    fn check_invariant_tri(&self, idx: TriIdx, msg: &str) -> Result<()> {
        let t = self.tri(idx);
        for i in 0..3 {
            let i = SubIdx(i);

            if let Some(idx_neighbor) = t.neighbor(i) {
                let n = self.tri(idx_neighbor);
                let violated = if let Some(j) = n.neighbor_idx(idx) {
                    t.vert(i) != n.vert(j.cw()) || t.vert(i.cw()) != n.vert(j)
                } else {
                    true
                };

                if violated {
                    anyhow::bail!(
                        "invariant violated: {}, {:?}={:?}, {:?}={:?}",
                        msg,
                        idx,
                        t,
                        idx_neighbor,
                        n
                    );
                }
            }

            let e = Edge::new(idx, i);
            if let Some(d) = self.edge_duel(&e) {
                assert_eq!(Some(e), self.edge_duel(&d));
            }
        }
        Ok(())
    }

    #[allow(unused)]
    #[cfg(not(debug_assertions))]
    fn check_invariant(&self, msg: &str) -> Result<()> {
        Ok(())
    }

    #[cfg(debug_assertions)]
    fn check_invariant(&self, msg: &str) -> Result<()> {
        for idx in 0..self.triangles.len() {
            self.check_invariant_tri(TriIdx(idx), msg)?;
        }
        Ok(())
    }

    fn debug_tri(&self, idx: TriIdx) -> String {
        let p0 = self.tri_vert(idx, SubIdx(0));
        let p1 = self.tri_vert(idx, SubIdx(1));
        let p2 = self.tri_vert(idx, SubIdx(2));
        format!(
            "{:?}={:?}, ({:?}, {:?}, {:?})",
            idx,
            self.tri(idx),
            p0,
            p1,
            p2
        )
    }

    fn maybe_swap(&mut self, idx0: TriIdx, reductions: &mut usize) -> Result<bool> {
        use Orientation::*;

        let t0_t1_idx = SubIdx(2);
        let idx1 = match self.tri(idx0).neighbor(t0_t1_idx) {
            Some(idx) => idx,
            None => return Ok(false),
        };

        if *reductions == 0 {
            return Ok(false);
        }
        *reductions -= 1;

        let t0 = self.tri(idx0);
        let t1 = self.tri(idx1);

        let t1_t0_idx = match t1.neighbor_idx(idx0) {
            Some(idx) => idx,
            None => {
                anyhow::bail!(
                    "invalid tri pair: \n{}\n{}\n{:#?}",
                    self.debug_tri(idx0),
                    self.debug_tri(idx1),
                    self,
                );
            }
        };

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

        if d0 == CoLinear || d1 == CoLinear || d2 == CoLinear || d3 == CoLinear {
            return Ok(false);
        }

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
            if !self.tri_mut(idx).update_neighbor(idx0, idx1) {
                anyhow::bail!(
                    "invalid tri pair: \nt0={}\nt1={}\nn1={}",
                    self.debug_tri(idx0),
                    self.debug_tri(idx1),
                    self.debug_tri(idx)
                );
            }
        }
        if let Some(idx) = n3 {
            if !self.tri_mut(idx).update_neighbor(idx1, idx0) {
                anyhow::bail!(
                    "invalid tri pair: \nt0={}\nt1={}\nn3={}",
                    self.debug_tri(idx0),
                    self.debug_tri(idx1),
                    self.debug_tri(idx)
                );
            }
        }

        self.check_invariant_tri(idx0, "pre-swap idx0")?;
        self.check_invariant_tri(idx1, "pre-swap idx1")?;

        self.check_invariant_tri_opt(n0, "pre-swap n0")?;
        self.check_invariant_tri_opt(n1, "pre-swap n1")?;
        self.check_invariant_tri_opt(n2, "pre-swap n2")?;
        self.check_invariant_tri_opt(n3, "pre-swap n3")?;

        self.maybe_swap(idx0, reductions)?;
        self.maybe_swap(idx1, reductions)?;

        self.check_invariant("post-swap")?;

        Ok(true)
    }

    pub fn insert(&mut self, p: &Point<T>, reductions: &mut usize) -> Result<VertIdx> {
        use TriangularNetworkLocation::*;

        if *reductions == 0 {
            anyhow::bail!("reduction");
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

                self.check_invariant_tri(idx_t0, "InTriangle(t0)")?;
                self.check_invariant_tri(idx_t1, "InTriangle(t1)")?;
                self.check_invariant_tri(idx_t2, "InTriangle(t2)")?;

                self.maybe_swap(idx_t0, reductions)?;
                self.maybe_swap(idx_t1, reductions)?;
                self.maybe_swap(idx_t2, reductions)?;

                self.check_invariant("post-InTriangle")?;
                Ok(idx_v)
            }
            OnVertex(tri, sub) => Ok(self.tri(tri).vert(sub)),

            OnEdge(Edge {
                tri: idx_t,
                sub: idx_neighbor,
            }) => {
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

                //       v2(v0)
                //     t3     t2
                // v1     idx_v    v1
                //     t1     t0
                //       v0(v2)

                *self.tri_mut(idx_t0) = Triangle {
                    vertices: [idx_v, v0, v1],
                    neighbors: [Some(idx_t2), idx_t1, t0.neighbor(idx_neighbor.ccw())],
                };

                let n = t0.neighbor(idx_neighbor.cw());
                *self.tri_mut(idx_t2) = Triangle {
                    vertices: [idx_v, v1, v2],
                    neighbors: [idx_t3, Some(idx_t0), n],
                };
                if let Some(n) = n {
                    self.tri_mut(n).update_neighbor(idx_t0, idx_t2);
                }

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

                    let n = t1.neighbor(idx_neighbor.ccw());
                    *self.tri_mut(idx_t3) = Triangle {
                        vertices: [idx_v, v0, v1],
                        neighbors: [Some(idx_t1), Some(idx_t2), n],
                    };
                    if let Some(n) = n {
                        self.tri_mut(n).update_neighbor(idx_t1, idx_t3);
                    }
                }

                self.check_invariant_tri(idx_t0, "Colinear(t0)")?;
                if let Some(idx_t1) = idx_t1 {
                    self.check_invariant_tri(idx_t1, "Colinear(t1)")?;
                }
                self.check_invariant_tri(idx_t2, "Colinear(t2)")?;
                if let Some(idx_t3) = idx_t3 {
                    self.check_invariant_tri(idx_t3, "Colinear(t3)")?;
                }

                self.maybe_swap(idx_t0, reductions)?;
                self.maybe_swap(idx_t2, reductions)?;
                if let Some(idx) = idx_t1 {
                    self.maybe_swap(idx, reductions)?;
                }
                if let Some(idx) = idx_t3 {
                    self.maybe_swap(idx, reductions)?;
                }

                self.check_invariant("post-Colinear")?;
                Ok(idx_v)
            }

            Outside(_e) => {
                todo!();
            }
        }
    }

    pub fn locate(&self, start: TriIdx, p: &Point<T>) -> TriangularNetworkLocation {
        use Orientation::*;
        use TriangularNetworkLocation::*;

        for i in 0..3 {
            if p == self.tri_vert(start, SubIdx(i)) {
                return OnVertex(start, SubIdx(i));
            }
        }

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
            (ClockWise, _, _) => Outside(Edge::new(start, SubIdx(1))),
            (_, ClockWise, _) => Outside(Edge::new(start, SubIdx(2))),
            (_, _, ClockWise) => Outside(Edge::new(start, SubIdx(0))),
            // TODO: (Colinear, Colinear, _) -> exact point?
            (CoLinear, CounterClockWise, CounterClockWise) => OnEdge(Edge::new(start, SubIdx(1))),
            (CounterClockWise, CoLinear, CounterClockWise) => OnEdge(Edge::new(start, SubIdx(2))),
            (CounterClockWise, CounterClockWise, CoLinear) => OnEdge(Edge::new(start, SubIdx(0))),
            _ => panic!("{:?}", (d0, d1, d2)),
        }
    }

    pub fn locate_recursive(&self, p: &Point<T>) -> TriangularNetworkLocation {
        let mut start = TriIdx(0);

        loop {
            use TriangularNetworkLocation::*;

            start = match self.locate(start, p) {
                Outside(e) => {
                    match self.tri(e.tri).neighbor(e.sub) {
                        Some(idx) => idx,
                        None => {
                            eprintln!("{:?}, {:?}", e, self.tri(e.tri));
                            todo!();
                            // return Outside(idx, idx_neighbor);
                        }
                    }
                }
                e => return e,
            };
        }
    }

    // https://arxiv.org/abs/1403.3905
    pub fn visibility(&self, sx: &SimplicalChain<T>, p: &Point<T>) -> Vec<(Point<T>, Point<T>)> {
        use TriangularNetworkLocation::*;

        if let Err(e) = self.check_invariant("visibility") {
            eprintln!("{:?}", e);
        }

        match self.locate_recursive(p) {
            InTriangle(idx) => {
                let t = self.tri(idx);
                let mut out = Vec::new();
                for i in 0..3 {
                    let sub = SubIdx(i);
                    let edge = Edge::new(idx, sub);
                    let cw = t.vert(sub.cw());
                    let ccw = t.vert(sub);
                    if let Some(duel) = self.edge_duel(&edge) {
                        self.visibility_tri(
                            sx,
                            duel,
                            VisibilityQuery {
                                src: p.clone(),
                                cw,
                                ccw,
                            },
                            &mut out,
                        );
                    }
                }

                let mut segments = out
                    .into_iter()
                    .map(|s| {
                        let t = self.tri(s.edge.tri);
                        let p_start = self.vert(t.vert(s.edge.sub));
                        let p_end = self.vert(t.vert(s.edge.sub.cw()));

                        let l = Line::new_through(p_start, p_end);

                        let l_ccw = Line::new_through(p, self.vert(s.ccw));
                        let l_cw = Line::new_through(p, self.vert(s.cw));

                        let p_ccw = l.intersection_point(&l_ccw).unwrap();
                        let p_cw = l.intersection_point(&l_cw).unwrap();

                        (p_cw, p_ccw)
                    })
                    .collect::<Vec<_>>();

                segments.sort_by(|a, b| p.ccw_cmp_around(&a.0, &b.0));
                segments
            }
            _ => todo!(),
        }
    }

    pub fn centroid(&self, tri: TriIdx) -> Point<T> {
        let t = self.tri(tri);

        let [v0, v1, v2] = t.vertices;
        let p0 = self.vert(v0);
        let p1 = self.vert(v1);
        let p2 = self.vert(v2);

        pt_mean(&[p0, p1, p2])
    }

    fn visibility_tri(
        &self,
        sx: &SimplicalChain<T>,
        e: Edge,
        q: VisibilityQuery<T>,
        out: &mut Vec<VisibilitySegment>,
    ) {
        use Orientation::*;

        let p = &q.src;
        let p_cw = self.vert(q.cw);
        let p_ccw = self.vert(q.ccw);

        let t = self.tri(e.tri);
        let center = self.centroid(e.tri);

        if sx.characteristic(&center) != 1.0 {
            out.push(VisibilitySegment {
                edge: e,
                ccw: q.ccw,
                cw: q.cw,
            });
            return;
        }

        let v1 = t.vert(e.sub.ccw());

        let d_ccw_1 = Point::orient_along_direction(p, Direction::Through(p_ccw), self.vert(v1));
        let d_cw_1 = Point::orient_along_direction(p, Direction::Through(p_cw), self.vert(v1));

        // ccw side
        if d_ccw_1 != CounterClockWise {
            let ccw = q.ccw;
            let cw = match d_cw_1 {
                CoLinear => v1,
                ClockWise => q.cw,
                CounterClockWise => v1,
            };

            let e = Edge::new(e.tri, e.sub.cw());
            if let Some(duel) = self.edge_duel(&e) {
                self.visibility_tri(
                    sx,
                    duel,
                    VisibilityQuery {
                        src: q.src.clone(),
                        ccw,
                        cw,
                    },
                    out,
                );
            }
        }

        // cw side
        if d_cw_1 != ClockWise {
            let cw = q.cw;
            let ccw = match d_ccw_1 {
                CoLinear => v1,
                CounterClockWise => q.ccw,
                ClockWise => v1,
            };

            let e = Edge::new(e.tri, e.sub.ccw());
            if let Some(duel) = self.edge_duel(&e) {
                self.visibility_tri(
                    sx,
                    duel,
                    VisibilityQuery {
                        src: q.src.clone(),
                        ccw,
                        cw,
                    },
                    out,
                );
            }
        }
    }
}

struct VisibilitySegment {
    edge: Edge,
    ccw: VertIdx,
    cw: VertIdx,
}

struct VisibilityQuery<T> {
    src: Point<T>,
    ccw: VertIdx,
    cw: VertIdx,
}

#[derive(Clone)]
pub struct Triangle {
    // counterclockwise
    pub vertices: [VertIdx; 3],
    //  - self, neighbors[0], neighbors[1] meets at vertices[0], counterclockwise
    //  - vertices[0], vertices[1] meets with neighbors[1], ...
    pub neighbors: [Option<TriIdx>; 3],
}

impl std::fmt::Debug for Triangle {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(fmt, "Tri{{v=(")?;
        for (idx, v) in self.vertices.iter().enumerate() {
            let prefix = if idx == 0 { "" } else { ", " };
            write!(fmt, "{}{}", prefix, v.0)?;
        }
        write!(fmt, "), n=(")?;
        for (idx, n) in self.neighbors.iter().enumerate() {
            let prefix = if idx == 0 { "" } else { ", " };
            match n {
                Some(idx) => write!(fmt, "{}{}", prefix, idx.0)?,
                None => write!(fmt, "{}_", prefix)?,
            }
        }
        write!(fmt, ")}}")
    }
}

impl Triangle {
    pub fn is_super(&self) -> bool {
        let [v0, v1, v2] = self.vertices;
        is_super(v0) || is_super(v1) || is_super(v2)
    }

    pub fn vert(&self, idx: SubIdx) -> VertIdx {
        self.vertices[idx.0]
    }

    fn neighbor(&self, idx: SubIdx) -> Option<TriIdx> {
        self.neighbors[idx.0]
    }

    fn neighbor_mut(&mut self, idx: SubIdx) -> &mut Option<TriIdx> {
        &mut self.neighbors[idx.0]
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

    #[allow(unused)]
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
            (0.5, 0.0, OnEdge(Edge::new(TriIdx(0), SubIdx(1)))),
            (0.5, -0.1, Outside(Edge::new(TriIdx(0), SubIdx(1)))),
            (1.0, 0.5, OnEdge(Edge::new(TriIdx(0), SubIdx(2)))),
            (1.1, 0.5, Outside(Edge::new(TriIdx(0), SubIdx(2)))),
            (0.5, 0.5, OnEdge(Edge::new(TriIdx(0), SubIdx(0)))),
            (0.4, 0.6, Outside(Edge::new(TriIdx(0), SubIdx(0)))),
            (0.5, 0.1, InTriangle(TriIdx(0))),
        ];

        for (x, y, expected) in cases {
            assert_eq!(net.locate(TriIdx(0), &Point::new([x, y])), expected);
        }
    }
}
