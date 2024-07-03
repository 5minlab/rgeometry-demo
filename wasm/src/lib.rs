use core::{
    aabb::AABB,
    boolean::SimplicalChain,
    delaunay::{TriangularNetwork, VertIdx},
    *,
};
use rgeometry::data::*;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

#[wasm_bindgen]
#[derive(Clone)]
pub struct Simplical {
    sx: SimplicalChain<f64>,
}

#[wasm_bindgen]
pub enum FillOp {
    Add,
    Subtract,
}

#[wasm_bindgen]
impl Simplical {
    pub fn new() -> Self {
        Self {
            sx: SimplicalChain::default(),
        }
    }

    pub fn dup(&self) -> Self {
        self.clone()
    }

    pub fn from_rect(
        x: f64,
        y: f64,
        extent_x: f64,
        extent_y: f64,
        rot: f64,
        subdivide: usize,
    ) -> Self {
        let r = Rect::new(extent_x, extent_y).pos(x, y).rot(rot);
        let p = r.polygon(subdivide);
        let sx = SimplicalChain::from_polygon(&p);

        Self { sx }
    }

    pub fn from_rect_seglen(
        x: f64,
        y: f64,
        extent_x: f64,
        extent_y: f64,
        rot: f64,
        seglen: f64,
    ) -> Self {
        let r = Rect::new(extent_x, extent_y).pos(x, y).rot(rot);
        let p = r.polygon_seglen(seglen);
        let sx = SimplicalChain::from_polygon(&p);

        Self { sx }
    }

    #[wasm_bindgen]
    pub fn simplices(&self) -> js_sys::Float32Array {
        let mut v = Vec::with_capacity(self.sx.simplices.len() * 4);
        for s in &self.sx.simplices {
            v.push(s.src.array[0] as f32);
            v.push(s.src.array[1] as f32);
            v.push(s.dst.array[0] as f32);
            v.push(s.dst.array[1] as f32);
        }
        js_sys::Float32Array::from(&v[..])
    }

    #[wasm_bindgen]
    pub fn union(&self, other: &Simplical) -> Self {
        Self {
            sx: self.sx.union(&other.sx),
        }
    }

    #[wasm_bindgen]
    pub fn intersect(&self, other: &Simplical) -> Self {
        Self {
            sx: self.sx.intersect(&other.sx),
        }
    }

    #[wasm_bindgen]
    pub fn subtract(&self, other: &Simplical) -> Self {
        Self {
            sx: self.sx.subtract(&other.sx),
        }
    }
}

#[wasm_bindgen]
pub struct Delaunay {
    net: TriangularNetwork<f64>,
}

#[wasm_bindgen]
impl Delaunay {
    pub fn from(coords: &[f64]) -> Self {
        let mut view = 0f64;
        let mut points = Vec::with_capacity(coords.len() / 2);
        for i in 0..coords.len() / 2 {
            let x = coords[i * 2];
            let y = coords[i * 2 + 1];
            view = view.max(x.max(y));
            points.push(Point::new([x, y]));
        }

        let v = view * 4.0;
        let mut net = TriangularNetwork::new(
            Point::new([-v, -v]),
            Point::new([v, -v]),
            Point::new([0.0, v]),
        );

        let mut reductions = std::usize::MAX;
        let mut indices = indexmap::IndexSet::new();
        for p in &points {
            let idx = net.insert(p, &mut reductions).unwrap();
            assert!(indices.insert(idx));
        }

        Self { net }
    }

    pub fn neighbors(&self) -> js_sys::Uint16Array {
        let mut v = Vec::with_capacity(self.net.triangles.len() * 6);
        for tri in &self.net.triangles {
            for i in 0..3 {
                let v0 = tri.vertices[i];
                let v1 = tri.vertices[(i + 1) % 3];
                if v0.is_super() || v1.is_super() {
                    continue;
                }

                v.push(v0.0 as u16 - 3);
                v.push(v1.0 as u16 - 3);
            }
        }
        js_sys::Uint16Array::from(&v[..])
    }
}

#[wasm_bindgen]
pub struct Triangulated {
    net: TriangularNetwork<f64>,
    constraints: Vec<(VertIdx, VertIdx)>,
}

#[wasm_bindgen]
impl Triangulated {
    pub fn from(sim: &Simplical) -> Self {
        let (net, constraints) = build_net(10000.0, &sim.sx, true);
        let triangulated = Triangulated { net, constraints };

        triangulated
    }

    pub fn visibility(&self, x: f64, y: f64, out_to_in: bool) -> Visibility {
        let origin = Point::new([x, y]);
        let vis = self
            .net
            .visibility_dir(&self.constraints, &origin, out_to_in)
            .unwrap_or(visibility::VisibilityResult::empty(origin.clone()));
        Visibility { vis }
    }

    pub fn connectivity(&self, coords: &[f64]) -> js_sys::Uint16Array {
        let mut v = Vec::new();

        let mut points = Vec::with_capacity(coords.len() / 2);
        for i in 0..coords.len() / 2 {
            let x = coords[i * 2];
            let y = coords[i * 2 + 1];
            points.push(Point::new([x, y]));
        }

        for i in 0..points.len() {
            let p0 = &points[i];
            let vis = self.net.visibility(&self.constraints, p0);
            let vis = match vis {
                Some(v) => v,
                None => continue,
            };

            for j in (i + 1)..points.len() {
                let p1 = &points[j];

                if vis.point_visible(p1) {
                    v.push(i as u16);
                    v.push(j as u16);
                }
            }
        }

        js_sys::Uint16Array::from(&v[..])
    }
}

#[wasm_bindgen]
pub struct Visibility {
    vis: visibility::VisibilityResult<f64>,
}

#[wasm_bindgen]
impl Visibility {
    pub fn limit(&mut self, limit: f64) {
        core::visibility_limit(&mut self.vis, limit);
    }

    pub fn clip(&mut self, points: &[f64]) {
        if let [p0x, p0y, p1x, p1y] = points {
            let p0 = Point::new([*p0x, *p0y]);
            let p1 = Point::new([*p1x, *p1y]);

            self.vis = self.vis.clip(p0, p1);
        }
    }

    pub fn serialize(&self) -> js_sys::Float32Array {
        let mut v = Vec::with_capacity(self.vis.pairs.len() * 4);
        for (from, to) in &self.vis.pairs {
            v.push(from.array[0] as f32);
            v.push(from.array[1] as f32);
            v.push(to.array[0] as f32);
            v.push(to.array[1] as f32);
        }
        js_sys::Float32Array::from(&v[..])
    }

    pub fn fill_op(&self, gridinfo: &[f64], grid: &mut [f32], op: FillOp, amount: f32) -> usize {
        let vis = &self.vis.pairs;
        let mut count = 0;

        if let [minx, miny, width, _, gridsize] = gridinfo {
            let grid_count_x = (width / gridsize) as usize;

            for (p0, p1) in vis {
                core::raster::raster(
                    *gridsize as usize,
                    &[self.vis.origin.clone(), p1.clone(), p0.clone()],
                    |x, y| {
                        let x = x - minx / gridsize;
                        let y = y - miny / gridsize;
                        if x < 0.0 || y < 0.0 {
                            return;
                        }
                        count += 1;
                        let x = x as usize;
                        let y = y as usize;
                        let idx = y * grid_count_x + x;
                        if idx >= grid.len() {
                            return;
                        }
                        let prev = grid[idx];

                        match op {
                            FillOp::Add => {
                                grid[idx] = (prev + amount).min(1.0);
                            }
                            FillOp::Subtract => {
                                grid[idx] = (prev - amount).max(0.0);
                            }
                        }
                    },
                )
            }
        }
        count
    }

    pub fn fill(&self, gridinfo: &[f64], grid: js_sys::Uint8Array) -> VisibilityResult {
        let vis = &self.vis.pairs;
        let mut count = 0;

        if let [minx, miny, width, height, gridsize] = gridinfo {
            let grid_count_x = (width / gridsize) as usize;
            let grid_count_y = (height / gridsize) as usize;

            // find bounding box
            let mut aabb: Option<AABB<f64>> = None;
            for (p0, p1) in vis {
                let aabb_other = core::raster::raster_bounds(
                    *gridsize as usize,
                    &[self.vis.origin.clone(), p1.clone(), p0.clone()],
                );
                aabb = match aabb {
                    Some(mut aabb) => {
                        aabb.extend_aabb(&aabb_other);
                        Some(aabb)
                    }
                    None => Some(aabb_other),
                }
            }

            let aabb = aabb.unwrap_or(AABB::new(&Point::new([*minx, *miny])));

            for (p0, p1) in vis {
                core::raster::raster(
                    *gridsize as usize,
                    &[self.vis.origin.clone(), p1.clone(), p0.clone()],
                    |x, y| {
                        let x = x - minx / gridsize;
                        let y = y - miny / gridsize;
                        if x < 0.0 || y < 0.0 {
                            return;
                        }
                        count += 1;
                        let x = x as usize;
                        let y = y as usize;
                        let idx = (y * grid_count_x + x) as u32;
                        grid.set_index(idx, 1);
                    },
                )
            }

            VisibilityResult {
                min_x: (aabb.min.array[0] - minx / gridsize).max(0.0),
                min_y: (aabb.min.array[1] - miny / gridsize).max(0.0),
                max_x: (aabb.max.array[0] - minx / gridsize).min(grid_count_x as f64),
                max_y: (aabb.max.array[1] - miny / gridsize).min(grid_count_y as f64),
                count,
            }
        } else {
            todo!();
        }
    }

    pub fn visible(&self, coords: &[f64]) -> js_sys::Uint8Array {
        let mut v = Vec::with_capacity(coords.len() * 2);

        for i in 0..coords.len() / 2 {
            if self.vis.pairs.is_empty() {
                v.push(0);
                continue;
            }

            let x = coords[i * 2];
            let y = coords[i * 2 + 1];
            let p1 = Point::new([x, y]);

            v.push(if self.vis.point_visible(&p1) { 1 } else { 0 });
        }

        js_sys::Uint8Array::from(&v[..])
    }

    pub fn raycast(&self, coords: &[f64]) -> js_sys::Float32Array {
        let mut v = Vec::with_capacity(coords.len() * 2);

        let p0 = &self.vis.origin;
        let [ox, oy] = p0.array;

        for i in 0..coords.len() / 2 {
            if self.vis.pairs.is_empty() {
                v.push(ox as f32);
                v.push(oy as f32);
                continue;
            }

            let x = coords[i * 2] + ox;
            let y = coords[i * 2 + 1] + oy;
            let p1 = Point::new([x, y]);

            let p2 = self.vis.raycast(&p1).unwrap_or_else(|| p0.clone());
            v.push(p2.array[0] as f32);
            v.push(p2.array[1] as f32);
        }

        js_sys::Float32Array::from(&v[..])
    }
}

#[wasm_bindgen]
pub struct Connectivity {
    points: Vec<Point<f64>>,
    connects: Vec<(usize, usize)>,
}

#[wasm_bindgen]
impl Connectivity {
    pub fn from(t: &Triangulated, coords: &[f64]) -> Self {
        let mut points = Vec::with_capacity(coords.len() / 2);
        for i in 0..coords.len() / 2 {
            let x = coords[i * 2];
            let y = coords[i * 2 + 1];
            points.push(Point::new([x, y]));
        }

        let mut connects = Vec::new();
        for i in 0..points.len() {
            let p0 = &points[i];
            let vis = t.net.visibility(&t.constraints, p0);
            let vis = match vis {
                Some(v) => v,
                None => continue,
            };

            for j in (i + 1)..points.len() {
                let p1 = &points[j];

                if vis.point_visible(p1) {
                    connects.push((i, j));
                    connects.push((j, i));
                }
            }
        }

        connects.sort();

        Self { points, connects }
    }

    pub fn connectivity(&self) -> js_sys::Uint16Array {
        let mut v = Vec::with_capacity(self.connects.len());

        for (i, j) in &self.connects {
            if j > i {
                continue;
            }

            v.push(*i as u16);
            v.push(*j as u16);
        }

        js_sys::Uint16Array::from(&v[..])
    }

    pub fn reachables(&self, t: &Triangulated, x: f64, y: f64) -> js_sys::Uint16Array {
        let mut v = Vec::new();
        let p0 = Point::new([x, y]);

        let vis = match t.net.visibility(&t.constraints, &p0) {
            Some(v) => v,
            None => {
                return js_sys::Uint16Array::from(&v[..]);
            }
        };

        if vis.pairs.is_empty() {
            return js_sys::Uint16Array::from(&v[..]);
        }

        for i in 0..self.points.len() {
            let p1 = &self.points[i];
            if vis.point_visible(p1) {
                v.push(i as u16);
            }
        }

        js_sys::Uint16Array::from(&v[..])
    }
}

#[wasm_bindgen]
pub struct VisibilityResult {
    pub min_x: f64,
    pub min_y: f64,
    pub max_x: f64,
    pub max_y: f64,

    pub count: usize,
}
