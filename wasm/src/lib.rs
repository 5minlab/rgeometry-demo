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

    pub fn visibility(&self, x: f64, y: f64) -> js_sys::Float32Array {
        let vis = self.net.visibility(&self.constraints, &Point::new([x, y]));

        visibility_serialize(&vis)
    }

    pub fn visibility_limit(&self, x: f64, y: f64, limit: f64) -> js_sys::Float32Array {
        let p = Point::new([x, y]);
        let mut vis = self.net.visibility(&self.constraints, &p);
        core::visibility_limit(&mut vis, &p, limit);
        visibility_serialize(&vis)
    }

    pub fn visibility_limit_fill(
        &self,
        x: f64,
        y: f64,
        limit: f64,
        gridinfo: &[f64],
        grid: js_sys::Uint8Array,
    ) -> VisibilityResult {
        let p = Point::new([x, y]);
        let mut vis = self.net.visibility(&self.constraints, &p);
        core::visibility_limit(&mut vis, &p, limit);
        let mut count = 0;

        if let [minx, miny, width, height, gridsize] = gridinfo {
            let grid_count_x = (width / gridsize) as usize;
            let grid_count_y = (height / gridsize) as usize;
            let l = grid_count_x * grid_count_y;
            let mut v = Vec::with_capacity(l);
            v.resize(l, 0.0f32);

            // find bounding box
            let mut aabb: Option<AABB<f64>> = None;
            for (p0, p1) in &vis {
                let aabb_other =
                    core::raster::raster_bounds(*gridsize as usize, &[p, p1.clone(), p0.clone()]);
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
                core::raster::raster(*gridsize as usize, &[p, p1, p0], |x, y| {
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
                })
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

    pub fn connectivity(&self, coords: &[f64]) -> js_sys::Uint16Array {
        let mut v = Vec::new();

        let mut points = Vec::new();
        for i in 0..coords.len() / 2 {
            let x = coords[i * 2];
            let y = coords[i * 2 + 1];
            points.push(Point::new([x, y]));
        }

        for i in 0..points.len() {
            let p0 = &points[i];
            let vis = self.net.visibility(&self.constraints, p0);
            if vis.is_empty() {
                continue;
            }

            for j in (i + 1)..points.len() {
                let p1 = &points[j];

                if core::visibility::point_visible(p0, &vis, p1) {
                    v.push(i as u16);
                    v.push(j as u16);
                }
                //
            }
        }

        js_sys::Uint16Array::from(&v[..])
    }

    pub fn raycast(&self, origin_x: f64, origin_y: f64, coords: &[f64]) -> js_sys::Float32Array {
        let mut v = Vec::with_capacity(coords.len() * 2);

        let p0 = Point::new([origin_x, origin_y]);

        let vis = self.net.visibility(&self.constraints, &p0);
        if vis.is_empty() {
            return js_sys::Float32Array::from(&v[..]);
        }

        for i in 0..coords.len() / 2 {
            let x = coords[i * 2] + origin_x;
            let y = coords[i * 2 + 1] + origin_y;
            let p1 = Point::new([x, y]);

            let p2 = core::visibility::raycast(&p0, &vis, &p1).unwrap_or(p0);
            v.push(p2.array[0] as f32);
            v.push(p2.array[1] as f32);
        }

        js_sys::Float32Array::from(&v[..])
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

fn visibility_serialize(vis: &[(Point<f64>, Point<f64>)]) -> js_sys::Float32Array {
    let mut v = Vec::with_capacity(vis.len() * 4);
    for (from, to) in vis {
        v.push(from.array[0] as f32);
        v.push(from.array[1] as f32);
        v.push(to.array[0] as f32);
        v.push(to.array[1] as f32);
    }
    js_sys::Float32Array::from(&v[..])
}
