use core::{
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

    #[wasm_bindgen]
    pub fn simplices(&self) -> js_sys::Float64Array {
        let mut v = Vec::with_capacity(self.sx.simplices.len() * 4);
        for s in &self.sx.simplices {
            v.push(s.src.array[0]);
            v.push(s.src.array[1]);
            v.push(s.dst.array[0]);
            v.push(s.dst.array[1]);
        }
        js_sys::Float64Array::from(&v[..])
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

    pub fn visibility(&self, x: f64, y: f64) -> js_sys::Float64Array {
        let vis = self.net.visibility(&self.constraints, &Point::new([x, y]));
        let mut v = Vec::with_capacity(vis.len() * 4);
        for (from, to) in &vis {
            v.push(from.array[0]);
            v.push(from.array[1]);
            v.push(to.array[0]);
            v.push(to.array[1]);
        }
        js_sys::Float64Array::from(&v[..])
    }
}
