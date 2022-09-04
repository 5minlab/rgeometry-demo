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
pub fn gen_points_cube(center: &[f64], extent: f64) -> js_sys::Float64Array {
    let points = points_cube(Point::new([center[0], center[1]]), extent);
    let mut v = Vec::with_capacity(points.len() * 2);
    for i in 0..points.len() {
        v.push(points[i].array[0]);
        v.push(points[i].array[1]);
    }
    js_sys::Float64Array::from(&v[..])
}

#[wasm_bindgen]
pub fn greet(name: &str) {
    log(&format!("Hello, {}!", name));
}

#[wasm_bindgen]
pub fn simplical_new() -> usize {
    let s = Box::new(SimplicalChain::<f64>::default());
    let ptr = Box::into_raw(s);
    ptr as usize
}

#[wasm_bindgen]
pub fn simplical_from_polygon(points: &[f64]) -> usize {
    let mut v = Vec::with_capacity(points.len() / 2);
    for i in 0..points.len() / 2 {
        let x = points[i * 2];
        let y = points[i * 2 + 1];
        let p = Point::new([x, y]);
        v.push(p);
    }
    let polygon = Polygon::new(v).unwrap();
    let s = Box::new(SimplicalChain::from_polygon(&polygon));

    Box::into_raw(s) as usize
}

#[wasm_bindgen]
pub fn simplical_simplices(ptr: usize) -> js_sys::Float64Array {
    let sx = from_ptr::<SimplicalChain<f64>>(ptr);

    let mut v = Vec::with_capacity(sx.simplices.len() * 4);
    for s in &sx.simplices {
        v.push(s.src.array[0]);
        v.push(s.src.array[1]);
    }
    std::mem::forget(sx);
    js_sys::Float64Array::from(&v[..])
}

#[wasm_bindgen]
pub fn simplical_free(ptr: usize) {
    let sx = from_ptr::<SimplicalChain<f64>>(ptr);
    std::mem::drop(sx);
}

#[wasm_bindgen]
pub fn simplical_union(ptr0: usize, ptr1: usize) -> usize {
    let sx0 = from_ptr::<SimplicalChain<f64>>(ptr0);
    let sx1 = from_ptr::<SimplicalChain<f64>>(ptr1);

    let sx = Box::new(sx0.bool_union(&sx1));

    std::mem::forget(sx0);
    std::mem::forget(sx1);

    Box::into_raw(sx) as usize
}

#[wasm_bindgen]
pub struct Triangulated {
    net: TriangularNetwork<f64>,
    constraints: Vec<(VertIdx, VertIdx)>,
}

#[wasm_bindgen]
impl Triangulated {
    pub fn from(ptr: usize) -> Self {
        let sx = from_ptr::<SimplicalChain<f64>>(ptr);
        let (net, constraints) = build_net(100.0, &sx, true);
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

fn from_ptr<T>(ptr: usize) -> Box<T> {
    unsafe { Box::from_raw(ptr as *mut T) }
}
