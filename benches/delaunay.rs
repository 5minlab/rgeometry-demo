use criterion::{criterion_group, criterion_main, Criterion};
use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use rgeometry::data::*;
use rgeometry_playground::delaunay::{TriIdx, TriangularNetwork, VertIdx};
use rgeometry_playground::demo::{points_grid, points_uniform};

fn rand_point<R: Rng>(rng: &mut R) -> Point<f64> {
    Point::new([rng.gen_range(-100.0..100.0), rng.gen_range(-100.0..100.0)])
}

fn gen_delaunay(view: f64, points: &[Point<f64>]) -> TriangularNetwork<f64> {
    let v = view * 4.0;
    let mut t = TriangularNetwork::new(
        Point::new([-v, -v]),
        Point::new([v, -v]),
        Point::new([0.0, v]),
    );

    let mut r = std::usize::MAX;
    for p in points {
        if let Err(e) = t.insert(&p, &mut r) {
            eprintln!("{:?}", e);
            break;
        }
    }

    t
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let seed: <ChaCha20Rng as SeedableRng>::Seed = Default::default();
    let mut rng = ChaCha20Rng::from_seed(seed);

    let p0 = rand_point(&mut rng);
    let p1 = rand_point(&mut rng);
    let p2 = rand_point(&mut rng);
    let p3 = rand_point(&mut rng);

    let view = 100.0f64;

    for i in [20, 100, 1000] {
        let points = points_uniform(&mut rng, view, i);
        c.bench_function(&format!("dalaunay {i}"), |b| {
            b.iter(|| gen_delaunay(view, &points))
        });
    }

    for size in [10, 100, 1000] {
        let points = points_uniform(&mut rng, view, size);
        let net = gen_delaunay(view, &points);
        c.bench_function(&format!("TriangularNetwork::locate {size}"), |b| {
            b.iter(|| net.locate(TriIdx(0), &p0))
        });

        c.bench_function(
            &format!("TriangularNetwork::locate_recursive {size}"),
            |b| b.iter(|| net.locate_recursive(&p0)),
        );

        let idx0 = rng.gen_range(0..size) + 3;
        let idx1 = rng.gen_range(0..size) + 3;

        c.bench_function(&format!("TriangularNetwork::cut {size}"), |b| {
            b.iter(|| net.cut(VertIdx(idx0), VertIdx(idx1)))
        });

        c.bench_function(&format!("TriangularNetwork::clone"), |b| {
            b.iter(|| net.clone())
        });
        let cut = net.cut(VertIdx(idx0), VertIdx(idx1));
        c.bench_function(&format!("TriangularNetwork::cut_apply {size}"), |b| {
            b.iter(|| net.clone().cut_apply(&cut))
        });
    }

    for i in [3, 10, 30] {
        let points = points_grid(view, i);
        c.bench_function(&format!("dalaunay rect {i}x{i}"), |b| {
            b.iter(|| gen_delaunay(view, &points))
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
