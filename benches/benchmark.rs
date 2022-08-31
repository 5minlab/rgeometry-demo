use criterion::{criterion_group, criterion_main, Criterion};
use rand::{thread_rng, Rng};
use rgeometry::data::*;
use rgeometry_playground::delaunay::TriangularNetwork;
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
    let mut rng = thread_rng();
    let p0 = rand_point(&mut rng);
    let p1 = rand_point(&mut rng);
    let p2 = rand_point(&mut rng);
    let p3 = rand_point(&mut rng);

    c.bench_function("Point::orient_along_direction", |b| {
        b.iter(|| Point::orient_along_direction(&p0, Direction::Through(&p1), &p2))
    });

    let l0 = Line::new_through(&p0, &p1);
    let l1 = Line::new_through(&p2, &p3);

    c.bench_function("Line::intersection_point", |b| {
        b.iter(|| l0.intersection_point(&l1))
    });

    let view = 100.0f64;

    for i in [20, 100, 1000] {
        let points = points_uniform(view, i);
        c.bench_function(&format!("dalaunay {i}"), |b| {
            b.iter(|| gen_delaunay(view, &points))
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
