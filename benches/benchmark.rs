use criterion::{criterion_group, criterion_main, Criterion};
use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use rgeometry::data::*;
use rgeometry_playground::raster::raster;

fn rand_point<R: Rng>(rng: &mut R) -> Point<f64> {
    Point::new([rng.gen_range(-100.0..100.0), rng.gen_range(-100.0..100.0)])
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let seed: <ChaCha20Rng as SeedableRng>::Seed = Default::default();
    let mut rng = ChaCha20Rng::from_seed(seed);

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

    for size in [1.0, 10.0, 100.0] {
        let p0 = Point::new([-size, -size]);
        let p1 = Point::new([size, -size]);
        let p2 = Point::new([0.0, size]);

        c.bench_function(&format!("raster {size}"), |b| {
            b.iter(|| {
                let mut sum = 0.0;
                raster(1, &[p0, p1, p2], |x, y| {
                    sum += x + y;
                });
            });
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
