use core::intersections::Intersections;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use rgeometry::data::*;

fn rand_point<R: Rng>(rng: &mut R) -> Point<f64> {
    Point::new([rng.gen_range(-100.0..100.0), rng.gen_range(-100.0..100.0)])
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let seed: <ChaCha20Rng as SeedableRng>::Seed = Default::default();
    let mut rng = ChaCha20Rng::from_seed(seed);

    let mut lines = Vec::new();
    for _ in 0..10000 {
        let p0 = rand_point(&mut rng);
        let p1 = rand_point(&mut rng);
        if p0.array[0] < p1.array[0] {
            lines.push((p0, p1));
        } else {
            lines.push((p1, p0));
        }
    }

    for i in [100, 1000, 10000] {
        c.bench_function(&format!("Intersections {i}"), |b| {
            b.iter(|| {
                let mut intersections = Intersections::new(&lines[0..i]);
                intersections.sweep();
            });
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
