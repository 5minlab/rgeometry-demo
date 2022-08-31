use criterion::{criterion_group, criterion_main, Criterion};
use rgeometry::data::*;
use rgeometry_playground::boolean::SimplicalChain;
use rgeometry_playground::demo::points_cube_subdivide;

pub fn criterion_benchmark(c: &mut Criterion) {
    for subdivide in [1, 10, 100] {
        let p0 = points_cube_subdivide(Point::new([0.0, 0.0]), 100.0, subdivide);
        let p1 = points_cube_subdivide(Point::new([50.0, 50.0]), 100.0, subdivide);
        let p0 = Polygon::new(p0).unwrap();
        let p1 = Polygon::new(p1).unwrap();

        let s0 = SimplicalChain::from_polygon(&p0);
        let s1 = SimplicalChain::from_polygon(&p1);

        c.bench_function(
            &format!("SimplicalChain::bool_intersect {subdivide}"),
            |b| b.iter(|| s0.bool_intersect(&s1)),
        );
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
