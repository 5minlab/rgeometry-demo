use core::{
    boolean::SimplicalChain, build_net, gen_rects, points_cube_subdivide, raster::raster, Rect,
};
use criterion::{criterion_group, criterion_main, Criterion};
use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use rgeometry::data::*;

pub fn criterion_benchmark(c: &mut Criterion) {
    let seed: <ChaCha20Rng as SeedableRng>::Seed = Default::default();
    let mut rng = ChaCha20Rng::from_seed(seed);

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

    // complex visibility
    {
        let view = 30.0;
        let mut rects = gen_rects(&mut rng, view, 100);
        rects.push(Rect::new(view / 5.0, view / 5.0));

        let mut sx = SimplicalChain::default();
        for r in rects {
            let p = r.polygon();
            let sx_r = SimplicalChain::from_polygon(&p);
            sx = sx.bool_union(&sx_r);
        }

        c.bench_function("build_net", |b| b.iter(|| build_net(view, &sx, true)));

        let (net, constraints) = build_net(view, &sx, true);

        let p0 = Point::new([0.0, 0.0]);
        c.bench_function("TriangularNetwork::visibility", |b| {
            b.iter(|| net.visibility(&constraints, &p0))
        });

        let vis = net.visibility(&constraints, &p0);
        c.bench_function("TriangularNetwork::visibility raster", |b| {
            b.iter(|| {
                let mut sum = 0.0;
                for (p1, p2) in &vis {
                    raster(1, &[p0, *p1, *p2], |x, y| {
                        sum += x + y;
                    });
                }
            })
        });

        c.bench_function("SimplicalChain::characteristic", |b| {
            b.iter(|| sx.characteristic(&p0))
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
