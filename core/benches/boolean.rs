use core::{
    boolean::SimplicalChain, build_net, gen_rects, points_circular, points_cube_subdivide,
    raster::raster, Rect,
};
use criterion::{criterion_group, criterion_main, Criterion};
use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use rgeometry::data::*;

pub fn criterion_benchmark(c: &mut Criterion) {
    let seed: <ChaCha20Rng as SeedableRng>::Seed = Default::default();
    let mut rng = ChaCha20Rng::from_seed(seed);

    // complex visibility
    {
        let view = 30.0;
        let mut rects = gen_rects(&mut rng, view, 100);
        rects.push(Rect::new(view / 5.0, view / 5.0));

        let mut sx = SimplicalChain::default();
        for r in rects {
            let p = r.polygon(1);
            let sx_r = SimplicalChain::from_polygon(&p);
            sx = sx.union(&sx_r);
        }

        c.bench_function("build_net", |b| b.iter(|| build_net(view, &sx, true)));

        c.bench_function("circle intersection", |b| {
            let points_circle = points_circular(10.0, 32);
            let p_circle = Polygon::new(points_circle).unwrap();
            let sx_circle = SimplicalChain::from_polygon(&p_circle);

            b.iter(|| sx.intersect(&sx_circle))
        });

        let (net, constraints) = build_net(view, &sx, true);

        let p0 = Point::new([0.0, 0.0]);
        c.bench_function("TriangularNetwork::visibility", |b| {
            b.iter(|| net.visibility(&constraints, &p0))
        });

        c.bench_function("TriangularNetwork::visibility circle", |b| {
            let points_circle = points_circular(10.0, 32);
            let p_circle = Polygon::new(points_circle).unwrap();
            let sx_circle = SimplicalChain::from_polygon(&p_circle);

            b.iter(|| {
                let sx = sx.intersect(&sx_circle);
                let (net, constraints) = build_net(view, &sx, true);
                net.visibility(&constraints, &p0)
            });
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

        c.bench_function("TriangularNetwork::visibility raster dist", |b| {
            b.iter(|| {
                let mut sum = 0.0;
                let limit = 200.0;
                for (p1, p2) in &vis {
                    raster(1, &[p0, *p1, *p2], |x, y| {
                        if x * x + y * y < limit {
                            sum += x + y;
                        }
                    });
                }
            })
        });

        c.bench_function("SimplicalChain::characteristic", |b| {
            b.iter(|| sx.characteristic(&p0))
        });
    }

    for subdivide in [1, 10, 100] {
        let p0 = points_cube_subdivide(Point::new([0.0, 0.0]), 100.0, subdivide);
        let p1 = points_cube_subdivide(Point::new([50.0, 50.0]), 100.0, subdivide);
        let p0 = Polygon::new(p0).unwrap();
        let p1 = Polygon::new(p1).unwrap();

        let s0 = SimplicalChain::from_polygon(&p0);
        let s1 = SimplicalChain::from_polygon(&p1);

        c.bench_function(
            &format!("SimplicalChain::bool_intersect {subdivide}"),
            |b| b.iter(|| s0.intersect(&s1)),
        );
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
