#[cfg(test)]
mod tests {
    use rand::prelude::*;
    use rand_chacha::*;
    use rgeometry::data::Point;
    use rgeometry_playground::boolean::SimplicalChain;
    use rgeometry_playground::demo::*;

    #[test]
    fn crash() {
        let view = 30f64;
        let len = 100;
        let seed: <ChaCha20Rng as SeedableRng>::Seed = [
            113, 193, 174, 249, 150, 225, 41, 118, 51, 61, 4, 220, 253, 36, 74, 134, 116, 13, 121,
            47, 138, 146, 27, 35, 230, 121, 125, 102, 189, 183, 119, 46,
        ];
        let t = unsafe { std::mem::transmute::<u64, f64>(4616977748265533440) };
        let mut rng = ChaCha20Rng::from_seed(seed);

        let mut rects = gen_rects(&mut rng, view, len);

        for r in &mut rects {
            r.rot = t;
        }

        let mut sx = SimplicalChain::default();
        for r in rects {
            let p = r.polygon();
            let sx_r = SimplicalChain::from_polygon(&p);
            sx = sx.bool_union(&sx_r);
        }

        let (net, c) = build_net(view, &sx, true);

        let _vis = net.visibility(&c, &Point::new([0.0, 0.0]));
    }
}
