use crate::delaunay::SubIdx;
use rgeometry::data::Point;

pub fn rester<F>(verts: &[Point<f64>; 3], mut f: F)
where
    F: FnMut(f64, f64),
{
    let [p0, p1, p2] = verts;

    let idx_topmost = if p0.array[1] < p1.array[1] {
        if p1.array[1] < p2.array[1] {
            2
        } else {
            1
        }
    } else {
        if p0.array[1] < p2.array[1] {
            2
        } else {
            0
        }
    };

    // ccw
    let mut left = SubIdx(idx_topmost).ccw();
    let mut right = SubIdx(idx_topmost).cw();

    fn midcoord(v: f64) -> f64 {
        (v + 0.5).floor() - 0.5
    }

    fn roundx(v: f64) -> f64 {
        (v + 0.5).ceil() - 0.5
    }

    fn interp(p0: &Point<f64>, p1: &Point<f64>, y: f64) -> f64 {
        let [p0x, p0y] = p0.array;
        let [p1x, p1y] = p1.array;
        if p0y == p1y {
            return p0x;
        }
        let t = (y - p0y) / (p1y - p0y);
        return (p1x - p0x) * t + p0x;
    }

    let mut y = midcoord(verts[idx_topmost].array[1]);
    loop {
        if y <= verts[left.0].array[1] {
            if verts[left.ccw().0].array[1] < y {
                left = left.ccw();
            } else {
                break;
            }
        }
        if y <= verts[right.0].array[1] {
            if verts[right.cw().0].array[1] < y {
                right = right.cw();
            } else {
                break;
            }
        }

        let b = verts[left.0].array[1].max(verts[right.0].array[1]);

        while y > b {
            let xl = interp(&verts[right.ccw().0], &verts[right.0], y);
            let xr = interp(&verts[left.cw().0], &verts[left.0], y);

            // round to next 0.5
            let mut x = roundx(xl);
            while x < xr {
                f(x, y);
                x += 1.0;
            }

            y -= 1.0;
        }
    }

    //
}
