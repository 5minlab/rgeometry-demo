use rgeometry::{data::*, Orientation, PolygonScalar};

pub fn visible_segment<'a, 'b, T: PolygonScalar>(
    origin: &'a Point<T>,
    vis: &'b [(Point<T>, Point<T>)],
    p: &'a Point<T>,
) -> Option<(&'b Point<T>, &'b Point<T>)> {
    use Orientation::*;

    let mut directions = Vec::with_capacity(vis.len());
    for (p0, _p1) in vis {
        let d = Point::orient_along_direction(&origin, Direction::Through(p0), p);
        directions.push(d);
    }
    directions.push(directions[0]);

    for i in 0..vis.len() {
        let (ref p0, ref p1) = &vis[i];

        match (directions[i], directions[i + 1]) {
            (CounterClockWise, ClockWise) => {
                return Some((p0, p1));
            }
            _ => continue,
        }
    }
    None
}

pub fn point_visible<T: PolygonScalar>(
    origin: &Point<T>,
    vis: &[(Point<T>, Point<T>)],
    p: &Point<T>,
) -> bool {
    use Orientation::*;

    match visible_segment(origin, vis, p) {
        Some((p0, p1)) => {
            Point::orient_along_direction(p0, Direction::Through(p1), p) == CounterClockWise
        }
        None => false,
    }
}

pub fn raycast<T: PolygonScalar>(
    origin: &Point<T>,
    vis: &[(Point<T>, Point<T>)],
    p: &Point<T>,
) -> Option<Point<T>> {
    match visible_segment(origin, vis, p) {
        Some((p0, p1)) => {
            let l0 = Line::new_through(origin, p);
            let l1 = Line::new_through(p0, p1);

            return l0.intersection_point(&l1);
        }
        None => None,
    }
}
