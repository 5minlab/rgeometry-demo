use rgeometry::{
    data::{Direction, Point},
    Orientation, PolygonScalar,
};

pub fn point_visible<T: PolygonScalar>(
    origin: &Point<T>,
    vis: &[(Point<T>, Point<T>)],
    p: &Point<T>,
) -> bool {
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
                return Point::orient_along_direction(p0, Direction::Through(p1), p)
                    == CounterClockWise
            }
            _ => continue,
        }
    }
    false
}
