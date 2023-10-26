use rgeometry::{data::*, Orientation, PolygonScalar};

pub struct VisibilityResult<T> {
    pub origin: Point<T>,
    pub pairs: Vec<(Point<T>, Point<T>)>,
}

impl<T> VisibilityResult<T> {
    pub fn empty(origin: Point<T>) -> Self {
        Self {
            origin,
            pairs: vec![],
        }
    }
}

impl<T: PolygonScalar> VisibilityResult<T> {
    pub fn visible_segment<'a, 'b>(
        &'b self,
        p: &'a Point<T>,
    ) -> Option<(&'b Point<T>, &'b Point<T>)> {
        use Orientation::*;

        let mut directions = Vec::with_capacity(self.pairs.len());
        for (p0, _p1) in &self.pairs {
            let d = Point::orient_along_direction(&self.origin, Direction::Through(p0), p);
            directions.push(d);
        }
        directions.push(directions[0]);

        for i in 0..self.pairs.len() {
            let (ref p0, ref p1) = &self.pairs[i];

            match (directions[i], directions[i + 1]) {
                (CounterClockWise, ClockWise) => {
                    return Some((p0, p1));
                }
                _ => continue,
            }
        }
        None
    }

    pub fn point_visible(&self, p: &Point<T>) -> bool {
        use Orientation::*;
        if self.pairs.len() == 0 {
            return false;
        }

        match self.visible_segment(p) {
            Some((p0, p1)) => {
                Point::orient_along_direction(p0, Direction::Through(p1), p) == CounterClockWise
            }
            None => false,
        }
    }

    pub fn raycast(&self, p: &Point<T>) -> Option<Point<T>> {
        match self.visible_segment(p) {
            Some((p0, p1)) => {
                let l0 = Line::new_through(&self.origin, p);
                let l1 = Line::new_through(p0, p1);

                return l0.intersection_point(&l1);
            }
            None => None,
        }
    }
}
