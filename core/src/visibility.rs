use rgeometry::{data::*, Orientation, PolygonScalar};

#[derive(Clone)]
pub struct VisibilityResult<T> {
    pub origin: Point<T>,
    pub pairs: Vec<(Point<T>, Point<T>)>,
    pub arc: bool,
}

impl<T> VisibilityResult<T> {
    pub fn empty(origin: Point<T>) -> Self {
        Self {
            origin,
            pairs: vec![],
            arc: false,
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

    pub fn clip(&self, d0: Direction<T, 2>, d1: Direction<T, 2>) -> Self {
        use Orientation::*;
        if self.pairs.is_empty() {
            return self.clone();
        }

        let mut pairs = vec![];
        let mut i = 0;
        while i < self.pairs.len() {
            let (ref p0, ref p1) = self.pairs[i];
            i += 1;
            let o0 = Point::orient_along_direction(&self.origin, d0, p0);
            let o1 = Point::orient_along_direction(&self.origin, d0, p1);

            match (o0, o1) {
                (ClockWise, CounterClockWise) => {
                    let l0 = Line::new(&self.origin, d0);
                    let l1 = Line::new_through(p0, p1);

                    if let Some(p0) = l0.intersection_point(&l1) {
                        pairs.push((p0, p1.clone()));
                    }
                    break;
                }
                (CoLinear, CounterClockWise) => {
                    pairs.push((p0.clone(), p1.clone()));
                }
                (ClockWise, CoLinear) => {
                    break;
                }
                _ => {
                    continue;
                }
            }
        }

        i = i % self.pairs.len();
        let end = (i + self.pairs.len() + 1) % self.pairs.len();
        while i != end {
            let pair = self.pairs[i].clone();
            pairs.push(pair);
            i = (i + 1) % self.pairs.len();
        }

        let mut pairs1 = vec![];

        for (p0, p1) in pairs {
            let o0 = Point::orient_along_direction(&self.origin, d1, &p0);
            let o1 = Point::orient_along_direction(&self.origin, d1, &p1);

            match (o0, o1) {
                (ClockWise, CounterClockWise) => {
                    let l0 = Line::new(&self.origin, d1);
                    let l1 = Line::new_through(&p0, &p1);

                    if let Some(p1) = l0.intersection_point(&l1) {
                        pairs1.push((p0, p1));
                    }
                    break;
                }
                (CoLinear, CounterClockWise) => {
                    break;
                }
                (ClockWise, CoLinear) => {
                    pairs1.push((p0, p1));
                    break;
                }
                _ => {
                    pairs1.push((p0, p1));
                }
            }
        }

        VisibilityResult {
            origin: self.origin.clone(),
            pairs: pairs1,
            arc: true,
        }
    }
}
