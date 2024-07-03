use rgeometry::{data::*, Orientation, PolygonScalar};

enum ClipTy<T> {
    CounterClockWise(Point<T>, Point<T>),
    In(Point<T>, Point<T>),
    ClockWise(Point<T>, Point<T>),
}

impl<T> std::fmt::Debug for ClipTy<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ClipTy::CounterClockWise(_p0, _p1) => write!(f, "CCW"),
            ClipTy::In(_p0, _p1) => write!(f, "IN"),
            ClipTy::ClockWise(_p0, _p1) => write!(f, "CW"),
        }
    }
}

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

    pub fn clip(&self, d0p: Point<T, 2>, d1p: Point<T, 2>) -> Self {
        use Orientation::*;
        if self.pairs.is_empty() {
            return self.clone();
        }

        let d0 = rgeometry::data::Direction::Through(&d0p);
        let d1 = rgeometry::data::Direction::Through(&d1p);

        let dir = Point::orient_along_direction(&self.origin, d0, &d1p) == CounterClockWise;

        let mut start = 0;
        let mut pairs = vec![];
        for (p0, p1) in &self.pairs {
            let d0p0 = Point::orient_along_direction(&self.origin, d0, p0);
            let d1p0 = Point::orient_along_direction(&self.origin, d1, p0);
            let d0p1 = Point::orient_along_direction(&self.origin, d0, p1);
            let d1p1 = Point::orient_along_direction(&self.origin, d1, p1);

            let p0_in = (d0p0 == CounterClockWise && d1p0 == ClockWise) ^ !dir;
            let p1_in = (d0p1 == CounterClockWise && d1p1 == ClockWise) ^ !dir;

            if !p0_in && !p1_in {
                continue;
            }

            if p0_in && p1_in {
                pairs.push(ClipTy::In(p0.clone(), p1.clone()));
            } else if p0_in {
                let l0 = Line::new(&self.origin, d1);
                let l1 = Line::new_through(p0, p1);

                if let Some(p1) = l0.intersection_point(&l1) {
                    pairs.push(ClipTy::ClockWise(p0.clone(), p1));
                } else {
                    pairs.push(ClipTy::ClockWise(p0.clone(), p1.clone()));
                }
            } else if p1_in {
                start = pairs.len();

                let l0 = Line::new(&self.origin, d0);
                let l1 = Line::new_through(p0, p1);

                if let Some(p0) = l0.intersection_point(&l1) {
                    pairs.push(ClipTy::CounterClockWise(p0, p1.clone()));
                } else {
                    pairs.push(ClipTy::CounterClockWise(p0.clone(), p1.clone()));
                }
            }
        }
        pairs.rotate_left(start);

        VisibilityResult {
            origin: self.origin.clone(),
            pairs: pairs
                .into_iter()
                .map(|p| match p {
                    ClipTy::CounterClockWise(p0, p1) => (p0, p1),
                    ClipTy::In(p0, p1) => (p0, p1),
                    ClipTy::ClockWise(p0, p1) => (p0, p1),
                })
                .collect(),
            arc: true,
        }
    }
}
