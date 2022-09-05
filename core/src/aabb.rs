use rgeometry::{data::Point, PolygonScalar};

pub struct AABB<T> {
    pub min: Point<T>,
    pub max: Point<T>,
}

fn segment_intersect<T: PolygonScalar>(a: &T, b: &T, c: &T, d: &T) -> bool {
    a < d && b > c
}

impl<T: PolygonScalar + Copy> AABB<T> {
    pub fn new(p: &Point<T>) -> Self {
        Self { min: *p, max: *p }
    }

    pub fn extend(&mut self, p: &Point<T>) {
        if p.array[0] < self.min.array[0] {
            self.min.array[0] = p.array[0];
        } else if p.array[0] > self.max.array[0] {
            self.max.array[0] = p.array[0];
        }

        if p.array[1] < self.min.array[1] {
            self.min.array[1] = p.array[1];
        } else if p.array[1] > self.max.array[1] {
            self.max.array[1] = p.array[1];
        }
    }

    pub fn interects(&self, other: &Self) -> bool {
        segment_intersect(
            &self.min.array[0],
            &self.max.array[0],
            &other.min.array[0],
            &other.max.array[0],
        ) || segment_intersect(
            &self.min.array[1],
            &self.max.array[1],
            &other.min.array[1],
            &other.max.array[1],
        )
    }
}
