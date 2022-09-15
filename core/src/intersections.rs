use rgeometry::data::{Direction, EndPoint, Line, LineSegment, Point};
use rgeometry::{Intersects, Orientation, PolygonScalar};
use std::collections::{BinaryHeap, HashMap};

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
enum SweepEventKind {
    Start(usize),
    End(usize),
    Intersection(usize, usize),
}

#[derive(Debug)]
struct SweepEvent<T: PolygonScalar> {
    p: Point<T>,
    kind: SweepEventKind,
}

impl<T: PolygonScalar> std::cmp::PartialEq for SweepEvent<T> {
    fn eq(&self, other: &Self) -> bool {
        self.p == other.p
    }
}

impl<T: PolygonScalar> std::cmp::Eq for SweepEvent<T> {}

impl<T: PolygonScalar> std::cmp::PartialOrd for SweepEvent<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(other.p.cmp(&self.p))
    }
}

impl<T: PolygonScalar> std::cmp::Ord for SweepEvent<T> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.p.cmp(&self.p).then(self.kind.cmp(&other.kind))
    }
}

pub struct Intersections<'a, T: PolygonScalar> {
    lines: &'a [(Point<T>, Point<T>)],

    // TODO: binary search tree
    sweepline: Vec<usize>,
    events: BinaryHeap<SweepEvent<T>>,
}

// https://www.ics.uci.edu/~goodrich/teach/geom/notes/LineSweep.pdf
impl<'a, T> Intersections<'a, T>
where
    T: PolygonScalar,
{
    pub fn new(lines: &'a [(Point<T>, Point<T>)]) -> Self {
        let mut events = BinaryHeap::new();

        for i in 0..lines.len() {
            events.push(SweepEvent {
                p: lines[i].0.clone(),
                kind: SweepEventKind::Start(i),
            });
            events.push(SweepEvent {
                p: lines[i].1.clone(),
                kind: SweepEventKind::End(i),
            });
        }

        Self {
            lines,
            sweepline: Vec::new(),
            events,
        }
    }

    pub fn find_pos(&self, line_idx: usize) -> usize {
        use ordslice::Ext;

        let l0 = &self.lines[line_idx];
        self.sweepline.upper_bound_by(|i| {
            let l1 = &self.lines[*i];
            match Point::orient_along_direction(&l1.0, Direction::Through(&l1.1), &l0.0) {
                Orientation::ClockWise => std::cmp::Ordering::Less,
                Orientation::CounterClockWise => std::cmp::Ordering::Greater,
                Orientation::CoLinear => {
                    match Point::orient_along_direction(&l1.0, Direction::Through(&l1.1), &l0.1) {
                        Orientation::ClockWise => std::cmp::Ordering::Less,
                        Orientation::CounterClockWise => std::cmp::Ordering::Greater,
                        _ => todo!(),
                    }
                }
            }
        })
    }

    fn intersects(&self, idx0: usize, idx1: usize, sweep: &T) -> Option<Point<T>> {
        let (ref p0, ref p1) = self.lines[idx0];
        let (ref p2, ref p3) = self.lines[idx1];
        let l0 = LineSegment::new(
            EndPoint::Exclusive(p0.clone()),
            EndPoint::Exclusive(p1.clone()),
        );
        let l1 = LineSegment::new(
            EndPoint::Exclusive(p2.clone()),
            EndPoint::Exclusive(p3.clone()),
        );

        if l0.intersect(&l1).is_none() {
            return None;
        }

        let l0 = Line::new_through(p0, p1);
        let l1 = Line::new_through(p2, p3);

        match l0.intersection_point(&l1) {
            Some(v) => {
                if &v[0] > sweep {
                    Some(v)
                } else {
                    None
                }
            }
            None => None,
        }
    }

    fn add_intersection_up(&mut self, pos: usize, sweep: &T) {
        let line_idx = self.sweepline[pos];
        if let Some(next_line_idx) = self.sweepline.get(pos + 1) {
            if let Some(p) = self.intersects(line_idx, *next_line_idx, sweep) {
                self.events.push(SweepEvent {
                    p,
                    kind: SweepEventKind::Intersection(line_idx, *next_line_idx),
                });
            }
        }
    }

    fn add_intersection_down(&mut self, pos: usize, sweep: &T) {
        let line_idx = self.sweepline[pos];
        if pos > 0 {
            let prev_line_idx = self.sweepline[pos - 1];
            if let Some(p) = self.intersects(line_idx, prev_line_idx, sweep) {
                self.events.push(SweepEvent {
                    p,
                    kind: SweepEventKind::Intersection(prev_line_idx, line_idx),
                });
            }
        }
    }

    fn add_intersections(&mut self, pos: usize, sweep: &T) {
        self.add_intersection_up(pos, sweep);
        self.add_intersection_down(pos, sweep);
    }

    fn insert(&mut self, line_idx: usize, pos: usize, sweep: &T) {
        self.sweepline.insert(pos, line_idx);
        self.add_intersections(pos, sweep);
    }

    pub fn sweep(&mut self) -> HashMap<(usize, usize), Point<T>> {
        let mut v = HashMap::new();
        while let Some(ev) = self.events.pop() {
            let sweep = &ev.p.array[0];
            match ev.kind {
                SweepEventKind::Start(i) => {
                    let pos = self.find_pos(i);
                    self.insert(i, pos, sweep);
                }
                SweepEventKind::End(i) => {
                    let pos = self.sweepline.iter().position(|&k| k == i).unwrap();
                    self.sweepline.remove(pos);
                    if pos > 0 {
                        self.add_intersection_up(pos - 1, sweep);
                    }
                }
                SweepEventKind::Intersection(i, j) => {
                    if v.contains_key(&(i, j)) {
                        continue;
                    }
                    let pos0 = self.sweepline.iter().position(|&k| k == i).unwrap();
                    let pos1 = pos0 + 1;
                    assert_eq!(self.sweepline[pos0], i);
                    assert_eq!(self.sweepline[pos1], j);

                    self.sweepline[pos0] = j;
                    self.sweepline[pos1] = i;

                    self.add_intersection_down(pos0, sweep);
                    self.add_intersection_up(pos1, sweep);

                    v.insert((i, j), ev.p);
                }
            }
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::points_uniform;

    #[test]
    fn it_works() {
        let mut rng = rand::thread_rng();
        let p0 = points_uniform(&mut rng, 100.0, 10);
        let p1 = points_uniform(&mut rng, 100.0, 10);
        let lines = (0..p0.len())
            .into_iter()
            .map(|i| match p0[i].cmp(&p1[i]) {
                std::cmp::Ordering::Less => (p0[i], p1[i]),
                _ => (p1[i], p0[i]),
            })
            .collect::<Vec<_>>();

        let mut s = Intersections::new(&lines);
        s.sweep();
    }

    #[test]
    fn intersections_touching_1() {
        let lines = vec![
            (Point::new([0.0, 0.0]), Point::new([1.0, 0.0])),
            (Point::new([1.0, 0.0]), Point::new([2.0, 0.0])),
        ];

        let mut s = Intersections::new(&lines);
        assert_eq!(s.sweep().len(), 0);
    }

    #[test]
    fn intersections_touching_2() {
        let lines = vec![
            (Point::new([1.0, 0.0]), Point::new([2.0, 0.0])),
            (Point::new([0.0, 0.0]), Point::new([1.0, 0.0])),
        ];

        let mut s = Intersections::new(&lines);
        assert_eq!(s.sweep().len(), 0);
    }
}
