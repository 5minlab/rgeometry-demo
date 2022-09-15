use rgeometry::data::{Direction, EndPoint, Line, LineSegment, Point};
use rgeometry::{Intersects, Orientation};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
enum SweepEventKind {
    Start(usize),
    End(usize),
    Intersection(usize, usize),
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SweepEvent {
    p: Point<f64>,
    kind: SweepEventKind,
}

pub struct Intersections<'a> {
    lines: &'a [(Point<f64>, Point<f64>)],

    // TODO: binary search tree
    sweepline: Vec<usize>,
    events: BinaryHeap<Reverse<SweepEvent>>,
}

// https://www.ics.uci.edu/~goodrich/teach/geom/notes/LineSweep.pdf
impl<'a> Intersections<'a> {
    pub fn new(lines: &'a [(Point<f64>, Point<f64>)]) -> Self {
        let mut events = BinaryHeap::new();

        for i in 0..lines.len() {
            events.push(Reverse(SweepEvent {
                p: lines[i].0,
                kind: SweepEventKind::Start(i),
            }));
            events.push(Reverse(SweepEvent {
                p: lines[i].1,
                kind: SweepEventKind::End(i),
            }));
        }

        Self {
            lines,
            sweepline: Vec::new(),
            events,
        }
    }

    pub fn find_pos(&self, line_idx: usize) -> usize {
        let l0 = self.lines[line_idx];
        for i in 0..self.sweepline.len() {
            let other_line_idx = self.sweepline[i];
            let l1 = self.lines[other_line_idx];
            match Point::orient_along_direction(&l1.0, Direction::Through(&l1.1), &l0.0) {
                Orientation::CounterClockWise => {
                    return i;
                }
                _ => continue,
            }
        }
        self.sweepline.len()
    }

    fn intersects(&self, idx0: usize, idx1: usize, sweep: f64) -> Option<Point<f64>> {
        let (p0, p1) = self.lines[idx0];
        let (p2, p3) = self.lines[idx1];
        let l0 = LineSegment::new(EndPoint::Exclusive(p0), EndPoint::Exclusive(p1));
        let l1 = LineSegment::new(EndPoint::Exclusive(p2), EndPoint::Exclusive(p3));

        if l0.intersect(&l1).is_none() {
            return None;
        }

        let l0 = Line::new_through(&p0, &p1);
        let l1 = Line::new_through(&p2, &p3);

        match l0.intersection_point(&l1) {
            Some(v) => {
                if v[0] > sweep {
                    Some(v)
                } else {
                    None
                }
            }
            None => None,
        }
    }

    fn add_intersection_up(&mut self, pos: usize, sweep: f64) {
        let line_idx = self.sweepline[pos];
        if let Some(next_line_idx) = self.sweepline.get(pos + 1) {
            if let Some(p) = self.intersects(line_idx, *next_line_idx, sweep) {
                self.events.push(Reverse(SweepEvent {
                    p,
                    kind: SweepEventKind::Intersection(line_idx, *next_line_idx),
                }));
            }
        }
    }

    fn add_intersection_down(&mut self, pos: usize, sweep: f64) {
        let line_idx = self.sweepline[pos];
        if pos > 0 {
            let prev_line_idx = self.sweepline[pos - 1];
            if let Some(p) = self.intersects(line_idx, prev_line_idx, sweep) {
                self.events.push(Reverse(SweepEvent {
                    p,
                    kind: SweepEventKind::Intersection(prev_line_idx, line_idx),
                }));
            }
        }
    }

    fn add_intersections(&mut self, pos: usize, sweep: f64) {
        self.add_intersection_up(pos, sweep);
        self.add_intersection_down(pos, sweep);
    }

    fn insert(&mut self, line_idx: usize, pos: usize, sweep: f64) {
        self.sweepline.insert(pos, line_idx);
        self.add_intersections(pos, sweep);
    }

    pub fn sweep(&mut self) -> HashMap<(usize, usize), Point<f64>> {
        let mut v = HashMap::new();
        while let Some(Reverse(ev)) = self.events.pop() {
            let sweep = ev.p.array[0];
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
}
