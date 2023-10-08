use super::{p_rg_to_egui, plot_line, pt_egui, Demo};
use core::boolean::*;
use core::{gen_rects, points_circular, Rect};
use eframe::egui::{self, epaint::Color32, Ui};
use egui_plot::{self, *};
use rgeometry::data::{Point, Polygon};

#[derive(Clone, Copy, PartialEq, Eq)]
enum CircleMode {
    Disabled,
    Union,
    Intersect,
    Subtract0,
    Subtract1,
}

pub struct DemoBoolean2 {
    opt_render_rect: bool,
    opt_render_union: bool,

    opt_circle_mode: CircleMode,

    circle: Vec<Point<f64>>,

    view: f64,

    count: usize,
    subdivide: usize,
    rects: Vec<Rect>,
    sx: SimplicalChain<f64>,

    rational: bool,
}

fn rect_union(rects: &[Rect], subdivide: usize, rational: bool) -> SimplicalChain<f64> {
    if !rational {
        let mut sx = SimplicalChain::default();
        for r in rects {
            let p = r.polygon(subdivide);
            let sx_r = SimplicalChain::from_polygon(&p);
            sx = sx.union(&sx_r);
        }
        sx
    } else {
        use num::{FromPrimitive, ToPrimitive};

        let mut sx = SimplicalChain::default();
        for r in rects {
            let p = r.polygon(subdivide);
            let points = p
                .iter_boundary()
                .map(|p| p.map(|v| num::BigRational::from_f64(v).unwrap()))
                .collect::<Vec<_>>();
            let p = Polygon::new_unchecked(points);
            let sx_r = SimplicalChain::from_polygon(&p);
            sx = sx.union(&sx_r);
        }

        SimplicalChain {
            simplices: sx
                .simplices
                .into_iter()
                .map(|s| Simplex {
                    src: s.src.map(|v| v.to_f64().unwrap()),
                    dst: s.dst.map(|v| v.to_f64().unwrap()),
                })
                .collect(),
        }
    }
}

impl DemoBoolean2 {
    #[allow(unused)]
    pub fn new(view: f64) -> Self {
        let mut rng = rand::thread_rng();
        let rects = {
            let mut rects = Vec::new();
            for [x, y, w, h, heading] in &[
                [
                    -123.54193216787064,
                    -52.91399622244809,
                    7.5,
                    2.5,
                    -0.4070386786719689,
                ],
                [
                    -17.937753007736486,
                    -98.44154111468106,
                    82.5,
                    2.5,
                    -0.4070386786719689,
                ],
                [
                    42.874641733514906,
                    93.13573479724495,
                    102.5,
                    2.5,
                    -0.4070386786719689,
                ],
                [
                    -88.54424874680433,
                    40.89518325993956,
                    2.5,
                    102.5,
                    -0.4070386786719689,
                ],
                [
                    95.11519327082028,
                    -38.283155683074284,
                    2.5,
                    102.5,
                    -0.4070386786719689,
                ],
                [
                    9.611163294907207,
                    -110.31829195613312,
                    52.5,
                    2.5,
                    -0.4070386786719689,
                ],
                [
                    49.200332766414135,
                    -18.488570947320827,
                    52.5,
                    2.5,
                    -0.4070386786719689,
                ],
                [
                    -34.3242387359236,
                    -85.93222116993906,
                    2.5,
                    7.5,
                    -0.4070386786719689,
                ],
                [
                    -8.591278579444097,
                    -26.242902514211057,
                    2.5,
                    32.5,
                    -0.4070386786719689,
                ],
                [
                    75.32060853506682,
                    -84.19801618748045,
                    2.5,
                    52.5,
                    -0.4070386786719689,
                ],
                [
                    49.200332766414135,
                    -18.488570947320827,
                    52.5,
                    2.5,
                    -0.4070386786719689,
                ],
                [
                    88.78950223792106,
                    73.34115006149148,
                    52.5,
                    2.5,
                    -0.4070386786719689,
                ],
                [
                    5.264930735583327,
                    5.89749983887325,
                    2.5,
                    7.5,
                    -0.4070386786719689,
                ],
                [
                    30.997890892062827,
                    65.58681849460125,
                    2.5,
                    32.5,
                    -0.4070386786719689,
                ],
                [
                    114.90977800657376,
                    7.631704821331866,
                    2.5,
                    52.5,
                    -0.4070386786719689,
                ],
                [
                    -105.17598796610817,
                    -60.83183011674947,
                    17.5,
                    2.5,
                    -0.4070386786719689,
                ],
                [
                    -26.406404841622212,
                    -67.56627696817658,
                    2.5,
                    17.5,
                    -0.4070386786719689,
                ],
                [
                    13.182764629884712,
                    24.263444040635708,
                    2.5,
                    17.5,
                    -0.4070386786719689,
                ],
                [141.74915755470903, 100.0, 5.0, 5.0, -0.616137379020819],
            ] {
                rects.push(Rect::new(*w, *h).pos(*x, *y).rot(*heading));
            }
            rects
        };
        let rational = false;
        let subdivide = 1;
        let sx = rect_union(&rects, subdivide, rational);

        Self {
            opt_render_rect: true,
            opt_render_union: true,

            opt_circle_mode: CircleMode::Disabled,

            circle: points_circular(view / 2.0, 32),

            view,

            count: rects.len(),
            subdivide,
            rects,
            sx,

            rational,
        }
    }
}

impl Demo for DemoBoolean2 {
    fn name(&self) -> &'static str {
        "boolean2"
    }

    fn ui(&mut self, _t: f64, _ctx: &egui::Context, ui: &mut Ui) {
        let mut rng = rand::thread_rng();
        ui.horizontal(|ui| {
            if ui.button("regenerate").clicked() {
                self.rects = gen_rects(&mut rng, self.view, 100);
            }
            ui.separator();

            for (mode, label) in &[
                (CircleMode::Disabled, "none"),
                (CircleMode::Union, "union"),
                (CircleMode::Intersect, "intersect"),
                (CircleMode::Subtract0, "subtract0"),
                (CircleMode::Subtract1, "subtract1"),
            ] {
                ui.radio_value(&mut self.opt_circle_mode, *mode, *label);
            }
        });

        ui.horizontal(|ui| {
            ui.checkbox(&mut self.opt_render_rect, "render rect");
            ui.checkbox(&mut self.opt_render_union, "render union");

            ui.checkbox(&mut self.rational, "rational");
            ui.separator();

            ui.add(egui::Slider::new(&mut self.count, 0..=self.rects.len()).text("counts"));
            ui.add(egui::Slider::new(&mut self.subdivide, 1..=10).text("subdivide"));
        });

        self.sx = rect_union(&self.rects[..self.count], self.subdivide, self.rational);

        if self.opt_circle_mode == CircleMode::Disabled {
            return;
        }

        let p = Polygon::new(self.circle.clone()).unwrap();
        let sx_circle = SimplicalChain::from_polygon(&p);

        match self.opt_circle_mode {
            CircleMode::Union => {
                self.sx = self.sx.union(&sx_circle);
            }
            CircleMode::Intersect => {
                self.sx = self.sx.intersect(&sx_circle);
            }
            CircleMode::Subtract0 => {
                self.sx = self.sx.subtract(&sx_circle);
            }
            CircleMode::Subtract1 => {
                self.sx = sx_circle.subtract(&self.sx);
            }
            _ => (),
        }
    }

    fn plot_ui(&self, plot_ui: &mut PlotUi) {
        if self.opt_render_rect {
            for r in &self.rects {
                let p = r.polygon(1);
                let pe = p_rg_to_egui(&p);

                plot_ui.polygon(pe.fill_color(Color32::RED));
            }
        }

        if self.opt_render_union {
            for s in &self.sx.simplices {
                plot_line(plot_ui, &[&s.src, &s.dst], Color32::GREEN);

                plot_ui.points(
                    egui_plot::Points::new(PlotPoints::Owned(vec![
                        pt_egui(&s.src),
                        pt_egui(&s.dst),
                    ]))
                    .color(Color32::WHITE)
                    .radius(2.0),
                );
            }
        }
    }
}
