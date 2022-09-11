use core::delaunay::{TriIdx, TriangularNetwork};

use eframe::{
    egui::{
        self,
        plot::{self, Legend, Plot, PlotPoint, PlotPoints, PlotUi},
        Key,
    },
    epaint::Color32,
};
use rgeometry::data::{Point, Polygon};

mod boolean;
mod boolean2;
mod boolean_tri;
mod cdt;
mod delaunay;
mod gridnet;
mod raster;

pub fn plot_line(plot_ui: &mut PlotUi, points: &[&Point<f64>], color: Color32) {
    let e_points = points.iter().map(|p| pt_egui(p)).collect();
    plot_ui.line(plot::Line::new(PlotPoints::Owned(e_points)).color(color));
}

pub fn plot_polygon(plot_ui: &mut PlotUi, p: &Polygon<f64>, color: Color32) {
    for edge in p.iter_boundary_edges() {
        let e_points = vec![pt_egui(edge.src), pt_egui(edge.dst)];
        plot_ui.line(plot::Line::new(PlotPoints::Owned(e_points)).color(color));
    }
}

fn plot_net(net: &TriangularNetwork<f64>, plot_ui: &mut PlotUi, render_supertri: bool) {
    for (_t_idx, t) in net.triangles.iter().enumerate() {
        let [v0, v1, v2] = t.vertices;
        let p0 = net.vert(v0);
        let p1 = net.vert(v1);
        let p2 = net.vert(v2);
        if !render_supertri && t.is_super() {
            continue;
        }

        plot_line(plot_ui, &[p0, p1, p2, p0], Color32::GREEN);

        if true {
            let center = net.centroid(TriIdx(_t_idx));
            let label = format!("{:?}={:?}", _t_idx, t);

            plot_ui.text(plot::Text::new(pt_egui(&center), label));
        }
    }
}

fn p_rg_to_egui(p: &Polygon<f64>) -> plot::Polygon {
    let points: Vec<[f64; 2]> = p
        .iter()
        .map(|p| [p.array[0] as f64, p.array[1] as f64])
        .collect::<Vec<_>>();
    plot::Polygon::new(points)
}

pub fn pt_egui(p: &Point<f64>) -> PlotPoint {
    let [x, y] = p.array;
    PlotPoint::new(x, y)
}

trait Demo {
    fn name(&self) -> &'static str;
    fn ui(&mut self, t: f64, ctx: &egui::Context, ui: &mut egui::Ui);
    fn plot_ui(&self, plot_ui: &mut PlotUi);
}

pub struct MyApp {
    pause: bool,
    t: f64,

    view: f64,

    selected: &'static str,
    demos: Vec<Box<dyn Demo>>,
}

impl Default for MyApp {
    fn default() -> Self {
        let view = 30f64;

        let demos: Vec<Box<dyn Demo>> = vec![
            Box::new(boolean_tri::DemoBooleanTri::new(view)),
            Box::new(boolean::DemoBoolean::new(view)),
            Box::new(boolean2::DemoBoolean2::new(view)),
            Box::new(cdt::DemoCDT::new(view)),
            Box::new(delaunay::DemoDelaunay::new(view)),
            Box::new(gridnet::DemoGridNet::new(view)),
            Box::new(raster::DemoRaster::new(view)),
        ];
        let selected = demos[0].name();

        Self {
            pause: false,
            t: 0.0,

            view,

            selected,
            demos,
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if ctx.input().key_pressed(Key::Space) {
            self.pause = !self.pause;
        }
        if ctx.input().key_pressed(Key::R) {
            self.t = 0.0;
        }
        if ctx.input().key_pressed(Key::A) {
            self.t -= 1.0 / 60.0;
        }
        if ctx.input().key_pressed(Key::S) {
            self.t += 1.0 / 60.0;
        }

        if !self.pause {
            self.t += ctx.input().predicted_dt as f64;
            ctx.request_repaint();
        }

        let view = self.view;

        let t = self.t;

        egui::TopBottomPanel::top("top").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("time={:.2}", t));
                ui.separator();
                ui.checkbox(&mut self.pause, "pause");
                ui.separator();
                ui.label("shortcuts: (Space) pause/resume | (R) Reset time | (S) single-step");
            });

            ui.horizontal(|ui| {
                ui.visuals_mut().button_frame = false;
                for d in &self.demos {
                    let name = d.name();
                    if ui.selectable_label(self.selected == name, name).clicked() {
                        self.selected = name;
                    }
                    ui.separator();
                }
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            for demo in &mut self.demos {
                if demo.name() == self.selected {
                    demo.ui(t, ctx, ui);
                }
            }

            let plot = Plot::new(0)
                .legend(Legend::default())
                .data_aspect(1.0)
                /*
                .center_x_axis(true)
                .center_y_axis(true)
                .include_x(-view)
                .include_x(view)
                .include_y(-view)
                .include_y(view)
                */
                .allow_zoom(true)
                .allow_drag(true);

            plot.show(ui, |plot_ui| {
                for demo in &self.demos {
                    if demo.name() == self.selected {
                        demo.plot_ui(plot_ui);
                    }
                }

                let area = view;
                plot_ui.points(
                    plot::Points::new(PlotPoints::Owned(vec![
                        PlotPoint::new(-area, -area),
                        PlotPoint::new(area, -area),
                        PlotPoint::new(-area, area),
                        PlotPoint::new(area, area),
                    ]))
                    .color(Color32::WHITE),
                );
            });
        });
    }
}

#[cfg(target_arch = "wasm32")]
use eframe::web::AppRunnerRef;

#[cfg(target_arch = "wasm32")]
use eframe::wasm_bindgen::{self, prelude::*};

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub struct WebHandle {
    handle: AppRunnerRef,
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl WebHandle {
    #[wasm_bindgen]
    #[cfg(target_arch = "wasm32")]
    pub fn stop_web(&self) -> Result<(), wasm_bindgen::JsValue> {
        let mut app = self.handle.lock();
        let res = app.destroy();

        // let numw = Arc::weak_count(&app);
        // let nums = Arc::strong_count(&app);
        // tracing::debug!("runner ref {:?}, {:?}", numw, nums);

        res
    }
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn init_wasm_hooks() {
    // Make sure panics are logged using `console.error`.
    console_error_panic_hook::set_once();

    // Redirect tracing to console.log and friends:
    tracing_wasm::set_as_global_default();
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn start_separate(canvas_id: &str) -> Result<WebHandle, wasm_bindgen::JsValue> {
    let web_options = eframe::WebOptions::default();
    let handle = eframe::start_web(
        canvas_id,
        web_options,
        Box::new(|_cc| Box::new(MyApp::default())),
    )
    .map(|handle| WebHandle { handle });

    handle
}

/// This is the entry-point for all the web-assembly.
/// This is called once from the HTML.
/// It loads the app, installs some callbacks, then returns.
/// You can add more callbacks like this if you want to call in to your code.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn start(canvas_id: &str) -> Result<WebHandle, wasm_bindgen::JsValue> {
    init_wasm_hooks();
    start_separate(canvas_id)
}
