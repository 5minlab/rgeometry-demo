//! Demo app for egui

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

// When compiling natively:
fn main() {
    let options = eframe::NativeOptions::default();
    let _ = eframe::run_native(
        "rgeometry playground",
        options,
        Box::new(|_cc| Ok(Box::new(gui::MyApp::default()))),
    );
}
