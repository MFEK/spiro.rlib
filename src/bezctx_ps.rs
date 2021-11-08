use crate::BezierContext;

pub fn move_to(_: &mut BezierContext<()>, x: f64, y: f64, _is_open: bool) {
    println!("{} {} moveto", x, y);
}

pub fn line_to(_: &mut BezierContext<()>, x: f64, y: f64) {
    println!("{} {} lineto", x, y);
}

pub fn curve_to(_: &mut BezierContext<()>, x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) {
    println!("{} {} {} {} {} {} curveto", x1, y1, x2, y2, x3, y3);
}

pub fn mark_knot(_: &mut BezierContext<()>, knot_idx: usize) {
    println!("% {} knot", knot_idx);
}

pub static mut POSTSCRIPT_BEZCTX: BezierContext<()> = BezierContext { move_fn: move_to, line_fn: line_to, curve_fn: curve_to, mark_knot_fn: mark_knot, data: None };
