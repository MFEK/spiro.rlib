use crate::BezierContext;
use std::io::{Result, Write};

pub type PostScriptBezierContext = BezierContext<PostScriptEmitter<Box<dyn Write>>>;
pub struct PostScriptEmitter<W: Write>(pub W);

impl<W> Write for PostScriptEmitter<W> where W: Write {
    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        self.0.write(buf)
    }
    fn flush(&mut self) -> Result<()> {
        self.0.flush()
    }
}

pub fn move_to(w: &mut PostScriptBezierContext, x: f64, y: f64, _is_open: bool) {
    writeln!(w.data, "{} {} moveto", x, y).unwrap();
}

pub fn line_to(w: &mut PostScriptBezierContext, x: f64, y: f64) {
    writeln!(w.data, "{} {} lineto", x, y).unwrap();
}

pub fn curve_to(w: &mut PostScriptBezierContext, x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) {
    writeln!(w.data, "{} {} {} {} {} {} curveto", x1, y1, x2, y2, x3, y3).unwrap();
}

pub fn mark_knot(w: &mut PostScriptBezierContext, knot_idx: usize) {
    writeln!(w.data, "% {} knot", knot_idx).unwrap();
}

impl PostScriptBezierContext {
    pub fn new<W: Write + 'static>(w: W) -> Self {
        Self {
            move_fn: move_to,
            line_fn: line_to,
            curve_fn: curve_to,
            mark_knot_fn: mark_knot,
            data: PostScriptEmitter(Box::new(w)),
        }
    }
}
