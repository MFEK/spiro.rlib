//! Provides a [`BezierContext`] that builds up a [`Vec<Operation>`].
//!
//! Note: If compiled with `default` (`log`) feature, [`log::trace`] is active.

use crate::BezierContext;
/// Create a list of Bézier operations suitable for SVG etc
#[cfg(feature = "log")]
use log::trace;
#[cfg(not(feature = "log"))]
macro_rules! trace {
    ($($arg:tt)*) => {} // noop
}

/// A Bézier pen operation.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Operation {
    /// `(x, y, is_open)`
    MoveTo(f64, f64, bool),
    /// `(x, y)`
    LineTo(f64, f64),
    /// `(x1, y1, x2, y2, x3, y3)`
    CurveTo(f64, f64, f64, f64, f64, f64),
    /// `(knot_idx)`
    MarkKnot(usize),
}

pub fn move_to(ctx: &mut BezierContext<Vec<Operation>>, x: f64, y: f64, is_open: bool) {
    ctx.data.as_mut().map(|v| v.push(Operation::MoveTo(x, y, is_open)));
    trace!("Spiro callback: M {}, {} ", x, y);
}
pub fn line_to(ctx: &mut BezierContext<Vec<Operation>>, x: f64, y: f64) {
    ctx.data.as_mut().map(|v| v.push(Operation::LineTo(x, y)));
    trace!("Spiro callback: L {}, {} ", x, y);
}
pub fn curve_to(ctx: &mut BezierContext<Vec<Operation>>, x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) {
    ctx.data.as_mut().map(|v| v.push(Operation::CurveTo(x1, y1, x2, y2, x3, y3)));
    trace!("Spiro callback: C {}, {}, {}, {}, {}, {} ", x1, y1, x2, y2, x3, y3);
}
pub fn mark_knot(ctx: &mut BezierContext<Vec<Operation>>, knot_idx: usize) {
    ctx.data.as_mut().map(|v| v.push(Operation::MarkKnot(knot_idx)));
    trace!("Spiro callback: KNOT {}", knot_idx);
}

impl BezierContext<Vec<Operation>> {
    pub fn new() -> Self {
        Self {
            move_fn: move_to,
            line_fn: line_to,
            curve_fn: curve_to,
            mark_knot_fn: mark_knot,
            data: Some(Vec::new()),
        }
    }
}
