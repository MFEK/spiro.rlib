//! Provides a [`BezierContext`] that builds up a [`glifparser::outline::PenOperationsPath`] (via
//! [`BezCtxGpPenOpsData`]).
//!
//! Note: If compiled with `default` (`log`) feature, [`log::trace`] is active.

use glifparser::outline::SplitPenOperations;
use glifparser::outline::{PenOperations, PenOperationsContour, PenOperationsPath};
use glifparser::point::{GlifPoint, PointType as GpPointType};

use crate::bezctx_oplist::trace;
use crate::BezierContext;

/// Data built up by this [`BezierContext`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BezCtxGpPenOpsData {
    ops: PenOperationsContour,
    pub ops_path: PenOperationsPath,
    pub knots: Vec<usize>,
    must_close: bool,
}

pub fn move_to(ctx: &mut BezierContext<BezCtxGpPenOpsData, ()>, x: f64, y: f64, is_open: bool) {
    ctx.check_end();
    ctx.data.ops.push(PenOperations::MoveTo(GlifPoint::from_x_y_type(
        (x as f32, y as f32),
        GpPointType::Move,
    )));
    ctx.data.must_close = !is_open;
    trace!("Spiro callback: M {}, {} ", x, y);
}
pub fn line_to(ctx: &mut BezierContext<BezCtxGpPenOpsData, ()>, x: f64, y: f64) {
    ctx.data.ops.push(PenOperations::LineTo(GlifPoint::from_x_y_type(
        (x as f32, y as f32),
        GpPointType::Line,
    )));
    trace!("Spiro callback: L {}, {} ", x, y);
}
pub fn curve_to(ctx: &mut BezierContext<BezCtxGpPenOpsData, ()>, x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) {
    ctx.data.ops.push(PenOperations::CurveTo(
        GlifPoint::from_x_y_type((x3 as f32, y3 as f32), GpPointType::Curve),
        GlifPoint::from_x_y_type((x2 as f32, y2 as f32), GpPointType::OffCurve),
        GlifPoint::from_x_y_type((x1 as f32, y1 as f32), GpPointType::OffCurve),
    ));
    trace!("Spiro callback: C {}, {}, {}, {}, {}, {} ", x1, y1, x2, y2, x3, y3);
}
pub fn mark_knot(ctx: &mut BezierContext<BezCtxGpPenOpsData, ()>, knot_idx: usize) {
    ctx.data.knots.push(knot_idx);
    trace!("Spiro callback: KNOT {}", knot_idx);
}
pub fn end(ctx: &mut BezierContext<BezCtxGpPenOpsData, ()>) {
    ctx.check_end();
    ctx.data.ops_path = ctx.data.ops.split_pen_operations();
    debug_assert!(!ctx.data.ops_path.is_empty());
}

impl BezierContext<BezCtxGpPenOpsData, ()> {
    pub fn new() -> Self {
        Self {
            move_fn: move_to,
            line_fn: line_to,
            curve_fn: curve_to,
            mark_knot_fn: mark_knot,
            start: |_| {},
            end,
            data: Default::default(),
        }
    }

    fn check_end(&mut self) {
        if self.data.must_close {
            self.data.ops.push(PenOperations::Close);
        }
    }
}
