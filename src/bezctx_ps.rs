//! Provides a [`BezierContext`] that writes PostScript code to an [`io::Write`].

use crate::BezierContext;
use io::{Result, Write};
use std::io;

pub type PostScriptBezierContext<'a, W, R> = BezierContext<PostScriptEmitter<&'a mut Box<W>>, R>;
pub struct PostScriptEmitter<W: Write>(pub W);

impl<W> Write for PostScriptEmitter<W>
where
    W: Write,
{
    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        self.0.write(buf)
    }
    fn flush(&mut self) -> Result<()> {
        self.0.flush()
    }
}

pub fn move_to<W: Write>(w: &mut PostScriptBezierContext<W, Result<()>>, x: f64, y: f64, _is_open: bool) -> Result<()> {
    writeln!(w.data, "{:.5} {:.5} moveto", x, y)
}

pub fn line_to<W: Write>(w: &mut PostScriptBezierContext<W, Result<()>>, x: f64, y: f64) -> Result<()> {
    writeln!(w.data, "{:.5} {:.5} lineto", x, y)
}

#[rustfmt::skip]
pub fn curve_to<W: Write>(w: &mut PostScriptBezierContext<W, Result<()>>, x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) -> Result<()> {
    writeln!(w.data, "{:.5} {:.5} {:.5} {:.5} {:.5} {:.5} curveto", x1, y1, x2, y2, x3, y3)
}

pub fn mark_knot<W: Write>(w: &mut PostScriptBezierContext<W, Result<()>>, knot_idx: usize) -> Result<()> {
    writeln!(w.data, "% {} knot", knot_idx)
}
pub fn start<W: Write>(w: &mut PostScriptBezierContext<W, Result<()>>) -> Result<()> {
    writeln!(w.data, "%!PS\n1 -1 scale\n0 -800 translate\n1 setlinewidth")
}
pub fn end<W: Write>(w: &mut PostScriptBezierContext<W, Result<()>>) -> Result<()> {
    writeln!(w.data, "closepath\nstroke\nshowpage")
}

impl<'a, W: Write> PostScriptBezierContext<'a, W, Result<()>> {
    pub fn new(w: &'a mut Box<W>) -> Self {
        Self {
            move_fn: move_to,
            line_fn: line_to,
            curve_fn: curve_to,
            mark_knot_fn: mark_knot,
            start,
            end,
            data: PostScriptEmitter(w),
        }
    }
}
