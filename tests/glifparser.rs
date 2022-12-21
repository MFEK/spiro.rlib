#![cfg(feature = "glifparser")]

#[macro_use]
mod data;

use glifparser::outline::ToOutline as _;

use spiro::*;

#[test]
fn test() {
    let mut ctx = BezierContext::<BezCtxGpPenOpsData, ()>::new();
    (ctx.start)(&mut ctx);
    let path = test_data!();
    let mut segs = setup_path(&path);
    solve_spiro(&mut segs);
    spiro_to_bpath(&mut segs, path.len(), &mut ctx);
    (ctx.end)(&mut ctx);
    let outline: glifparser::Outline<()> = ctx.data.ops_path.to_outline();
    let mut glif = glifparser::Glif::<()>::new();
    glif.outline = Some(outline);
    let expected_glif = include_str!("data/spiro.glif");
    let out_glif = glifparser::glif::write(&glif).unwrap();
    eprintln!("{}", &out_glif);
    assert_eq!(out_glif, expected_glif);
}
