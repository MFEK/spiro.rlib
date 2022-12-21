#![cfg(feature = "glifparser")]

#[macro_use]
mod data;

use glifparser::outline::ToOutline as _;

use spiro::*;

#[test]
fn test() {
    let mut ctx = BezierContext::<BezCtxGpPenOpsData, ()>::new();
    let path = test_data!();
    ctx.run_spiro(&path);
    let outline: glifparser::Outline<()> = ctx.data.ops_path.to_outline();
    let mut glif = glifparser::Glif::<()>::new();
    glif.outline = Some(outline);
    let expected_glif = include_str!("data/spiro.glif");
    let out_glif = glifparser::glif::write(&glif).unwrap();
    eprintln!("{}", &out_glif);
    assert_eq!(out_glif, expected_glif);
}
