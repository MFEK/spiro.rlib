use spiro_inner::*;
use std::convert::TryInto as _;

#[test]
fn test_curve() {
    let mut path: Vec<SpiroCP> = vec![
        SpiroCP { x: 334., y: 117., ty: 'v' },
        SpiroCP { x: 305., y: 176., ty: 'v' },
        SpiroCP { x: 212., y: 142., ty: 'c' },
        SpiroCP { x: 159., y: 171., ty: 'c' },
        SpiroCP { x: 224., y: 237., ty: 'c' },
        SpiroCP { x: 347., y: 335., ty: 'c' },
        SpiroCP { x: 202., y: 467., ty: 'c' },
        SpiroCP { x: 81., y: 429., ty: 'v' },
        SpiroCP { x: 114., y: 368., ty: 'v' },
        SpiroCP { x: 201., y: 402., ty: 'c' },
        SpiroCP { x: 276., y: 369., ty: 'c' },
        SpiroCP { x: 218., y: 308., ty: 'c' },
        SpiroCP { x: 91., y: 211., ty: 'c' },
        SpiroCP { x: 124., y: 111., ty: 'c' },
        SpiroCP { x: 229., y: 82., ty: 'c' },
    ];
    let mut segs = Vec::new();
    let mut i = 0;

    const TEST_ITERATIONS: usize = 1_000;
    let path_len = path.len().try_into().unwrap();
    while i < TEST_ITERATIONS {
        segs = setup_path(&path);
        solve_spiro(&mut segs, path.len().try_into().unwrap());
        i += 1;
    }

    let oplist = run_spiro(&mut path);

    println!("100 800 translate 1 -1 scale 1 setlinewidth");
    spiro_to_bpath(&mut segs, path_len, &mut bezctx_ps::PostScriptBezierContext::new(bezctx_ps::PostScriptEmitter(std::io::stdout())));
    println!("stroke");
    println!("showpage");

    eprintln!("{:?}", oplist);
}
