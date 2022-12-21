#[macro_use]
mod data;

use spiro_inner::*;

const TEST_ITERATIONS: usize = 1_000;

#[test]
fn test() {
    let path = test_data!();
    for _ in 0..TEST_ITERATIONS {
        let mut segs = setup_path(&path);
        solve_spiro(&mut segs);
    }
}
