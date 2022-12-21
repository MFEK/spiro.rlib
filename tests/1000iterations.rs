#[macro_use]
mod data;

use spiro_inner::*;

const TEST_ITERATIONS: usize = 1_000;

#[test]
fn test() {
    let path = test_data!();
    let mut i = 0;
    while i < TEST_ITERATIONS {
        let mut segs = setup_path(&path);
        solve_spiro(&mut segs, path.len().try_into().unwrap());
        i += 1;
    }
}
