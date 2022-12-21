#[macro_use]
mod data;

use spiro::*;

#[test]
fn test() {
    let mut path = test_data!();
    let oplist = run_spiro(&mut path);
    eprintln!("{:?}", oplist);
}
