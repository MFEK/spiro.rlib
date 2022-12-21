#[macro_use]
mod data;

use spiro_inner::*;
use std::str;

#[test]
fn test() {
    let path = test_data!();

    let mut buffer = ['\0' as u8; 10240];
    {
        let mut writer = std::io::BufWriter::new(buffer.as_mut());
        let mut boxed_writer = Box::new(&mut writer);
        let mut ctx = bezctx_ps::PostScriptBezierContext::new(&mut boxed_writer);
        ctx.run_spiro(&path);
    }

    let bufstr: &str = str::from_utf8(&buffer).unwrap();
    eprintln!("{}", bufstr);
    assert_eq!(bufstr.trim_matches('\0'), include_str!("data/spiro.ps"));
}
