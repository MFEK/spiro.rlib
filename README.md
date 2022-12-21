# `spiro.rlib` v1.0.0 (⏫︎2022-12-20)

This is Raph Levien's C implementation of Spiro in pure Rust, transpiled by C2Rust and then heavily edited by Fredrick R. Brennan (@ctrlcctrlv).

## Features

* Full support for all Spiro point types.
* No more unsafe code as of v1.0.0 (thanks Seth).
* Optional [`glifparser`](https://github.com/MFEK/glifparser.rlib) support, through which (`ToKurbo` on an `Outline` convertible from a `PenOperationsPath`) one can also get Raph Levien's [`kurbo`](https://docs.rs/kurbo) working.
  * This provides `.glif` output. See § «Usage example».
* Easy PostScript output context.

## Usage example

```rust
use spiro::BezierContext;
use glifparser::{Glif, outline::ToOutline as _};

let path = test_data!();
ctx.run_spiro(&path);
let mut glif = Glif::<()>::new();
glif.outline = Some(ctx.data.ops_path.to_outline());
eprintln!("{}", glifparser::glif::write(&glif).unwrap());
```

## License (this edit of C2Rust output)
```
Copyright (C) 2020–2022 Fredrick R. Brennan, Seth Erfurt and MFEK Authors
Copyright (C) 2007 Raph Levien

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this software or any of the provided source code files except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied.  See the License for the
specific language governing permissions and limitations under the License.
```

## License (original C code)
```
Copyright (C) 2007 Raph Levien

Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
<LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
option. This file may not be copied, modified, or distributed
except according to those terms.
```
