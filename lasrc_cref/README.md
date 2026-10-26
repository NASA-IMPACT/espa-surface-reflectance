# lasrc_cref - C vs Rust differential test harness

Links the original C LaSRC leaf functions (compiled from
`../lasrc/c_version/src`) against the `lasrc_core` Rust port and fuzzes both with
identical inputs, asserting they agree. The goal is to localize where the port
diverges from C -- in particular the documented f32/f64 precision gap in the
aerosol retrieval -- instead of reverse-engineering it from whole-scene diffs.

This is a standalone Cargo workspace: building it never touches the production
`lasrc_py` / `lasrc-rs` build or its `Cargo.lock`.

## Staged by dependency weight

The C leaves split into groups by how much of the ESPA/HDF chain they pull in:

| Stage | Functions | C deps | Where it runs |
|-------|-----------|--------|---------------|
| 0 | `quick_select` | none (own header only) | any Rust toolchain |
| 1 | `atmcorlamb2_new`, `poly_coeff` | common.h -> hdf, lut_subr.h -> espa_* | conda/pixi env (Docker on macOS) |
| 2 | `subaeroret_new` | same as stage 1 | conda/pixi env |

Bottom-up on purpose: prove the FFI + proptest plumbing on the trivial pure
function first, then the per-band correction kernel, then the aerosol retrieval
that composes it (the documented precision-gap source).

## Running

Stage 0 (no conda needed):

```bash
cd lasrc_cref
cargo test
```

Stage 1+ (needs the conda headers/libs). On Linux:

```bash
cd lasrc_cref
pixi run test-lut          # sets PREFIX=$CONDA_PREFIX, builds --features lut
```

On macOS (osx-arm64 has no C packages), run stage 1+ in the linux-64 container.
Build context is the repo root:

```bash
docker build -f lasrc_cref/Dockerfile --platform linux/amd64 \
  --secret id=aws_access_key_id,env=AWS_ACCESS_KEY_ID \
  --secret id=aws_secret_access_key,env=AWS_SECRET_ACCESS_KEY \
  --secret id=aws_session_token,env=AWS_SESSION_TOKEN \
  -t lasrc-cref .
docker run --rm lasrc-cref            # runs `pixi run test-lut`
```

The AWS secrets let pixi fetch espa-product-formatter / espa-surface-reflectance
from the private `s3://hls-conda-channels` channel.

## Fidelity notes

- The C is compiled with `-ffp-contract=off` so the reference does not fuse
  multiply-adds. Run inside the activated conda env so `cc` inherits the same
  `$CC`/`$CFLAGS` used for the shipped binary; otherwise differences may be
  compiler codegen artifacts rather than real port bugs.
- Inputs are generated once and fed to both sides as identical values. For
  functions C takes as `float`, generate `f32` and widen, so both start from the
  same number.
- Compare on the value that matters: bit-exact for pure-arithmetic leaves;
  tight ULP tolerance where transcendentals (exp/log/pow/cos) are involved; and
  always the final scaled-int output, which is the real product contract.

## First finding

`quick_select` agrees bit-for-bit with C for all odd-length arrays but selects a
different element for even lengths (C uses index `(n-1)/2`, the Rust port uses
`n/2`). LaSRC only ever uses odd windows, so it is harmless in practice -- but
it is exactly the off-by-one the harness is meant to catch. See
`tests/quick_select.rs`.
