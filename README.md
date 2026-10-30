# openNMR

Free, open-source NMR processing for the desktop.

openNMR takes a Bruker FID or processed spectrum through phasing, baseline
correction, referencing and annotation, then produces a publication figure and a
multiplet report. It is built with Tauri 2, Svelte 5, TypeScript and a Rust
numerical core.

> **Status: early pre-release (0.1.1).** openNMR is under active development.
> Its processing and fitted areas have not been validated for general research
> use. Check the phase, baseline, referencing and integrations before relying on
> any result. Release builds are currently unsigned; macOS builds are not
> notarized by Apple.

## Features

- **Bruker import** — opens experiment folders containing `fid`, `acqus` and
  `pdata/N`, including processed `1r` data and 2D experiments.
- **Reproducible processing** — processing is a list of removable operations
  replayed from the raw data, not a destructive edit.
- **Phasing and baseline** — automatic and manual phase correction with a
  visible pivot, automatic and point-based baseline correction, smoothing and
  window functions.
- **Peaks and integrals** — automatic or manual peak picking, Sum integration,
  and fitted areas from bounded generalized-Lorentzian deconvolution.
- **Multiplet analysis** — automatic multiplet grouping, coupling constants and
  a copyable report.
- **Figures and export** — page layout in millimetres, SVG and PDF export with
  the project attached, plus CSV peak and integral tables.
- **Projects** — `.opennmr` files carry the data and settings needed to
  reproduce a result on another computer.
- **Keyboard-first** — every action is reachable from the keyboard, with an
  in-app list of keys and mouse gestures.

The browser preview under `apps/desktop/preview.html` is a development aid on
synthetic data, not a supported product surface.

## Download

Prebuilt installers for Windows, macOS and Linux are attached to
[GitHub Releases](../../releases).

The builds are not code-signed. On macOS, Gatekeeper may warn that the developer
cannot be verified; on Windows, SmartScreen may show a similar warning. Only
download builds from this repository's Releases page.

## Try the example

`examples/synthetic-ethanol/` contains a simulated proton spectrum:

1. Choose **Open folder**.
2. Select `examples/synthetic-ethanol/1` — the folder containing `fid`, `acqus`
   and `pdata`.
3. Do not select the `fid` file on its own.

Expected signals are at approximately 1.22, 2.61, 3.69, 7.26 and 0 ppm. This is
synthetic data, not an experimental measurement.

## Build from source

### Prerequisites

- Rust stable (edition 2024, so Rust 1.85 or newer)
- Node.js 24 and npm
- The [Tauri 2 system dependencies](https://v2.tauri.app/start/prerequisites/)
  for your platform

On Debian/Ubuntu:

```sh
sudo apt-get install -y libwebkit2gtk-4.1-dev libxdo-dev libssl-dev \
  libayatana-appindicator3-dev librsvg2-dev build-essential
```

### Build and run

```sh
cd apps/desktop
npm ci
npm run tauri dev      # development
npm run tauri build    # release bundle in src-tauri/target/release/bundle/
```

The Tauri configuration points `frontendDist` at `apps/desktop/dist`, so build
the frontend before running Rust checks from a clean checkout:

```sh
cd apps/desktop
npm ci
npm run build
cd ../..
cargo test --workspace
```

### Tests and checks

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace

cd apps/desktop
npm run check
npm run lint
npm test
```

## Repository layout

```text
apps/desktop/        Tauri desktop app (Svelte frontend, Rust backend)
crates/nmr-core/     Reading, processing and analysis of NMR spectra
examples/            Synthetic example data
testdata/            Placeholder for real Bruker data (not committed)
tools/reference/     Comparison and synthetic-data scripts
```

## Scientific notes

Fitted peak areas use a generalized Lorentzian line shape with a kurtosis
parameter. The model is used for high-precision qNMR by:

> T. Schoenberger, S. Menges, M. A. Bernstein, M. Perez, F. Seoane, S. Sykora
> and C. Cobas, “Improving the Performance of High-Precision qNMR Measurements
> by a Double Integration Procedure in Practical Cases”, *Analytical Chemistry*
> **88**, 3836–3843 (2016).
> DOI: [10.1021/acs.analchem.5b04911](https://doi.org/10.1021/acs.analchem.5b04911)

The same line-shape family is used for multiplet deconvolution by:

> D. Jeannerat and C. Cobas, “Application of multiplet structure deconvolution
> to extract scalar coupling constants from 1D nuclear magnetic resonance
> spectra”, *Magnetic Resonance* **2**, 545–555 (2021).
> DOI: [10.5194/mr-2-545-2021](https://doi.org/10.5194/mr-2-545-2021)

The area is the analytic integral of the fitted line. Detection, optimization
and model selection are local choices.

## Licence

openNMR is dual-licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in this project by you, as defined in the Apache-2.0 licence,
shall be dual-licensed as above, without any additional terms or conditions.

## Disclaimer

This software is provided “as is”, without warranty of any kind. It is a
research and productivity tool, not a validated analytical instrument. Verify
all results independently before using them in publications, regulatory work or
any other setting where an incorrect result could cause harm.
