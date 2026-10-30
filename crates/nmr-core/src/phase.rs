//! Phase correction: applying a phase, and finding one automatically.

use crate::{Complex64, Metadata};

/// Applies phase in degrees, with `p0` measured at `pivot` (a fraction of the
/// length). nmrglue `proc_base.ps`: `exp(i*(p0 + p1*k/n))`.
pub fn apply(data: &mut [Complex64], p0: f64, p1: f64, pivot: f64) {
    let n = data.len() as f64;
    let (p0, p1) = ((p0 - p1 * pivot).to_radians(), p1.to_radians());
    for (k, z) in data.iter_mut().enumerate() {
        *z *= Complex64::from_polar(1.0, p0 + p1 * k as f64 / n);
    }
}

/// Blocks of baseline at least this wide in Hz stand between two regions, and
/// half of it is added to each side (Zorin et al.: 50 Hz, 0.05 ppm at 500 MHz).
const REGION_GAP_HZ: f64 = 50.0;
/// Rounds of the per-region search, each 20 steps across two steps of the last.
const REGION_ROUNDS: usize = 3;
const REGION_STEPS: usize = 20;
/// A region whose phase is further than this from the fitted line, in radians,
/// is an outlier (Zorin et al.).
const OUTLIER_RAD: f64 = 0.6;
/// Wider than the noise stays alike over, narrower than a baseline rolls.
const ROLL_POINTS: f64 = 32.0;
/// Width at half height of the smoothing each region gets before it is phased,
/// in ppm of the nucleus, for anything but ¹H and ¹⁹F (Zorin et al.: 0.056 ppm,
/// 7 Hz for ¹³C at 125 MHz): lines a point or two wide phase badly on their own.
const SMOOTH_PPM: f64 = 0.056;
/// The first order searched for either way, in degrees.
const MAX_P1: i32 = 720;
/// A maximum counts as a peak of its region when it stands this many times the
/// noise above the region's edges.
const REGION_PEAK_SIGMA: f64 = 5.0;

/// Phase `(p0, p1)` in degrees, pivot at the first point, by the "regions"
/// algorithm of Zorin, Bernstein & Cobas (Magn. Reson. Chem. 55 (2017) 738):
/// the regions that hold peaks are found on the magnitude spectrum, each is
/// given the zero-order phase that leaves least area under the line between its
/// edges, and `p0 + p1 k/n` is fitted through those phases weighted by the
/// peaks in each, dropping outliers. `meta` gives the spacing of the points,
/// the nucleus, and whether the experiment has lines of both signs. The fitted
/// first order is refined within each region before its local phase is measured
/// again, so the phase really refers to the region midpoint.
///
/// Weighted by the square root of each region's peaks rather than the peaks
/// themselves: on the datasets in `testdata` the water and solvent regions
/// otherwise decided the first order between them alone, and the phase across
/// the signals missed the one the operator recorded by up to 24°, where it now
/// misses by 6° at most (2.5° on average), misphased on purpose or not. The
/// peak-argument fit this replaces was 180° out on one of them. In an APT or a
/// DEPT phases are fitted modulo 180° within ±90° of first order, with the
/// tallest line taken to point up. The tallest line, not the heaviest region:
/// two CH2 lines side by side can outweigh the tallest CH.
pub fn auto_phase(data: &[Complex64], meta: &Metadata) -> (f64, f64) {
    let n = data.len();
    let hz_per_point = meta.sw_hz / n as f64;
    if n < 64 || hz_per_point <= 0.0 {
        return (0.0, 0.0);
    }
    let magnitude: Vec<f64> = data.iter().map(|z| z.norm()).collect();
    let mut quiet = still(&derivative(&magnitude));
    // Zero filling and line broadening make neighbouring points of the noise
    // alike, and a spike spread over several points reads as a peak to the
    // derivative; the paper smooths over the distance the noise stays alike.
    let reach = correlation_length(data, &quiet);
    let smoothed = if reach > 1 {
        blur(&magnitude, reach as f64)
    } else {
        magnitude.clone()
    };
    if reach > 1 {
        quiet = still(&derivative(&smoothed));
    }
    // The top of a broad line or the floor between two close ones is flat, and
    // the derivative takes it for baseline: what stands high above the rest of
    // the baseline is not.
    flag_high(&smoothed, &mut quiet);
    let spike = 2 * reach.max(1);
    let gap = (REGION_GAP_HZ / hz_per_point).round() as usize;
    let blocks = merge(&blocks_of(&quiet, reach.max(2), spike), gap);

    let level = noise(&magnitude);
    let smooth = if matches!(meta.nucleus.as_str(), "1H" | "19F") {
        None
    } else {
        // Width at half height to a standard deviation.
        let sigma = SMOOTH_PPM * meta.obs_mhz / hz_per_point / 2.355;
        let re = blur(&data.iter().map(|z| z.re).collect::<Vec<_>>(), sigma);
        let im = blur(&data.iter().map(|z| z.im).collect::<Vec<_>>(), sigma);
        Some(
            re.into_iter()
                .zip(im)
                .map(|(a, b)| Complex64::new(a, b))
                .collect::<Vec<_>>(),
        )
    };
    let phased = smooth.as_deref().unwrap_or(data);
    let edge = (4 * reach).max(4);
    // Remove the estimated slope within each region before measuring its phase
    // at the midpoint; a constant rotation otherwise biases asymmetric regions.
    let mut slope = 0.0;
    let mut result = (0.0, 0.0);
    for _ in 0..5 {
        let mut found: Vec<(f64, f64, f64)> = Vec::new();
        // The phase of the region holding the tallest line, and that line's height.
        let mut tallest = (0.0, 0.0, f64::NEG_INFINITY);
        for &(from, to) in &blocks {
            // Its edges have to be baseline for the line under it to be right.
            if from < gap || to + gap > n {
                continue;
            }
            let (from, to) = (from - gap / 2, (to + gap / 2).min(n));
            let floor = (magnitude[from].max(magnitude[to - 1])).max(0.0);
            let weight: f64 = (from + 3..to.saturating_sub(3))
                .filter(|&i| (i - 3..=i + 3).all(|j| magnitude[j] <= magnitude[i]))
                .map(|i| magnitude[i] - floor)
                .filter(|h| *h > REGION_PEAK_SIGMA * level)
                .sum();
            if weight <= 0.0 || to - from <= 8 {
                continue;
            }
            // A long noise correlation estimate must not consume the whole region.
            let mut region = phased[from..to].to_vec();
            apply(&mut region, 0.0, slope * (to - from) as f64 / n as f64, 0.5);
            let phase = region_phase(&region, edge.min(((to - from) / 4).max(4)));
            let top = magnitude[from..to].iter().copied().fold(0.0, f64::max);
            if top > tallest.2 {
                tallest = ((from + to) as f64 / 2.0 / n as f64, phase, top);
            }
            found.push(((from + to) as f64 / 2.0 / n as f64, phase, weight.sqrt()));
        }
        if found.is_empty() {
            return (0.0, 0.0);
        }
        // Unwrapped against a line fitted round the circle, each phase within half a
        // turn of it. Unwrapping each within half a turn of the last, as the paper
        // does, gained a turn at every jump of more than 180° between neighbours: on
        // a ¹³C whose regions scatter by ±80° that summed to a first order of -627°
        // where -85° was recorded, and the spectrum came out upside down.
        let (mut c0, c1) = circular_fit(&found, meta.signed);
        let period = if meta.signed { 180.0 } else { 360.0 };
        if meta.signed {
            // Negative lines differ by half a turn; the tallest fixes the overall sign.
            c0 += 180.0 * ((tallest.1 - c0 - c1 * tallest.0) / 180.0).round();
        }
        for f in &mut found {
            f.1 += period * ((c0 + c1 * f.0 - f.1) / period).round();
        }
        let (p0, p1) = regress(&mut found, if meta.signed { 90.0 } else { f64::INFINITY });
        // A slope beyond the range the circle was searched over is the regions
        // agreeing by accident, not a phase: the ⁸⁷Rb in `testdata`, one broad line
        // and three specks, fitted PH1 = 1292°, right at the line and rolling
        // everywhere else. The heaviest region alone then gives the zero order.
        if p1.abs() > f64::from(MAX_P1) {
            let heaviest = found.iter().max_by(|a, b| a.2.total_cmp(&b.2));
            let p0 = heaviest.map_or(0.0, |f| f.1);
            return ((p0 + 180.0).rem_euclid(360.0) - 180.0, 0.0);
        }
        result = ((p0 + 180.0).rem_euclid(360.0) - 180.0, p1);
        if (p1 - slope).abs() < 0.05 {
            break;
        }
        slope = p1;
    }
    result
}

/// Holoborodko's smooth noise-robust differentiator for eleven points, the one
/// Zorin et al. use (their Eq. 5).
fn derivative(y: &[f64]) -> Vec<f64> {
    const WEIGHTS: [f64; 5] = [42.0, 48.0, 27.0, 8.0, 1.0];
    let last = y.len() as isize - 1;
    let at = |i: isize| y[i.clamp(0, last) as usize];
    (0..=last)
        .map(|k| {
            WEIGHTS
                .iter()
                .zip(1..)
                .map(|(w, j)| w * (at(k + j) - at(k - j)))
                .sum::<f64>()
                / 512.0
        })
        .collect()
}

/// Which points are baseline: those within three standard deviations of the
/// rest, removing points until none is left to remove (Dietrich et al.).
fn still(values: &[f64]) -> Vec<bool> {
    let mut quiet = vec![true; values.len()];
    loop {
        let (mean, sd) = spread(values, &quiet);
        let mut changed = false;
        for (q, v) in quiet.iter_mut().zip(values) {
            if *q && (v - mean).abs() > 3.0 * sd {
                *q = false;
                changed = true;
            }
        }
        if !changed {
            return quiet;
        }
    }
}

/// `still` once more, on the spectrum itself and upwards only.
fn flag_high(values: &[f64], quiet: &mut [bool]) {
    loop {
        let (mean, sd) = spread(values, quiet);
        let mut changed = false;
        for (q, v) in quiet.iter_mut().zip(values) {
            if *q && v - mean > 3.0 * sd {
                *q = false;
                changed = true;
            }
        }
        if !changed {
            return;
        }
    }
}

fn spread(values: &[f64], quiet: &[bool]) -> (f64, f64) {
    let kept = values
        .iter()
        .zip(quiet)
        .filter(|(_, q)| **q)
        .map(|(v, _)| *v);
    let (count, sum): (f64, f64) = kept.clone().fold((0.0, 0.0), |(c, s), v| (c + 1.0, s + v));
    let mean = sum / count.max(1.0);
    let var = kept.map(|v| (v - mean).powi(2)).sum::<f64>() / count.max(1.0);
    (mean, var.sqrt())
}

/// How many neighbouring points the noise stays alike over: twice the shift at
/// which its autocorrelation falls below 0.6 of its value at none (Zorin et al.
/// Eqs. 6-7), measured on both channels where both points are baseline.
/// Using both keeps the estimate unchanged by a zero-order phase rotation.
fn correlation_length(data: &[Complex64], quiet: &[bool]) -> usize {
    // Less the baseline's slow roll, which is alike over hundreds of points and
    // on its own kept the correlation above 0.6 as far as it was measured.
    let real: Vec<f64> = data.iter().map(|z| z.re).collect();
    let imag: Vec<f64> = data.iter().map(|z| z.im).collect();
    let slow_re = blur(&real, ROLL_POINTS);
    let slow_im = blur(&imag, ROLL_POINTS);
    let values: Vec<Complex64> = data
        .iter()
        .zip(slow_re.iter().zip(&slow_im))
        .map(|(z, (re, im))| z - Complex64::new(*re, *im))
        .collect();
    let (count, sum) = values
        .iter()
        .zip(quiet)
        .filter(|(_, q)| **q)
        .fold((0usize, Complex64::ZERO), |(count, sum), (z, _)| {
            (count + 1, sum + z)
        });
    let mean = sum / count.max(1) as f64;
    let g = |shift: usize| {
        (0..values.len() - shift)
            .filter(|&x| quiet[x] && quiet[x + shift])
            .map(|x| ((values[x] - mean) * (values[x + shift] - mean).conj()).re)
            .sum::<f64>()
    };
    let at_zero = g(0);
    if at_zero <= 0.0 {
        return 1;
    }
    (1..64)
        .find(|&s| g(s) < 0.6 * at_zero)
        .map_or(64, |s| 2 * s)
}

/// Gaussian smoothing of width `sigma` points, as three running means.
fn blur(y: &[f64], sigma: f64) -> Vec<f64> {
    let radius = (((4.0 * sigma * sigma + 1.0).sqrt() - 1.0) / 2.0)
        .round()
        .max(1.0) as usize;
    let mut out = y.to_vec();
    for _ in 0..3 {
        let mut sums = Vec::with_capacity(out.len() + 1);
        sums.push(0.0);
        for v in &out {
            sums.push(sums[sums.len() - 1] + v);
        }
        out = (0..out.len())
            .map(|i| {
                let (a, b) = (i.saturating_sub(radius), (i + radius + 1).min(out.len()));
                (sums[b] - sums[a]) / (b - a) as f64
            })
            .collect();
    }
    out
}

/// Runs of non-baseline points, those `join` or fewer apart made one, and
/// those `spike` or fewer wide dropped as noise.
fn blocks_of(quiet: &[bool], join: usize, spike: usize) -> Vec<(usize, usize)> {
    let mut blocks: Vec<(usize, usize)> = Vec::new();
    let mut start = None;
    for (i, q) in quiet.iter().enumerate() {
        match (q, start) {
            (false, None) => start = Some(i),
            (true, Some(s)) => {
                blocks.push((s, i));
                start = None;
            }
            _ => {}
        }
    }
    if let Some(s) = start {
        blocks.push((s, quiet.len()));
    }
    let blocks = merge(&blocks, join);
    blocks.into_iter().filter(|(a, b)| b - a > spike).collect()
}

/// Blocks closer than `gap` made one.
fn merge(blocks: &[(usize, usize)], gap: usize) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    for &(a, b) in blocks {
        match out.last_mut() {
            Some(last) if a < last.1 + gap => last.1 = b,
            _ => out.push((a, b)),
        }
    }
    out
}

/// The zero-order phase in degrees that leaves least area of the region's real
/// part under the straight line between its two edges, each the mean of `edge`
/// points. The area is not smooth in the phase, so it is mapped rather than
/// descended: 20 steps round the circle, then 20 across the best two, three times.
fn region_phase(region: &[Complex64], edge: usize) -> f64 {
    let n = region.len();
    let below = |phase: f64| {
        let turn = Complex64::from_polar(1.0, phase.to_radians());
        let real: Vec<f64> = region.iter().map(|z| (z * turn).re).collect();
        let left = real[..edge].iter().sum::<f64>() / edge as f64;
        let right = real[n - edge..].iter().sum::<f64>() / edge as f64;
        real.iter()
            .enumerate()
            .map(|(i, v)| {
                let line = left + (right - left) * i as f64 / (n - 1) as f64;
                (line - v).max(0.0)
            })
            .sum::<f64>()
    };
    let (mut centre, mut span) = (0.0, 360.0);
    for _ in 0..REGION_ROUNDS {
        let step = span / REGION_STEPS as f64;
        centre = (0..=REGION_STEPS)
            .map(|k| centre - span / 2.0 + k as f64 * step)
            .min_by(|a, b| below(*a).total_cmp(&below(*b)))
            .unwrap_or(centre);
        span = 2.0 * step;
    }
    centre
}

/// The line `p0 + p1 x` (degrees, `x` a fraction of the spectrum) the phases
/// of `points` lie closest to on the circle: every whole degree of `p1` within
/// two turns, each with `p0` the weighted circular mean of what is left, keeping
/// the one whose phases agree best. Signed spectra use half turns and ±90°.
fn circular_fit(points: &[(f64, f64, f64)], signed: bool) -> (f64, f64) {
    let limit = if signed { 90 } else { MAX_P1 };
    let factor = if signed { 2.0 } else { 1.0 };
    let mut best = (f64::NEG_INFINITY, 0.0, 0.0);
    for p1 in -limit..=limit {
        let p1 = f64::from(p1);
        let sum: Complex64 = points
            .iter()
            .map(|&(x, phase, w)| {
                Complex64::from_polar(w, (factor * (phase - p1 * x)).to_radians())
            })
            .sum();
        if sum.norm() > best.0 {
            best = (sum.norm(), sum.arg().to_degrees() / factor, p1);
        }
    }
    (best.1, best.2)
}

/// Weighted line through region phases, dropping outliers. A lone position
/// determines only zero order; otherwise retain the measured slope, even near
/// zero, so an added first-order phase is not mistaken for a change in confidence.
/// Bound signed fits here, so the intercept is fitted to the bounded slope too.
fn regress(points: &mut Vec<(f64, f64, f64)>, max_slope: f64) -> (f64, f64) {
    loop {
        let total: f64 = points.iter().map(|p| p.2).sum();
        let mx = points.iter().map(|p| p.0 * p.2).sum::<f64>() / total;
        let my = points.iter().map(|p| p.1 * p.2).sum::<f64>() / total;
        let sxx: f64 = points.iter().map(|p| p.2 * (p.0 - mx).powi(2)).sum();
        let sxy: f64 = points.iter().map(|p| p.2 * (p.0 - mx) * (p.1 - my)).sum();
        let slope = if sxx > 0.0 { sxy / sxx } else { 0.0 }.clamp(-max_slope, max_slope);
        let fit = |x: f64| my + slope * (x - mx);
        let worst = points
            .iter()
            .enumerate()
            .map(|(k, p)| (k, (p.1 - fit(p.0)).abs()))
            .max_by(|a, b| a.1.total_cmp(&b.1));
        if let Some((k, miss)) = worst
            && miss.to_radians() > OUTLIER_RAD
            && points.len() > 2
        {
            points.remove(k);
            continue;
        }
        return if sxx == 0.0 {
            (my, 0.0)
        } else {
            (my - slope * mx, slope)
        };
    }
}

/// Position of the largest magnitude as a fraction of the spectrum length: the
/// default pivot for manual phasing.
pub fn largest_peak(data: &[Complex64]) -> f64 {
    let Some((index, _)) = data
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.norm_sqr().total_cmp(&b.1.norm_sqr()))
    else {
        return 0.0;
    };
    index as f64 / data.len() as f64
}

/// Noise from the median absolute difference between neighbours, which a broad
/// baseline does not inflate the way the deviation from the median does.
fn noise(magnitude: &[f64]) -> f64 {
    let mut differences: Vec<f64> = magnitude.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    if differences.is_empty() {
        return 0.0;
    }
    let middle = differences.len() / 2;
    let (_, median, _) = differences.select_nth_unstable_by(middle, f64::total_cmp);
    1.4826 * *median / std::f64::consts::SQRT_2
}

#[cfg(test)]
mod tests {
    use super::*;

    const LINES: [(f64, f64); 5] = [(0.15, 1.0), (0.4, 0.6), (0.42, 0.6), (0.7, 0.3), (0.9, 0.8)];

    /// Absorption-mode Lorentzians with a known phase error and small noise.
    fn dephased(n: usize, lines: &[(f64, f64)], p0: f64, p1: f64) -> Vec<Complex64> {
        (0..n)
            .map(|k| {
                let x = k as f64;
                let line: Complex64 = lines
                    .iter()
                    .map(|&(c, h)| h * 4.0 / Complex64::new(4.0, x - c * n as f64))
                    .sum();
                let noise = ((k * 7919 % 997) as f64 / 997.0 - 0.5) * 1e-3;
                let angle = (p0 + p1 * x / n as f64).to_radians();
                (line + noise) * Complex64::from_polar(1.0, angle)
            })
            .collect()
    }

    /// One Hz a point over `n` points.
    fn meta(nucleus: &str, n: usize, signed: bool) -> Metadata {
        Metadata {
            nucleus: nucleus.into(),
            sw_hz: n as f64,
            obs_mhz: 100.0,
            carrier_hz: 0.0,
            group_delay: 0.0,
            solvent: None,
            phase: None,
            conjugate: false,
            signed,
            size: None,
            line_broadening: None,
            delay_at_ft: false,
        }
    }

    fn angle_difference(a: f64, b: f64) -> f64 {
        ((a - b + 180.0).rem_euclid(360.0) - 180.0).abs()
    }

    #[test]
    fn noise_correlation_does_not_depend_on_zero_order_phase() {
        let data: Vec<_> = (0..1024)
            .map(|i| {
                Complex64::new(
                    (i * 7919 % 997) as f64 / 997.0 - 0.5,
                    (i as f64 * 0.37).sin(),
                )
            })
            .collect();
        let quiet = vec![true; data.len()];
        let expected = correlation_length(&data, &quiet);
        for angle in [45.0, 90.0, 135.0] {
            let mut rotated = data.clone();
            apply(&mut rotated, angle, 0.0, 0.0);
            assert_eq!(correlation_length(&rotated, &quiet), expected);
        }
    }

    #[test]
    fn small_first_order_is_not_dropped_when_region_phases_scatter() {
        let mut points = vec![
            (0.1, 1.0, 1.0),
            (0.3, -1.0, 1.0),
            (0.5, 0.0, 1.0),
            (0.9, 1.0, 1.0),
        ];
        let mut shifted = points
            .iter()
            .map(|&(x, phase, weight)| (x, phase + 50.0 - 80.0 * x, weight))
            .collect();
        let (p0, p1) = regress(&mut points, f64::INFINITY);
        let (q0, q1) = regress(&mut shifted, f64::INFINITY);
        assert!((q0 - p0 - 50.0).abs() < 1e-10);
        assert!((q1 - p1 + 80.0).abs() < 1e-10);
    }

    #[test]
    fn asymmetric_regions_recover_sparse_spectra_without_reference_phases() {
        // Close lines make wide regions whose midpoint is not a peak position.
        let lines = [
            (0.15, 1.0),
            (0.38, 0.6),
            (0.395, 0.5),
            (0.425, 0.8),
            (0.433, 0.8),
            (0.6264, 1.0),
            (0.6272, 1.0),
            (0.628, 1.0),
            (0.751, 0.8),
            (0.808, 0.6),
            (0.838, 0.5),
            (0.855, 0.5),
        ];
        let n = 16384;
        for nucleus in ["1H", "13C"] {
            for width in [2.0, 4.0, 8.0] {
                let expected: Vec<Complex64> = (0..n)
                    .map(|i| {
                        let line: Complex64 = lines
                            .iter()
                            .map(|(c, h)| {
                                h * width / Complex64::new(width, i as f64 - c * n as f64)
                            })
                            .sum();
                        line + ((i * 7919 % 997) as f64 / 997.0 - 0.5) * 1e-3
                    })
                    .collect();
                let mut metadata = meta(nucleus, n, false);
                metadata.obs_mhz = 150.0;
                for (e0, e1) in [(40.0, 60.0), (40.0, 250.0), (-100.0, -260.0)] {
                    let mut data = expected.clone();
                    apply(&mut data, e0, e1, 0.0);
                    let (p0, p1) = auto_phase(&data, &metadata);
                    apply(&mut data, p0, p1, 0.0);
                    let error = data
                        .iter()
                        .zip(&expected)
                        .map(|(a, b)| (a.re - b.re).powi(2))
                        .sum::<f64>();
                    let signal = expected.iter().map(|z| z.re.powi(2)).sum::<f64>();
                    let rms = (error / signal).sqrt();
                    let tolerance = if nucleus == "13C" { 0.005 } else { 0.01 };
                    assert!(
                        rms < tolerance,
                        "{nucleus} width {width}, applied {e0}/{e1}, auto {p0}/{p1}, relative RMS {rms}"
                    );
                }
            }
        }
    }

    #[test]
    fn recovers_applied_phase() {
        // 13C is smoothed before its regions are phased; the answer is the same.
        for nucleus in ["1H", "13C"] {
            for (e0, e1) in [(70.0, -120.0), (-150.0, 45.0), (10.0, 250.0)] {
                let data = dephased(8192, &LINES, e0, e1);
                let (p0, p1) = auto_phase(&data, &meta(nucleus, 8192, false));
                assert!(
                    angle_difference(p0, -e0) < 2.0,
                    "{nucleus} p0 {p0} for {e0}"
                );
                assert!((p1 + e1).abs() < 3.0, "{nucleus} p1 {p1} for {e1}");
            }
        }
    }

    #[test]
    fn asymmetric_peak_regions_recover_phase_at_signal_positions() {
        let lines = [
            (0.15, 0.8),
            (0.4, 1.0),
            (0.405, 0.12),
            (0.7, 0.6),
            (0.9, 0.7),
        ];
        let data = dephased(8192, &lines, 40.0, 600.0);
        let (p0, p1) = auto_phase(&data, &meta("1H", 8192, false));
        for (position, _) in lines {
            let error = angle_difference(p0 + p1 * position, -40.0 - 600.0 * position);
            assert!(
                error < 0.1,
                "position {position}, phase error {error}, p0 {p0}, p1 {p1}"
            );
        }
    }

    #[test]
    fn noise_free_spectrum_has_a_determinable_phase() {
        let n = 8192;
        let lines = [
            (0.16, 1.0),
            (0.4, 1.0),
            (0.405, 0.05),
            (0.70, 0.5),
            (0.90, 0.8),
        ];
        let data: Vec<_> = (0..n)
            .map(|k| {
                let x = k as f64;
                let line: Complex64 = lines
                    .iter()
                    .map(|(c, h)| h * 4.0 / Complex64::new(4.0, x - c * n as f64))
                    .sum();
                line * Complex64::from_polar(1.0, (40.0 + 250.0 * x / n as f64).to_radians())
            })
            .collect();
        let (p0, p1) = auto_phase(&data, &meta("1H", n, false));
        assert!(angle_difference(p0, -40.0) < 3.0, "p0 {p0}, p1 {p1}");
        assert!((p1 + 250.0).abs() < 4.0, "p0 {p0}, p1 {p1}");
    }

    #[test]
    fn lines_pointing_down_are_kept_down_where_they_may() {
        // A DEPT-135: CH2 lines point down, and are not a phase error.
        let lines = [
            (0.15, 1.0),
            (0.4, -0.6),
            (0.42, -0.6),
            (0.7, 0.3),
            (0.9, -0.8),
        ];
        for (e0, e1) in [
            (40.0, 60.0),
            (-150.0, 60.0),
            (40.0, -60.0),
            (40.0, 90.0),
            (0.0, 89.9),
            (0.0, -89.9),
            (90.0, 89.9),
            (90.0, -89.9),
        ] {
            let mut data = dephased(8192, &lines, e0, e1);
            let (p0, p1) = auto_phase(&data, &meta("13C", 8192, true));
            assert!(
                angle_difference(p0, -e0) < 3.0,
                "applied {e0}/{e1}, p0 {p0}"
            );
            assert!((p1 + e1).abs() < 4.0, "applied {e0}/{e1}, p1 {p1}");
            apply(&mut data, p0, p1, 0.0);
            let expected = dephased(8192, &lines, 0.0, 0.0);
            let error = data
                .iter()
                .zip(&expected)
                .map(|(a, b)| (a.re - b.re).powi(2))
                .sum::<f64>();
            let signal = expected.iter().map(|z| z.re.powi(2)).sum::<f64>();
            assert!(
                (error / signal).sqrt() < 0.01,
                "applied {e0}/{e1}: signed RMS"
            );
        }
    }

    #[test]
    fn largest_peak_is_a_fraction() {
        let data = dephased(1000, &LINES, 90.0, 0.0);
        assert!((largest_peak(&data) - 0.15).abs() < 1e-3);
        assert_eq!(largest_peak(&[]), 0.0);
    }

    #[test]
    fn empty_or_flat_input() {
        assert_eq!(auto_phase(&[], &meta("1H", 0, false)), (0.0, 0.0));
        let flat = [Complex64::ZERO; 256];
        assert_eq!(auto_phase(&flat, &meta("1H", 256, false)), (0.0, 0.0));
    }
}
