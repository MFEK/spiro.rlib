//! Inner functions of libspiro, translated essentially verbatim by C2Rust.
//!
//! For more information please see: [Levien, Raph Linus (2009). «From Spiral to Spline: Optimal
//! Techniques in Interactive Curve Design». A dissertation submitted in partial satisfaction of
//! the requirements for the degree of Doctor of Philosophy in Engineering–Electrical Engineering
//! and Computer Sciences in the Graduate Division. Berkeley: University of
//! California.](https://levien.com/phd/thesis.pdf), [Chapter 8: Numerical
//! Toolbox](https://levien.com/phd/toolbox_frame.pdf).

use super::{BezierContext, SpiroSegment};

#[derive(Copy, Clone, Debug, Default)]
pub struct BandMath {
    pub(crate) a: [f64; 11],
    pub(crate) al: [f64; 5],
}

/// Integrate polynomial spiral curve over range `-.5`…`.5`.
pub fn integrate_spiro(ks: [f64; 4], xy: &mut [f64; 2]) {
    let th1: f64 = ks[0];
    let th2: f64 = 0.5 * ks[1];
    let th3: f64 = 1.0 / 6. * ks[2];
    let th4: f64 = 1.0 / 24. * ks[3];
    let mut x: f64 = 0.;
    let mut y: f64 = 0.;
    let ds: f64 = 1.0 / 4 as f64;
    let ds2: f64 = ds * ds;
    let ds3: f64 = ds2 * ds;
    let k0: f64 = ks[0] * ds;
    let k1: f64 = ks[1] * ds;
    let k2: f64 = ks[2] * ds;
    let k3: f64 = ks[3] * ds;
    let mut i: isize = 0;
    let mut s: f64 = 0.5 * ds - 0.5;
    while i < 4 {
        let mut u;
        let mut v;
        let km0;
        let km1;
        let km2;
        let km3;

        km0 = ((1.0 / 6. * k3 * s + 0.5 * k2) * s + k1) * s + k0;
        km1 = ((0.5 * k3 * s + k2) * s + k1) * ds;
        km2 = (k3 * s + k2) * ds2;

        km3 = k3 * ds3;
        let t1_1: f64 = km0;
        let t1_2: f64 = 0.5 * km1;
        let t1_3: f64 = 1.0 / 6. * km2;
        let t1_4: f64 = 1.0 / 24. * km3;
        let t2_2: f64 = t1_1 * t1_1;
        let t2_3: f64 = 2. * (t1_1 * t1_2);
        let t2_4: f64 = 2. * (t1_1 * t1_3) + t1_2 * t1_2;
        let t2_5: f64 = 2. * (t1_1 * t1_4 + t1_2 * t1_3);
        let t2_6: f64 = 2. * (t1_2 * t1_4) + t1_3 * t1_3;
        let t2_7: f64 = 2. * (t1_3 * t1_4);
        let t2_8: f64 = t1_4 * t1_4;
        let t3_4: f64 = t2_2 * t1_2 + t2_3 * t1_1;
        let t3_6: f64 = t2_2 * t1_4 + t2_3 * t1_3 + t2_4 * t1_2 + t2_5 * t1_1;
        let t3_8: f64 = t2_4 * t1_4 + t2_5 * t1_3 + t2_6 * t1_2 + t2_7 * t1_1;
        let t3_10: f64 = t2_6 * t1_4 + t2_7 * t1_3 + t2_8 * t1_2;
        let t4_4: f64 = t2_2 * t2_2;
        let t4_5: f64 = 2. * (t2_2 * t2_3);
        let t4_6: f64 = 2. * (t2_2 * t2_4) + t2_3 * t2_3;
        let t4_7: f64 = 2. * (t2_2 * t2_5 + t2_3 * t2_4);
        let t4_8: f64 = 2. * (t2_2 * t2_6 + t2_3 * t2_5) + t2_4 * t2_4;
        let t4_9: f64 = 2. * (t2_2 * t2_7 + t2_3 * t2_6 + t2_4 * t2_5);
        let t4_10: f64 = 2. * (t2_2 * t2_8 + t2_3 * t2_7 + t2_4 * t2_6) + t2_5 * t2_5;
        let t5_6: f64 = t4_4 * t1_2 + t4_5 * t1_1;
        let t5_8: f64 = t4_4 * t1_4 + t4_5 * t1_3 + t4_6 * t1_2 + t4_7 * t1_1;
        let t5_10: f64 = t4_6 * t1_4 + t4_7 * t1_3 + t4_8 * t1_2 + t4_9 * t1_1;
        let t6_6: f64 = t4_4 * t2_2;
        let t6_7: f64 = t4_4 * t2_3 + t4_5 * t2_2;
        let t6_8: f64 = t4_4 * t2_4 + t4_5 * t2_3 + t4_6 * t2_2;
        let t6_9: f64 = t4_4 * t2_5 + t4_5 * t2_4 + t4_6 * t2_3 + t4_7 * t2_2;
        let t6_10: f64 = t4_4 * t2_6 + t4_5 * t2_5 + t4_6 * t2_4 + t4_7 * t2_3 + t4_8 * t2_2;
        let t7_8: f64 = t6_6 * t1_2 + t6_7 * t1_1;
        let t7_10: f64 = t6_6 * t1_4 + t6_7 * t1_3 + t6_8 * t1_2 + t6_9 * t1_1;
        let t8_8: f64 = t6_6 * t2_2;
        let t8_9: f64 = t6_6 * t2_3 + t6_7 * t2_2;
        let t8_10: f64 = t6_6 * t2_4 + t6_7 * t2_3 + t6_8 * t2_2;
        let t9_10: f64 = t8_8 * t1_2 + t8_9 * t1_1;
        let t10_10: f64 = t8_8 * t2_2;
        u = 1.;
        v = 0.;
        v += 1.0 / 12. * t1_2 + 1.0 / 80. * t1_4;
        u -= 1.0 / 24. * t2_2 + 1.0 / 160. * t2_4 + 1.0 / 896. * t2_6 + 1.0 / 4608. * t2_8;
        v -= 1.0 / 480. * t3_4 + 1.0 / 2688. * t3_6 + 1.0 / 13824. * t3_8 + 1.0 / 67584. * t3_10;
        u += 1.0 / 1920. * t4_4 + 1.0 / 10752. * t4_6 + 1.0 / 55296. * t4_8 + 1.0 / 270336. * t4_10;
        v += 1.0 / 53760. * t5_6 + 1.0 / 276480. * t5_8 + 1.0 / 1.35168e+06 * t5_10;
        u -= 1.0 / 322560. * t6_6 + 1.0 / 1.65888e+06 * t6_8 + 1.0 / 8.11008e+06 * t6_10;
        v -= 1.0 / 1.16122e+07 * t7_8 + 1.0 / 5.67706e+07 * t7_10;
        u += 1.0 / 9.28973e+07 * t8_8 + 1.0 / 4.54164e+08 * t8_10;
        v += 1.0 / 4.08748e+09 * t9_10;
        u -= 1.0 / 4.08748e+10 * t10_10;

        let th: f64 = (((th4 * s + th3) * s + th2) * s + th1) * s;
        let cth: f64 = (th).cos();
        let sth: f64 = (th).sin();
        x += cth * u - sth * v;
        y += cth * v + sth * u;
        s += ds;

        i += 1
    }
    xy[0] = x * ds;
    xy[1] = y * ds;
}

pub(crate) fn compute_ends(ks: [f64; 4], ends: &mut [[f64; 4]; 2], seg_ch: f64) -> f64 {
    let mut xy: [f64; 2] = [0.; 2];
    let ch;
    let th;
    let l;
    let l2;
    let l3;
    let th_even;
    let th_odd;
    let k0_even;
    let k0_odd;
    let k1_even;
    let k1_odd;
    let k2_even;
    let k2_odd;
    integrate_spiro(ks, &mut xy);
    ch = xy[0].hypot(xy[1]);
    th = xy[1].atan2(xy[0]);
    l = ch / seg_ch;
    th_even = 0.5 * ks[0] + 1.0 / 48. * ks[2];
    th_odd = 0.125 * ks[1] + 1.0 / 384. * ks[3] - th;
    (ends[0])[0] = th_even - th_odd;
    (ends[1])[0] = th_even + th_odd;
    k0_even = l * (ks[0] + 0.125 * ks[2]);
    k0_odd = l * (0.5 * ks[1] + 1.0 / 48. * ks[3]);
    (ends[0])[1] = k0_even - k0_odd;
    (ends[1])[1] = k0_even + k0_odd;
    l2 = l * l;
    k1_even = l2 * (ks[1] + 0.125 * ks[3]);
    k1_odd = l2 * 0.5 * ks[2];
    (ends[0])[2] = k1_even - k1_odd;
    (ends[1])[2] = k1_even + k1_odd;
    l3 = l2 * l;
    k2_even = l3 * ks[2];
    k2_odd = l3 * 0.5 * ks[3];
    (ends[0])[3] = k2_even - k2_odd;
    (ends[1])[3] = k2_even + k2_odd;
    return l;
}

/// Computes Spiro’s Jacobian matrix (matrix of all the first-order partial derivatives of a vector-valued function of several variables)
pub(crate) fn compute_pderivs(s: &SpiroSegment, ends: &mut [[f64; 4]; 2], derivs: &mut [[[f64; 4]; 2]; 4], jinc: usize) {
    let recip_d: f64 = 2e6;
    let delta: f64 = 1.0 / recip_d;
    let mut try_ks: [f64; 4] = [0.; 4];
    let mut try_ends: [[f64; 4]; 2] = [[0.; 4]; 2];
    let mut i: usize = 0;
    let mut j;
    let mut k;
    compute_ends(s.ks, ends, s.seg_ch);
    while i < jinc + 1 {
        j = 0;
        while j < 4 {
            try_ks[j] = s.ks[j];
            j += 1
        }
        try_ks[i] += delta;
        compute_ends(try_ks, &mut try_ends, s.seg_ch);
        k = 0;
        while k < 2 {
            j = 0;
            while j < 4 {
                (derivs[j])[k][i] = recip_d * (try_ends[k][j] - (ends[k])[j]);
                j += 1
            }
            k += 1
        }
        i += 1
    }
}

use std::f64::consts::PI;
pub(crate) fn mod2pi(th: f64) -> f64 {
    let u = th / (2. * PI);
    return 2. * PI * (u - (u + 0.5).floor());
}

pub(crate) fn bandec11(m: &mut [BandMath], perm: &mut [usize], n: usize) {
    let mut i: usize = 0;
    let mut j;
    let mut k: usize = 0;
    let mut l: usize = 5;
    /* pack top triangle to the left. */
    while i < 5 {
        j = 0;
        while j < i + 6 {
            m[i as usize].a[j as usize] = m[i as usize].a[(j + 5 - i) as usize];
            j += 1
        }
        while j < 11 {
            m[i as usize].a[j as usize] = 0.0;
            j += 1
        }
        i += 1
    }
    while k < n {
        let mut pivot: usize = k;
        let mut pivot_val: f64 = m[k as usize].a[0];
        if l < n {
            l += 1
        }
        j = k + 1;
        while j < l {
            if m[j as usize].a[0].abs() > pivot_val.abs() {
                pivot_val = m[j as usize].a[0];
                pivot = j;
            }
            j += 1
        }
        perm[k as usize] = pivot;
        if pivot != k {
            j = 0;
            while j < 11 {
                let tmp: f64 = m[k as usize].a[j as usize];
                m[k as usize].a[j as usize] = m[pivot as usize].a[j as usize];
                m[pivot as usize].a[j as usize] = tmp;
                j += 1
            }
        }
        if pivot_val.abs() < 1e-12 {
            pivot_val = 1e-12
        }
        let pivot_scale = 1.0 / pivot_val;
        i = k + 1;
        while i < l {
            let x: f64 = m[i as usize].a[0] * pivot_scale;
            m[k as usize].al[(i - k - 1) as usize] = x;
            j = 1;
            while j < 11 {
                m[i as usize].a[(j - 1) as usize] = m[i as usize].a[j as usize] - x * m[k as usize].a[j as usize];
                j += 1
            }
            m[i as usize].a[10] = 0.0;
            i += 1
        }
        k += 1
    }
}

pub(crate) fn banbks11(m: &[BandMath], perm: &[usize], v: &mut [f64], n: usize) {
    let mut i: usize;
    let mut k: usize = 0;
    let mut l: usize = 5;
    /* forward substitution */
    while k < n {
        i = perm[k as usize];
        if i != k {
            let tmp = v[k as usize];
            v[k as usize] = v[i as usize];
            v[i as usize] = tmp;
        }
        if l < n {
            l += 1
        }
        i = k + 1;
        while i < l {
            v[i as usize] -= m[k as usize].al[(i - k - 1) as usize] * v[k as usize];
            i += 1
        }
        k += 1
    }
    /* back substitution */
    l = 1;
    i = n - 1;
    loop {
        let mut x: f64 = v[i as usize];
        k = 1;
        while k < l {
            x -= m[i as usize].a[k as usize] * v[(k + i) as usize];
            k += 1
        }
        v[i as usize] = x / m[i as usize].a[0];
        if l < 11 {
            l += 1
        }
        if i == 0 {
            break;
        }
        i -= 1
    }
}

pub(crate) fn compute_jinc(ty0: char, ty1: char) -> usize {
    if ty0 == 'o' || ty1 == 'o' || ty0 == ']' || ty1 == '[' {
        return 4;
    } else if ty0 == 'c' && ty1 == 'c' {
        return 2;
    } else if (ty0 == '{' || ty0 == 'v' || ty0 == '[') && ty1 == 'c' || ty0 == 'c' && (ty1 == '}' || ty1 == 'v' || ty1 == ']') {
        return 1;
    } else {
        return 0;
    };
}

pub fn count_vec(s: &[SpiroSegment]) -> usize {
    let mut n = 0;
    let mut iter = s.iter().peekable();
    while let Some(seg) = iter.next() {
        if let Some(next) = iter.peek() {
            n += compute_jinc(seg.ty, next.ty);
        }
    }
    n
}

fn add_mat_line(
    m: &mut [BandMath],
    v: &mut [f64],
    derivs: &mut [f64],
    x: f64,
    y: f64,
    j: usize,
    jj: Option<usize>,
    jinc: usize,
    nmat: usize,
) {
    let mut k: usize = 0;
    if let (Some(jj),) = (jj,) {
        let joff: usize = (j + 5 - jj + nmat) % nmat;
        v[jj as usize] += x;
        while k <= jinc {
            m[jj as usize].a[(joff + k) as usize] += y * derivs[k as usize];
            k += 1
        }
    };
}

/// Computes the outer loop of the Newton iteration
pub fn spiro_iter(s: &mut [SpiroSegment], m: &mut [BandMath], perm: &mut [usize], v: &mut [f64]) -> f64 {
    let cyclic = s[0].ty != '{' && s[0].ty != 'v';
    let mut i: usize = 0;
    let mut j: usize;
    let mut jj: usize;
    let nmat: usize = count_vec(s);
    let mut norm: f64;
    let n_invert: usize;
    while i < nmat {
        v[i as usize] = 0.0;
        j = 0;
        while j < 11 {
            m[i as usize].a[j as usize] = 0.0;
            j += 1
        }
        j = 0;
        while j < 5 {
            m[i as usize].al[j as usize] = 0.0;
            j += 1
        }
        i += 1
    }
    j = 0;
    if s[0].ty == 'o' {
        jj = nmat - 2;
    } else if s[0].ty == 'c' || s[0].ty == '[' || s[0].ty == ']' {
        jj = nmat - 1;
    } else {
        jj = 0;
    }
    i = 0;
    while i < s.len() - 1 {
        let ty0: char = s[i as usize].ty;
        let ty1: char = s[i as usize + 1].ty;
        let jinc: usize = compute_jinc(ty0, ty1);
        let th: f64 = s[i as usize].bend_th;
        let mut ends: [[f64; 4]; 2] = [[0.; 4]; 2];
        let mut derivs: [[[f64; 4]; 2]; 4] = [[[0.; 4]; 2]; 4];
        let mut jthl: Option<usize> = None;
        let mut jk0l: Option<usize> = None;
        let mut jk1l: Option<usize> = None;
        let mut jk2l: Option<usize> = None;
        let mut jthr: Option<usize> = None;
        let mut jk0r: Option<usize> = None;
        let mut jk1r: Option<usize> = None;
        let mut jk2r: Option<usize> = None;
        compute_pderivs(&mut s[i as usize], &mut ends, &mut derivs, jinc);
        /* constraints crossing left */
        if ty0 == 'o' || ty0 == 'c' || ty0 == '[' || ty0 == ']' {
            let fresh0 = jj;
            jj = jj + 1;
            jthl = Some(fresh0);
            jj %= nmat;
            let fresh1 = jj;
            jj = jj + 1;
            jk0l = Some(fresh1);
        }
        if ty0 == 'o' {
            jj %= nmat;
            let fresh2 = jj;
            jj = jj + 1;
            jk1l = Some(fresh2);
            let fresh3 = jj;
            jj = jj + 1;
            jk2l = Some(fresh3);
        }
        /* constraints on left */
        if (ty0 == '[' || ty0 == 'v' || ty0 == '{' || ty0 == 'c') && jinc == 4 {
            if ty0 != 'c' {
                let fresh4 = jj;
                jj = jj + 1;
                jk1l = Some(fresh4);
            }
            let fresh5 = jj;
            jj = jj + 1;
            jk2l = Some(fresh5);
        }
        /* constraints on right */
        if (ty1 == ']' || ty1 == 'v' || ty1 == '}' || ty1 == 'c') && jinc == 4 {
            if ty1 != 'c' {
                let fresh6 = jj;
                jj = jj + 1;
                jk1r = Some(fresh6);
            }
            let fresh7 = jj;
            jj = jj + 1;
            jk2r = Some(fresh7);
        }
        /* constraints crossing right */
        if ty1 == 'o' || ty1 == 'c' || ty1 == '[' || ty1 == ']' {
            jthr = Some(jj);
            jk0r = Some((jj + 1) % nmat);
        }
        if ty1 == 'o' {
            jk1r = Some((jj + 2) % nmat);
            jk2r = Some((jj + 3) % nmat);
        }

        add_mat_line(m, v, &mut derivs[0][0], th - ends[0][0], 1., j, jthl, jinc, nmat);
        add_mat_line(m, v, &mut derivs[1][0], ends[0][1], -1., j, jk0l, jinc, nmat);
        add_mat_line(m, v, &mut derivs[2][0], ends[0][2], -1., j, jk1l, jinc, nmat);
        add_mat_line(m, v, &mut derivs[3][0], ends[0][3], -1., j, jk2l, jinc, nmat);
        add_mat_line(m, v, &mut derivs[0][1], -ends[1][0], 1., j, jthr, jinc, nmat);
        add_mat_line(m, v, &mut derivs[1][1], -ends[1][1], 1., j, jk0r, jinc, nmat);
        add_mat_line(m, v, &mut derivs[2][1], -ends[1][2], 1., j, jk1r, jinc, nmat);
        add_mat_line(m, v, &mut derivs[3][1], -ends[1][3], 1., j, jk2r, jinc, nmat);

        j += jinc;
        i += 1
    }
    if cyclic {
        let u_nmat: usize = nmat.try_into().unwrap();
        m.copy_within(..u_nmat, u_nmat);
        m.copy_within(..u_nmat, 2 * u_nmat);
        v.copy_within(..u_nmat, u_nmat);
        v.copy_within(..u_nmat, 2 * u_nmat);
        n_invert = 3 * nmat;
        j = nmat
    } else {
        n_invert = nmat;
        j = 0
    }
    bandec11(m, perm, n_invert);
    banbks11(m, perm, v, n_invert);
    norm = 0.0;
    i = 0;
    while i < s.len() - 1 {
        let ty0_0: char = s[i as usize].ty;
        let ty1_0: char = s[(i + 1) as usize].ty;
        let jinc_0: usize = compute_jinc(ty0_0, ty1_0);
        let mut k: usize = 0;
        while k < jinc_0 {
            let fresh8 = j;
            j = j + 1;
            let dk: f64 = v[fresh8 as usize];
            s[i as usize].ks[k as usize] += dk;
            norm += dk * dk;
            k += 1
        }
        i += 1
    }
    return norm;
}

pub(crate) fn spiro_seg_to_bpath<T, A>(ks: [f64; 4], x0: f64, y0: f64, x1: f64, y1: f64, bc: &mut BezierContext<T, A>, depth: isize) {
    let bend: f64 = (ks[0]).abs() + (0.5 * ks[1]).abs() + (0.125 * ks[2]).abs() + (1.0 / 48. * ks[3]).abs();
    if bend == 0. {
        bc.line_to(x1, y1);
    } else {
        let seg_ch: f64 = (x1 - x0).hypot(y1 - y0);
        let seg_th: f64 = (y1 - y0).atan2(x1 - x0);
        let mut xy: [f64; 2] = [0.; 2];
        let ch: f64;
        let th: f64;
        let scale: f64;
        let rot: f64;
        let th_even: f64;
        let th_odd: f64;
        let ul: f64;
        let vl: f64;
        let ur: f64;
        let vr: f64;
        integrate_spiro(ks, &mut xy);
        ch = xy[0].hypot(xy[1]);
        th = xy[1].atan2(xy[0]);
        scale = seg_ch / ch;
        rot = seg_th - th;
        if depth > 5 || bend < 1.0 {
            th_even = 1.0 / 384. * ks[3] + 1.0 / 8. * ks[1] + rot;
            th_odd = 1.0 / 48. * ks[2] + 0.5 * ks[0];
            ul = scale * (1.0 / 3.) * (th_even - th_odd).cos();
            vl = scale * (1.0 / 3.) * (th_even - th_odd).sin();
            ur = scale * (1.0 / 3.) * (th_even + th_odd).cos();
            vr = scale * (1.0 / 3.) * (th_even + th_odd).sin();
            bc.curve_to(x0 + ul, y0 + vl, x1 - ur, y1 - vr, x1, y1);
        } else {
            /* subdivide */
            let mut ksub: [f64; 4] = [0.; 4];
            let thsub: f64;
            let mut xysub: [f64; 2] = [0.; 2];
            let xmid: f64;
            let ymid: f64;
            let cth: f64;
            let sth: f64;
            ksub[0] = 0.5 * ks[0] - 0.125 * ks[1] + 1.0 / 64. * ks[2] - 1.0 / 768. * ks[3];
            ksub[1] = 0.25 * ks[1] - 1.0 / 16. * ks[2] + 1.0 / 128. * ks[3];
            ksub[2] = 0.125 * ks[2] - 1.0 / 32. * ks[3];
            ksub[3] = 1.0 / 16. * ks[3];
            thsub = rot - 0.25 * ks[0] + 1.0 / 32. * ks[1] - 1.0 / 384. * ks[2] + 1.0 / 6144. * ks[3];
            cth = 0.5 * scale * (thsub).cos();
            sth = 0.5 * scale * (thsub).sin();
            integrate_spiro(ksub, &mut xysub);
            xmid = x0 + cth * xysub[0] - sth * xysub[1];
            ymid = y0 + cth * xysub[1] + sth * xysub[0];
            spiro_seg_to_bpath(ksub, x0, y0, xmid, ymid, bc, depth + 1);
            ksub[0] += 0.25 * ks[1] + 1.0 / 384. * ks[3];
            ksub[1] += 0.125 * ks[2];
            ksub[2] += 1.0 / 16. * ks[3];
            spiro_seg_to_bpath(ksub, xmid, ymid, x1, y1, bc, depth + 1);
        }
    };
}

/// Get knot theta (ϑ) at index
pub fn get_knot_theta(s: &[SpiroSegment], i: usize) -> f64 {
    let mut ends: [[f64; 4]; 2] = [[0.; 4]; 2];
    if i == 0 {
        compute_ends(s[i as usize].ks, &mut ends, s[i].seg_ch);
        return s[i].seg_th - ends[0][0];
    } else {
        compute_ends(s[(i - 1)].ks, &mut ends, s[(i - 1)].seg_ch);
        return s[(i - 1)].seg_th + ends[1][0];
    };
}
