use subtle::Choice;

use crate::{fp::Fp, fp12::Fp12, fp2::Fp2, G1Affine, G2Affine, Scalar, BLS_X, BLS_X_IS_NEGATIVE};

// pub(crate) const BLS_BETA: Fp2 = Fp2 {
//     c0: Fp::from_raw_unchecked([1, 0, 0, 0, 0, 0]),
//     c1: Fp::from_raw_unchecked([1, 0, 0, 0, 0, 0]),
// };

const FROBENIUS_COEFFS: [Fp2; 3] = [
    Fp2 {
        c0: Fp([
            0x14fe_c701_e8fb_0ce9,
            0xed5e_6427_3c4f_538b,
            0x1797_ab14_58a8_8de9,
            0x343e_a979_1495_6dc8,
            0x7fe1_1274_d898_fafb,
            0xf4d3_8259_380b_4820,
        ]),
        c1: Fp([
            0x0502_4ae8_5084_d9b0,
            0x5dbd_438f_06fc_594c,
            0x4cdf_a070_9adc_84d6,
            0x32f2_2927_e21b_885b,
            0x9eca_ed89_d8bb_0503,
            0xc52b_7da6_c7f4_628b,
        ]),
    },
    Fp2 {
        c0: Fp([
            0x14fe_c701_e8fb_0ce9,
            0xed5e_6427_3c4f_538b,
            0x1797_ab14_58a8_8de9,
            0x343e_a979_1495_6dc8,
            0x7fe1_1274_d898_fafb,
            0xf4d3_8259_380b_4820,
        ]),
        c1: Fp([
            0x14fe_c701_e8fb_0ce9,
            0xed5e_6427_3c4f_538b,
            0x1797_ab14_58a8_8de9,
            0x343e_a979_1495_6dc8,
            0x7fe1_1274_d898_fafb,
            0xf4d3_8259_380b_4820,
        ]),
    },
    Fp2 {
        c0: Fp([
            0x14fe_c701_e8fb_0ce9,
            0xed5e_6427_3c4f_538b,
            0x1797_ab14_58a8_8de9,
            0x343e_a979_1495_6dc8,
            0x7fe1_1274_d898_fafb,
            0xf4d3_8259_380b_4820,
        ]),
        c1: Fp([
            0x0502_4ae8_5084_d9b0,
            0x5dbd_438f_06fc_594c,
            0x4cdf_a070_9adc_84d6,
            0x32f2_2927_e21b_885b,
            0x9eca_ed89_d8bb_0503,
            0xc52b_7da6_c7f4_628b,
        ]),
    },
];

/// 6U+2 for in NAF form
pub const SIX_U_PLUS_2_NAF: [i8; 65] = [
    0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 0,
    1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 0,
    0, 1, 0, 1, 1,
];

fn line_evaluation(f: &mut Fp12, lambda: Fp2, p1: &G1Affine, t_x: Fp2, t_y: Fp2) {
    let c1 = Fp2::new(lambda.c0 * &p1.x, lambda.c1 * &p1.x);

    let t = lambda * t_x - t_y;
    let c3 = Fp2::new(t.c0 * &p1.y, t.c1 * &p1.y);

    f.mul_by_034(&Fp2::one(), &c1, &c3);
}

fn add_eval(f: &mut Fp12, acc: &mut G2Affine, p2: &G2Affine, p1: &G1Affine, neg: bool) {
    let t_x = (*acc).x;
    let t_y = (*acc).y;
    let lambda = add(acc, p2, neg);

    line_evaluation(f, lambda, p1, t_x, t_y);
}

pub(crate) fn add(acc: &mut G2Affine, p2: &G2Affine, neg: bool) -> Fp2 {
    let t0 = if neg { acc.y + p2.y } else { acc.y - p2.y };
    let t1 = (acc.x - p2.x).invert().unwrap();
    let lambda = t0 * t1;
    let x3 = lambda.square() - acc.x - p2.x;
    let y3 = lambda * (acc.x - x3) - acc.y;

    acc.x = x3;
    acc.y = y3;
    lambda
}

fn double_eval(f: &mut Fp12, acc: &mut G2Affine, p1: &G1Affine) {
    let t_x = (*acc).x;
    let t_y = (*acc).y;
    let lambda = double(acc);

    line_evaluation(f, lambda, p1, t_x, t_y);
}

pub(crate) fn double(acc: &mut G2Affine) -> Fp2 {
    let x2 = acc.x.square();
    let t0 = x2.double() + x2;
    let t1 = acc.y.double().invert().unwrap();
    let lambda = t0 * t1;
    let x3 = lambda.square() - acc.x.double();
    let y3 = lambda * (acc.x - x3) - acc.y;

    acc.x = x3;
    acc.y = y3;
    lambda
}

/// Perform the Miller loop.
fn miller_loop(p: &G1Affine, q: &G2Affine) -> Fp12 {
    let mut f = Fp12::one();
    let mut acc = (*q).clone();

    SIX_U_PLUS_2_NAF.iter().enumerate().for_each(|(i, naf)| {
        (i != 0).then(|| f = f.square());
        double_eval(&mut f, &mut acc, p);

        if naf.pow(2) == 1 {
            add_eval(&mut f, &mut acc, q, p, *naf == -1);
        }
    });

    let x = if BLS_X_IS_NEGATIVE {
        -Scalar::from(BLS_X)
    } else {
        Scalar::from(BLS_X)
    };

    assert_eq!(
        acc,
        (q * Scalar::from(Scalar::from(6) * x + Scalar::from(2))).into()
    );

    FROBENIUS_COEFFS.iter().for_each(|frob| {
        let _q = G2Affine {
            x: q.x * frob,
            y: q.y * frob,
            infinity: Choice::from(0),
        };
        add_eval(&mut f, &mut acc, &_q, &p, false)
    });

    f
}

/// a^((p - 1) / 2)
fn legendre(a: Fp12) -> Fp12 {
    a.pow_vartime(&[
        0x0d00_88f5_1cbf_f34d,
        0x258d_d3db_21a5_d66b,
        0xb23b_a5c2_79c2_895f,
        0xb398_6950_7b58_7b12,
        0x0f55_ffff_58a9_ffff,
        0xdcff_7fff_ffff_d555,
    ])
}

/// Toneelli-Shanks algorithm to compute cubic roots in a finite field.
/// p - 1 = 2^s * t where s == 0 in this case
/// if s == 1:
///    return pow(n, (p + 1) // 4, p)
pub(crate) fn tonelli_shanks(a: Fp12) -> Fp12 {
    assert_eq!(
        legendre(a),
        Fp12::one(),
        "not a cubic residue, cannot compute cubic root",
    );
    // let s = 2; // s = 2^s * t
    let t = Fp([
        0x02e3_aca8_3f47_199f,
        0x5dad_bd4d_23eb_f6c2,
        0x9962_969c_fe9d_0215,
        0x445a_c211_e285_70ae,
        0xae13_1c71_a1ec_e38e,
        0x311c_5555_5555_4bda,
    ]);

    // t = 3k +/- 1
    let k = Fp([
        0x00f6_8ee2_bfc2_5ddf,
        0xc9e4_946f_0bf9_5240,
        0xddcb_8789_aa34_5607,
        0x16c8_eb5b_4b81_d03a,
        0x3a06_5ed0_8b4e_f684,
        0xbb09_71c7_1c71_c3f3,
    ]);

    let three_times_k_plus_one = Fp([
        0x02e3_aca8_3f47_199f,
        0x5dad_bd4d_23eb_f6c2,
        0x9962_969c_fe9d_0215,
        0x445a_c211_e285_70ae,
        0xae13_1c71_a1ec_e38e,
        0x311c_5555_5555_4bda,
    ]);
    let c = Fp12::random(&mut rand::thread_rng()).pow_vartime(&t.0); // c = b^t
    let mut r = a.pow_vartime(&t.0); // r = a^t

    let mut h = Fp12::one();
    let mut c = c.invert().unwrap();
    let c_3 = c.square() * c;

    let d = r.pow_vartime(&[1, 0, 0, 0, 0, 0]); // d = r^(3^(s - i - 1))

    // in this case, c' = c^(3^(s - i - 1)) = c^3
    if d == c_3 {
        h = h * c;
        r = r * c_3;
    } else if d != Fp12::one() {
        h = h * c.square();
        r = r * c_3.square();
    }

    r = a.pow_vartime(&k.0) * h;
    if t == three_times_k_plus_one {
        r = r.invert().unwrap()
    }

    r
}

#[test]
fn test_tonelli_shanks() {
    let a = Fp12::random(&mut rand::thread_rng());
    let cubic_root = tonelli_shanks(a);
    assert_eq!(cubic_root.square() * cubic_root, a);
}
