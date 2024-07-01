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
pub fn miller_loop(p: &G1Affine, q: &G2Affine) -> Fp12 {
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
