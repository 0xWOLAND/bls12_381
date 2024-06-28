use crate::{fp12::Fp12, fp2::Fp2, G1Affine, G2Affine};

fn eval(f: &mut Fp12, lambda: Fp2, p1: &G1Affine, t_x: Fp2, t_y: Fp2) {
    let c1 = Fp2::new(lambda.c0 * &p1.x, lambda.c1 * &p1.x);

    let t = lambda * t_x - t_y;
    let c3 = Fp2::new(t.c0 * &p1.y, t.c1 * &p1.y);

    f.mul_by_034(&Fp2::one(), &c1, &c3);
}

fn add_eval(f: &mut Fp12, acc: &mut G2Affine, p2: &G2Affine, p1: &G1Affine, neg: bool) {
    let t_x = (*acc).x;
    let t_y = (*acc).y;
    let lambda = add(acc, p2, neg);

    eval(f, lambda, p1, t_x, t_y);
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

    eval(f, lambda, p1, t_x, t_y);
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
