use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use curves::weierstrass::{WeierstrassCurve, Point as WeierstrassPoint};
use curves::montgomery::{MontgomeryCurve, Point as MontgomeryPoint};
use curves::twisted_edwards::{TwistedEdwardsCurve, Point as TwistedEdwardsPoint};
use curves::traits::{CurvePoint, EllipticCurve};

// Benchmark configurations
const BENCHMARK_PRIME: i32 = 13; // Use smaller prime for faster point enumeration 
const SCALAR_VALUES: &[i32] = &[3, 5, 7, 11]; // Different scalar multiplication values

fn setup_weierstrass_curve() -> (WeierstrassCurve, WeierstrassPoint, WeierstrassPoint) {
    let curve = WeierstrassCurve::new(1, 1, BENCHMARK_PRIME);
    let all_points = curve.get_all_points();
    let finite_points: Vec<_> = all_points.iter()
        .filter(|p| !p.is_infinity_point())
        .cloned()
        .collect();
    
    // Return the curve and first two finite points
    (curve, finite_points[0].clone(), finite_points[1].clone())
}

fn setup_montgomery_curve() -> (MontgomeryCurve, MontgomeryPoint, MontgomeryPoint) {
    let curve = MontgomeryCurve::new(1, 1, BENCHMARK_PRIME).unwrap();
    let all_points = curve.get_all_points();
    let finite_points: Vec<_> = all_points.iter()
        .filter(|p| !p.is_infinity_point())
        .cloned()
        .collect();
    
    // Return the curve and first two finite points
    (curve, finite_points[0].clone(), finite_points[1].clone())
}

fn setup_twisted_edwards_curve() -> (TwistedEdwardsCurve, TwistedEdwardsPoint, TwistedEdwardsPoint) {
    let curve = TwistedEdwardsCurve::new(1, 2, BENCHMARK_PRIME).unwrap();
    let all_points = curve.get_all_points();
    let finite_points: Vec<_> = all_points.iter()
        .filter(|p| !p.is_infinity_point())
        .cloned()
        .collect();
    
    // Return the curve and first two finite points
    (curve, finite_points[0].clone(), finite_points[1].clone())
}

fn bench_point_addition(c: &mut Criterion) {
    let mut group = c.benchmark_group("point_addition");
    
    // Weierstrass
    let (weierstrass_curve, w_p1, w_p2) = setup_weierstrass_curve();
    group.bench_with_input(
        BenchmarkId::new("Weierstrass", "point_add"),
        &(&weierstrass_curve, &w_p1, &w_p2),
        |b, (curve, p1, p2)| {
            b.iter(|| curve.point_add(black_box(p1), black_box(p2)))
        },
    );
    
    // Montgomery
    let (montgomery_curve, m_p1, m_p2) = setup_montgomery_curve();
    group.bench_with_input(
        BenchmarkId::new("Montgomery", "point_add"),
        &(&montgomery_curve, &m_p1, &m_p2),
        |b, (curve, p1, p2)| {
            b.iter(|| curve.point_add(black_box(p1), black_box(p2)))
        },
    );
    
    // Twisted Edwards
    let (twisted_edwards_curve, t_p1, t_p2) = setup_twisted_edwards_curve();
    group.bench_with_input(
        BenchmarkId::new("TwistedEdwards", "point_add"),
        &(&twisted_edwards_curve, &t_p1, &t_p2),
        |b, (curve, p1, p2)| {
            b.iter(|| curve.point_add(black_box(p1), black_box(p2)))
        },
    );
    
    group.finish();
}

fn bench_point_doubling(c: &mut Criterion) {
    let mut group = c.benchmark_group("point_doubling");
    
    // Weierstrass
    let (weierstrass_curve, w_p1, _) = setup_weierstrass_curve();
    group.bench_with_input(
        BenchmarkId::new("Weierstrass", "point_double"),
        &(&weierstrass_curve, &w_p1),
        |b, (curve, p)| {
            b.iter(|| curve.point_double(black_box(p)))
        },
    );
    
    // Montgomery
    let (montgomery_curve, m_p1, _) = setup_montgomery_curve();
    group.bench_with_input(
        BenchmarkId::new("Montgomery", "point_double"),
        &(&montgomery_curve, &m_p1),
        |b, (curve, p)| {
            b.iter(|| curve.point_double(black_box(p)))
        },
    );
    
    // Twisted Edwards
    let (twisted_edwards_curve, t_p1, _) = setup_twisted_edwards_curve();
    group.bench_with_input(
        BenchmarkId::new("TwistedEdwards", "point_double"),
        &(&twisted_edwards_curve, &t_p1),
        |b, (curve, p)| {
            b.iter(|| curve.point_double(black_box(p)))
        },
    );
    
    group.finish();
}

fn bench_scalar_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("scalar_multiplication");
    
    // Setup curves and points
    let (weierstrass_curve, w_p1, _) = setup_weierstrass_curve();
    let (montgomery_curve, m_p1, _) = setup_montgomery_curve();
    let (twisted_edwards_curve, t_p1, _) = setup_twisted_edwards_curve();
    
    for &scalar in SCALAR_VALUES {
        // Weierstrass
        group.bench_with_input(
            BenchmarkId::new("Weierstrass", scalar),
            &(&weierstrass_curve, scalar, &w_p1),
            |b, (curve, k, p)| {
                b.iter(|| curve.scalar_multiply(black_box(*k), black_box(p)))
            },
        );
        
        // Montgomery
        group.bench_with_input(
            BenchmarkId::new("Montgomery", scalar),
            &(&montgomery_curve, scalar, &m_p1),
            |b, (curve, k, p)| {
                b.iter(|| curve.scalar_multiply(black_box(*k), black_box(p)))
            },
        );
        
        // Twisted Edwards
        group.bench_with_input(
            BenchmarkId::new("TwistedEdwards", scalar),
            &(&twisted_edwards_curve, scalar, &t_p1),
            |b, (curve, k, p)| {
                b.iter(|| curve.scalar_multiply(black_box(*k), black_box(p)))
            },
        );
    }
    
    group.finish();
}

fn bench_comprehensive_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("comprehensive_operations");
    
    // Setup curves and points
    let (weierstrass_curve, w_p1, w_p2) = setup_weierstrass_curve();
    let (montgomery_curve, m_p1, m_p2) = setup_montgomery_curve();
    let (twisted_edwards_curve, t_p1, t_p2) = setup_twisted_edwards_curve();
    
    // Mixed operations: point addition + doubling + small scalar multiplication
    group.bench_with_input(
        BenchmarkId::new("Weierstrass", "mixed_ops"),
        &(&weierstrass_curve, &w_p1, &w_p2),
        |b, (curve, p1, p2)| {
            b.iter(|| {
                let sum = curve.point_add(black_box(p1), black_box(p2)).unwrap();
                let doubled = curve.point_double(black_box(p1)).unwrap();
                let scalar_result = curve.scalar_multiply(black_box(7), black_box(p1)).unwrap();
                (sum, doubled, scalar_result)
            })
        },
    );
    
    group.bench_with_input(
        BenchmarkId::new("Montgomery", "mixed_ops"),
        &(&montgomery_curve, &m_p1, &m_p2),
        |b, (curve, p1, p2)| {
            b.iter(|| {
                let sum = curve.point_add(black_box(p1), black_box(p2)).unwrap();
                let doubled = curve.point_double(black_box(p1)).unwrap();
                let scalar_result = curve.scalar_multiply(black_box(7), black_box(p1)).unwrap();
                (sum, doubled, scalar_result)
            })
        },
    );
    
    group.bench_with_input(
        BenchmarkId::new("TwistedEdwards", "mixed_ops"),
        &(&twisted_edwards_curve, &t_p1, &t_p2),
        |b, (curve, p1, p2)| {
            b.iter(|| {
                let sum = curve.point_add(black_box(p1), black_box(p2)).unwrap();
                let doubled = curve.point_double(black_box(p1)).unwrap();
                let scalar_result = curve.scalar_multiply(black_box(7), black_box(p1)).unwrap();
                (sum, doubled, scalar_result)
            })
        },
    );
    
    group.finish();
}

criterion_group!(
    benches,
    bench_point_addition,
    bench_point_doubling,
    bench_scalar_multiplication,
    bench_comprehensive_operations
);
criterion_main!(benches);