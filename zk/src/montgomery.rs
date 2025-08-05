use crate::traits::{CurvePoint, EllipticCurve};

/// Represents a point on a Montgomery curve
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Point {
    pub x: i32,
    pub y: i32,
    pub is_infinity: bool,
}

impl CurvePoint for Point {
    fn new(x: i32, y: i32) -> Self {
        Point { x, y, is_infinity: false }
    }

    fn infinity() -> Self {
        Point { x: 0, y: 0, is_infinity: true }
    }

    fn is_infinity_point(&self) -> bool {
        self.is_infinity
    }
    
    fn x(&self) -> i32 {
        self.x
    }
    
    fn y(&self) -> i32 {
        self.y
    }
}

/// Represents a Montgomery curve B·y² = x³ + A·x² + x over a prime field
/// Constraints: B ≠ 0, A² ≠ 4 (mod p), and p > 3
#[derive(Debug)]
pub struct MontgomeryCurve {
    pub a: i32,
    pub b: i32,
    pub prime: i32,
}

impl MontgomeryCurve {
    /// Creates a new Montgomery curve with the given parameters
    /// Validates the constraints: B ≠ 0, A² ≠ 4 (mod p), and p > 3
    pub fn new(a: i32, b: i32, prime: i32) -> Option<Self> {
        // Check constraints
        if prime <= 3 {
            return None; // p must be > 3
        }
        
        if b == 0 {
            return None; // B ≠ 0
        }
        
        let a_squared = (a * a) % prime;
        if a_squared == 4 {
            return None; // A² ≠ 4 (mod p)
        }
        
        Some(MontgomeryCurve { a, b, prime })
    }



    /// Checks if a point lies on the curve
    pub fn is_on_curve(&self, point: &Point) -> bool {
        // The point at infinity is always considered to be on the curve
        if point.is_infinity {
            return true;
        }
        
        let x = point.x;
        let y = point.y;
        
        // Calculate B·y² mod prime
        let mut by_squared = (self.b * y * y) % self.prime;
        if by_squared < 0 {
            by_squared += self.prime;
        }
        
        // Calculate x³ + A·x² + x mod prime
        let x_cubed = (x * x * x) % self.prime;
        let ax_squared = (self.a * x * x) % self.prime;
        let mut right_side = (x_cubed + ax_squared + x) % self.prime;
        
        // Handle negative values
        if right_side < 0 {
            right_side += self.prime;
        }
        
        by_squared == right_side
    }

    /// Point doubling for Montgomery curves
    /// For P = (x, y) with y ≠ 0:
    /// x' = ((3x² + 2Ax + 1)/(2By))² - 2x - A
    /// y' = (3x² + 2Ax + 1)/(2By) * (x - x') - y
    pub fn point_double(&self, point: &Point) -> Option<Point> {
        if point.is_infinity {
            return Some(Point::infinity());
        }
        
        let x = point.x;
        let y = point.y;
        
        // Check if y = 0 (vertical tangent)
        if y == 0 {
            return Some(Point::infinity());
        }
        
        // Calculate 3x² + 2Ax + 1
        let three_x_squared = (3 * x * x) % self.prime;
        let two_ax = (2 * self.a * x) % self.prime;
        let numerator = (three_x_squared + two_ax + 1) % self.prime;
        
        // Calculate 2By
        let two_by = (2 * self.b * y) % self.prime;
        
        // Calculate λ = (3x² + 2Ax + 1)/(2By)
        let lambda = self.modular_divide(numerator, two_by)?;
        
        // Calculate x' = λ² - 2x - A
        let lambda_squared = (lambda * lambda) % self.prime;
        let x_prime = (lambda_squared - 2 * x - self.a) % self.prime;
        
        // Handle negative values
        let x_prime = if x_prime < 0 { x_prime + self.prime } else { x_prime };
        
        // Calculate y' = λ(x - x') - y
        let x_diff = (x - x_prime) % self.prime;
        let x_diff = if x_diff < 0 { x_diff + self.prime } else { x_diff };
        let y_prime = (lambda * x_diff - y) % self.prime;
        let y_prime = if y_prime < 0 { y_prime + self.prime } else { y_prime };
        
        Some(Point::new(x_prime, y_prime))
    }

    /// Point addition for Montgomery curves
    /// For P = (x1, y1) and Q = (x2, y2) with x1 ≠ x2:
    /// x' = ((y2 - y1)/(x2 - x1))² * B - (x1 + x2) - A
    /// y' = (y2 - y1)/(x2 - x1) * (x1 - x') - y1
    pub fn point_add(&self, p: &Point, q: &Point) -> Option<Point> {
        // Handle point at infinity cases
        if p.is_infinity {
            return Some(q.clone());
        }
        if q.is_infinity {
            return Some(p.clone());
        }
        
        let x1 = p.x;
        let y1 = p.y;
        let x2 = q.x;
        let y2 = q.y;
        
        // Check if points are the same (use doubling)
        if x1 == x2 {
            if y1 == y2 {
                return self.point_double(p);
            } else {
                return Some(Point::infinity()); // Points are inverses
            }
        }
        
        // Calculate y2 - y1
        let y_diff = (y2 - y1) % self.prime;
        let y_diff = if y_diff < 0 { y_diff + self.prime } else { y_diff };
        
        // Calculate x2 - x1
        let x_diff = (x2 - x1) % self.prime;
        let x_diff = if x_diff < 0 { x_diff + self.prime } else { x_diff };
        
        // Calculate λ = (y2 - y1)/(x2 - x1)
        let lambda = self.modular_divide(y_diff, x_diff)?;
        
        // Calculate x' = λ² * B - (x1 + x2) - A
        let lambda_squared = (lambda * lambda) % self.prime;
        let lambda_squared_b = (lambda_squared * self.b) % self.prime;
        let x_sum = (x1 + x2) % self.prime;
        let x_prime = (lambda_squared_b - x_sum - self.a) % self.prime;
        let x_prime = if x_prime < 0 { x_prime + self.prime } else { x_prime };
        
        // Calculate y' = λ(x1 - x') - y1
        let x1_minus_x_prime = (x1 - x_prime) % self.prime;
        let x1_minus_x_prime = if x1_minus_x_prime < 0 { x1_minus_x_prime + self.prime } else { x1_minus_x_prime };
        let y_prime = (lambda * x1_minus_x_prime - y1) % self.prime;
        let y_prime = if y_prime < 0 { y_prime + self.prime } else { y_prime };
        
        Some(Point::new(x_prime, y_prime))
    }

    /// Returns the inverse of a point: -P = (x, -y)
    pub fn point_inverse(&self, point: &Point) -> Point {
        if point.is_infinity {
            return Point::infinity();
        }
        
        let inverse_y = (self.prime - point.y) % self.prime;
        Point::new(point.x, inverse_y)
    }

    /// Returns all possible points on the curve
    pub fn get_all_points(&self) -> Vec<Point> {
        let mut points = Vec::new();
        
        // Add the point at infinity (neutral element)
        points.push(Point::infinity());
        
        // Check all possible x and y values in the field
        for x in 0..self.prime {
            for y in 0..self.prime {
                let point = Point::new(x, y);
                if self.is_on_curve(&point) {
                    points.push(point);
                }
            }
        }
        
        points
    }

    /// Returns the number of points on the curve
    pub fn point_count(&self) -> usize {
        self.get_all_points().len()
    }

    /// Modular division: a/b mod p
    fn modular_divide(&self, a: i32, b: i32) -> Option<i32> {
        let b_inv = self.modular_inverse(b)?;
        let result = (a * b_inv) % self.prime;
        let result = if result < 0 { result + self.prime } else { result };
        Some(result)
    }

    /// Calculates modular multiplicative inverse using extended Euclidean algorithm
    fn modular_inverse(&self, a: i32) -> Option<i32> {
        let mut t = 0;
        let mut new_t = 1;
        let mut r = self.prime;
        let mut new_r = a;
        
        while new_r != 0 {
            let quotient = r / new_r;
            let temp_t = t;
            t = new_t;
            new_t = temp_t - quotient * new_t;
            let temp_r = r;
            r = new_r;
            new_r = temp_r - quotient * new_r;
        }
        
        if r > 1 {
            return None; // No inverse exists
        }
        
        if t < 0 {
            t += self.prime;
        }
        
        Some(t)
    }
}

impl EllipticCurve for MontgomeryCurve {
    type Point = Point;
    
    fn new(a: i32, b: i32, prime: i32) -> Self {
        MontgomeryCurve::new(a, b, prime).unwrap()
    }
    
    fn a(&self) -> i32 {
        self.a
    }
    
    fn b(&self) -> i32 {
        self.b
    }
    
    fn prime(&self) -> i32 {
        self.prime
    }
    
    fn is_on_curve(&self, point: &Self::Point) -> bool {
        self.is_on_curve(point)
    }
    
    fn point_double(&self, point: &Self::Point) -> Option<Self::Point> {
        self.point_double(point)
    }
    
    fn point_add(&self, p: &Self::Point, q: &Self::Point) -> Option<Self::Point> {
        self.point_add(p, q)
    }
    
    fn point_inverse(&self, point: &Self::Point) -> Self::Point {
        self.point_inverse(point)
    }
    
    fn scalar_multiply(&self, k: i32, point: &Self::Point) -> Option<Self::Point> {
        // Montgomery curves don't have a direct scalar multiplication implementation
        // We'll implement it using double-and-add like Weierstrass
        if k == 0 {
            return Some(Point::infinity());
        }
        
        let mut result = Point::infinity();
        let mut current = point.clone();
        let mut scalar = k.abs();
        
        while scalar > 0 {
            if scalar & 1 == 1 {
                result = self.point_add(&result, &current)?;
            }
            current = self.point_double(&current)?;
            scalar >>= 1;
        }
        
        // Handle negative scalar
        if k < 0 {
            result = self.point_inverse(&result);
        }
        
        Some(result)
    }
    
    fn get_all_points(&self) -> Vec<Self::Point> {
        self.get_all_points()
    }
    

    
    fn equation(&self) -> String {
        format!("{}y² = x³ + {}x² + x (mod {})", self.b, self.a, self.prime)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_montgomery_curve_creation() {
        // Valid curve
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        assert_eq!(curve.a, 1);
        assert_eq!(curve.b, 1);
        assert_eq!(curve.prime, 13);
        
        // Invalid curve: A² = 4 mod p
        let invalid_curve = MontgomeryCurve::new(2, 1, 13);
        assert!(invalid_curve.is_none());
        
        // Invalid curve: B = 0
        let invalid_curve2 = MontgomeryCurve::new(1, 0, 13);
        assert!(invalid_curve2.is_none());
        
        // Invalid curve: p ≤ 3
        let invalid_curve3 = MontgomeryCurve::new(1, 1, 3);
        assert!(invalid_curve3.is_none());
    }

    #[test]
    fn test_montgomery_constraints() {
        // Test valid curve
        let curve = MontgomeryCurve::new(1, 1, 13);
        assert!(curve.is_some());
        
        // Test invalid curves
        let invalid_curve1 = MontgomeryCurve::new(2, 1, 13);
        assert!(invalid_curve1.is_none()); // A² = 4 mod 13
        
        let invalid_curve2 = MontgomeryCurve::new(1, 0, 13);
        assert!(invalid_curve2.is_none()); // B = 0
        
        let invalid_curve3 = MontgomeryCurve::new(1, 1, 3);
        assert!(invalid_curve3.is_none()); // p ≤ 3
    }

    #[test]
    fn test_point_on_curve() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        
        // Test some points that should be on the curve y² = x³ + x² + x (mod 13)
        let point1 = Point::new(0, 0); // 0² = 0, 0³ + 0² + 0 = 0 ✓
        let point2 = Point::new(1, 1); // 1² = 1, 1³ + 1² + 1 = 1 + 1 + 1 = 3 ✗
        let point3 = Point::new(2, 2); // 2² = 4, 2³ + 2² + 2 = 8 + 4 + 2 = 14 ≡ 1 (mod 13) ✗
        let point4 = Point::new(3, 0); // 0² = 0, 3³ + 3² + 3 = 27 + 9 + 3 = 39 ≡ 0 (mod 13) ✓
        
        assert!(curve.is_on_curve(&point1));
        assert!(!curve.is_on_curve(&point2));
        assert!(!curve.is_on_curve(&point3));
        assert!(curve.is_on_curve(&point4));
    }

    #[test]
    fn test_all_points() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        let points = curve.get_all_points();
        
        // For the curve y² = x³ + x² + x (mod 13), some valid points are:
        // x=0: y² = 0, so y = 0
        // x=3: y² = 0, so y = 0
        // x=6: y² = 0, so y = 0
        // x=9: y² = 0, so y = 0
        // Plus the point at infinity
        
        assert!(points.len() > 0);
        assert!(points.iter().any(|p| p.is_infinity_point()));
        
        // Check that all result points are on the curve
        for point in &points {
            assert!(curve.is_on_curve(point), 
                   "Point {:?} is not on the curve", point);
        }
        
        // Check that there's exactly one point at infinity
        let infinity_count = points.iter().filter(|p| p.is_infinity_point()).count();
        assert_eq!(infinity_count, 1);
    }

    #[test]
    fn test_point_count() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        let count = curve.point_count();
        assert!(count > 0);
        println!("Montgomery curve has {} points", count);
    }

    #[test]
    fn test_point_at_infinity() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        
        // Test point at infinity creation
        let infinity = Point::infinity();
        assert!(infinity.is_infinity_point());
        assert!(infinity.is_infinity);
        
        // Test that point at infinity is on the curve
        assert!(curve.is_on_curve(&infinity));
        
        // Test that point at infinity is included in all points
        let all_points = curve.get_all_points();
        assert!(all_points.iter().any(|p| p.is_infinity_point()));
        
        // Test that there's exactly one point at infinity
        let infinity_count = all_points.iter().filter(|p| p.is_infinity_point()).count();
        assert_eq!(infinity_count, 1);
    }

    #[test]
    fn test_point_doubling() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        
        // Find a valid point on the curve
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if !finite_points.is_empty() {
            let p = finite_points[0].clone();
            let doubled = curve.point_double(&p).unwrap();
            
            // Verify the doubled point is on the curve
            assert!(curve.is_on_curve(&doubled));
        }
        
        // Test doubling the point at infinity
        let infinity = Point::infinity();
        let doubled_infinity = curve.point_double(&infinity).unwrap();
        assert!(doubled_infinity.is_infinity_point());
    }

    #[test]
    fn test_point_addition() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        
        // Find two different valid points on the curve
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if finite_points.len() >= 2 {
            let p = finite_points[0].clone();
            let q = finite_points[1].clone();
            let sum = curve.point_add(&p, &q).unwrap();
            
            // Verify the sum is on the curve
            assert!(curve.is_on_curve(&sum));
            
            // Test adding a point to itself (should use doubling)
            let doubled = curve.point_double(&p).unwrap();
            let added = curve.point_add(&p, &p).unwrap();
            assert_eq!(doubled, added);
        }
        
        // Test adding a point to infinity
        let infinity = Point::infinity();
        if !finite_points.is_empty() {
            let p = finite_points[0].clone();
            let result = curve.point_add(&p, &infinity).unwrap();
            assert_eq!(result, p);
        }
    }

    #[test]
    fn test_group_law_properties() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if finite_points.len() >= 3 {
            let p = finite_points[0].clone();
            let q = finite_points[1].clone();
            let r = finite_points[2].clone();
            
            // Test associativity: (P + Q) + R = P + (Q + R)
            let left = curve.point_add(&curve.point_add(&p, &q).unwrap(), &r).unwrap();
            let right = curve.point_add(&p, &curve.point_add(&q, &r).unwrap()).unwrap();
            assert_eq!(left, right);
            
            // Test commutativity: P + Q = Q + P
            let sum1 = curve.point_add(&p, &q).unwrap();
            let sum2 = curve.point_add(&q, &p).unwrap();
            assert_eq!(sum1, sum2);
            
            // Test identity: P + ∞ = P
            let infinity = Point::infinity();
            let identity_test = curve.point_add(&p, &infinity).unwrap();
            assert_eq!(identity_test, p);
        }
    }

    #[test]
    fn test_infinity_rules() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if !finite_points.is_empty() {
            let p = finite_points[0].clone();
            let infinity = Point::infinity();
            
            // Test neutral element rule: P ⊕ O = P
            let result1 = curve.point_add(&p, &infinity).unwrap();
            assert_eq!(result1, p);
            
            // Test neutral element rule: O ⊕ P = P
            let result2 = curve.point_add(&infinity, &p).unwrap();
            assert_eq!(result2, p);
            
            // Test neutral element rule: O ⊕ O = O
            let result3 = curve.point_add(&infinity, &infinity).unwrap();
            assert!(result3.is_infinity_point());
            
            // Test inverse elements rule: P ⊕ (-P) = O
            let p_inverse = curve.point_inverse(&p);
            let result4 = curve.point_add(&p, &p_inverse).unwrap();
            assert!(result4.is_infinity_point());
        }
    }

    #[test]
    fn test_scalar_multiplication() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if !finite_points.is_empty() {
            let p = finite_points[0].clone();
            
            // Test scalar multiplication with various values
            let zero_p = curve.scalar_multiply(0, &p).unwrap();
            assert!(zero_p.is_infinity_point());
            
            let one_p = curve.scalar_multiply(1, &p).unwrap();
            assert_eq!(one_p, p);
            
            let two_p = curve.scalar_multiply(2, &p).unwrap();
            let doubled = curve.point_double(&p).unwrap();
            assert_eq!(two_p, doubled);
            
            // Test negative scalar
            let neg_one_p = curve.scalar_multiply(-1, &p).unwrap();
            let p_inverse = curve.point_inverse(&p);
            assert_eq!(neg_one_p, p_inverse);
        }
    }

    #[test]
    fn test_point_inverse() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if !finite_points.is_empty() {
            let p = finite_points[0].clone();
            let p_inverse = curve.point_inverse(&p);
            
            // Test that inverse has same x-coordinate
            assert_eq!(p_inverse.x, p.x);
            
            // For Montgomery curves, if y = 0, then the inverse is the same point
            if p.y == 0 {
                assert_eq!(p_inverse.y, p.y);
            } else {
                // Otherwise, y-coordinates should be different
                assert_ne!(p_inverse.y, p.y);
            }
            
            // Test that adding a point to its inverse gives infinity
            let sum = curve.point_add(&p, &p_inverse).unwrap();
            assert!(sum.is_infinity_point());
            
            // Test that inverse of inverse is the original point
            let p_inverse_inverse = curve.point_inverse(&p_inverse);
            assert_eq!(p_inverse_inverse, p);
        }
    }

    #[test]
    fn test_montgomery_equation() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        let equation = curve.equation();
        assert_eq!(equation, "1y² = x³ + 1x² + x (mod 13)");
        
        let curve2 = MontgomeryCurve::new(3, 2, 17).unwrap();
        let equation2 = curve2.equation();
        assert_eq!(equation2, "2y² = x³ + 3x² + x (mod 17)");
    }

    #[test]
    fn test_montgomery_laws() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if finite_points.len() >= 2 {
            let p = finite_points[0].clone();
            let q = finite_points[1].clone();
            let infinity = Point::infinity();
            
            // Test Law 1: Neutral element - P ⊕ O = P
            let result1 = curve.point_add(&p, &infinity).unwrap();
            assert_eq!(result1, p, "Law 1 failed: P ⊕ O should equal P");
            
            let result2 = curve.point_add(&infinity, &p).unwrap();
            assert_eq!(result2, p, "Law 1 failed: O ⊕ P should equal P");
            
            // Test Law 2: Inverse elements - P ⊕ (-P) = O
            let p_inverse = curve.point_inverse(&p);
            let result3 = curve.point_add(&p, &p_inverse).unwrap();
            assert!(result3.is_infinity_point(), "Law 2 failed: P ⊕ (-P) should equal O");
            
            // Test Law 3: Tangent rule (doubling) - P ⊕ P = 2P
            if p.y != 0 {
                let doubled = curve.point_double(&p).unwrap();
                let added = curve.point_add(&p, &p).unwrap();
                assert_eq!(doubled, added, "Law 3 failed: P ⊕ P should equal 2P");
                
                // Verify the doubled point is on the curve
                assert!(curve.is_on_curve(&doubled), "Doubled point should be on the curve");
            }
            
            // Test Law 4: Chord rule (addition) - P ⊕ Q = R
            if p.x != q.x {
                let sum = curve.point_add(&p, &q).unwrap();
                assert!(curve.is_on_curve(&sum), "Sum point should be on the curve");
                
                // Test commutativity: P ⊕ Q = Q ⊕ P
                let sum2 = curve.point_add(&q, &p).unwrap();
                assert_eq!(sum, sum2, "Addition should be commutative");
            }
        }
    }

    #[test]
    fn test_montgomery_formulas() {
        let curve = MontgomeryCurve::new(1, 1, 13).unwrap();
        
        // Test with specific points that we know are on the curve
        let p = Point::new(0, 0); // This should be on the curve y² = x³ + x² + x (mod 13)
        let q = Point::new(3, 0); // This should also be on the curve
        
        if curve.is_on_curve(&p) && curve.is_on_curve(&q) {
            // Test doubling formula: x' = ((3x² + 2Ax + 1)/(2By))² - 2x - A
            let doubled = curve.point_double(&p).unwrap();
            assert!(curve.is_on_curve(&doubled), "Doubled point should be on the curve");
            
            // Test addition formula: x' = ((y2 - y1)/(x2 - x1))² * B - (x1 + x2) - A
            let sum = curve.point_add(&p, &q).unwrap();
            assert!(curve.is_on_curve(&sum), "Sum point should be on the curve");
            
            // Verify that the formulas produce correct results
            println!("P = ({}, {})", p.x, p.y);
            println!("Q = ({}, {})", q.x, q.y);
            println!("2P = ({}, {})", doubled.x, doubled.y);
            println!("P + Q = ({}, {})", sum.x, sum.y);
        }
    }
}