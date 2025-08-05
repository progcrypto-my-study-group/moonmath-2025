use crate::traits::{CurvePoint, EllipticCurve};

/// Represents a point on a Weierstrass curve
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

/// Represents a short Weierstrass curve y² = x³ + ax + b over a prime field
#[derive(Debug)]
pub struct WeierstrassCurve {
    pub a: i32,
    pub b: i32,
    pub prime: i32,
}

impl WeierstrassCurve {
    /// Creates a new Weierstrass curve with the given parameters
    pub fn new(a: i32, b: i32, prime: i32) -> Self {
        WeierstrassCurve { a, b, prime }
    }

    /// Point doubling: P ⊕ P = 2P (Tangent Rule)
    /// For P = (x, y) with y ≠ 0:
    /// x' = (3x² + a)/(2y)² - 2x
    /// y' = (3x² + a)/(2y) * (x - x') - y
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
        
        // Calculate 3x² + a
        let three_x_squared = (3 * x * x) % self.prime;
        let numerator = (three_x_squared + self.a) % self.prime;
        
        // Calculate 2y
        let two_y = (2 * y) % self.prime;
        
        // Calculate λ = (3x² + a)/(2y)
        let lambda = self.modular_divide(numerator, two_y)?;
        
        // Calculate x' = λ² - 2x
        let lambda_squared = (lambda * lambda) % self.prime;
        let x_prime = (lambda_squared - 2 * x) % self.prime;
        
        // Handle negative values
        let x_prime = if x_prime < 0 { x_prime + self.prime } else { x_prime };
        
        // Calculate y' = λ(x - x') - y
        let x_diff = (x - x_prime) % self.prime;
        let x_diff = if x_diff < 0 { x_diff + self.prime } else { x_diff };
        let y_prime = (lambda * x_diff - y) % self.prime;
        let y_prime = if y_prime < 0 { y_prime + self.prime } else { y_prime };
        
        Some(Point::new(x_prime, y_prime))
    }

    /// Point addition: P ⊕ Q = R (Chord Rule)
    /// For P = (x1, y1) and Q = (x2, y2) with x1 ≠ x2:
    /// x3 = ((y2 - y1)/(x2 - x1))² - x1 - x2
    /// y3 = ((y2 - y1)/(x2 - x1)) * (x1 - x3) - y1
    /// 
    /// Special cases:
    /// (Neutral element) If Q = O, then P ⊕ Q = P
    /// (Inverse elements) If P = (x, y) and Q = (x, -y), then P ⊕ Q = O
    pub fn point_add(&self, p: &Point, q: &Point) -> Option<Point> {
        // Handle point at infinity cases (Neutral element rule)
        if p.is_infinity {
            return Some(q.clone()); // O ⊕ Q = Q
        }
        if q.is_infinity {
            return Some(p.clone()); // P ⊕ O = P
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
                // Inverse elements rule: P = (x, y) and Q = (x, -y) → P ⊕ Q = O
                return Some(Point::infinity());
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
        
        // Calculate x3 = λ² - x1 - x2
        let lambda_squared = (lambda * lambda) % self.prime;
        let x3 = (lambda_squared - x1 - x2) % self.prime;
        let x3 = if x3 < 0 { x3 + self.prime } else { x3 };
        
        // Calculate y3 = λ(x1 - x3) - y1
        let x1_minus_x3 = (x1 - x3) % self.prime;
        let x1_minus_x3 = if x1_minus_x3 < 0 { x1_minus_x3 + self.prime } else { x1_minus_x3 };
        let y3 = (lambda * x1_minus_x3 - y1) % self.prime;
        let y3 = if y3 < 0 { y3 + self.prime } else { y3 };
        
        Some(Point::new(x3, y3))
    }

    /// Scalar multiplication: k·P using double-and-add algorithm
    /// Computes the sum of P added to itself k times
    pub fn scalar_multiply(&self, k: i32, point: &Point) -> Option<Point> {
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

    /// Returns the inverse of a point: -P = (x, -y)
    pub fn point_inverse(&self, point: &Point) -> Point {
        if point.is_infinity {
            return Point::infinity(); // O is its own inverse
        }
        
        let inverse_y = (self.prime - point.y) % self.prime;
        Point::new(point.x, inverse_y)
    }

    /// Checks if the curve parameters satisfy the discriminant condition 4a³ + 27b² ≠ 0
    pub fn is_nonsingular(&self) -> bool {
        let a_cubed = (self.a * self.a * self.a) % self.prime;
        let b_squared = (self.b * self.b) % self.prime;
        let term1 = (4 * a_cubed) % self.prime;
        let term2 = (27 * b_squared) % self.prime;
        let discriminant = (term1 + term2) % self.prime;
        
        // Handle negative values
        let discriminant = if discriminant < 0 { discriminant + self.prime } else { discriminant };
        
        discriminant != 0
    }

    /// Checks if a point lies on the curve
    pub fn is_on_curve(&self, point: &Point) -> bool {
        // The point at infinity is always considered to be on the curve
        if point.is_infinity {
            return true;
        }
        
        let x = point.x;
        let y = point.y;
        
        // Calculate y² mod prime
        let mut y_squared = (y * y) % self.prime;
        if y_squared < 0 {
            y_squared += self.prime;
        }
        
        // Calculate x³ + ax + b mod prime
        let x_cubed = (x * x * x) % self.prime;
        let ax = (self.a * x) % self.prime;
        let mut right_side = (x_cubed + ax + self.b) % self.prime;
        
        // Handle negative values
        if right_side < 0 {
            right_side += self.prime;
        }
        
        y_squared == right_side
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

    /// Isomorphism function I = (c²·x, c³·y) that maps points from this curve to another
    /// 
    /// This function implements the isomorphism transformation that maps points
    /// from the current curve to a target curve. For the transformation I = (c²·x, c³·y),
    /// if the original curve is y² = x³ + ax + b, then points mapped via this
    /// transformation will lie on the curve y'² = x'³ + a'x' + b' where:
    /// - a' = a/c⁴ (mod p) 
    /// - b' = b/c⁶ (mod p)
    /// 
    /// Example: For tinyjubjub13 (a=8, b=8) with c=6:
    /// - c² = 36 ≡ 10 (mod 13)
    /// - c³ = 216 ≡ 8 (mod 13) 
    /// - Maps to curve with a' = 8/6⁴ ≡ 8/10 ≡ 7 (mod 13), b' = 8/6⁶ ≡ 8/8 ≡ 5 (mod 13)
    pub fn isomorphism_map(&self, point: &Point, c: i32) -> Option<Point> {
        if point.is_infinity {
            return Some(Point::infinity());
        }
        
        // Calculate c² mod p
        let c_squared = (c * c) % self.prime;
        let c_squared = if c_squared < 0 { c_squared + self.prime } else { c_squared };
        
        // Calculate c³ mod p
        let c_cubed = (c * c * c) % self.prime;
        let c_cubed = if c_cubed < 0 { c_cubed + self.prime } else { c_cubed };
        
        // Apply transformation: (x', y') = (c²·x, c³·y)
        let x_prime = (c_squared * point.x) % self.prime;
        let x_prime = if x_prime < 0 { x_prime + self.prime } else { x_prime };
        
        let y_prime = (c_cubed * point.y) % self.prime;
        let y_prime = if y_prime < 0 { y_prime + self.prime } else { y_prime };
        
        Some(Point::new(x_prime, y_prime))
    }

    /// Inverse isomorphism function I⁻¹ = (x/c², y/c³) that maps points back
    pub fn isomorphism_inverse_map(&self, point: &Point, c: i32) -> Option<Point> {
        if point.is_infinity {
            return Some(Point::infinity());
        }
        
        // Calculate c⁻² = (c²)⁻¹ mod p
        let c_squared = (c * c) % self.prime;
        let c_squared = if c_squared < 0 { c_squared + self.prime } else { c_squared };
        let c_squared_inv = self.modular_inverse(c_squared)?;
        
        // Calculate c⁻³ = (c³)⁻¹ mod p
        let c_cubed = (c * c * c) % self.prime;
        let c_cubed = if c_cubed < 0 { c_cubed + self.prime } else { c_cubed };
        let c_cubed_inv = self.modular_inverse(c_cubed)?;
        
        // Apply inverse transformation: (x', y') = (x/c², y/c³)
        let x_prime = (point.x * c_squared_inv) % self.prime;
        let x_prime = if x_prime < 0 { x_prime + self.prime } else { x_prime };
        
        let y_prime = (point.y * c_cubed_inv) % self.prime;
        let y_prime = if y_prime < 0 { y_prime + self.prime } else { y_prime };
        
        Some(Point::new(x_prime, y_prime))
    }

    /// Creates an isomorphic curve with the given transformation parameter c
    /// The new curve will have parameters a' = a·c⁴ and b' = b·c⁶
    pub fn create_isomorphic_curve(&self, c: i32) -> Option<WeierstrassCurve> {
        // Calculate c⁴ = c² * c²
        let c_squared = (c * c) % self.prime;
        let c_squared = if c_squared < 0 { c_squared + self.prime } else { c_squared };
        let c_fourth = (c_squared * c_squared) % self.prime;
        let c_fourth = if c_fourth < 0 { c_fourth + self.prime } else { c_fourth };
        
        // Calculate c⁶ = c² * c⁴
        let c_sixth = (c_squared * c_fourth) % self.prime;
        let c_sixth = if c_sixth < 0 { c_sixth + self.prime } else { c_sixth };
        
        // Calculate a' = a·c⁴
        let a_prime = (self.a * c_fourth) % self.prime;
        let a_prime = if a_prime < 0 { a_prime + self.prime } else { a_prime };
        
        // Calculate b' = b·c⁶
        let b_prime = (self.b * c_sixth) % self.prime;
        let b_prime = if b_prime < 0 { b_prime + self.prime } else { b_prime };
        
        Some(WeierstrassCurve::new(a_prime, b_prime, self.prime))
    }

    /// Finds the transformation parameter c such that the target curve maps to this curve
    /// via the inverse transformation (x, y) = (x'/c², y'/c³)
    /// Returns None if no such transformation exists
    pub fn find_isomorphism_parameter(&self, target_a: i32, target_b: i32) -> Option<i32> {
        // We need to find c such that applying the inverse transformation
        // (x, y) = (x'/c², y'/c³) to the target curve y'² = x'³ + target_a·x' + target_b
        // gives us this curve y² = x³ + self.a·x + self.b
        //
        // This means: self.a = target_a / c⁴ (mod p)
        //            self.b = target_b / c⁶ (mod p)
        //
        // So: c⁴ = target_a / self.a (mod p)
        //     c⁶ = target_b / self.b (mod p)
        
        for c in 1..self.prime {
            // Calculate c⁴ and c⁶
            let c_squared = (c * c) % self.prime;
            let c_fourth = (c_squared * c_squared) % self.prime;
            let c_sixth = (c_squared * c_fourth) % self.prime;
            
            // Check if target_a / c⁴ = self.a
            if let Some(c4_inv) = self.modular_inverse(c_fourth) {
                let computed_a = (target_a * c4_inv) % self.prime;
                let computed_a = if computed_a < 0 { computed_a + self.prime } else { computed_a };
                
                // Check if target_b / c⁶ = self.b  
                if let Some(c6_inv) = self.modular_inverse(c_sixth) {
                    let computed_b = (target_b * c6_inv) % self.prime;
                    let computed_b = if computed_b < 0 { computed_b + self.prime } else { computed_b };
                    
                    if computed_a == self.a && computed_b == self.b {
                        return Some(c);
                    }
                }
            }
        }
        None
    }
}

impl EllipticCurve for WeierstrassCurve {
    type Point = Point;
    
    fn new(a: i32, b: i32, prime: i32) -> Self {
        WeierstrassCurve { a, b, prime }
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
        self.scalar_multiply(k, point)
    }
    
    fn get_all_points(&self) -> Vec<Self::Point> {
        self.get_all_points()
    }
    

    
    fn equation(&self) -> String {
        format!("y² = x³ + {}x + {} (mod {})", self.a, self.b, self.prime)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weierstrass_curve_creation() {
        let curve = WeierstrassCurve::new(1, 1, 5);
        assert_eq!(curve.a, 1);
        assert_eq!(curve.b, 1);
        assert_eq!(curve.prime, 5);
    }

    #[test]
    fn test_point_on_curve() {
        let curve = WeierstrassCurve::new(1, 1, 5);
        
        // Test some points that should be on the curve y² = x³ + x + 1 (mod 5)
        let point1 = Point::new(0, 1); // 1² = 0³ + 0 + 1 = 1 ✓
        let point2 = Point::new(0, 4); // 4² = 16 ≡ 1 (mod 5), 0³ + 0 + 1 = 1 ✓
        let point3 = Point::new(1, 0); // 0² = 0, 1³ + 1 + 1 = 3 ✗
        let point4 = Point::new(2, 0); // 0² = 0, 2³ + 2 + 1 = 8 + 2 + 1 = 11 ≡ 1 (mod 5) ✗
        
        assert!(curve.is_on_curve(&point1));
        assert!(curve.is_on_curve(&point2));
        assert!(!curve.is_on_curve(&point3));
        assert!(!curve.is_on_curve(&point4));
    }

    #[test]
    fn test_all_points() {
        let curve = WeierstrassCurve::new(1, 1, 5);
        let points = curve.get_all_points();
        
        // For the curve y² = x³ + x + 1 (mod 5), the correct points are:
        // x=0: y² = 1, so y = 1, 4
        // x=1: y² = 3, but 3 is not a quadratic residue mod 5, so no points
        // x=2: y² = 1, so y = 1, 4  
        // x=3: y² = 1, so y = 1, 4
        // x=4: y² = 4, so y = 2, 3
        // Plus the point at infinity
        let expected_points = vec![
            Point::infinity(),                    // Point at infinity
            Point::new(0, 1), Point::new(0, 4),  // x=0: y² = 1, so y = 1, 4
            Point::new(2, 1), Point::new(2, 4),  // x=2: y² = 1, so y = 1, 4
            Point::new(3, 1), Point::new(3, 4),  // x=3: y² = 1, so y = 1, 4
            Point::new(4, 2), Point::new(4, 3),  // x=4: y² = 4, so y = 2, 3
        ];
        
        assert_eq!(points.len(), expected_points.len());
        
        // Check that all expected points are in the result
        for expected_point in &expected_points {
            assert!(points.contains(expected_point), 
                   "Expected point {:?} not found in result", expected_point);
        }
        
        // Check that all result points are on the curve
        for point in &points {
            assert!(curve.is_on_curve(point), 
                   "Point {:?} is not on the curve", point);
        }
    }

    #[test]
    fn test_point_count() {
        let curve = WeierstrassCurve::new(1, 1, 5);
        assert_eq!(curve.point_count(), 9); // 8 finite points + 1 point at infinity
    }

    #[test]
    fn test_point_at_infinity() {
        let curve = WeierstrassCurve::new(1, 1, 5);
        
        // Test point at infinity creation
        let infinity = Point::infinity();
        assert!(infinity.is_infinity_point());
        assert!(infinity.is_infinity); // The field should be true
        
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
        let curve = WeierstrassCurve::new(1, 1, 5);
        
        // Test doubling a point
        let p = Point::new(0, 1);
        let doubled = curve.point_double(&p).unwrap();
        
        // Verify the doubled point is on the curve
        assert!(curve.is_on_curve(&doubled));
        
        // Test doubling the point at infinity
        let infinity = Point::infinity();
        let doubled_infinity = curve.point_double(&infinity).unwrap();
        assert!(doubled_infinity.is_infinity_point());
        
        // Test doubling a point with y = 0 (should return infinity)
        let p_y_zero = Point::new(3, 0); // This point has y = 0
        let doubled_y_zero = curve.point_double(&p_y_zero).unwrap();
        assert!(doubled_y_zero.is_infinity_point());
    }

    #[test]
    fn test_point_addition() {
        let curve = WeierstrassCurve::new(1, 1, 5);
        
        // Test adding two different points
        let p = Point::new(0, 1);
        let q = Point::new(2, 1);
        let sum = curve.point_add(&p, &q).unwrap();
        
        // Verify the sum is on the curve
        assert!(curve.is_on_curve(&sum));
        
        // Test adding a point to itself (should use doubling)
        let doubled = curve.point_double(&p).unwrap();
        let added = curve.point_add(&p, &p).unwrap();
        assert_eq!(doubled, added);
        
        // Test adding a point to infinity
        let infinity = Point::infinity();
        let result = curve.point_add(&p, &infinity).unwrap();
        assert_eq!(result, p);
        
        // Test adding a point to its inverse (should return infinity)
        let p_inverse = Point::new(0, 4); // (0,4) is the inverse of (0,1) since 4 ≡ -1 (mod 5)
        let sum_inverse = curve.point_add(&p, &p_inverse).unwrap();
        assert!(sum_inverse.is_infinity_point());
    }

    #[test]
    fn test_group_law_properties() {
        let curve = WeierstrassCurve::new(1, 1, 5);
        
        // Test associativity: (P + Q) + R = P + (Q + R)
        let p = Point::new(0, 1);
        let q = Point::new(2, 1);
        let r = Point::new(4, 2);
        
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

    #[test]
    fn test_infinity_rules() {
        let curve = WeierstrassCurve::new(1, 1, 5);
        let p = Point::new(0, 1);
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
        let p_inverse = Point::new(0, 4); // (0,4) is the inverse of (0,1) since 4 ≡ -1 (mod 5)
        let result4 = curve.point_add(&p, &p_inverse).unwrap();
        assert!(result4.is_infinity_point());
        
        // Test another inverse pair
        let q = Point::new(2, 1);
        let q_inverse = Point::new(2, 4); // (2,4) is the inverse of (2,1) since 4 ≡ -1 (mod 5)
        let result5 = curve.point_add(&q, &q_inverse).unwrap();
        assert!(result5.is_infinity_point());
        
        // Verify that the inverse points are actually inverses
        assert_eq!(p_inverse.y, (curve.prime - p.y) % curve.prime);
        assert_eq!(q_inverse.y, (curve.prime - q.y) % curve.prime);
    }

    #[test]
    fn test_scalar_multiplication() {
        let curve = WeierstrassCurve::new(1, 1, 5);
        let p = Point::new(0, 1);
        
        // Test scalar multiplication with various values
        let zero_p = curve.scalar_multiply(0, &p).unwrap();
        assert!(zero_p.is_infinity_point());
        
        let one_p = curve.scalar_multiply(1, &p).unwrap();
        assert_eq!(one_p, p);
        
        let two_p = curve.scalar_multiply(2, &p).unwrap();
        let doubled = curve.point_double(&p).unwrap();
        assert_eq!(two_p, doubled);
        
        let three_p = curve.scalar_multiply(3, &p).unwrap();
        let expected_three_p = curve.point_add(&p, &curve.point_add(&p, &p).unwrap()).unwrap();
        assert_eq!(three_p, expected_three_p);
        
        // Test negative scalar
        let neg_one_p = curve.scalar_multiply(-1, &p).unwrap();
        let p_inverse = curve.point_inverse(&p);
        assert_eq!(neg_one_p, p_inverse);
        
        // Test that scalar multiplication is consistent
        let five_p = curve.scalar_multiply(5, &p).unwrap();
        let manual_five_p = curve.point_add(&p, &curve.point_add(&p, &curve.point_add(&p, &curve.point_add(&p, &p).unwrap()).unwrap()).unwrap()).unwrap();
        assert_eq!(five_p, manual_five_p);
    }

    #[test]
    fn test_isomorphic_curves() {
        // Original Tiny-jubjub curve: y² = x³ + 8x + 8 (mod 13)
        let tiny_jubjub = WeierstrassCurve::new(8, 8, 13);
        
        // Create isomorphic curve with a=7, b=5
        // We need to find u such that: 7 = 8/u² and 5 = 8/u³
        // Let's calculate this directly
        let isomorphic_curve = WeierstrassCurve::new(7, 5, 13);
        
        // Verify both curves are nonsingular
        assert!(tiny_jubjub.is_nonsingular());
        assert!(isomorphic_curve.is_nonsingular());
        
        // Find a valid point on the original curve
        let original_points = tiny_jubjub.get_all_points();
        let finite_points: Vec<_> = original_points.iter().filter(|p| !p.is_infinity_point()).collect();
        let _original_point = finite_points[1].clone(); // Use the second finite point
        
        // For this example, we'll demonstrate the concept by showing both curves
        // have the same number of points (they should be isomorphic)
        let isomorphic_points = isomorphic_curve.get_all_points();
        assert_eq!(original_points.len(), isomorphic_points.len());
    }

    #[test]
    fn test_discriminant_condition() {
        // Test that the discriminant condition works correctly
        
        // Tiny-jubjub curve (should be nonsingular)
        let tiny_jubjub = WeierstrassCurve::new(8, 8, 13);
        assert!(tiny_jubjub.is_nonsingular(), "Tiny-jubjub should be nonsingular");
        
        // Original test curve (should be nonsingular)
        let original_curve = WeierstrassCurve::new(1, 1, 5);
        assert!(original_curve.is_nonsingular(), "Original curve should be nonsingular");
        
        // Test a curve that would be singular: a=0, b=0 (y² = x³)
        // This would have discriminant = 0, making it singular
        let singular_curve = WeierstrassCurve::new(0, 0, 17);
        assert!(!singular_curve.is_nonsingular(), "Curve with a=0, b=0 should be singular");
        
        // Test another singular case: a=0, b=0 over different prime
        let singular_curve2 = WeierstrassCurve::new(0, 0, 7);
        assert!(!singular_curve2.is_nonsingular(), "Curve with a=0, b=0 should be singular");
    }

    #[test]
    fn test_isomorphism_tinyjubjub_to_7_5() {
        // Original Tiny-jubjub curve: y² = x³ + 8x + 8 (mod 13)
        let tiny_jubjub = WeierstrassCurve::new(8, 8, 13);
        
        // Target curve: y² = x³ + 7x + 5 (mod 13)
        let target_curve = WeierstrassCurve::new(7, 5, 13);
        
        // Find the transformation parameter c
        let c = tiny_jubjub.find_isomorphism_parameter(7, 5);
        assert!(c.is_some(), "Should find a transformation parameter c");
        let c = c.unwrap();
        
        println!("Found transformation parameter c = {}", c);
        
        // Verify that the isomorphic curve created with this c matches the target
        let isomorphic_curve = tiny_jubjub.create_isomorphic_curve(c).unwrap();
        assert_eq!(isomorphic_curve.a, 7, "Transformed curve should have a = 7");
        assert_eq!(isomorphic_curve.b, 5, "Transformed curve should have b = 5");
        
        // Get all points on the original curve
        let original_points = tiny_jubjub.get_all_points();
        let target_points = target_curve.get_all_points();
        
        // Verify same number of points (isomorphic property)
        assert_eq!(original_points.len(), target_points.len(), 
                  "Isomorphic curves should have the same number of points");
        
        // Test the isomorphism mapping on several points
        for original_point in &original_points {
            // Map the point using I = (c²·x, c³·y)
            let mapped_point = tiny_jubjub.isomorphism_map(original_point, c).unwrap();
            
            // Verify the mapped point is on the target curve
            assert!(target_curve.is_on_curve(&mapped_point), 
                   "Mapped point {:?} should be on target curve", mapped_point);
            
            // Test the inverse mapping
            let inverse_mapped = target_curve.isomorphism_inverse_map(&mapped_point, c).unwrap();
            assert_eq!(inverse_mapped, *original_point, 
                      "Inverse mapping should return original point");
        }
        
        // Test specific point mappings to demonstrate the isomorphism
        let finite_points: Vec<_> = original_points.iter()
            .filter(|p| !p.is_infinity_point())
            .take(3)
            .collect();
            
        if finite_points.len() >= 2 {
            let p1 = finite_points[0];
            let p2 = finite_points[1];
            
            // Map individual points
            let mapped_p1 = tiny_jubjub.isomorphism_map(p1, c).unwrap();
            let mapped_p2 = tiny_jubjub.isomorphism_map(p2, c).unwrap();
            
            // Add points on original curve
            if let Some(sum_original) = tiny_jubjub.point_add(p1, p2) {
                // Map the sum
                let mapped_sum = tiny_jubjub.isomorphism_map(&sum_original, c).unwrap();
                
                // Add mapped points on target curve
                if let Some(sum_mapped) = target_curve.point_add(&mapped_p1, &mapped_p2) {
                    // The isomorphism should preserve the group operation
                    assert_eq!(mapped_sum, sum_mapped, 
                              "Isomorphism should preserve group operations: I(P + Q) = I(P) + I(Q)");
                    
                    println!("Verified group operation preservation:");
                    println!("  Original: ({},{}) + ({},{}) = ({},{})", 
                            p1.x(), p1.y(), p2.x(), p2.y(), sum_original.x(), sum_original.y());
                    println!("  Mapped: ({},{}) + ({},{}) = ({},{})", 
                            mapped_p1.x(), mapped_p1.y(), mapped_p2.x(), mapped_p2.y(), sum_mapped.x(), sum_mapped.y());
                }
            }
        }
        
        println!("Successfully demonstrated isomorphism I = ({}²·x, {}³·y)", c, c);
        println!("Maps tinyjubjub13 (a=8, b=8) to target curve (a=7, b=5) over F₁₃");
    }

    #[test]
    fn test_isomorphism_properties() {
        // Test basic properties of the isomorphism functions
        let curve = WeierstrassCurve::new(8, 8, 13);
        let c = 2; // Use c = 2 as an example
        
        // Get a test point
        let points = curve.get_all_points();
        let test_point = &points[1]; // Use a non-infinity point
        
        // Test that mapping and inverse mapping are inverses
        let mapped = curve.isomorphism_map(test_point, c).unwrap();
        let inverse_mapped = curve.isomorphism_inverse_map(&mapped, c).unwrap();
        assert_eq!(inverse_mapped, *test_point, "I⁻¹(I(P)) should equal P");
        
        // Test that infinity maps to infinity
        let infinity = Point::infinity();
        let mapped_infinity = curve.isomorphism_map(&infinity, c).unwrap();
        assert!(mapped_infinity.is_infinity_point(), "Infinity should map to infinity");
        
        // Test the transformation formulas directly
        if !test_point.is_infinity_point() {
            let c_squared = (c * c) % 13;
            let c_cubed = (c * c * c) % 13;
            
            let expected_x = (c_squared * test_point.x()) % 13;
            let expected_y = (c_cubed * test_point.y()) % 13;
            
            assert_eq!(mapped.x(), expected_x, "x-coordinate should be c²·x");
            assert_eq!(mapped.y(), expected_y, "y-coordinate should be c³·y");
        }
    }
}