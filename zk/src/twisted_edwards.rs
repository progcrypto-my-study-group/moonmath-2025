use crate::traits::{CurvePoint, EllipticCurve};

/// Represents a point on a Twisted Edwards curve
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Point {
    pub x: i32,
    pub y: i32,
}

impl CurvePoint for Point {
    fn new(x: i32, y: i32) -> Self {
        Point { x, y }
    }

    fn infinity() -> Self {
        Point { x: 0, y: 1 } // Neutral element for Twisted Edwards curves
    }

    fn is_infinity_point(&self) -> bool {
        self.x == 0 && self.y == 1 // Check if it's the neutral element
    }
    
    fn x(&self) -> i32 {
        self.x
    }
    
    fn y(&self) -> i32 {
        self.y
    }
}

/// Represents a Twisted Edwards curve a·x² + y² = 1 + d·x²y² over a prime field
/// Constraints: a ≠ 0, d ≠ 0, a ≠ d, and p > 3
#[derive(Debug)]
pub struct TwistedEdwardsCurve {
    pub a: i32,
    pub d: i32,
    pub prime: i32,
}

impl TwistedEdwardsCurve {
    /// Creates a new Twisted Edwards curve with the given parameters
    /// Validates the constraints: a ≠ 0, d ≠ 0, a ≠ d, and p > 3
    pub fn new(a: i32, d: i32, prime: i32) -> Option<Self> {
        // Check constraints
        if prime <= 3 {
            return None; // p must be > 3
        }
        
        if a == 0 {
            return None; // a ≠ 0
        }
        
        if d == 0 {
            return None; // d ≠ 0
        }
        
        if a == d {
            return None; // a ≠ d
        }
        
        Some(TwistedEdwardsCurve { a, d, prime })
    }

    /// Checks if a point lies on the curve
    pub fn is_on_curve(&self, point: &Point) -> bool {
        // The point at infinity is always considered to be on the curve
        if point.is_infinity_point() {
            return true;
        }
        
        let x = point.x;
        let y = point.y;
        
        // Calculate a·x² mod prime
        let ax_squared = (self.a * x * x) % self.prime;
        let ax_squared = if ax_squared < 0 { ax_squared + self.prime } else { ax_squared };
        
        // Calculate y² mod prime
        let y_squared = (y * y) % self.prime;
        let y_squared = if y_squared < 0 { y_squared + self.prime } else { y_squared };
        
        // Calculate left side: a·x² + y²
        let left_side = (ax_squared + y_squared) % self.prime;
        let left_side = if left_side < 0 { left_side + self.prime } else { left_side };
        
        // Calculate d·x²y² mod prime
        let x_squared = (x * x) % self.prime;
        let x_squared = if x_squared < 0 { x_squared + self.prime } else { x_squared };
        let dx_squared_y_squared = (self.d * x_squared * y_squared) % self.prime;
        let dx_squared_y_squared = if dx_squared_y_squared < 0 { dx_squared_y_squared + self.prime } else { dx_squared_y_squared };
        
        // Calculate right side: 1 + d·x²y²
        let right_side = (1 + dx_squared_y_squared) % self.prime;
        let right_side = if right_side < 0 { right_side + self.prime } else { right_side };
        
        left_side == right_side
    }

    /// Point addition for Twisted Edwards curves
    /// For P = (x1, y1) and Q = (x2, y2):
    /// x3 = (x1*y2 + y1*x2)/(1 + d*x1*x2*y1*y2)
    /// y3 = (y1*y2 - a*x1*x2)/(1 - d*x1*x2*y1*y2)
    pub fn point_add(&self, p: &Point, q: &Point) -> Option<Point> {
        // Handle neutral element cases: (0,1) is the neutral element
        if p.is_infinity_point() {
            return Some(q.clone());
        }
        if q.is_infinity_point() {
            return Some(p.clone());
        }
        
        let x1 = p.x;
        let y1 = p.y;
        let x2 = q.x;
        let y2 = q.y;
        
        // Check if points are inverses (P + (-P) = (0,1))
        let p_inverse = self.point_inverse(p);
        if x2 == p_inverse.x && y2 == p_inverse.y {
            return Some(Point::infinity()); // Return (0,1)
        }
        
        // Check if points are the same (use doubling)
        if x1 == x2 && y1 == y2 {
            return self.point_double(p);
        }
        
        // Calculate x1*y2 + y1*x2
        let x1_y2 = (x1 * y2) % self.prime;
        let y1_x2 = (y1 * x2) % self.prime;
        let numerator_x = (x1_y2 + y1_x2) % self.prime;
        let numerator_x = if numerator_x < 0 { numerator_x + self.prime } else { numerator_x };
        
        // Calculate y1*y2 - a*x1*x2
        let y1_y2 = (y1 * y2) % self.prime;
        let ax1_x2 = (self.a * x1 * x2) % self.prime;
        let numerator_y = (y1_y2 - ax1_x2) % self.prime;
        let numerator_y = if numerator_y < 0 { numerator_y + self.prime } else { numerator_y };
        
        // Calculate d*x1*x2*y1*y2
        let dx1_x2_y1_y2 = (self.d * x1 * x2 * y1 * y2) % self.prime;
        let dx1_x2_y1_y2 = if dx1_x2_y1_y2 < 0 { dx1_x2_y1_y2 + self.prime } else { dx1_x2_y1_y2 };
        
        // Calculate denominators
        let denom_x = (1 + dx1_x2_y1_y2) % self.prime;
        let denom_x = if denom_x < 0 { denom_x + self.prime } else { denom_x };
        
        let denom_y = (1 - dx1_x2_y1_y2) % self.prime;
        let denom_y = if denom_y < 0 { denom_y + self.prime } else { denom_y };
        
        // Check for division by zero
        if denom_x == 0 || denom_y == 0 {
            return Some(Point::infinity()); // Return (0,1)
        }
        
        // Calculate x3 and y3
        let x3 = self.modular_divide(numerator_x, denom_x)?;
        let y3 = self.modular_divide(numerator_y, denom_y)?;
        
        Some(Point::new(x3, y3))
    }

    /// Point doubling for Twisted Edwards curves
    /// For P = (x, y):
    /// x2 = 2*x*y/(a*x² + y²)
    /// y2 = (y² - a*x²)/(2 - a*x² - y²)
    pub fn point_double(&self, point: &Point) -> Option<Point> {
        if point.is_infinity_point() {
            return Some(Point::infinity()); // Return (0,1)
        }
        
        let x = point.x;
        let y = point.y;
        
        // Calculate a*x² + y²
        let ax_squared = (self.a * x * x) % self.prime;
        let ax_squared = if ax_squared < 0 { ax_squared + self.prime } else { ax_squared };
        let y_squared = (y * y) % self.prime;
        let y_squared = if y_squared < 0 { y_squared + self.prime } else { y_squared };
        let ax_squared_plus_y_squared = (ax_squared + y_squared) % self.prime;
        let ax_squared_plus_y_squared = if ax_squared_plus_y_squared < 0 { ax_squared_plus_y_squared + self.prime } else { ax_squared_plus_y_squared };
        
        // Check for division by zero
        if ax_squared_plus_y_squared == 0 {
            return Some(Point::infinity()); // Return (0,1)
        }
        
        // Calculate 2*x*y
        let two_xy = (2 * x * y) % self.prime;
        let two_xy = if two_xy < 0 { two_xy + self.prime } else { two_xy };
        
        // Calculate x2 = 2*x*y/(a*x² + y²)
        let x2 = self.modular_divide(two_xy, ax_squared_plus_y_squared)?;
        
        // Calculate y² - a*x²
        let y_squared_minus_ax_squared = (y_squared - ax_squared) % self.prime;
        let y_squared_minus_ax_squared = if y_squared_minus_ax_squared < 0 { y_squared_minus_ax_squared + self.prime } else { y_squared_minus_ax_squared };
        
        // Calculate 2 - a*x² - y²
        let two_minus_ax_squared_minus_y_squared = (2 - ax_squared - y_squared) % self.prime;
        let two_minus_ax_squared_minus_y_squared = if two_minus_ax_squared_minus_y_squared < 0 { two_minus_ax_squared_minus_y_squared + self.prime } else { two_minus_ax_squared_minus_y_squared };
        
        // Check for division by zero
        if two_minus_ax_squared_minus_y_squared == 0 {
            return Some(Point::infinity()); // Return (0,1)
        }
        
        // Calculate y2 = (y² - a*x²)/(2 - a*x² - y²)
        let y2 = self.modular_divide(y_squared_minus_ax_squared, two_minus_ax_squared_minus_y_squared)?;
        
        Some(Point::new(x2, y2))
    }

    /// Returns the inverse of a point: -P = (-x, y)
    pub fn point_inverse(&self, point: &Point) -> Point {
        if point.is_infinity_point() {
            return Point::infinity(); // Return (0,1)
        }
        
        let inverse_x = (self.prime - point.x) % self.prime;
        Point::new(inverse_x, point.y)
    }

    /// Returns all possible points on the curve
    pub fn get_all_points(&self) -> Vec<Point> {
        let mut points = Vec::new();
        
        // Check all possible (x, y) pairs
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

    /// Modular division: a / b mod prime
    fn modular_divide(&self, a: i32, b: i32) -> Option<i32> {
        let b_inverse = self.modular_inverse(b)?;
        Some((a * b_inverse) % self.prime)
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

impl EllipticCurve for TwistedEdwardsCurve {
    type Point = Point;
    
    fn new(a: i32, b: i32, prime: i32) -> Self {
        // For Twisted Edwards curves, we use 'a' and 'd' parameters
        // The 'b' parameter from the trait is ignored, we use 'd' instead
        TwistedEdwardsCurve::new(a, b, prime).unwrap()
    }
    
    fn a(&self) -> i32 {
        self.a
    }
    
    fn b(&self) -> i32 {
        self.d // Return d as b for trait compatibility
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
        format!("{}x² + y² = 1 + {}x²y² (mod {})", self.a, self.d, self.prime)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_twisted_edwards_curve_creation() {
        // Valid curve
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        assert_eq!(curve.a, 1);
        assert_eq!(curve.d, 2);
        assert_eq!(curve.prime, 13);
        
        // Invalid curve: a = 0
        let invalid_curve = TwistedEdwardsCurve::new(0, 2, 13);
        assert!(invalid_curve.is_none());
        
        // Invalid curve: d = 0
        let invalid_curve2 = TwistedEdwardsCurve::new(1, 0, 13);
        assert!(invalid_curve2.is_none());
        
        // Invalid curve: a = d
        let invalid_curve3 = TwistedEdwardsCurve::new(1, 1, 13);
        assert!(invalid_curve3.is_none());
        
        // Invalid curve: p ≤ 3
        let invalid_curve4 = TwistedEdwardsCurve::new(1, 2, 3);
        assert!(invalid_curve4.is_none());
    }

    #[test]
    fn test_point_on_curve() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        
        // Test some points that should be on the curve x² + y² = 1 + 2x²y² (mod 13)
        let point1 = Point::new(0, 1); // 0² + 1² = 1, 1 + 2*0²*1² = 1 ✓
        let point2 = Point::new(0, 12); // 0² + 12² = 144 ≡ 1 (mod 13), 1 + 2*0²*12² = 1 ✓
        let point3 = Point::new(1, 0); // 1² + 0² = 1, 1 + 2*1²*0² = 1 ✓
        let point4 = Point::new(1, 1); // 1² + 1² = 2, 1 + 2*1²*1² = 3 ✗
        
        assert!(curve.is_on_curve(&point1));
        assert!(curve.is_on_curve(&point2));
        assert!(curve.is_on_curve(&point3));
        assert!(!curve.is_on_curve(&point4));
    }

    #[test]
    fn test_all_points() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        let points = curve.get_all_points();
        
        assert!(points.len() > 0);
        // Check that (0,1) is included as the neutral element
        assert!(points.iter().any(|p| p.x == 0 && p.y == 1));
        
        // Check that all result points are on the curve
        for point in &points {
            assert!(curve.is_on_curve(point), 
                   "Point {:?} is not on the curve", point);
        }
        
        // Check that there's exactly one neutral element (0,1)
        let neutral_count = points.iter().filter(|p| p.x == 0 && p.y == 1).count();
        assert_eq!(neutral_count, 1);
    }

    #[test]
    fn test_point_count() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        let count = curve.point_count();
        assert!(count > 0);
        println!("Twisted Edwards curve has {} points", count);
    }

    #[test]
    fn test_point_at_infinity() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        
        // Test neutral element creation
        let neutral = Point::infinity();
        assert!(neutral.is_infinity_point());
        assert!(neutral.x == 0 && neutral.y == 1);
        
        // Test that neutral element is on the curve
        assert!(curve.is_on_curve(&neutral));
        
        // Test that neutral element is included in all points
        let all_points = curve.get_all_points();
        assert!(all_points.iter().any(|p| p.x == 0 && p.y == 1));
        
        // Test that there's exactly one neutral element
        let neutral_count = all_points.iter().filter(|p| p.x == 0 && p.y == 1).count();
        assert_eq!(neutral_count, 1);
    }

    #[test]
    fn test_point_doubling() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        
        // Find a valid point on the curve
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if !finite_points.is_empty() {
            let p = finite_points[0].clone();
            let doubled = curve.point_double(&p).unwrap();
            
            // Verify the doubled point is on the curve
            assert!(curve.is_on_curve(&doubled));
        }
        
        // Test doubling the neutral element
        let neutral = Point::infinity();
        let doubled_neutral = curve.point_double(&neutral).unwrap();
        assert!(doubled_neutral.is_infinity_point());
    }

    #[test]
    fn test_point_addition() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        
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
        
        // Test adding a point to neutral element
        let neutral = Point::infinity();
        if !finite_points.is_empty() {
            let p = finite_points[0].clone();
            let result = curve.point_add(&p, &neutral).unwrap();
            assert_eq!(result, p);
        }
    }

    #[test]
    fn test_group_law_properties() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
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
            
            // Test identity: P + (0,1) = P
            let neutral = Point::infinity();
            let identity_test = curve.point_add(&p, &neutral).unwrap();
            assert_eq!(identity_test, p);
        }
    }

    #[test]
    fn test_infinity_rules() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if !finite_points.is_empty() {
            let p = finite_points[0].clone();
            let neutral = Point::infinity();
            
            // Test neutral element rule: P ⊕ (0,1) = P
            let result1 = curve.point_add(&p, &neutral).unwrap();
            assert_eq!(result1, p);
            
            // Test neutral element rule: (0,1) ⊕ P = P
            let result2 = curve.point_add(&neutral, &p).unwrap();
            assert_eq!(result2, p);
            
            // Test neutral element rule: (0,1) ⊕ (0,1) = (0,1)
            let result3 = curve.point_add(&neutral, &neutral).unwrap();
            assert!(result3.is_infinity_point());
            
            // Test inverse elements rule: P ⊕ (-P) = (0,1)
            let p_inverse = curve.point_inverse(&p);
            let result4 = curve.point_add(&p, &p_inverse).unwrap();
            assert!(result4.is_infinity_point());
        }
    }

    #[test]
    fn test_scalar_multiplication() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
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
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        let all_points = curve.get_all_points();
        let finite_points: Vec<_> = all_points.iter().filter(|p| !p.is_infinity_point()).collect();
        
        if !finite_points.is_empty() {
            let p = finite_points[0].clone();
            let p_inverse = curve.point_inverse(&p);
            
            // Test that inverse has different x-coordinate but same y-coordinate
            // For Twisted Edwards curves, if x = 0, then the inverse is the same point
            if p.x == 0 {
                assert_eq!(p_inverse.x, p.x);
            } else {
                assert_ne!(p_inverse.x, p.x);
            }
            assert_eq!(p_inverse.y, p.y);
            
            // Test that adding a point to its inverse gives neutral element
            let sum = curve.point_add(&p, &p_inverse).unwrap();
            assert!(sum.is_infinity_point());
            
            // Test that inverse of inverse is the original point
            let p_inverse_inverse = curve.point_inverse(&p_inverse);
            assert_eq!(p_inverse_inverse, p);
        }
    }

    #[test]
    fn test_twisted_edwards_equation() {
        let curve = TwistedEdwardsCurve::new(1, 2, 13).unwrap();
        let equation = curve.equation();
        assert_eq!(equation, "1x² + y² = 1 + 2x²y² (mod 13)");
        
        let curve2 = TwistedEdwardsCurve::new(3, 4, 17).unwrap();
        let equation2 = curve2.equation();
        assert_eq!(equation2, "3x² + y² = 1 + 4x²y² (mod 17)");
    }
}