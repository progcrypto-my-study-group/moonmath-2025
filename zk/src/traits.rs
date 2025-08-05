use std::fmt::Debug;

/// Common trait for points on elliptic curves
pub trait CurvePoint: Clone + Debug + PartialEq + Eq {
    /// Creates a new point with the given coordinates
    fn new(x: i32, y: i32) -> Self;
    
    /// Creates the point at infinity (neutral element)
    fn infinity() -> Self;
    
    /// Checks if this is the point at infinity
    fn is_infinity_point(&self) -> bool;
    
    /// Returns the x-coordinate of the point
    fn x(&self) -> i32;
    
    /// Returns the y-coordinate of the point
    fn y(&self) -> i32;
}

/// Common trait for elliptic curves
pub trait EllipticCurve: Debug {
    type Point: CurvePoint;
    
    /// Creates a new curve with the given parameters
    fn new(a: i32, b: i32, prime: i32) -> Self;
    
    /// Returns the curve parameters
    fn a(&self) -> i32;
    fn b(&self) -> i32;
    fn prime(&self) -> i32;
    
    /// Checks if a point lies on the curve
    fn is_on_curve(&self, point: &Self::Point) -> bool;
    
    /// Point doubling: P ⊕ P = 2P
    fn point_double(&self, point: &Self::Point) -> Option<Self::Point>;
    
    /// Point addition: P ⊕ Q = R
    fn point_add(&self, p: &Self::Point, q: &Self::Point) -> Option<Self::Point>;
    
    /// Returns the inverse of a point: -P
    fn point_inverse(&self, point: &Self::Point) -> Self::Point;
    
    /// Scalar multiplication: k·P
    fn scalar_multiply(&self, k: i32, point: &Self::Point) -> Option<Self::Point>;
    
    /// Returns all possible points on the curve
    fn get_all_points(&self) -> Vec<Self::Point>;
    
    /// Returns the number of points on the curve
    fn point_count(&self) -> usize {
        self.get_all_points().len()
    }
    

    
    /// Returns the curve equation as a string
    fn equation(&self) -> String;
}