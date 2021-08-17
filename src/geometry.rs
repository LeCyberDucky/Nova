use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use rand::{prelude::ThreadRng, Rng};

use crate::image::{self, Color};

pub trait Obstacle {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit>;
}

pub struct ObstacleCollection {
    obstacles: Vec<Box<dyn Obstacle>>,
}

impl ObstacleCollection {
    pub fn new(obstacles: Vec<Box<dyn Obstacle>>) -> Self {
        Self { obstacles }
    }
    fn clear(&mut self) {
        self.obstacles.clear();
    }

    fn add(&mut self, obstacle: Box<dyn Obstacle>) {
        self.obstacles.push(obstacle);
    }
}

impl Obstacle for ObstacleCollection {
    fn hit(&self, ray: &Ray, t_min: f64, mut t_max: f64) -> Option<Hit> {
        let mut closest_hit = None;

        for obstacle in &self.obstacles {
            if let Some(hit) = obstacle.hit(ray, t_min, t_max) {
                t_max = hit.t;
                closest_hit = Some(hit);
            }
        }
        closest_hit
    }
}

pub struct Hit {
    p: Point3D,
    normal: Vec3D,
    t: f64,
    front_face: bool,
}

impl Hit {
    fn get_face_normal(ray: &Ray, outward_normal: Vec3D) -> (bool, Vec3D) {
        let front_face = (ray.direction * outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        (front_face, normal)
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Sphere {
    center: Point3D,
    radius: f64,
}

impl Sphere {
    pub fn new<T: Into<f64>>(center: Point3D, radius: T) -> Self {
        Self {
            center,
            radius: radius.into(),
        }
    }

    pub fn outward_normal(&self, surface_point: Point3D) -> Vec3D {
        (surface_point - self.center) / self.radius
    }
}

impl Obstacle for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
        let co = ray.origin - self.center;
        let a = ray.direction.magnitude_squared();
        let half_b = co * ray.direction;
        let c = co.magnitude_squared() - self.radius.powi(2);

        let discriminant = half_b.powi(2) - a * c;
        if discriminant < 0.0 {
            return None;
        };
        let discriminant_sqrt = discriminant.sqrt();

        // Find the nearest root that lies in the acceptable range
        let mut root = (-half_b - discriminant_sqrt) / a;
        if root < t_min || t_max < root {
            root = (discriminant_sqrt - half_b) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }

        let p = ray.at(root);
        let outward_normal = self.outward_normal(p);
        let (front_face, normal) = Hit::get_face_normal(ray, outward_normal);
        Some(Hit {
            p,
            normal,
            t: root,
            front_face,
        })
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Ray {
    origin: Point3D,
    direction: Vec3D,
}

impl Ray {
    pub fn new(origin: Point3D, direction: Vec3D) -> Self {
        Self { origin, direction }
    }

    pub fn at<T: Into<f64>>(&self, t: T) -> Point3D {
        self.origin + t.into() * self.direction
    }

    pub fn color(
        &self,
        obstacles: &ObstacleCollection,
        rng: &mut ThreadRng,
        depth: usize,
    ) -> image::Color {
        if depth == 0 {
            return Color::BLACK;
        }

        // If the ray hits an obstacle, let the obstacle absorb some light (reduce color), and let the ray bounce back in a random direction
        if let Some(hit) = obstacles.hit(self, f64::EPSILON, f64::INFINITY) {
            let target_direction = hit.normal + Vec3D::random_in_unit_sphere(rng);
            let ray = Ray::new(hit.p, target_direction);
            return 0.5 * ray.color(obstacles, rng, depth - 1);
        }

        let direction = self.direction.normalized();
        let t = 0.5 * (direction.y + 1.0);
        (1.0 - t) * Color::WHITE + t * Color::SKY_BLUE
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Vec3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

pub type Point3D = Vec3D;

impl Vec3D {
    pub fn new<A: Into<f64>, B: Into<f64>, C: Into<f64>>(x: A, y: B, z: C) -> Self {
        Self {
            x: x.into(),
            y: y.into(),
            z: z.into(),
        }
    }

    pub fn magnitude_squared(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }

    pub fn magnitude(&self) -> f64 {
        self.magnitude_squared().sqrt()
    }

    pub fn cross_product(&self, rhs: &Vec3D) -> Vec3D {
        let x = self.y * rhs.z - self.z * rhs.y;
        let y = self.z * rhs.x - self.x * rhs.z;
        let z = self.x * rhs.y - self.y * rhs.x;

        Vec3D::new(x, y, z)
    }

    pub fn normalize(&mut self) {
        *self /= self.magnitude()
    }

    pub fn normalized(&self) -> Self {
        *self / self.magnitude()
    }

    /// Creates a Vec3D with random components in the half-open range [0; 1[
    pub fn random(rng: &mut ThreadRng) -> Vec3D {
        Self::new(rng.gen::<f64>(), rng.gen::<f64>(), rng.gen::<f64>())
    }

    /// Creates a Vec3D with random components in the half-open range [min; max[
    pub fn random_in_range<A: Into<f64>, B: Into<f64>>(
        rng: &mut ThreadRng,
        min: A,
        max: B,
    ) -> Vec3D {
        let (min, max) = (min.into(), max.into());
        min + (max - min) * Self::random(rng)
    }

    pub fn random_in_unit_sphere(rng: &mut ThreadRng) -> Vec3D {
        let mut vector = Self::random_in_range(rng, -1.0, 1.0); // Random vector in unit cube
        let magnitude_squared = vector.magnitude_squared();
        if magnitude_squared >= 1.0 {
            // Vector needs scaling to fit inside unit sphere
            let max_scale = 1.0 / magnitude_squared.sqrt();
            let scale = rng.gen_range(0.0..max_scale) * 2.0 - max_scale; // Obtain scaling factor in range ]-max_scale, max_scale[. Directly using max_scale would lead to a weird distribution with a too many vectors directly at the sphere border
            vector *= scale;
        }
        vector
    }
}

impl<T: Into<f64>> Add<T> for Vec3D {
    type Output = Vec3D;

    fn add(self, rhs: T) -> Self::Output {
        let rhs = rhs.into();
        Self::new(self.x + rhs, self.y + rhs, self.z + rhs)
    }
}

impl Add<Vec3D> for Vec3D {
    type Output = Vec3D;

    fn add(self, rhs: Vec3D) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Add<Vec3D> for i32 {
    type Output = Vec3D;

    fn add(self, rhs: Vec3D) -> Self::Output {
        rhs + self
    }
}

impl Add<Vec3D> for f64 {
    type Output = Vec3D;

    fn add(self, rhs: Vec3D) -> Self::Output {
        rhs + self
    }
}

impl<T: Into<f64>> AddAssign<T> for Vec3D {
    fn add_assign(&mut self, rhs: T) {
        let rhs = rhs.into();
        self.x += rhs;
        self.y += rhs;
        self.z += rhs;
    }
}

impl AddAssign<Vec3D> for Vec3D {
    fn add_assign(&mut self, rhs: Vec3D) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T: Into<f64>> Sub<T> for Vec3D {
    type Output = Vec3D;

    fn sub(self, rhs: T) -> Self::Output {
        self.add(-rhs.into())
    }
}

impl Sub<Vec3D> for Vec3D {
    type Output = Vec3D;

    fn sub(self, rhs: Vec3D) -> Self::Output {
        self.add(-rhs)
    }
}

// Does subtracting a vector from a scalar actually make sense?
// Eh, this non-symmetric operator stuff has mushed my brain 🧠
impl Sub<Vec3D> for i32 {
    type Output = Vec3D;

    fn sub(self, rhs: Vec3D) -> Self::Output {
        -(rhs - self)
    }
}

impl Sub<Vec3D> for f64 {
    type Output = Vec3D;

    fn sub(self, rhs: Vec3D) -> Self::Output {
        -(rhs - self)
    }
}

impl<T: Into<f64>> SubAssign<T> for Vec3D {
    fn sub_assign(&mut self, rhs: T) {
        self.add_assign(-rhs.into())
    }
}

impl SubAssign<Vec3D> for Vec3D {
    fn sub_assign(&mut self, rhs: Vec3D) {
        self.add_assign(-rhs)
    }
}

impl<T: Into<f64>> Mul<T> for Vec3D {
    type Output = Vec3D;

    fn mul(self, rhs: T) -> Self::Output {
        let rhs = rhs.into();
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl<T: Into<f64>> MulAssign<T> for Vec3D {
    fn mul_assign(&mut self, rhs: T) {
        let rhs = rhs.into();
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl Mul<&Vec3D> for Vec3D {
    type Output = f64;

    fn mul(self, rhs: &Vec3D) -> Self::Output {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }
}

impl Mul<Vec3D> for Vec3D {
    type Output = f64;

    fn mul(self, rhs: Vec3D) -> Self::Output {
        self * &rhs
    }
}

impl Mul<Vec3D> for i32 {
    type Output = Vec3D;

    fn mul(self, rhs: Vec3D) -> Self::Output {
        rhs * self
    }
}

impl Mul<Vec3D> for f64 {
    type Output = Vec3D;

    fn mul(self, rhs: Vec3D) -> Self::Output {
        rhs * self
    }
}

impl<T: Into<f64>> Div<T> for Vec3D {
    type Output = Vec3D;

    fn div(self, rhs: T) -> Self::Output {
        let rhs = rhs.into();
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl<T: Into<f64>> DivAssign<T> for Vec3D {
    fn div_assign(&mut self, rhs: T) {
        let rhs = rhs.into();
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

impl Neg for &Vec3D {
    type Output = Vec3D;

    fn neg(self) -> Self::Output {
        Self::Output::new(-self.x, -self.y, -self.z)
    }
}

impl Neg for Vec3D {
    type Output = Vec3D;

    fn neg(self) -> Self::Output {
        -&self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let a = Vec3D::new(1, 2, 3);
        let b = Vec3D::new(4, 5, 6);
        let c = Vec3D::new(5, 7, 9);
        let d = Vec3D::new(2, 3, 4);
        assert_eq!(a + b, c);
        assert_eq!(a + 1, d);
        assert_eq!(1 + a, d);
        assert_eq!(1.0 + a, d);
    }

    #[test]
    fn test_sub() {
        let a = Vec3D::new(1, 2, 3);
        let b = Vec3D::new(4, 5, 6);
        let c = Vec3D::new(-3, -3, -3);
        let d = Vec3D::new(0, 1, 2);
        assert_eq!(a - b, c);
        assert_eq!(a - 1, d);
        assert_eq!(-(1 - a), d);
        assert_eq!(-(1.0 - a), d);
    }

    #[test]
    fn test_add_assign() {
        let mut a = Vec3D::new(1, 2, 3);
        let b = Vec3D::new(4, 5, 6);
        let c = Vec3D::new(5, 7, 9);
        let d = Vec3D::new(6, 8, 10);
        a += b;
        assert_eq!(a, c);
        a += 1;
        assert_eq!(a, d);
    }

    #[test]
    fn test_sub_assign() {
        let mut a = Vec3D::new(1, 2, 3);
        let b = Vec3D::new(4, 5, 6);
        let c = Vec3D::new(-3, -3, -3);
        let d = Vec3D::new(-4, -4, -4);
        a -= b;
        assert_eq!(a, c);
        a -= 1;
        assert_eq!(a, d);
    }

    #[test]
    fn test_mul() {
        let a = Vec3D::new(1, 2, 3);
        let b = Vec3D::new(-16, 5, 2);
        assert_eq!(a * b, 0.0);
        assert_eq!(a * (-2), -a * 2);
        assert_eq!(-2 * a, -a * 2);
    }

    #[test]
    fn test_mul_assign() {
        let a = Vec3D::new(1, 2, 3);
        let mut b = a;
        b *= -2;
        assert_eq!(b, -a * 2);
    }

    #[test]
    fn test_div() {
        let a = Vec3D::new(1, 2, 3);
        let b = Vec3D::new(0.5, 1, 1.5);
        assert_eq!(a / 2, b);
    }

    #[test]
    fn test_div_assign() {
        let mut a = Vec3D::new(1, 2, 3);
        let b = Vec3D::new(0.5, 1, 1.5);
        a /= 2;
        assert_eq!(a, b);
    }

    #[test]
    fn test_neg() {
        let a = -Vec3D::new(1, 2, 3);
        let b = Vec3D::new(-1, -2, -3);
        assert_eq!(a, b);
    }

    #[test]
    fn test_magnitude() {
        assert_eq!(Vec3D::new(1, 2, 2).magnitude(), 3.0);
        assert_eq!(Vec3D::new(-2, -3, -6).magnitude(), 7.0);
        assert_eq!(Vec3D::new(1, -4, 8).magnitude(), 9.0);
    }

    #[test]
    fn test_cross_product() {
        let a = Vec3D::new(1, 2, 3);
        let b = Vec3D::new(4, 5, 6);
        let c = Vec3D::new(-12, 3, 2);

        assert_eq!(a.cross_product(&b), Vec3D::new(-3, 6, -3));
        assert_eq!(b.cross_product(&a), Vec3D::new(3, -6, 3));
        assert_eq!(a * c, 0.0); // a and c need to be perpendicular for the next tests
        assert_eq!(a.cross_product(&c) * a, 0.0);
        assert_eq!(a.cross_product(&c) * c, 0.0);
        assert_eq!(c.cross_product(&a) * a, 0.0);
        assert_eq!(c.cross_product(&a) * c, 0.0);
    }

    #[test]
    fn test_normalize() {
        let mut a = Vec3D::new(1, -1, 1);
        a.normalize();
        assert_eq!(a.magnitude(), 1.0);
        assert_eq!(
            a,
            Vec3D::new(
                1.0 / 3.0_f64.sqrt(),
                -1.0 / 3.0_f64.sqrt(),
                1.0 / 3.0_f64.sqrt()
            )
        );
    }
}
