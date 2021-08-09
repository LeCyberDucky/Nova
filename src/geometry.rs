use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use num::Float;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct Vec3D<T: Float + Copy> {
    x: T,
    y: T,
    z: T,
}

impl<T: Float + Copy> Vec3D<T> {
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }

    fn x(&self) -> &T {
        &self.x
    }

    fn y(&self) -> &T {
        &self.y
    }

    fn z(&self) -> &T {
        &self.z
    }

    fn x_mut(&mut self) -> &mut T {
        &mut self.x
    }

    fn y_mut(&mut self) -> &mut T {
        &mut self.y
    }

    fn z_mut(&mut self) -> &mut T {
        &mut self.z
    }

    pub fn magnitude(&self) -> T {
        (self.x().powi(2) + self.y().powi(2) + self.z().powi(2)).sqrt()
    }
}

impl<T: Add<Output = T> + Copy + Float> Add<T> for Vec3D<T> {
    type Output = Vec3D<T>;

    fn add(self, rhs: T) -> Self::Output {
        Self::new(self.x + rhs, self.y + rhs, self.z + rhs)
    }
}

impl<T: Add<Output = T> + Copy + Float> Add<Vec3D<T>> for Vec3D<T> {
    type Output = Vec3D<T>;

    fn add(self, rhs: Vec3D<T>) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl<T: AddAssign + Copy + Float> AddAssign<T> for Vec3D<T> {
    fn add_assign(&mut self, rhs: T) {
        *self.x_mut() += rhs;
        *self.y_mut() += rhs;
        *self.z_mut() += rhs;
    }
}

impl<T: AddAssign + Copy + Float> AddAssign<Vec3D<T>> for Vec3D<T> {
    fn add_assign(&mut self, rhs: Vec3D<T>) {
        *self.x_mut() += *rhs.x();
        *self.y_mut() += *rhs.y();
        *self.z_mut() += *rhs.z();
    }
}

impl<T: Sub<Output = T> + Copy + Float> Sub<T> for Vec3D<T> {
    type Output = Vec3D<T>;

    fn sub(self, rhs: T) -> Self::Output {
        Self::new(self.x - rhs, self.y - rhs, self.z - rhs)
    }
}

impl<T: Sub<Output = T> + Copy + Float> Sub<Vec3D<T>> for Vec3D<T> {
    type Output = Vec3D<T>;

    fn sub(self, rhs: Vec3D<T>) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<T: SubAssign + Copy + Float> SubAssign<T> for Vec3D<T> {
    fn sub_assign(&mut self, rhs: T) {
        *self.x_mut() -= rhs;
        *self.y_mut() -= rhs;
        *self.z_mut() -= rhs;
    }
}

impl<T: SubAssign + Copy + Float> SubAssign<Vec3D<T>> for Vec3D<T> {
    fn sub_assign(&mut self, rhs: Vec3D<T>) {
        *self.x_mut() -= *rhs.x();
        *self.y_mut() -= *rhs.y();
        *self.z_mut() -= *rhs.z();
    }
}

impl<T: Mul<Output = T> + Copy + Float> Mul<T> for Vec3D<T> {
    type Output = Vec3D<T>;

    fn mul(self, rhs: T) -> Self::Output {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl<T: MulAssign + Copy + Float> MulAssign<T> for Vec3D<T> {
    fn mul_assign(&mut self, rhs: T) {
        *self.x_mut() *= rhs;
        *self.y_mut() *= rhs;
        *self.z_mut() *= rhs;
    }
}

impl<T: Mul<Output = T> + Add<Output = T> + Copy + Float> Mul<Vec3D<T>> for Vec3D<T> {
    type Output = T;

    fn mul(self, rhs: Vec3D<T>) -> Self::Output {
        *self.x() * *rhs.x() + *self.y() * *rhs.y() + *self.z() * *rhs.z()
    }
}

impl<T: Div<Output = T> + Copy + Float> Div<T> for Vec3D<T> {
    type Output = Vec3D<T>;

    fn div(self, rhs: T) -> Self::Output {
        Self::new(self.x/rhs, self.y/rhs, self.z/rhs)
    }
}

impl<T: Neg<Output = T> + Float + Copy> Neg for Vec3D<T> {
    type Output = Vec3D<T>;

    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}