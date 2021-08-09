use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Vec3D {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3D {
    pub fn new<T: Into<f64>>(x: T, y: T, z: T) -> Self {
        Self {
            x: x.into(),
            y: y.into(),
            z: z.into(),
        }
    }

    fn x(&self) -> f64 {
        self.x
    }

    fn y(&self) -> f64 {
        self.y
    }

    fn z(&self) -> f64 {
        self.z
    }

    fn x_mut(&mut self) -> &mut f64 {
        &mut self.x
    }

    fn y_mut(&mut self) -> &mut f64 {
        &mut self.y
    }

    fn z_mut(&mut self) -> &mut f64 {
        &mut self.z
    }

    pub fn magnitude(&self) -> f64 {
        (self.x().powi(2) + self.y().powi(2) + self.z().powi(2)).sqrt()
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

impl<T: Into<f64>> AddAssign<T> for Vec3D {
    fn add_assign(&mut self, rhs: T) {
        let rhs = rhs.into();
        *self.x_mut() += rhs;
        *self.y_mut() += rhs;
        *self.z_mut() += rhs;
    }
}

impl AddAssign<Vec3D> for Vec3D {
    fn add_assign(&mut self, rhs: Vec3D) {
        *self.x_mut() += rhs.x();
        *self.y_mut() += rhs.y();
        *self.z_mut() += rhs.z();
    }
}

impl<T: Into<f64>> Sub<T> for Vec3D {
    type Output = Vec3D;

    fn sub(self, rhs: T) -> Self::Output {
        let rhs = rhs.into();
        Self::new(self.x - rhs, self.y - rhs, self.z - rhs)
    }
}

impl Sub<Vec3D> for Vec3D {
    type Output = Vec3D;

    fn sub(self, rhs: Vec3D) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<T: Into<f64>> SubAssign<T> for Vec3D {
    fn sub_assign(&mut self, rhs: T) {
        let rhs = rhs.into();
        *self.x_mut() -= rhs;
        *self.y_mut() -= rhs;
        *self.z_mut() -= rhs;
    }
}

impl SubAssign<Vec3D> for Vec3D {
    fn sub_assign(&mut self, rhs: Vec3D) {
        *self.x_mut() -= rhs.x();
        *self.y_mut() -= rhs.y();
        *self.z_mut() -= rhs.z();
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
        *self.x_mut() *= rhs;
        *self.y_mut() *= rhs;
        *self.z_mut() *= rhs;
    }
}

impl Mul<Vec3D> for Vec3D {
    type Output = f64;

    fn mul(self, rhs: Vec3D) -> Self::Output {
        self.x() * rhs.x() + self.y() * rhs.y() + self.z() * rhs.z()
    }
}

impl<T: Into<f64>> Div<T> for Vec3D {
    type Output = Vec3D;

    fn div(self, rhs: T) -> Self::Output {
        let rhs = rhs.into();
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl Neg for Vec3D {
    type Output = Vec3D;

    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}
