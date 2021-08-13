use std::{
    convert::From,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use anyhow::Result;
use ndarray::{self, Array2};

#[derive(Debug, Default, Clone, Copy)]
pub struct Color {
    // Floating point components, as this struct will be used for a lot of math that would probably suffer from rounding errors with integer math. The components should lie in the range [0; 1]
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl Color {
    pub fn new(r: f64, g: f64, b: f64) -> Self {
        Self { r, g, b }
    }
    pub fn clamped(&self, min: f64, max: f64) -> Self {
        Self::new(
            self.r.clamp(min, max),
            self.g.clamp(min, max),
            self.b.clamp(min, max),
        )
    }

    pub const RED: Self = Self {
        r: 1.0,
        g: 0.0,
        b: 0.0
    };

    pub const SKY_BLUE: Self = Self {
        r: 135.0 / 255.0,
        g: 206.0 / 255.0,
        b: 235.0 / 255.0,
    };

    pub const WHITE: Self = Self {
        r: 1.0,
        g: 1.0,
        b: 1.0,
    };
}

impl<T: Into<f64>> Add<T> for Color {
    type Output = Color;

    fn add(self, rhs: T) -> Self::Output {
        let rhs = rhs.into();
        Self::new(self.r + rhs, self.g + rhs, self.b + rhs)
    }
}

impl Add<Color> for Color {
    type Output = Color;

    fn add(self, rhs: Color) -> Self::Output {
        Self::new(self.r + rhs.r, self.g + rhs.g, self.b + rhs.b)
    }
}

impl Add<Color> for i32 {
    type Output = Color;

    fn add(self, rhs: Color) -> Self::Output {
        rhs + self
    }
}

impl Add<Color> for f64 {
    type Output = Color;

    fn add(self, rhs: Color) -> Self::Output {
        rhs + self
    }
}

impl<T: Into<f64>> AddAssign<T> for Color {
    fn add_assign(&mut self, rhs: T) {
        let rhs = rhs.into();
        self.r += rhs;
        self.g += rhs;
        self.b += rhs;
    }
}

impl AddAssign<Color> for Color {
    fn add_assign(&mut self, rhs: Color) {
        self.r += rhs.r;
        self.g += rhs.g;
        self.b += rhs.b;
    }
}

impl<T: Into<f64>> Sub<T> for Color {
    type Output = Color;

    fn sub(self, rhs: T) -> Self::Output {
        self.add(-rhs.into())
    }
}

impl Sub<Color> for Color {
    type Output = Color;

    fn sub(self, rhs: Color) -> Self::Output {
        self.add(-rhs)
    }
}

// Does subtracting a vector from a scalar actually make sense?
// Eh, this non-symmetric operator stuff has mushed my brain 🧠
impl Sub<Color> for i32 {
    type Output = Color;

    fn sub(self, rhs: Color) -> Self::Output {
        -(rhs - self)
    }
}

impl Sub<Color> for f64 {
    type Output = Color;

    fn sub(self, rhs: Color) -> Self::Output {
        -(rhs - self)
    }
}

impl<T: Into<f64>> SubAssign<T> for Color {
    fn sub_assign(&mut self, rhs: T) {
        self.add_assign(-rhs.into())
    }
}

impl SubAssign<Color> for Color {
    fn sub_assign(&mut self, rhs: Color) {
        self.add_assign(-rhs)
    }
}

impl<T: Into<f64>> Mul<T> for Color {
    type Output = Color;

    fn mul(self, rhs: T) -> Self::Output {
        let rhs = rhs.into();
        Self::new(self.r * rhs, self.g * rhs, self.b * rhs)
    }
}

impl<T: Into<f64>> MulAssign<T> for Color {
    fn mul_assign(&mut self, rhs: T) {
        let rhs = rhs.into();
        self.r *= rhs;
        self.g *= rhs;
        self.b *= rhs;
    }
}

impl Mul<Color> for Color {
    type Output = f64;

    fn mul(self, rhs: Color) -> Self::Output {
        self.r * rhs.r + self.g * rhs.g + self.b * rhs.b
    }
}

impl Mul<Color> for i32 {
    type Output = Color;

    fn mul(self, rhs: Color) -> Self::Output {
        rhs * self
    }
}

impl Mul<Color> for f64 {
    type Output = Color;

    fn mul(self, rhs: Color) -> Self::Output {
        rhs * self
    }
}

impl<T: Into<f64>> Div<T> for Color {
    type Output = Color;

    fn div(self, rhs: T) -> Self::Output {
        let rhs = rhs.into();
        Self::new(self.r / rhs, self.g / rhs, self.b / rhs)
    }
}

impl<T: Into<f64>> DivAssign<T> for Color {
    fn div_assign(&mut self, rhs: T) {
        let rhs = rhs.into();
        self.r /= rhs;
        self.g /= rhs;
        self.b /= rhs;
    }
}

impl Neg for Color {
    type Output = Color;

    fn neg(self) -> Self::Output {
        Self::new(-self.r, -self.g, -self.b)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Pixel {
    r: u8,
    g: u8,
    b: u8,
    a: u8,
}

impl Pixel {
    pub fn new(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self { r, g, b, a }
    }

    pub fn set(&mut self, r: u8, g: u8, b: u8, a: u8) {
        self.r = r;
        self.g = g;
        self.b = b;
        self.a = a;
    }

    /// Get a reference to the pixel's a.
    pub fn a(&self) -> &u8 {
        &self.a
    }

    /// Get a mutable reference to the pixel's a.
    pub fn a_mut(&mut self) -> &mut u8 {
        &mut self.a
    }
}

impl Default for Pixel {
    fn default() -> Self {
        Self::new(0, 0, 0, 255)
    }
}

impl From<Color> for Pixel {
    fn from(color: Color) -> Self {
        let color = color.clamped(0.0, 1.0);
        let pixel_max = u8::MAX as f64;
        Pixel {
            r: (color.r * pixel_max).round() as u8,
            g: (color.g * pixel_max).round() as u8,
            b: (color.b * pixel_max).round() as u8,
            ..Default::default()
        }
    }
}

impl From<Color> for image::Rgba<u8> {
    fn from(color: Color) -> Self {
        Pixel::from(color).into()
    }
}

impl From<Pixel> for image::Rgba<u8> {
    fn from(pixel: Pixel) -> Self {
        image::Rgba([pixel.r, pixel.g, pixel.b, pixel.a])
    }
}

pub struct Image {
    data: Array2<Pixel>,
    width: usize,
    height: usize,
}

impl Image {
    fn from_data(data: Array2<Pixel>, width: usize, height: usize) -> Self {
        Self {
            data,
            width,
            height,
        }
    }

    pub fn new(width: usize, height: usize) -> Self {
        let data = Array2::<Pixel>::default((height, width));
        Self::from_data(data, width, height)
    }

    pub fn save<P>(&self, file_path: P) -> Result<()>
    where
        P: AsRef<std::path::Path>,
    {
        image::ImageBuffer::from_fn(self.width as u32, self.height as u32, |x, y| {
            image::Rgba::from(self.data[(y as usize, x as usize)])
        })
        .save(file_path)?;
        Ok(())
    }

    /// Get a reference to the image's width.
    pub fn width(&self) -> &usize {
        &self.width
    }

    /// Get a reference to the image's height.
    pub fn height(&self) -> &usize {
        &self.height
    }

    /// Get a mutable reference to the image's data.
    pub fn pixels(&mut self) -> &mut Array2<Pixel> {
        &mut self.data
    }
}
