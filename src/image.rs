use std::convert::From;

use anyhow::Result;
use ndarray::{self, Array2};

#[derive(Clone, Copy)]
pub struct Pixel {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
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
}

impl From<Pixel> for image::Rgba<u8> {
    fn from(pixel: Pixel) -> Self {
        image::Rgba([pixel.r, pixel.g, pixel.b, pixel.a])
    }
}

impl Default for Pixel {
    fn default() -> Self {
        Self::new(0, 0, 0, 255)
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
