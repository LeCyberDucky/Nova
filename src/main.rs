use std::convert::From;

use anyhow::Result;
use ndarray::{self, Array2};

#[derive(Clone, Copy)]
struct Pixel {
    r: u8,
    g: u8,
    b: u8,
    a: u8,
}

impl Pixel {
    fn new(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self { r, g, b, a }
    }

    fn r(&self) -> u8 {
        self.r
    }

    fn g(&self) -> u8 {
        self.g
    }

    fn b(&self) -> u8 {
        self.b
    }

    fn a(&self) -> u8 {
        self.a
    }

    fn set(&mut self, r: u8, g: u8, b: u8, a: u8) {
        self.r = r;
        self.g = g;
        self.b = b;
        self.a = a;
    }
}

impl From<Pixel> for image::Rgba<u8> {
    fn from(pixel: Pixel) -> Self {
        image::Rgba([pixel.r(), pixel.g(), pixel.b(), pixel.a()])
    }
}

impl Default for Pixel {
    fn default() -> Self {
        Self::new(0, 0, 0, 255)
    }
}

struct Image {
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

    fn new(width: usize, height: usize) -> Self {
        let data = Array2::<Pixel>::default((height, width));
        Self::from_data(data, width, height)
    }

    fn save<P>(&self, file_path: P) -> Result<()>
    where
        P: AsRef<std::path::Path>,
    {
        image::ImageBuffer::from_fn(self.width as u32, self.height as u32, |x, y| {
            image::Rgba::from(self.data[(y as usize, x as usize)])
        })
        .save(file_path)?;
        Ok(())
    }
}

fn main() -> Result<()> {
    let mut image = Image::new(1024, 1024);
    for y in 0..image.height {
        for x in 0..image.width {
            let j = image.height - 1 - y;
            let i = x;

            let r = i as f64 / (image.width - 1) as f64;
            let g = j as f64 / (image.height - 1) as f64;
            let b = 0.25;

            let r = (255.999 * r) as u8;
            let g = (255.999 * g) as u8;
            let b = (255.999 * b) as u8;

            image.data[(y, x)].set(r, g, b, 255);
        }
    }

    image.save("Nova.png")
}
