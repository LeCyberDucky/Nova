use std::io::Write;

use anyhow::Result;
use ndarray::{self, Array2, array};

struct Pixel {
    data: ndarray::Array1<u8>,
}

impl Pixel {
    fn new(r: u8, g: u8, b: u8) -> Self {
        Self {
            data: array![r, g, b]
        }
    }

    fn r(&self) -> u8 {
        self.data[0]
    }

    fn g(&self) -> u8 {
        self.data[1]
    }

    fn b(&self) -> u8 {
        self.data[2]
    }

    fn set(&mut self, r: u8, g: u8, b: u8) {
        self.data[0] = r;
        self.data[1] = g;
        self.data[2] = b;
    }
}

impl Default for Pixel {
    fn default() -> Self {
        Self::new(0, 0, 0)
    }
}

struct Image {
    data: Array2<Pixel>,
    width: usize,
    height: usize
}

impl Image {
    fn from_data(data: Array2<Pixel>, width: usize, height: usize) -> Self { Self { data, width, height } }

    fn new(width: usize, height: usize) -> Self {
        let data = Array2::<Pixel>::default((height, width));
        Self::from_data(data, width, height)
    }

    fn save<P>(&self, file_path: P) -> Result<()> 
    where P: AsRef<std::path::Path>
    {
        let mut file = std::fs::File::create(file_path)?;
        let header = format!("P3\n\
        {} {}\n\
        {}\n", self.width, self.height, u8::MAX);
        file.write_all(header.as_bytes())?;

        for y in 0..self.height {
            for x in 0..self.width {
                let pixel = &self.data[(y, x)];
                let content = format!("{} {} {}\n", pixel.r(), pixel.g(), pixel.b());
                file.write_all(content.as_bytes())?;
            }
        }

        Ok(())
    }
}


fn main() -> Result<()> {
    let mut image = Image::new(256, 256);
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

            image.data[(y, x)].set(r, g, b);
        }
    }

    image.save("Nova.ppm")?;

    Ok(())
}
