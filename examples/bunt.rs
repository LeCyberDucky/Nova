use anyhow::Result;
use indicatif::{ProgressBar, ProgressIterator, ProgressStyle};
use nova::image::Image;

fn main() -> Result<()> {
    let mut image = Image::new(1024, 1024);
    let progress_bar = ProgressBar::new(*image.height() as u64).with_style(
        ProgressStyle::default_bar().progress_chars("🌞🌈☔").template("{wide_bar:.125} {percent}/100%\n{spinner:.cyan} Elapsed time: {elapsed_precise} | Estimated total time: {duration_precise}")
    );
    for y in (0..*image.height()).progress_with(progress_bar) {
        for x in 0..*image.width() {
            let j = image.height() - 1 - y;
            let i = x;

            let r = i as f64 / (image.width() - 1) as f64;
            let g = j as f64 / (image.height() - 1) as f64;
            // let b = 0.25;
            let b = (i + j) as f64 / (image.width() + image.height() - 2) as f64;
            // let b = (image.width() - i - 1 + j) as f64 / (image.width() + image.height() - 2) as f64;

            let r = (255.999 * r) as u8;
            let g = (255.999 * g) as u8;
            let b = (255.999 * b) as u8;

            image.pixels()[(y, x)].set(r, g, b, 255);
        }
    }
    image.save("Bunt.png")
}
