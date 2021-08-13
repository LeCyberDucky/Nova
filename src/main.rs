use anyhow::Result;
// use indicatif::{ProgressBar, ProgressIterator, ProgressStyle};

use indicatif::{ProgressIterator, ProgressStyle};
use nova::{
    geometry::{Point3D, Ray, Vec3D},
    image::Image,
};

fn main() -> Result<()> {
    // Image
    let aspect_ratio = (16, 9);
    let image_scale = 80;
    let mut image = Image::new(image_scale * aspect_ratio.0, image_scale * aspect_ratio.1);

    // Camera
    let view_port = (aspect_ratio.0 as i32, aspect_ratio.1 as i32);
    let focal_length = 1;

    let origin = Point3D::new(0, 0, 0);
    let horizontal = Vec3D::new(view_port.0, 0, 0);
    let vertical = Vec3D::new(0, view_port.1, 0);
    let lower_left_corner = origin - horizontal / 2 - vertical / 2 - Vec3D::new(0, 0, focal_length);

    // Render
    let progress_style = ProgressStyle::default_bar().progress_chars("🌞🌈☔").template("{wide_bar:.125} {percent}/100%\n{spinner:.cyan} Elapsed time: {elapsed_precise} | Estimated total time: {duration_precise}");

    for y in (0..*image.height()).progress_with_style(progress_style) {
        for x in 0..*image.width() {
            let j = image.height() - 1 - y;
            let i = x;

            let u = i as f64 / (image.width() - 1) as f64;
            let v = j as f64 / (image.height() - 1) as f64;

            let direction = lower_left_corner + u * horizontal + v * vertical - origin;
            let ray = Ray::new(origin, direction);
            let pixel = ray.color().into();

            image.pixels()[(y, x)] = pixel;
        }
    }
    image.save("Nova.png")?;

    Ok(())
}
