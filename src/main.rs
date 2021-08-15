use anyhow::Result;
use indicatif::{ProgressIterator, ProgressStyle};

use nova::{
    camera::Camera,
    geometry::{ObstacleCollection, Point3D, Sphere, Vec3D},
    image::{Color, Image},
};
use rand::Rng;

fn main() -> Result<()> {
    // Image
    let aspect_ratio = (16, 9);
    let image_scale = 80;
    let mut image = Image::new(image_scale * aspect_ratio.0, image_scale * aspect_ratio.1);
    let samples_per_pixel = 100;

    // World
    let world = ObstacleCollection::new(vec![
        Box::new(Sphere::new(Vec3D::new(0, 0, -1), 0.5)),
        Box::new(Sphere::new(Vec3D::new(0, -100.5, -1), 100)),
    ]);

    // Camera
    let camera = Camera::new(
        (aspect_ratio.0 as i32, aspect_ratio.1 as i32),
        1,
        Point3D::new(0, 0, 0),
    );

    // Render
    let mut rng = rand::thread_rng();
    let progress_style = ProgressStyle::default_bar().progress_chars("🌞🌈☔").template("{wide_bar:.125} {percent}/100%\n{spinner:.cyan} Elapsed time: {elapsed_precise} | Estimated total time: {duration_precise}");

    for y in (0..*image.height()).progress_with_style(progress_style) {
        for x in 0..*image.width() {
            let mut color = Color::BLACK;
            for _s in 0..samples_per_pixel {
                // Sampling pixles with random offset to simulate antialiasing
                let j = image.height() - 1 - y;
                let i = x;

                let u = (i as f64 + rng.gen::<f64>()) / (image.width() - 1) as f64;
                let v = (j as f64 + rng.gen::<f64>()) / (image.height() - 1) as f64;

                let ray = camera.get_ray(u, v);
                color += ray.color(&world) / samples_per_pixel as f64;
            }

            image.pixels()[(y, x)] = color.into();
        }
    }
    image.save("Nova.png")?;

    Ok(())
}
