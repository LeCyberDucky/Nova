use anyhow::Result;
use indicatif::{ProgressIterator, ProgressStyle};

use nova::{
    camera::Camera,
    geometry::Point3D,
    image::{Color, Image},
    obstacle::{material, ObstacleCollection, Sphere},
};
use rand::Rng;

fn main() -> Result<()> {
    // Image
    let aspect_ratio = (16, 9);
    let image_scale = 80;
    let mut image = Image::new(image_scale * aspect_ratio.0, image_scale * aspect_ratio.1);
    let samples_per_pixel = 100;
    let max_ray_recursion_depth = 50;

    // World
    let world = ObstacleCollection::new(vec![
        Box::new(Sphere::new(
            Point3D::new(0, -100.5, -1),
            100,
            material::Lambertian::new(0.8, Color::new(0.8, 0.8, 0.0)),
        )),
        Box::new(Sphere::new(
            Point3D::new(0, 0, -1),
            0.5,
            material::Dielectric::new(1.5, Color::WHITE),
        )),
        Box::new(Sphere::new(
            Point3D::new(-1, 0, -1),
            0.5,
            material::Dielectric::new(1.5, Color::WHITE),
        )),
        Box::new(Sphere::new(
            Point3D::new(1, 0, -1),
            0.5,
            material::Metal::new(Color::new(0.8, 0.6, 0.2), 1.0),
        )),
    ]);

    // Camera
    let camera = Camera::new(
        (aspect_ratio.0 as i32, aspect_ratio.1 as i32),
        2.0,
        1.0,
        Point3D::new(0, 0, 0),
    );

    // Render
    let mut rng = rand::thread_rng();
    let progress_style = ProgressStyle::default_bar().progress_chars("🌞🌈☔").template("{wide_bar:.125} {percent}/100%\n{spinner:.cyan} Elapsed time: {elapsed_precise} | Estimated total time: {duration_precise}");

    for y in (0..*image.height()).progress_with_style(progress_style) {
        for x in 0..*image.width() {
            let j = image.height() - 1 - y;
            let i = x;
            let mut color = Color::BLACK;
            for _ in 0..samples_per_pixel {
                // Sampling pixles with random offset to simulate antialiasing
                let u = (i as f64 + rng.gen::<f64>()) / (image.width() - 1) as f64;
                let v = (j as f64 + rng.gen::<f64>()) / (image.height() - 1) as f64;

                let ray = camera.get_ray(u, v);
                color += ray.color(&world, &mut rng, max_ray_recursion_depth);
            }
            image.pixels()[(y, x)] = (color / samples_per_pixel).into();
        }
    }
    image.save("Nova.png")?;

    Ok(())
}
