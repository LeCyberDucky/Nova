use anyhow::Result;
use indicatif::{ProgressIterator, ProgressStyle};

use nova::{
    camera::Camera,
    geometry::{Point3D, Vec3D},
    image::{Color, Image},
    obstacle::{material, ObstacleCollection, Sphere},
};
use rand::Rng;

fn main() -> Result<()> {
    let mut rng = rand::thread_rng();

    // Image
    let aspect_ratio = (16, 9);
    let image_scale = 80;
    let mut image = Image::new(image_scale * aspect_ratio.0, image_scale * aspect_ratio.1);
    let samples_per_pixel = 100;
    let max_ray_recursion_depth = 50;

    // World
    // let R = (PI/4.0).cos();
    let mut world = ObstacleCollection::new(vec![
        Box::new(Sphere::new(Point3D::new(0, -1000, 0), 1000, material::Lambertian::new(0.8, Color::new(0.5, 0.5, 0.5)))),
        Box::new(Sphere::new(
            Point3D::new(0, 1, 0), 1, material::Dielectric::new(1.5, Color::WHITE)
        )),
        Box::new(Sphere::new(Point3D::new(-4, 1, 0), 1, material::Lambertian::new(0.8, Color::new(0.4, 0.2, 0.1)))),
        Box::new(Sphere::new(Point3D::new(4, 1, 0), 1, material::Metal::new(Color::new(0.7, 0.6, 0.5), 0.0)))
    ]);

    for a in -11..11 {
        for b in -11..11 {
            let material = rng.gen::<f64>();
            let center = Point3D::new(a as f64 + 0.9*rng.gen::<f64>(), 0.2, b as f64 + 0.9*rng.gen::<f64>());

            if (center - Point3D::new(4, 0.2, 0)).magnitude() > 0.9 {
                if material < 0.8 {
                    // Diffuse
                    let albedo = rng.gen();
                    let attenuation = Vec3D::random_in_unit_sphere(&mut rng);
                    let attenuation = Color::new(attenuation.x.abs(), attenuation.y.abs(), attenuation.z.abs());
                    world.add(Box::new(
                        Sphere::new(center, 0.2, material::Lambertian::new(albedo, attenuation))
                    ));
                } else if material < 0.95 {
                    // Metal
                    let attenuation = Vec3D::random_in_range(&mut rng, 0.5, 1.0);
                    let attenuation = Color::new(attenuation.x, attenuation.y, attenuation.z);
                    let fuzz = rng.gen_range(0.0..=0.5);
                    world.add(Box::new(
                        Sphere::new(center, 0.2, material::Metal::new(attenuation, fuzz))
                    ));
                } else {
                    // Glass
                    world.add(Box::new(
                        Sphere::new(center, 0.2, material::Dielectric::new(1.5, Color::WHITE))
                    ));
                }
            }
        }
    }

    // Camera
    let look_from = Point3D::new(13, 2, 3);
    let look_at = Point3D::new(0, 0, 0);
    let camera = Camera::new(
        look_from,
        look_at,
        Vec3D::new(0, 1, 0),
        (aspect_ratio.0 as i32, aspect_ratio.1 as i32),
        20.0,
        0.1,
        10.0,
    );

    // Render
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

                let ray = camera.get_ray(u, v, &mut rng);
                color += ray.color(&world, &mut rng, max_ray_recursion_depth);
            }
            image.pixels()[(y, x)] = (color / samples_per_pixel).into();
        }
    }
    image.save("Nova.png")?;

    Ok(())
}
