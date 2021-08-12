use anyhow::Result;
// use indicatif::{ProgressBar, ProgressIterator, ProgressStyle};

// use nova::{geometry::Vec3D, image::Image};
use nova::geometry::Vec3D;

fn main() -> Result<()> {
    let a = Vec3D::new(1, 2, 3);
    let b = Vec3D::new(4, 5, 6);
    println!("{:#?} * {:#?} = {:#?}", a, b, a * b);
    println!("Magnitude(a) = {}", a.magnitude());
    println!("{:#?}", (a + 1) == ((a + 2) - 1));
    println!("{:#?}", (a - a) == Vec3D::new(0, 0, 0));

    Ok(())
}
