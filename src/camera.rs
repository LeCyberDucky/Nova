use crate::geometry::{Point3D, Ray, Vec3D};

pub struct Camera {
    look_from: Point3D,
    look_at: Point3D,
    view_up: Vec3D,
    origin: Point3D,
    lower_left_corner: Point3D,
    horizontal: Vec3D,
    vertical: Vec3D,
    vertical_fov: f64, // Degrees
    aspect_ratio: (i32, i32),
}

impl Camera {
    pub fn new(
        look_from: Point3D,
        look_at: Point3D,
        view_up: Vec3D,
        aspect_ratio: (i32, i32),
        vertical_fov: f64,
    ) -> Self {
        let theta = vertical_fov.to_radians();
        let h = (theta / 2.0).tan();

        let viewport_height = 2.0 * h;
        let viewport_width = viewport_height * (aspect_ratio.0 as f64 / aspect_ratio.1 as f64);

        let w = (look_from - look_at).normalized();
        let u = view_up.cross_product(&w).normalized();
        let v = w.cross_product(&u);

        let origin = look_from;
        let horizontal = viewport_width * u;
        let vertical = viewport_height * v;
        let lower_left_corner = origin - horizontal / 2 - vertical / 2 - w;

        Self {
            look_from,
            look_at,
            view_up,
            origin,
            lower_left_corner,
            horizontal,
            vertical,
            vertical_fov,
            aspect_ratio,
        }
    }

    pub fn get_ray(&self, s: f64, t: f64) -> Ray {
        let direction =
            self.lower_left_corner + s * self.horizontal + t * self.vertical - self.origin;
        Ray::new(self.origin, direction)
    }
}
