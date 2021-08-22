use rand::prelude::ThreadRng;

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
    aperture: f64,
    focus_distance: f64,
    lens_radius: f64,
    u: Vec3D,
    v: Vec3D,
    w: Vec3D,
}

impl Camera {
    pub fn new(
        look_from: Point3D,
        look_at: Point3D,
        view_up: Vec3D,
        aspect_ratio: (i32, i32),
        vertical_fov: f64,
        aperture: f64,
        focus_distance: f64,
    ) -> Self {
        let theta = vertical_fov.to_radians();
        let h = (theta / 2.0).tan();

        let viewport_height = 2.0 * h;
        let viewport_width = viewport_height * (aspect_ratio.0 as f64 / aspect_ratio.1 as f64);

        let w = (look_from - look_at).normalized();
        let u = view_up.cross_product(&w).normalized();
        let v = w.cross_product(&u);

        let origin = look_from;
        let horizontal = focus_distance * viewport_width * u;
        let vertical = focus_distance * viewport_height * v;
        let lower_left_corner = origin - horizontal / 2 - vertical / 2 - focus_distance * w;

        let lens_radius = aperture / 2.0;

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
            aperture,
            focus_distance,
            lens_radius,
            u,
            v,
            w,
        }
    }

    pub fn get_ray(&self, s: f64, t: f64, rng: &mut ThreadRng) -> Ray {
        let rd = self.lens_radius * Vec3D::random_in_unit_disk(rng);
        let offset = self.u * rd.x + self.v * rd.y;

        let direction =
            self.lower_left_corner + s * self.horizontal + t * self.vertical - self.origin - offset;
        Ray::new(self.origin + offset, direction)
    }
}
