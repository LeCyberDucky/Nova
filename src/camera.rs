use crate::geometry::{Point3D, Ray, Vec3D};

pub struct Camera {
    origin: Point3D,
    lower_left_corner: Point3D,
    horizontal: Vec3D,
    vertical: Vec3D,
}

impl Camera {
    pub fn new(
        aspect_ratio: (i32, i32),
        viewport_height: f64,
        focal_length: f64,
        origin: Vec3D,
    ) -> Self {
        // let viewport_height = 2.0;
        let viewport_width = viewport_height * (aspect_ratio.0 as f64 / aspect_ratio.1 as f64);

        let horizontal = Vec3D::new(viewport_width, 0, 0);
        let vertical = Vec3D::new(0, viewport_height, 0);

        // let view_port = aspect_ratio;
        // let horizontal = Vec3D::new(view_port.0, 0, 0);
        // let vertical = Vec3D::new(0, view_port.1, 0);
        let lower_left_corner =
            origin - horizontal / 2 - vertical / 2 - Vec3D::new(0, 0, focal_length);

        Self {
            origin,
            lower_left_corner,
            horizontal,
            vertical,
        }
    }

    pub fn get_ray(&self, u: f64, v: f64) -> Ray {
        let direction =
            self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin;
        Ray::new(self.origin, direction)
    }
}
