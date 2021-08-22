use rand::prelude::ThreadRng;

use crate::{
    geometry::{Point3D, Ray, Vec3D},
    image::Color,
};

use self::material::Material;

pub trait Obstacle: Material {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit>;
}

pub struct ObstacleCollection {
    obstacles: Vec<Box<dyn Obstacle>>,
}

impl ObstacleCollection {
    pub fn new(obstacles: Vec<Box<dyn Obstacle>>) -> Self {
        Self { obstacles }
    }
    fn clear(&mut self) {
        self.obstacles.clear();
    }

    fn add(&mut self, obstacle: Box<dyn Obstacle>) {
        self.obstacles.push(obstacle);
    }

    pub fn hit(&self, ray: &Ray, t_min: f64, mut t_max: f64) -> Option<(usize, Hit)> {
        let mut closest_hit = None;
        let mut hit_id = 0;

        for (id, obstacle) in self.obstacles.iter().enumerate() {
            if let Some(hit) = obstacle.hit(ray, t_min, t_max) {
                t_max = hit.t;
                closest_hit = Some(hit);
                hit_id = id;
            }
        }
        Some((hit_id, closest_hit?))
    }

    pub fn scatter(
        &self,
        incident_ray: &Ray,
        hit: &Hit,
        obstacle_id: usize,
        rng: &mut ThreadRng,
    ) -> Option<(Ray, Color)> {
        self.obstacles[obstacle_id].scatter(incident_ray, hit, rng)
    }
}

pub struct Hit {
    p: Point3D,
    normal: Vec3D,
    t: f64,
    front_face: bool,
}

impl Hit {
    fn get_face_normal(ray: &Ray, outward_normal: Vec3D) -> (bool, Vec3D) {
        let front_face = (ray.direction() * outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        (front_face, normal)
    }

    /// Get a reference to the hit's normal.
    pub fn normal(&self) -> &Vec3D {
        &self.normal
    }

    /// Get a reference to the hit's p.
    pub fn p(&self) -> &Point3D {
        &self.p
    }

    /// Minimum threshold for hits. Choosing this value too low will cause shadow acne, since determining the intersection of rays and obstacles is sensitive to floatin point inaccuracies
    pub const THRESHOLD: f64 = 1e-4;
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Sphere<A: Material> {
    center: Point3D,
    radius: f64,
    material: A,
}

impl<A: Material> Sphere<A> {
    pub fn new<T: Into<f64>>(center: Point3D, radius: T, material: A) -> Self {
        Self {
            center,
            radius: radius.into(),
            material,
        }
    }

    pub fn outward_normal(&self, surface_point: Point3D) -> Vec3D {
        (surface_point - self.center) / self.radius
    }
}

impl<A: Material> Obstacle for Sphere<A> {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<Hit> {
        let co = ray.origin() - self.center;
        let a = ray.direction().magnitude_squared();
        let half_b = co * ray.direction();
        let c = co.magnitude_squared() - self.radius.powi(2);

        let discriminant = half_b.powi(2) - a * c;
        if discriminant < 0.0 {
            return None;
        };
        let discriminant_sqrt = discriminant.sqrt();

        // Find the nearest root that lies in the acceptable range
        let mut root = (-half_b - discriminant_sqrt) / a;
        // if root < t_min || t_max < root {
        if !(t_min..=t_max).contains(&root) {
            root = (discriminant_sqrt - half_b) / a;
            // if root < t_min || t_max < root {
            if !(t_min..=t_max).contains(&root) {
                return None;
            }
        }

        let p = ray.at(root);
        let outward_normal = self.outward_normal(p);
        let (front_face, normal) = Hit::get_face_normal(ray, outward_normal);
        Some(Hit {
            p,
            normal,
            t: root,
            front_face,
        })
    }
}

impl<A: Material> Material for Sphere<A> {
    fn scatter(&self, incident_ray: &Ray, hit: &Hit, rng: &mut ThreadRng) -> Option<(Ray, Color)> {
        self.material.scatter(incident_ray, hit, rng)
    }
}

pub mod material {
    use rand::Rng;

    use super::*;
    //     use rand::prelude::ThreadRng;

    // use crate::{geometry::{Hit, Ray, Vec3D}, image::Color};

    pub trait Material {
        /// Rays are absorbed with probability (1 - reflectance). Reflected rays lose energy given by the attenuation. This is a color instead of a floating point value, since the material might especially absorb a specific color (From a physics stand point, it would probably be more correct to give each ray its own energy level (color) and work with that, but eh...)
        fn scatter(
            &self,
            incident_ray: &Ray,
            hit: &Hit,
            rng: &mut ThreadRng,
        ) -> Option<(Ray, Color)>;
    }

    #[derive(Clone, Copy)]
    pub struct Lambertian {
        albedo: f64,      // Fraction of incident energy that is reflected  [0; 1]
        reflectance: f64, // Fraction of rays being reflected [0; 1]
        attenuation: Color,
        // Incident rays are completely absorbed with probability (1 - reflectance). Reflected rays lose an amount of energy energy (attenuation) that balances albedo and reflectance. I.e. if 50% of rays are reflected, but albedo is only 10%, the reflected rays need to lose 80% of their energy, in order to bring the total reflected energy (albedo) to 10%
        // It must hold, that albedo <= reflectance. Otherwise, the material would have to increase the energy of the reflected rays, or produce additional rays
        // Attenuation = albedo/reflectance
    }

    impl Lambertian {
        pub fn new(albedo: f64, attenuation: Color) -> Self {
            assert!(
                (0.0..=1.0).contains(&albedo),
                "Albedo should be between 0.0 and 1.0"
            );
            assert!(attenuation >= Color::BLACK && attenuation <= Color::WHITE);
            let max_energy = Color::WHITE.r + Color::WHITE.g + Color::WHITE.b;

            // Assuming the red, green, and blue components all weigh equally in terms of energy
            let attenuated_energy = max_energy / (attenuation.r + attenuation.g + attenuation.b);
            assert!(attenuated_energy >= albedo); // If the energy left after attenuating is lower than the albedo, we need more than 100% reflectance to reach the desired albedo. I.e. the material would need to actively produce light, instead of only reflecting

            let reflectance = albedo / attenuated_energy;

            Self {
                albedo,
                reflectance,
                attenuation,
            }
        }
    }

    impl Default for Lambertian {
        fn default() -> Self {
            Self::new(0.85, Color::new(0.8, 0.8, 0.8))
        }
    }

    impl Material for Lambertian {
        fn scatter(
            &self,
            _incident_ray: &Ray,
            hit: &Hit,
            rng: &mut ThreadRng,
        ) -> Option<(Ray, Color)> {
            if !rng.gen_bool(self.reflectance) {
                // No reflection. Ray absorbed
                return None;
            }
            let mut scatter_direction = hit.normal + Vec3D::random_on_unit_sphere(rng);
            if scatter_direction.near_zero() {
                // Would lead to NaNs and infinities later on
                scatter_direction = hit.normal
            }
            let scattered_ray = Ray::new(hit.p, scatter_direction);
            Some((scattered_ray, self.attenuation))
        }
    }

    #[derive(Clone, Copy)]
    pub struct Metal {
        attenuation: Color,
        fuzz: f64, // Level of fuzziness [0; 1]
    }

    impl Metal {
        pub fn new(attenuation: Color, fuzz: f64) -> Self {
            // Not sure about all these asserts here. Should this return a result instead?
            // Eh, I guess panicking is okay, because this means that there is a logic error that needs to be taken care of
            assert!(
                (0.0..=1.0).contains(&fuzz),
                "Level of fuzziness should be in range [0; 1]"
            );

            assert!(attenuation >= Color::BLACK && attenuation <= Color::WHITE);

            Self { attenuation, fuzz }
        }
    }

    impl Default for Metal {
        fn default() -> Self {
            Self::new(Color::new(0.8, 0.8, 0.8), 0.0)
        }
    }

    impl Material for Metal {
        fn scatter(
            &self,
            incident_ray: &Ray,
            hit: &Hit,
            rng: &mut ThreadRng,
        ) -> Option<(Ray, Color)> {
            let mut reflection_direction =
                incident_ray.direction().normalized().reflect(&hit.normal);
            if self.fuzz > 0.0 {
                reflection_direction += self.fuzz * Vec3D::random_in_unit_sphere(rng);
            }
            let scattered_ray = Ray::new(hit.p, reflection_direction);
            ((scattered_ray.direction() * hit.normal) > 0.0)
                .then_some((scattered_ray, self.attenuation))
        }
    }

    pub struct Dielectric {
        refraction_index: f64,
        attenuation: Color,
    }

    impl Dielectric {
        pub fn new(refraction_index: f64, attenuation: Color) -> Self {
            Self {
                refraction_index,
                attenuation,
            }
        }

        pub fn reflectance(cosine: f64, refraction_index: f64) -> f64 {
            // Use Schlick's approximation for reflectance
            let r0 = ((1.0 - refraction_index) / (1.0 + refraction_index)).powi(2);
            r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
        }
    }

    impl Default for Dielectric {
        fn default() -> Self {
            Self {
                refraction_index: 1.0,
                attenuation: Color::BLACK,
            }
        }
    }

    impl Material for Dielectric {
        fn scatter(
            &self,
            incident_ray: &Ray,
            hit: &Hit,
            _rng: &mut ThreadRng,
        ) -> Option<(Ray, Color)> {
            let refraction_ratio = if hit.front_face {
                1.0 / self.refraction_index
            } else {
                self.refraction_index
            };

            let unit_direction = incident_ray.direction().normalized();

            let cos_theta = (-unit_direction * hit.normal()).min(1.0);
            let sin_theta = (1.0 - cos_theta.powi(2)).sqrt();

            let cannot_refract = (refraction_ratio * sin_theta) > 1.0;

            let direction =
                if cannot_refract || Dielectric::reflectance(cos_theta, refraction_ratio) > 1.0 {
                    unit_direction.reflect(hit.normal())
                } else {
                    unit_direction.refract(hit.normal(), refraction_ratio)
                };

            Some((Ray::new(hit.p, direction), self.attenuation))
        }
    }
}
