import numpy as np
from raytracer.physics import refract

class OpticalElement:
    """
    Base class for optical elements. Not meant to be instantiated directly.

    Methods:
    intercept(ray):
        Determine where the given ray intersects the element.
    
    propagate_ray(ray):
        Propagate the ray through this element.
    """

    def intercept(self, ray):
        """
        Calculate where the ray intersects the element.
        
        Parameters:
        ray : Ray
            The ray to test for intersection.

        Returns:
        intercept_point : numpy array or None
            The point of intersection, or None if there is no intersection.
        """
        raise NotImplementedError("intercept() needs to be implemented in a subclass.")

    def propagate_ray(self, ray):
        """
        Propagate a ray through the element.

        Parameters:
        ray : Ray
            The ray to be propagated.

        1. Check for intersection,
        2. Apply changes (e.g. reflection/refraction),
        3. Update the ray’s internal state.
        """
        raise NotImplementedError("propagate_ray() needs to be implemented in a subclass.")

class SphericalRefraction(OpticalElement):
    """
    Represents a spherical refracting surface.

    Parameters:
    z_0 : float
        Intercept of the surface with the z-axis.
    aperture : float
        Aperture radius (maximum extent from optical axis).
    curvature : float
        Surface curvature (1 / radius of curvature).
    n_1 : float
        Refractive index before the surface (z < z_0).
    n_2 : float
        Refractive index after the surface (z > z_0).
    """

    def __init__(self, z_0, aperture, curvature, n_1, n_2):
        self._z_0 = z_0
        self._aperture = aperture
        self._curvature = curvature
        self._n_1 = n_1
        self._n_2 = n_2

    def z_0(self):
        #Return the intercept of the surface with the z-axis
        return self._z_0

    def aperture(self):
        #Return the aperture radius
        return self._aperture

    def curvature(self):
        #Return the curvature of the surface
        return self._curvature

    def n_1(self):
        #Return refractive index on the left (z < z_0)
        return self._n_1

    def n_2(self):
        #Return refractive index on the right (z > z_0)
        return self._n_2

    def intercept(self, ray):
        #Ray position and direction
        P = ray.pos()
        k = ray.direc()

        #Sphere centre (on z-axis at z_0)
        O = np.array([0.0, 0.0, self._z_0])
        R = abs(1.0 / self._curvature) if self._curvature != 0 else np.inf

        r = P - O
        k_dot_r = np.dot(k, r)
        r_squared = np.dot(r, r)

        discriminant = k_dot_r**2 - (r_squared - R**2)
        if discriminant < 0:
            return None  #No real intersection

        #Compute the two roots
        sqrt_disc = np.sqrt(discriminant)
        l1 = -k_dot_r + sqrt_disc
        l2 = -k_dot_r - sqrt_disc

        #Choose smallest positive distance
        l_candidates = [l for l in [l1, l2] if l > 1e-6]
        if not l_candidates:
            return None

        l_min = min(l_candidates)
        Q = P + l_min * k  #intersection point

        #Check if it's within the aperture
        radial_distance = np.linalg.norm(Q[:2])  #x-y plane distance
        if radial_distance <= self._aperture:
            return Q
        return None

    def propagate_ray(self, ray):
        #Find intersection
        intercept_point = self.intercept(ray)
        if intercept_point is None:
            return  #No valid intersection

        #Compute surface normal at the intercept point
        #Sphere centre
        centre = np.array([0.0, 0.0, self._z_0])
        normal = intercept_point - centre
        normal = normal / np.linalg.norm(normal)

        #Flip normal if ray is entering from outside (dot product < 0)
        if np.dot(ray.direc(), normal) > 0:
            normal = -normal

        #Apply Snell’s law 
        new_direc = refract(ray.direc(), normal, self._n_1, self._n_2)
        if new_direc is None:
            return  #Total internal reflection – stop propagation

        #Append new point and direction to the ray
        ray.append(intercept_point, new_direc)

    def focal_point(self):
        #Return the z-position of the paraxial focal point (in mm)
        f = self._n_2 / ((self._n_2 - self._n_1) * self._curvature)
        return self._z_0 + f

class OutputPlane(OpticalElement):
    def __init__(self, z_0):
        self._z_0 = z_0

    def z_0(self):
        #Return the z-position of the output plane
        return self._z_0

    def intercept(self, ray):
        #Find intersection of ray with the plane z = z_0.
        #Returns the intersection point, or None if direction is parallel to the plane.

        pos = ray.pos()
        direc = ray.direc()

        if direc[2] == 0:
            return None  #Ray is parallel to plane, no intersection

        l = (self._z_0 - pos[2]) / direc[2]
        if l <= 0:
            return None  #Intersection is behind the ray

        intercept_point = pos + l * direc
        return intercept_point

    def propagate_ray(self, ray):
        #Append the intercept point to the ray without changing its direction.
        
        point = self.intercept(ray)
        if point is not None:
            ray.append(point, ray.direc())  #Direction stays the same