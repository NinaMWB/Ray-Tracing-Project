from raytracer.elements import SphericalRefraction

class PlanoConvex:
    def __init__(self, z_0, curvature1, curvature2, n_inside, n_outside, thickness):
        """
        Constructs a plano-convex (or convex-plano) lens as two spherical surfaces.
        Parameters:
            z_0         : z-position of the first surface
            curvature1  : curvature of the first surface (positive for convex)
            curvature2  : curvature of the second surface
            n_inside    : refractive index inside the lens
            n_outside   : refractive index of surrounding medium
            thickness   : axial distance between surfaces
        """
        self.z_0 = z_0
        self.z_1 = z_0 + thickness

        self.surface1 = SphericalRefraction(
            z_0=z_0,
            aperture=10,
            curvature=curvature1,
            n_1=n_outside,
            n_2=n_inside
        )

        self.surface2 = SphericalRefraction(
            z_0=self.z_1,
            aperture=10,
            curvature=curvature2,
            n_1=n_inside,
            n_2=n_outside
        )

    def elements(self):
        #Return the list of two optical surfaces (in order)
        return [self.surface1, self.surface2]
