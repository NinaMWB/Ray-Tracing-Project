from matplotlib import pyplot as plt
import numpy as np

class Ray:
    """
    A class to represent an optical ray in 3D space.
    
    Attributes:
    _points : list of np.ndarray
        List of 3D points along the ray path.
    _direction : np.ndarray
        The current (normalised) direction of the ray.
    """

    def __init__(self, pos=None, direc=None):
        """
        Initialise a new Ray object.

        Parameters:
        pos : list or np.ndarray, optional
            Initial position of the ray. Default is [0, 0, 0].
        direc : list or np.ndarray, optional
            Initial direction of the ray. Default is [0, 0, 1].
        """
        self._points = [np.array(pos) if pos is not None else np.array([0.0, 0.0, 0.0])]
        init_dir = np.array(direc) if direc is not None else np.array([0.0, 0.0, 1.0])
        self._direction = self._normalise(init_dir)

    def _normalise(self, v):
        #Return the normalised version of a vector
        norm = np.linalg.norm(v)
        if norm == 0:
            raise ValueError("Direction vector cannot be zero.")
        return v / norm

    def pos(self):
        #Return the current (last) point of the ray
        return self._points[-1]

    def direc(self):
        #Return the current direction of the ray
        return self._direction

    def append(self, pos, direc):
        #Append a new point and direction to the ray.
        self._points.append(np.array(pos))
        self._direction = self._normalise(np.array(direc))

    def vertices(self):
        #Return all the points along the ray path.
        return self._points

class RayBundle:
    def __init__(self, rmax=5.0, nrings=5, multi=6):
        self._rays = []
        # Generate rings
        for i in range(nrings):
            r = (i + 1) * rmax / nrings
            for j in range(multi):
                theta = 2 * np.pi * j / multi
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                self._rays.append(Ray(pos=[x, y, 0], direc=[0, 0, 1]))
        # Add one central (paraxial) ray
        self._rays.append(Ray(pos=[0, 0, 0], direc=[0, 0, 1]))

    def rays(self):
        #Return the list of Ray objects
        return self._rays

    def propagate_bundle(self, elements):
        #Propagate all rays through each optical element in order
        for ray in self._rays:
            for elem in elements:
                elem.propagate_ray(ray)

    def track_plot(self):
        #Plot the path of each ray in the bundle. Returns matplotlib figure
        fig, ax = plt.subplots()
        for ray in self._rays:
            path = np.array(ray.vertices())
            ax.plot(path[:, 2], path[:, 0], alpha=0.7)  # z vs x
        ax.set_title("Ray Bundle Propagation")
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("x (mm)")
        ax.grid(True)
        return fig
    
    def rms(self):
        #Return the RMS (root mean square) spot size from optical axis
        positions = [ray.pos()[:2] for ray in self._rays]  #final (x, y) only
        squared_radii = [np.sum(pos**2) for pos in positions]
        rms = np.sqrt(np.mean(squared_radii))
        return rms
    
    def spot_plot(self):
        #Plot the x-y spot diagram of the ray bundle (at final z).
        x_vals = [ray.pos()[0] for ray in self._rays]
        y_vals = [ray.pos()[1] for ray in self._rays]

        fig, ax = plt.subplots()
        ax.scatter(x_vals, y_vals, alpha=0.7)
        ax.set_aspect('equal')
        ax.set_title("Spot Diagram at Focal Plane")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.grid(True)
        return fig
