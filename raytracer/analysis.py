"""Analysis module."""
import matplotlib.pyplot as plt
import numpy as np
from rays import Ray, RayBundle
from elements import SphericalRefraction, OutputPlane
from physics import refract
from lenses import PlanoConvex

def find_focal_point(lens, rmax=5.0):
    """
    Estimate focal point by propagating paraxial rays and finding z where rays cross.
    """
    #Try a range of z values to find minimal RMS
    zs = np.linspace(120, 200, 100)
    best_rms = float('inf')
    best_z = None

    for z in zs:
        bundle = RayBundle(rmax=rmax, nrings=5, multi=6)
        output = OutputPlane(z_0=z)
        bundle.propagate_bundle(lens.elements() + [output])
        rms = bundle.rms()
        if rms < best_rms:
            best_rms = rms
            best_z = z

    return best_z

def task8():
    """
    Task 8.

    In this function you should check your propagate_ray function properly
    finds the correct intercept and correctly refracts a ray. Don't forget
    to check that the correct values are appended to your Ray object.
    """
    #Create a spherical surface
    surface = SphericalRefraction(z_0=10, aperture=5, curvature=0.2, n_1=1.0, n_2=1.5)

    #Test ray 1: Straight on-axis
    ray1 = Ray(pos=[0, 0, 0], direc=[0, 0, 1])
    surface.propagate_ray(ray1)
    print("Ray 1 Vertices:", ray1.vertices())
    print("Ray 1 Direction:", ray1.direc())

    #Test ray 2: Slightly angled, still within aperture
    ray2 = Ray(pos=[0, 0, 0], direc=[0.1, 0, 0.995])  # slightly off-axis
    surface.propagate_ray(ray2)
    print("Ray 2 Vertices:", ray2.vertices())
    print("Ray 2 Direction:", ray2.direc())

    #Test ray 3: Way outside aperture
    ray3 = Ray(pos=[10, 0, 0], direc=[0, 0, 1])
    surface.propagate_ray(ray3)
    print("Ray 3 Vertices:", ray3.vertices())
    print("Ray 3 Direction:", ray3.direc())

def task10():
    """
    Task 10.

    In this function you should create Ray objects with the given initial positions.
    These rays should be propagated through the surface, up to the output plane.
    You should then plot the tracks of these rays.
    This function should return the matplotlib figure of the ray paths.

    Returns:
        Figure: the ray path plot.
    """
    #Define optical elements
    sph_surface = SphericalRefraction(z_0=100, aperture=10, curvature=0.03, n_1=1.0, n_2=1.5)
    output = OutputPlane(z_0=250)
    elements = [sph_surface, output]

    #Initial y-positions for input rays
    y_positions = [4, 1, 0.2, 0, -0.2, -1, -4]
    
    rays = []
    for y in y_positions:
        ray = Ray(pos=[0, y, 0], direc=[0, 0, 1])
        for elem in elements:
            elem.propagate_ray(ray)
        rays.append(ray)

    #Plotting
    fig, ax = plt.subplots()
    for ray in rays:
        path = np.array(ray.vertices())
        ax.plot(path[:, 2], path[:, 1])  #z vs y

    ax.set_xlabel("z (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Ray Paths through Spherical Refracting Surface")
    ax.axvline(250, color='grey', linestyle='--', label='Output Plane')
    ax.legend()
    ax.grid(True)
    
    return fig

def task11():
    """
    Task 11.

    In this function you should propagate the three given paraxial rays through the system
    to the output plane and the tracks of these rays should then be plotted.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for ray paths
    2. the calculated focal point.

    Returns:
        tuple[Figure, float]: the ray path plot and the focal point
    """
    #Create the spherical surface
    surface = SphericalRefraction(z_0=100, aperture=10, curvature=0.03, n_1=1.0, n_2=1.5)

    #Calculate paraxial focal point
    focal_z = surface.focal_point()

    #Create output plane at focal point
    output = OutputPlane(z_0=focal_z)

    #Near-axis input rays
    input_positions = [[0.1, 0.1, 0],[0, 0, 0],[-0.1, -0.1, 0]]

    rays = []
    for pos in input_positions:
        ray = Ray(pos=pos, direc=[0, 0, 1])
        for elem in [surface, output]:
            elem.propagate_ray(ray)
        rays.append(ray)

    #Plot rays
    fig, ax = plt.subplots()
    for ray in rays:
        path = np.array(ray.vertices())
        ax.plot(path[:, 2], path[:, 0], label=f'y={path[0,1]:.1f} mm')

    ax.axvline(focal_z, color='gray', linestyle='--', label='Focal Plane')
    ax.set_title('Paraxial Rays Converging at Focal Point')
    ax.set_xlabel('z (mm)')
    ax.set_ylabel('x (mm)')
    ax.grid(True)
    ax.legend()

    return fig, focal_z

def task12():
    """
    Task 12.

    In this function you should create a RayBunble and propagate it to the output plane
    before plotting the tracks of the rays.
    This function should return the matplotlib figure of the track plot.

    Returns:
        Figure: the track plot.
    """
    #Set up optical elements
    surface = SphericalRefraction(z_0=100, aperture=10, curvature=0.03, n_1=1.0, n_2=1.5)
    output = OutputPlane(z_0=surface.focal_point())

    #Create and propagate bundle
    bundle = RayBundle()
    bundle.propagate_bundle([surface, output])

    #Plot and return figure
    return bundle.track_plot()

def task13():
    """
    Task 13.

    In this function you should again create and propagate a RayBundle to the output plane
    before plotting the spot plot.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot
    2. the simulation RMS

    Returns:
        tuple[Figure, float]: the spot plot and rms
    """
    #Optical system setup
    surface = SphericalRefraction(z_0=100, curvature=0.03, aperture=10, n_1=1.0, n_2=1.5)
    output = OutputPlane(z_0=surface.focal_point())

    #Ray bundle propagation
    bundle = RayBundle()
    bundle.propagate_bundle([surface, output])

    #Generate plot and RMS
    fig = bundle.spot_plot()
    spot_rms = bundle.rms()

    return fig, spot_rms

def task14():
    """
    Task 14.

    In this function you will trace a number of RayBundles through the optical system and
    plot the RMS and diffraction scale dependence on input beam radii.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the simulation RMS for input beam radius 2.5
    3. the diffraction scale for input beam radius 2.5

    Returns:
        tuple[Figure, float, float]: the plot, the simulation RMS value, the diffraction scale.
    """
    lam = 588e-6  #mm
    radii = np.linspace(0.1, 10, 30)  #from 0.1 mm to 10 mm
    rms_values = []
    diff_values = []

    #Setup surface once
    surface = SphericalRefraction(z_0=100, aperture=10, curvature=0.03, n_1=1.0, n_2=1.5)
    focal_length = surface.focal_point() - surface.z_0()

    for r in radii:
        bundle = RayBundle(rmax=r, nrings=5, multi=6)
        output = OutputPlane(z_0=surface.focal_point())
        bundle.propagate_bundle([surface, output])
        rms_values.append(bundle.rms())

        D = 2 * r
        diff_scale = lam * focal_length / D
        diff_values.append(diff_scale)

    #Plot results
    fig, ax = plt.subplots()
    ax.plot(radii, rms_values, label='RMS Spot Size (mm)')
    ax.plot(radii, diff_values, label='Diffraction Scale (mm)', linestyle='--')
    ax.set_xlabel("Bundle Radius (mm)")
    ax.set_ylabel("Spot Size / Diffraction Scale (mm)")
    ax.set_title("RMS vs Diffraction Limit vs Bundle Radius")
    ax.grid(True)
    ax.legend()

    #Also return value at r = 2.5 mm
    idx_2_5 = np.argmin(np.abs(radii - 2.5))
    return fig, rms_values[idx_2_5], diff_values[idx_2_5]

def task15():
    """
    Task 15.

    In this function you will create plano-convex lenses in each orientation and propagate a RayBundle
    through each to their respective focal point. You should then plot the spot plot for each orientation.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot for the plano-convex system
    2. the focal point for the plano-convex lens
    3. the matplotlib figure object for the spot plot for the convex-plano system
    4  the focal point for the convex-plano lens


    Returns:
        tuple[Figure, float, Figure, float]: the spot plots and rms for plano-convex and convex-plano.
    """
    
    #Common lens properties
    z0 = 100
    thickness = 5
    n_air = 1.0
    n_glass = 1.5168
    rmax = 5

    #Plane-Convex
    pc = PlanoConvex(z_0=z0, curvature1=0.0, curvature2=-0.02, n_inside=n_glass, n_outside=n_air, thickness=thickness)
    pc_focal = find_focal_point(pc, rmax)
    bundle1 = RayBundle(rmax=rmax)
    bundle1.propagate_bundle(pc.elements() + [OutputPlane(z_0=pc_focal)])
    fig1 = bundle1.spot_plot()

    #Convex-Plane
    cp = PlanoConvex(z_0=z0, curvature1=0.02, curvature2=0.0, n_inside=n_glass, n_outside=n_air, thickness=thickness)
    cp_focal = find_focal_point(cp, rmax)
    bundle2 = RayBundle(rmax=rmax)
    bundle2.propagate_bundle(cp.elements() + [OutputPlane(z_0=cp_focal)])
    fig2 = bundle2.spot_plot()

    return fig1, pc_focal, fig2, cp_focal

def task16():
    """
    Task 16.

    In this function you will be again plotting the radial dependence of the RMS and diffraction values
    for each orientation of your lens.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the RMS for input beam radius 3.5 for the plano-convex system
    3. the RMS for input beam radius 3.5 for the convex-plano system
    4  the diffraction scale for input beam radius 3.5

    Returns:
        tuple[Figure, float, float, float]: the plot, RMS for plano-convex, RMS for convex-plano, diffraction scale.
    """
    lam = 588e-6 #mm 
    z0 = 100 #location of first surface
    thickness = 5
    n_in = 1.5168
    n_out = 1.0
    radii = np.linspace(0.1, 10, 30)

    rms_pc, rms_cp, diff_vals = [], [], []

    #Fix focal points once at r = 3.5 mm to avoid unstable estimates
    pc = PlanoConvex(z_0=z0, curvature1=0.0, curvature2=-0.02, n_inside=n_in, n_outside=n_out, thickness=thickness)
    cp = PlanoConvex(z_0=z0, curvature1=0.02, curvature2=0.0, n_inside=n_in, n_outside=n_out, thickness=thickness)
    f_pc = find_focal_point(pc, 3.5)
    f_cp = find_focal_point(cp, 3.5)

    for r in radii:
        #Plane–Convex
        bundle_pc = RayBundle(rmax=r)
        bundle_pc.propagate_bundle(pc.elements() + [OutputPlane(z_0=f_pc)])
        rms_pc.append(bundle_pc.rms())

        #Convex–Plane
        bundle_cp = RayBundle(rmax=r)
        bundle_cp.propagate_bundle(cp.elements() + [OutputPlane(z_0=f_cp)])

        # Defensive check for broken rays
        if all(len(ray.vertices()) > 1 for ray in bundle_cp.rays()):
            rms_cp.append(bundle_cp.rms())
        else:
            rms_cp.append(np.nan)  # skip broken rays

        #Diffraction scale
        diff = lam * f_pc / (2 * r) # type: ignore
        diff_vals.append(diff)

    #Convert for plotting
    rms_pc = np.array(rms_pc)
    rms_cp = np.array(rms_cp)
    diff_vals = np.array(diff_vals)
    radii = np.array(radii)

    #Plot
    fig, ax = plt.subplots()
    ax.plot(radii, rms_pc, label="Plane–Convex RMS", lw=1.5)
    ax.plot(radii[~np.isnan(rms_cp)], rms_cp[~np.isnan(rms_cp)], label="Convex–Plane RMS", lw=1.5)
    ax.plot(radii, diff_vals, '--', label="Diffraction Scale", lw=1.5)
    ax.set_xlabel("Input Beam Radius (mm)")
    ax.set_ylabel("Spot Size (mm)")
    ax.set_title("RMS vs Diffraction Limit for Plano-Convex Lens")
    ax.grid(True)
    ax.legend()

    #Output values at r = 3.5 mm
    idx_3_5 = np.argmin(np.abs(radii - 3.5))
    return fig, rms_pc[idx_3_5], rms_cp[idx_3_5], diff_vals[idx_3_5]

if __name__ == "__main__":

    # Run task 8 function
    task8()

    # Run task 10 function
    FIG10 = task10()

    # Run task 11 function
    FIG11, FOCAL_POINT = task11()

    # Run task 12 function
    FIG12 = task12()

    # Run task 13 function
    FIG13, TASK13_RMS = task13()

    # Run task 14 function
    FIG14, TASK14_RMS, TASK14_DIFF_SCALE = task14()

    # Run task 15 function
    FIG15_PC, FOCAL_POINT_PC, FIG15_CP, FOCAL_POINT_CP = task15()

    # Run task 16 function
    FIG16, PC_RMS, CP_RMS, TASK16_DIFF_SCALE = task16()

    plt.show()
