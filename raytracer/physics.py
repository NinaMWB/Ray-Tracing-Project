import numpy as np

def refract(direc, normal, n_1, n_2):
    """
    Apply Snell's law to compute the refracted ray direction.

    Parameters:
    direc : np.ndarray
        Incident direction (normalised 3D vector).
    normal : np.ndarray
        Surface normal vector (normalised).
    n_1 : float
        Refractive index of the current medium.
    n_2 : float
        Refractive index of the next medium.

    Returns:
    np.ndarray or None
        Refracted direction (normalised), or None if total internal reflection occurs.
    """
    #Ensure all vectors are unit vectors
    direc = direc / np.linalg.norm(direc)
    normal = normal / np.linalg.norm(normal)

    cos_theta_i = -np.dot(direc, normal)
    sin2_theta_i = 1.0 - cos_theta_i**2
    eta = n_1 / n_2
    sin2_theta_t = eta**2 * sin2_theta_i

    #Check for total internal reflection
    if sin2_theta_t > 1.0:
        return None

    cos_theta_t = np.sqrt(1.0 - sin2_theta_t)

    #Compute refracted direction
    refracted = eta * direc + (eta * cos_theta_i - cos_theta_t) * normal
    return refracted / np.linalg.norm(refracted)
