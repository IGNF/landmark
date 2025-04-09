cdef extern from "math.h":
    double sin(double x)
    double sqrt(double x)
    double fabs(double x)
    double atan2(double y, double x)
    double M_PI

cdef double EPSILON = 1e-7  # Approx equivalent to np.finfo(np.float32).eps

cpdef tuple facet(double e0, double e1, double e2, double delta_x, double delta_y):
    """
    Cython-accelerated version of the facet function.

    Parameters
    ----------
    e0 : float
        Elevation at the central point.
    e1 : float
        Elevation at the first neighbor point.
    e2 : float
        Elevation at the second neighbor point.
    delta_x : float
        Grid spacing in the x-direction.
    delta_y : float
        Grid spacing in the y-direction.

    Returns
    -------
    tuple of float
        r : float
            Computed aspect in radians.
        s_max_facet : float
            Maximum slope of the facet (positive downward).
    """
    cdef double s1, s2, r_val, sp, sd, s_max_facet

    s1 = (e0 - e1) / delta_x
    s2 = (e1 - e2) / delta_x

    if fabs(s1) < EPSILON:
        if s2 >= 0.0:
            r_val = M_PI / 2.0
        else:
            r_val = -M_PI / 2.0
    else:
        r_val = atan2(s2, s1)

    sp = sqrt(s1 * s1 + s2 * s2)
    sd = (e0 - e2) / sqrt(delta_x * delta_x + delta_y * delta_y)

    if 0.0 <= r_val <= M_PI / 4.0 and s1 >= 0.0:
        s_max_facet = sp
    elif s1 > sd:
        s_max_facet = s1
        r_val = 0.0
    else:
        s_max_facet = sd
        r_val = M_PI / 4.0

    return r_val, s_max_facet
