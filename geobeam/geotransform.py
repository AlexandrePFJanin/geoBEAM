# -*- coding: utf-8 -*-
"""
@author: Alexandre JANIN
@aim:    Geographic transformation package
"""

# External dependencies:
import numpy as np


# ----------------- FUNCTIONS -----------------


def normal2fault(strike, dip):
    """
    Compute the unit normal vector to a fault plane
    from strike and dip angles.

    Args:
        strike (float): Strike angle in degrees (clockwise from North)
        dip (float): Dip angle in degrees (0 = horizontal, 90 = vertical),
                     dipping to the right relative to strike

    Returns:
        n (ndarray, shape (3,)): Unit normal vector to the fault plane (x, y, z)
    """
    strike_rad = np.deg2rad(strike)
    dip_rad = np.deg2rad(dip)

    # Unit vector along strike
    s = np.stack([
        np.sin(strike_rad),
        np.cos(strike_rad),
        np.zeros_like(strike_rad)
    ], axis=-1)

    # Horizontal dip azimuth (to the right of strike)
    dip_az = strike_rad + np.pi/2

    # Unit vector along dip (POSITIVE UPWARD)
    d = np.stack([
        -np.cos(dip_rad) * np.sin(dip_az),
        -np.cos(dip_rad) * np.cos(dip_az),
         np.sin(dip_rad)
    ], axis=-1)

    # Right-hand normal
    n = np.cross(s, d)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    return n




def displacement_sdt_to_xyz(u_s, u_d, u_t, strike, dip):
    """
    Convert displacement(s) from (strike, dip, tensile)
    to (x, y, z) coordinates.

    Args:
        u_s, u_d, u_t (float or array-like, shape (N,)):
            Displacement components in strike, dip (upward), tensile directions
        strike, dip (float or array-like, shape (N,)):
            Strike (deg, clockwise from North) and dip (deg, right-dipping)

    Returns:
        u_xyz (ndarray, shape (3,) for scalar input or (N, 3) for array input):
            Output x,y,z displacement.
    """

    u_s = np.asarray(u_s)
    u_d = np.asarray(u_d)
    u_t = np.asarray(u_t)
    strike = np.asarray(strike)
    dip = np.asarray(dip)

    strike_rad = np.deg2rad(strike)
    dip_rad = np.deg2rad(dip)

    # Unit vector along strike
    s = np.stack([
        np.sin(strike_rad),
        np.cos(strike_rad),
        np.zeros_like(strike_rad)
    ], axis=-1)

    # Horizontal dip azimuth (to the right of strike)
    dip_az = strike_rad + np.pi/2

    # Unit vector along dip (POSITIVE UPWARD)
    d = np.stack([
        -np.cos(dip_rad) * np.sin(dip_az),
        -np.cos(dip_rad) * np.cos(dip_az),
         np.sin(dip_rad)
    ], axis=-1)

    # Right-hand normal
    n = np.cross(s, d)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)

    # Transformation
    u_xyz = (
        u_s[..., None] * s +
        u_d[..., None] * d +
        u_t[..., None] * n
    )

    return u_xyz