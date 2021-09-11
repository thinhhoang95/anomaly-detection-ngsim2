import numpy as np
from scipy.interpolate import UnivariateSpline


def lengthen_trajectory_by_spline(t_vec, x_vec, new_t_vec, k=2):
    extrapolator = UnivariateSpline(t_vec, x_vec, k=k)
    return extrapolator(new_t_vec)


def trim_trajectory(t_vec, x_vec, new_t_vec):
    start_index = np.where(t_vec == new_t_vec[0])[0][0]
    end_index = np.where(t_vec == new_t_vec[-1])[0][0]
    if start_index >= end_index:
        raise Exception("Invalid extrapolation time vector. Matching start index and end index caused logic error")
    return x_vec[start_index:end_index + 1]
