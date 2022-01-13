# coding: utf-8

'''
This module implements core functions for processing raw IRIS measures (in Volts) to SSA and optical radius.

Authors : Paul Billecocq, Céline Vargel
'''

import numpy as np

#CONSTANTS
DENSITY_OF_ICE = 916.7
RAYONNEMENT_I = 1.3e-6
IM_REFRACTION = 1.3e-5
b = 4.53
k01 = 1.26
k02 = 1.205
k03 = 1.258

def calibration_polynom_fit(spectralon, calibration_values):
    """
    Computes a cubic calibration polynomial function from the calibration measured voltages and the spectralon values
    :param spectralon: percentage of reflectance of the spectralon used for calibration
    :type spectralon: numpy array of floats. For example, reflectance are logged as 80 for 80% in the spreadsheet
    :param calibration_values: voltage measured by IRIS
    :param calibration_values: numpy array of floats
    """

    polynom = np.polyfit(calibration_values, spectralon, deg=3)

    return polynom


def voltage_to_reflectance(voltage, polynom):
    """
    Computes reflectance based on a voltage measurement and the cubic polynomial calibration function
    :param voltage: voltage (Volt) measured by IRIS
    :type voltage: float
    :param polynom: calibration polynom computed for the considered IRIS measurements session
    :type polynom: numpy polyfit object
    """

    reflectance = polynom[0] * voltage**3 + polynom[1] * voltage**2 + polynom[2] * voltage + polynom[3]

    return np.around(reflectance, 2)


def reflectance_to_ssa(reflectance, version):
    """
        Computes SSA based on a reflectance computation and the IRIS version used calibration constant
        :param reflectance: reflectance value computed by voltage_to_reflectance
        :type reflectance: float
        :param version: IRIS version number used for the measurement session
        :type version: int
    """
    if version == 'IRIS_1':
        k0 = k01
    elif version == 'IRIS_2':
        k0 = k02
    elif version == 'IRIS_3':
        k0 = k03
    ssa = 4 * np.pi * IM_REFRACTION / RAYONNEMENT_I * 6 / DENSITY_OF_ICE * (k0 * b / np.log(reflectance / 100)) ** 2


    return np.around(ssa,2)

def ssa_to_optical_radius(ssa):
    """
        Computes optical radius (mm) based on ssa value computated by reflectance_to_ssa
        :param reflectance: SSA value in m :sup:`2`.m :sup:`-3`
        :type reflectance: float
    """

    optical_radius = (3.0 / (DENSITY_OF_ICE * ssa)) * 10**3

    return optical_radius
