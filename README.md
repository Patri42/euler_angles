# Euler Angles Rotation Library

A Python library for working with 3D rotation matrices and Euler angles conversions.

## Description

This library provides functions to create rotation matrices, convert between different Euler angle conventions (XYZ and ZXZ sequences), and extract Euler angles from rotation matrices. It's useful for robotics, computer graphics, aerospace applications, and any field requiring 3D orientation calculations.

## Features
Angle Conversion: Convert between radians and degrees
Rotation Matrices: Generate rotation matrices for X, Y, and Z axes
Euler Angle Sequences: Support for XYZ and ZXZ Euler angle conventions
Matrix to Angles: Extract Euler angles from rotation matrices
Sequence Conversion: Convert between XYZ and ZXZ Euler angle sequences

## Usage Example
import numpy as np

# Define angles in degrees
phi_deg = -30.0
theta_deg = 65.0
psi_deg = -45.0

# Convert to radians
phi = degToRad(phi_deg)
theta = degToRad(theta_deg)
psi = degToRad(psi_deg)

# Create rotation matrix using XYZ Euler angles
R = rotationEulerXYZ(phi, theta, psi)

# Extract ZXZ Euler angles from the rotation matrix
angles_zxz = eulerAnglesFromRzxz(R)

# Convert back to degrees
phi_zxz = radToDeg(angles_zxz[0])
theta_zxz = radToDeg(angles_zxz[1])
psi_zxz = radToDeg(angles_zxz[2])

print(f"ZXZ Euler Angles: [{phi_zxz}, {theta_zxz}, {psi_zxz}]")
