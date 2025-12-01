import numpy as np

def radToDeg(rad):
    """Convert radians to degrees."""
    return rad * (180.0 / np.pi)

def degToRad(deg):
    """Convert degrees to radians."""
    return deg * (np.pi / 180.0)

def rotationX(theta):
    """
    Create a rotation matrix for rotation around the X-axis.
    
    Args:
        theta: Rotation angle in radians
    
    Returns:
        3x3 rotation matrix
    """
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)], 
                     [0, np.sin(theta), np.cos(theta)]])

def rotationY(theta):
    """
    Create a rotation matrix for rotation around the Y-axis.
    
    Args:
        theta: Rotation angle in radians
    
    Returns:
        3x3 rotation matrix
    """
    return np.array([[np.cos(theta), 0, np.sin(theta)],
                     [0, 1, 0],
                     [-np.sin(theta), 0, np.cos(theta)]])

def rotationZ(theta):
    """
    Create a rotation matrix for rotation around the Z-axis.
    
    Args:
        theta: Rotation angle in radians
    
    Returns:
        3x3 rotation matrix
    """
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta), np.cos(theta), 0],
                     [0, 0, 1]])

def rotationEulerXYZ(phi, theta, psi):
    """
    Create a rotation matrix from Euler angles using XYZ sequence.
    
    Args:
        phi: Rotation around X-axis in radians
        theta: Rotation around Y-axis in radians
        psi: Rotation around Z-axis in radians
    
    Returns:
        3x3 rotation matrix
    """
    Rx = rotationX(phi)
    Ry = rotationY(theta)
    Rz = rotationZ(psi)
    R = np.matmul(Rx, np.matmul(Ry, Rz))
    return R

def rotationEulerZXZ(phi, theta, psi):
    """
    Create a rotation matrix from Euler angles using ZXZ sequence.
    
    Args:
        phi: First rotation around Z-axis in radians
        theta: Rotation around X-axis in radians
        psi: Second rotation around Z-axis in radians
    
    Returns:
        3x3 rotation matrix
    """
    Rz2 = rotationZ(phi)
    Rx = rotationX(theta)
    Rz1 = rotationZ(psi)
    R = np.matmul(Rz2, np.matmul(Rx, Rz1))
    return R

def eulerAnglesFromRzxz(Rzxz):
    """
    Extract Euler angles (ZXZ sequence) from a rotation matrix.
    
    Args:
        Rzxz: 3x3 rotation matrix
    
    Returns:
        Tuple of (phi, theta, psi) in radians
    """
    phi = np.arctan2(Rzxz[0][2], -Rzxz[1][2]) 
    theta = np.arccos(Rzxz[2][2])
    psi = np.arctan2(Rzxz[2][0], -Rzxz[2][1])
    return (phi, theta, psi)

def eulerAnglesFromRxyz(Rxyz):
    """
    Extract Euler angles (XYZ sequence) from a rotation matrix.
    
    Args:
        Rxyz: 3x3 rotation matrix
    
    Returns:
        Tuple of (phi, theta, psi) in radians
    """
    phi = np.arctan2(Rxyz[1][2], Rxyz[2][2])
    theta = -np.arcsin(Rxyz[0][2]) 
    psi = np.arctan2(Rxyz[0][1], Rxyz[0][0])
    return (phi, theta, psi)

def eulerAngleSequenceXYZToZXZ(phi_xyz, theta_xyz, psi_xyz):
    """
    Convert Euler angles from XYZ sequence to ZXZ sequence.
    
    Args:
        phi_xyz: Rotation around X-axis in radians
        theta_xyz: Rotation around Y-axis in radians
        psi_xyz: Rotation around Z-axis in radians
    
    Returns:
        Tuple of (phi_zxz, theta_zxz, psi_zxz) in radians
    """
    phi_zxz = np.arctan2(-np.sin(theta_xyz), np.sin(phi_xyz) * np.cos(theta_xyz))
    theta_zxz = np.arccos(np.cos(phi_xyz) * np.cos(theta_xyz))
    a = np.cos(psi_xyz) * np.sin(theta_xyz) + np.sin(psi_xyz) * np.sin(phi_xyz) * np.cos(theta_xyz)
    b = np.sin(psi_xyz) * np.sin(theta_xyz) - np.cos(psi_xyz) * np.sin(phi_xyz) * np.cos(theta_xyz)
    psi_zxz = np.arctan2(a, b)
    return (phi_zxz, theta_zxz, psi_zxz)

# Test angles in degrees
phi_xyz_deg = -30.0
theta_xyz_deg = 65.0
psi_xyz_deg = -45.0

print("Euler Angles XYZ")

# Convert to radians
phi_xyz = degToRad(phi_xyz_deg)
theta_xyz = degToRad(theta_xyz_deg)
psi_xyz = degToRad(psi_xyz_deg)

# Create rotation matrix from XYZ Euler angles
Rxyz = rotationEulerXYZ(phi_xyz, theta_xyz, psi_xyz)

# Extract ZXZ Euler angles from rotation matrix
attitude_zxz = eulerAnglesFromRzxz(Rxyz)

# Convert back to degrees and print
phi_zxz_deg = radToDeg(attitude_zxz[0])
theta_zxz_deg = radToDeg(attitude_zxz[1])
psi_zxz_deg = radToDeg(attitude_zxz[2])
print("Euler Angles ZXZ [{}, {}, {}]".format(phi_zxz_deg, theta_zxz_deg, psi_zxz_deg))

# Direct conversion from XYZ to ZXZ
attitude_zxz2 = eulerAngleSequenceXYZToZXZ(phi_xyz, theta_xyz, psi_xyz)

# Convert back to degrees and print
phi_zxz_deg2 = radToDeg(attitude_zxz2[0])
theta_zxz_deg2 = radToDeg(attitude_zxz2[1])
psi_zxz_deg2 = radToDeg(attitude_zxz2[2])
print("Euler Angles ZXZ [{}, {}, {}]".format(phi_zxz_deg2, theta_zxz_deg2, psi_zxz_deg2))
