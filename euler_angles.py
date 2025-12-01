import numpy as np

def radToDeg (rad):
  return rad * (180.0/np.pi)

def degToRad(deg):
  return deg * (np.pi/180.0)

def rotationX(theta):
  return np.array([1,0,0],
                  [0,np.cos(theta), np.sin(theta)],
                  [0, -np.sin(theta), np.cos(theta)]])
  
def rotationY(theta):
  return np.array ([np.cos(theta), 0, -np.sin(theta)],
                   [0, 1, 0],
                   [np.sin (theta), 0, np.cos(theta)]])

def rotationZ(theta):
  return np.array ([[np.cos(theta), np.sin(theta),0],
                    [-np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])

def rotationEulerXYZ (phi, theta, psi):
  Rx = rotationX(phi)
  Ry = rotationY(theta)
  Rz = rotation(psi)
  R = np.matmul (Rx, np.matmul(Ry, Rz))
  return R

def rotationEulerZXZ (phi, theta, psi):
  Rz2 = rotationZ(phi)
  Rx = rotationX(theta)
  Rz1 = rotation(psi)
  R = np.matmul (Rz2, np.matmul (Rx, Rz1))
  return R
  
