from calendar import c
import math
from re import A
import rlcompleter
from symtable import Symbol
from tkinter import BaseWidget
from sympy import symbols, simplify, pprint,  expand_trig, trigsimp, pi, atan2, cos, sin, tan, sqrt, atan, atanh
from sympy.matrices import Matrix

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
import numpy as np

import transformations

from test1 import R_corr

global robot_d1, robot_d2, robot_d3, robot_d4, robot_d5, robot_d6, robot_d7
global robot_a0, robot_a1, robot_a2, robot_a3, robot_a4, robot_a5, robot_a6
global robot_alpha0, robot_alpha1, robot_alpha2, robot_alpha3, robot_alpha4, robot_alpha5, robot_alpha6
global robot_theta1, robot_theta2, robot_theta3, robot_theta4, robot_theta5, robot_theta6, robot_theta7

# Doosan
robot_d1 = 0.135
robot_d2 = 0.00625
robot_d3 = 0.0
robot_d4 = 0.368
robot_d5 = 0.0
robot_d6 = 0.121
robot_d7 = 0.0

robot_a0 = 0.0
robot_a1 = 0.0
robot_a2 = 0.411
robot_a3 = 0.0
robot_a4 = 0.0
robot_a5 = 0.0
robot_a6 = 0.0

robot_alpha0 = 0.0
robot_alpha1 = -pi/2
robot_alpha2 = 0.0
robot_alpha3 = pi/2
robot_alpha4 = -pi/2
robot_alpha5 = pi/2
robot_alpha6 = 0.0

robot_theta1 = 0.0 
robot_theta2 = -pi/2 
robot_theta3 = pi/2 
robot_theta4 = 0.0
robot_theta5 = 0.0
robot_theta6 = 0.0
robot_theta7 = 0.0



def printMatrix(a):
  print("Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]")
  rows = a.shape[0]
  cols = a.shape[1]
  for i in range(0,rows):
    for j in range(0,cols):
        print ("%6.f\t" %a[i,j], end = '')
    print()
  print()      

def printMatrixE(a):
  print("Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]")
  rows = a.shape[0]
  cols = a.shape[1]
  for i in range(0,rows):
    for j in range(0,cols):
        print("%6.15f\t" %a[i,j], end = '')
    print()
  print()    

# Transformation matrix Classic i-1Ti (https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters))
def transformationMatrixClassic(theta, alpha, a, d):
  # returns the pose T of one joint frame i with respect to the previous joint frame (i - 1)
  # given the parameters:
  # theta: theta[i]
  # alpha: alpha[i-1]
  # a: a[i-1]
  # d: d[i]

  r11 = math.cos(theta)
  r12 = -math.sin(theta) * math.cos(alpha)
  r13 = math.sin(theta) * math.sin(alpha) 
  r14 = a * math.cos(theta)

  r21 = math.sin(theta)
  r22 = math.cos(theta) * math.cos(alpha)
  r23 = -math.cos(theta) * math.sin(alpha)
  r24 = a * math.sin(theta)

  r31 = 0.0
  r32 = math.sin(alpha)
  r33 = math.cos(alpha)
  r34 = d
    
  T = Matrix([[r11, r12, r13, r14],
              [r21, r22, r23, r24],
              [r31, r32, r33, r34],
              [0.0, 0.0, 0.0, 1.0]])

  #T = simplify(T)

  # print("Transformation for each link")
  # print("T : ")
  # pprint(T)

  return T

# Transformation matrix Classic i-1Ti (https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters)
def transformationMatrixClassicSymbol(theta, alpha, a, d):
  # returns the pose T of one joint frame i with respect to the previous joint frame (i - 1)
  # given the parameters:
  # theta: theta[i]
  # alpha: alpha[i-1]
  # a: a[i-1]
  # d: d[i]

  r11 = cos(theta)
  r12 = -sin(theta) * cos(alpha)
  r13 = sin(theta) * sin(alpha) 
  r14 = a * cos(theta)

  r21 = sin(theta)
  r22 = cos(theta) * cos(alpha)
  r23 = -cos(theta) * sin(alpha)
  r24 = a * sin(theta)

  r31 = 0.0
  r32 = sin(alpha)
  r33 = cos(alpha)
  r34 = d
    
  T = Matrix([[r11, r12, r13, r14],
              [r21, r22, r23, r24],
              [r31, r32, r33, r34],
              [0.0, 0.0, 0.0, 1.0]])

  #T = simplify(T)

  # print("Transformation for each link")
  # print("T : ")
  # pprint(T)

  return T


# Transformation matrix Modified i-1Ti (https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters)
def transformationMatrixModified(theta, alpha, a, d):
  # returns the pose T of one joint frame i with respect to the previous joint frame (i - 1)
  # given the parameters:
  # theta: theta[i]
  # alpha: alpha[i-1]
  # a: a[i-1]
  # d: d[i]

  r11 = math.cos(theta)
  r12 = -math.sin(theta)
  r13 = 0.0
  r14 = a

  r21 = math.sin(theta) * math.cos(alpha)
  r22 = math.cos(theta) * math.cos(alpha)
  r23 = -math.sin(alpha)
  r24 = -d * math.sin(alpha)

  r31 = math.sin(theta) * math.sin(alpha)
  r32 = math.cos(theta) * math.sin(alpha)
  r33 = math.cos(alpha)
  r34 = d * math.cos(alpha)
    
  T = Matrix([[r11, r12, r13, r14],
              [r21, r22, r23, r24],
              [r31, r32, r33, r34],
              [0.0, 0.0, 0.0, 1.0]])

  #T = simplify(T)

  # print("Transformation for each link")
  # print("T : ")
  # pprint(T)

  return T

# Transformation matrix Modified i-1Ti (https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters)
def transformationMatrixModifiedSymbol(theta, alpha, a, d):
  # returns the pose T of one joint frame i with respect to the previous joint frame (i - 1)
  # given the parameters:
  # theta: theta[i]
  # alpha: alpha[i-1]
  # a: a[i-1]
  # d: d[i]

  r11 = cos(theta)
  r12 = -sin(theta)
  r13 = 0
  r14 = a

  r21 = sin(theta) * cos(alpha)
  r22 = cos(theta) * cos(alpha)
  r23 = -sin(alpha)
  r24 = -d * sin(alpha)

  r31 = sin(theta) * sin(alpha)
  r32 = cos(theta) * sin(alpha)
  r33 = cos(alpha)
  r34 = d * cos(alpha)
    
  T = Matrix([[r11, r12, r13, r14],
              [r21, r22, r23, r24],
              [r31, r32, r33, r34],
              [0.0, 0.0, 0.0, 1.0]])

  #T = simplify(T)

  # print("Transformation for each link")
  # print("T : ")
  # pprint(T)

  return T


# Rotation matrix Rotxn(alphan) (https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters)
def rotationMatrixSymbol(theta, alpha, a, d):
  # returns the pose T of one joint frame i with respect to the previous joint frame (i - 1)
  # given the parameters:
  # theta: theta[i]
  # alpha: alpha[i-1]
  # a: a[i-1]
  # d: d[i]

  r11 = cos(theta)
  r12 = -sin(theta)
  r13 = 0.0
  r14 = 0.0

  r21 = sin(theta)
  r22 = cos(theta)
  r23 = 0.0
  r24 = 0.0

  r31 = 0.0
  r32 = 0.0
  r33 = 1.0
  r34 = 0.0
    
  T = Matrix([[r11, r12, r13, r14],
              [r21, r22, r23, r24],
              [r31, r32, r33, r34],
              [0.0, 0.0, 0., 1,0]])
  
  #T = simplify(T)

  # print("Transformation for each link")
  # print("T : ")
  # pprint(T)

  return T






def transformationMatrixModifiedNewSymbol(theta, alpha, a, d):
  # returns the pose T of one joint frame i with respect to the previous joint frame (i - 1)
  # given the parameters:
  # theta: theta[i]
  # alpha: alpha[i-1]
  # a: a[i-1]
  # d: d[i]

  r21 = a * cos(theta)
  r22 = cos(theta)
  r23 = -cos(alpha) * sin(theta)
  r24 = sin(alpha) * sin(theta)

  r31 = a * sin(theta)
  r32 = sin(theta)
  r33 = cos(alpha) * cos(theta)
  r34 = -cos(theta) * sin(alpha)

  r41 = d
  r42 = 0.0
  r43 = sin(alpha)
  r44 = cos(alpha)

  T = Matrix([[1.0, 0.0, 0.0, 0.0],
              [r21, r22, r23, r24],
              [r31, r32, r33, r34],
              [r41, r42, r43, r44]])

  #T = simplify(T)

  # print("Transformation for each link")
  # print("T : ")
  # pprint(T)

  return T






# Calculates Rotation Matrix given euler angles (https://mathworld.wolfram.com/RotationMatrix.html, https://mathworld.wolfram.com/EulerAngles.html, https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix)
def eulerAnglesToRotationMatrix(angle) :     
    R_x = np.array([[1,               0,              	  0               ],
                    [0,               cos(angle[0]),      -sin(angle[0])  ],
                    [0,               sin(angle[0]),      cos(angle[0])   ]])

    R_z1 = np.array([[cos(angle[0]),   -sin(angle[0]),    0               ],
                     [sin(angle[0]),   cos(angle[0]),     0               ],
                     [0,               0,                 1               ]])    

    R_y = np.array([[cos(angle[1]),    0,                 sin(angle[1])  ],
                    [0,                1,                 0               ],
                    [-sin(angle[1]),    0,                cos(angle[1])   ]])
                 
    R_z2 = np.array([[cos(angle[2]),   -sin(angle[2]),    0               ],
                     [sin(angle[2]),   cos(angle[2]),     0               ],
                     [0,               0,                 1               ]])                     
                     
    R = np.dot(R_z1, np.dot(R_y, R_z2))
    #R = R_z1 * R_y * R_z2
    # R = np.dot(R_z, np.dot( R_y, R_x ))
 
    return R










def rot_x(q):
    r_x = Matrix([[ 1,      0,       0,       0],
                  [ 0, cos(q), -sin(q),       0],
                  [ 0, sin(q),  cos(q),       0],
                  [ 0,      0,       0,       1]])
    return r_x


def rot_y(q):
    r_y = Matrix([[  cos(q), 0, sin(q),       0],
                  [       0, 1,      0,       0],
                  [ -sin(q), 0, cos(q),       0],
                  [       0, 0,      0,       1]])
    return r_y


def rot_z(q):
    r_z = Matrix([[ cos(q), -sin(q), 0,       0],
                  [ sin(q),  cos(q), 0,       0],
                  [      0,       0, 1,       0],
                  [      0,       0, 0,       1]])
    return r_z













def Rx(angle):
  R_x = Matrix([[1,               0,                0,              0],
                [0,               cos(angle),       sin(angle),     0],
                [0,               -sin(angle),      cos(angle),     0],
                [0,               0,                0,              1]])
  return R_x

def Ry(angle):
  R_y = Matrix([[cos(angle),      0,                -sin(angle),    0],
                [0,               1,                0,              0],
                [sin(angle),      0,                cos(angle),     0],
                [0,               0,                0,              1]])
  return R_y

def Rz(angle):
  R_z = Matrix([[cos(angle),      sin(angle),       0,              0],
                [-sin(angle),     cos(angle),       0,              0],
                [0,               0,                1,              0],
                [0,               0,                0,              1]])    
  return R_z



# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6

# Calculates rotation matrix to euler angles
# the result is the same as MATLAB except the order
# of the euler angles ( x and z are swapped ).
def rotationMatrixToEulerAngles(R) : 
    assert(isRotationMatrix(R))
     
    sy = sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
     
    singular = sy < 1e-6
 
    if not singular :
        x = atan2(R[2,1] , R[2,2])
        y = atan2(-R[2,0], sy)
        z = atan2(R[1,0], R[0,0])
    else :
        x = atan2(-R[1,2], R[1,1])
        y = atan2(-R[2,0], sy)
        z = 0
 
    return np.array([x, y, z])

def get_hypotenuse(a, b):
  # calculate the longest side given the two shorter sides 
  # of a right triangle using pythagorean theorem
  return sqrt(a*a + b*b)

def get_cosine_law_angle(a, b, c):    
  # given all sides of a triangle a, b, c
  # calculate angle gamma between sides a and b using cosine law
  cos_gamma = (a*a + b*b - c*c) / (2*a*b)
  sin_gamma = sqrt(1 - cos_gamma * cos_gamma)
  gamma = atan2(sin_gamma, cos_gamma)

  return gamma

def get_wrist_center(gripper_point, R0g, d6, d7):
  # get the coordinates of the wrist center wrt to the base frame (xw, yw, zw)
  # given the following info:
  # the coordinates of the gripper (end effector) (x, y, z)
  # the rotation of the gripper in gripper frame wrt to the base frame (R0u)
  # the distance between gripper and wrist center dg which is along common z axis
  # check WRITEUP.pdf for more info
  xu, yu, zu = gripper_point 
    
  nx, ny, nz = R0g[0, 2], R0g[1, 2], R0g[2, 2]
  # xw = xu - dg * nx
  # yw = yu - dg * ny
  # zw = zu - dg * nz 

  xw = xu - (d6 + d7) * nx
  yw = yu - (d6 + d7) * ny
  zw = zu - (d6 + d7) * nz


  return xw, yw, zw

def get_first_three_angles(wrist_center):
  # given the wrist center which a tuple of 3 numbers x, y, z
  # (x, y, z) is the wrist center point wrt base frame
  # return the angles q1, q2, q3 for each respective joint
  # given geometry of the kuka kr210
  # check WRITEUP.pdf for more info
  x, y, z  = wrist_center
    
  l = get_hypotenuse(robot_d4, abs(robot_a3))
  phi = atan2(robot_d4, abs(robot_a3))
  
  x_prime = get_hypotenuse(x, y)
  mx = x_prime -  robot_a1
  my = z - robot_d1 
  m = get_hypotenuse(mx, my)
  alpha = atan2(my, mx)
  
  gamma = get_cosine_law_angle(l, robot_a2, m)
  beta = get_cosine_law_angle(m, robot_a2, l)
  
  q1 = atan2(y, x)
  q2 = pi/2 - beta - alpha 
  q3 = -(gamma - phi)
    
  return q1, q2, q3 

def get_last_three_angles(R):
  #Recall that from our simplification, R36 (R) equals the following:
  #From trigonometry we can get q4, q5, q6 if we know numerical values of all cells of matrix R36 (R)
  #check WRITEUP.pdf for more info    
  sin_q4 = R[2, 2]
  cos_q4 =  -R[0, 2]
    
  sin_q5 = sqrt(R[0, 2]**2 + R[2, 2]**2) 
  cos_q5 = R[1, 2]
    
  sin_q6 = -R[1, 1]
  cos_q6 = R[1, 0] 
  
  q4 = atan2(sin_q4, cos_q4)
  q5 = atan2(sin_q5, cos_q5)
  q6 = atan2(sin_q6, cos_q6)
    
  return q4, q5, q6






def forwardKinematic(q11, q22, q33, q44, q55, q66):
  print("Forward kinematic start ...")

  X = []
  Y = []
  Z = []

  # q1 = q1 * (math.pi / 180)
  q1 = math.radians(q11)
  print("q1: ", q1, "rad; ", q11, "°") 
  # q2 = q2 * (math.pi / 180)
  q2 = math.radians(q22)
  print("q2: ", q2, "rad; ", q22, "°") 
  #q3 = q3 * (math.pi / 180)
  q3 = math.radians(q33)
  print("q3: ", q3, "rad; ", q33, "°") 
  # q4 = q4 * (math.pi / 180)
  q4 = math.radians(q44)
  print("q4: ", q4, "rad; ", q44, "°") 
  # q5 = q5 * (math.pi / 180)
  q5 = math.radians(q55)
  print("q5: ", q5, "rad; ", q55, "°") 
  # q6 = q6 * (math.pi / 180)
  q6 = math.radians(q66)
  print("q6: ", q6, "rad; ", q66, "°") 
  print()

  # T01
  T01 = transformationMatrixModified(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  # print("T01: ", Matrix(T01))
  T0g = T01
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

   # T12
  T12 = transformationMatrixModified(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  # print("T12: ", Matrix(T12))
  T0g = T0g * T12
  px,py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T23
  T23 = transformationMatrixModified(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  # print("T23: ", Matrix(T23))
  T0g = T0g * T23
  px, py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T34
  T34 = transformationMatrixModified(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  # print("T34: ", Matrix(T34))
  T0g = T0g * T34
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T45
  T45 = transformationMatrixModified(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  # print("T45: ", Matrix(T45))
  T0g = T0g * T45
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T56
  T56 = transformationMatrixModified(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  # print("T56: ", Matrix(T56))
  T0g = T0g * T56
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)
  
  # T6g (final position and rotation)
  T6g = transformationMatrixModified(0 + robot_theta7, robot_alpha6, robot_a6, robot_d7)
  # print("T6g:")
  # printMatrixE(T6g)

  T0g = T0g * T6g
  print("T0g:")
  printMatrixE(T0g)
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)




  # Correction Needd to Account of Orientation Difference BEtween Deinition of
  # Gripper_link in URDF versus DH Convention
  # R_corr = simplify(R_z * R_y)
  #R_corr = Rz(math.pi) * Ry(-math.pi/2)
  #T0g = T0g * R_corr




  X = np.reshape(X, (1,7))
  Y = np.reshape(Y, (1,7))
  Z = np.reshape(Z, (1,7))



  # eps = 0.0001

  # bw = math.atan2(-T0g[0,2], math.sqrt(math.pow(T0g[0,0], 2) + math.pow(T0g[0,1], 2)))
  # aw = 0.0
  # cw = 0.0

  # if(abs(abs(bw) - math.pi / 2.0) > eps):
  #   aw = math.atan2(T0g[1,2], T0g[2,2])
  #   cw = math.atan2(T0g[0,1], T0g[0,0])
  # elif(bw > 0.0):
  #   aw = math.atan2(T0g[1,0], T0g[1,1])
  #   bw = math.pi / 2.0
  #   cw = 0.0
  # else:
  #   aw = -math.atan2(T0g[1,0], T0g[1,1])
  #   bw = -math.pi / 2.0
  #   cw = 0.0

  # if(aw < 0.0):
  #   aw = aw + 2.0 * math.pi


  cw = math.atan2(math.sqrt(math.pow(T0g[0,2], 2) + math.pow(T0g[1,2], 2)), T0g[2,2],)
  bw = math.atan2(T0g[2,1], -T0g[2,0])
  aw = math.atan2(T0g[1,2], T0g[0,2])




  R1 = eulerAnglesToRotationMatrix(np.array([aw, cw, bw]))
  print("R1:")
  printMatrixE(R1)




  print("x: ", X[0,5])
  print("y: ", Y[0,5])
  print("z: ", Z[0,5])
  print("Roll : ", aw, "rad; ", math.degrees(aw), "°")
  print("Pitch : ", cw, "rad; ", math.degrees(cw), "°")
  print("Yaw : ", bw, "rad; ", math.degrees(bw), "°")
  print("Forward kinematic end ...")
  print("########################################")
  return X, Y, Z, aw, cw, bw


def inversKinematic(x, y, z, roll, pitch, yaw):
  # input: given position and orientation of the gripper_URDF wrt base frame
  # output: angles q1, q2, q3, q4, q5, q6

  print("Invers kinematic start ...")
  print("Input: ")
  print("x: ", x)
  print("y: ", y)
  print("z: ", z)
  print("Roll : ", roll, "°; ", math.radians(roll), "rad")
  print("Pitch : ", pitch, "°; ", math.radians(pitch), "rad")
  print("Yaw : ", yaw, "°; ", math.radians(yaw), "rad")

  roll = math.radians(roll)
  pitch = math.radians(pitch)
  yaw = math.radians(yaw)
    
  gripper_point = x, y, z

  ################################################################################
  # All important symbolic transformations matrices are declared below 
  ################################################################################

  q1, q2, q3, q4, q5, q6 = symbols('q1:7')
  alpha, beta, gamma = symbols('alpha beta gamma', real = True)
  px, py, pz = symbols('px py pz', real = True)

  # X = []
  # Y = []
  # Z = []

  ############################## Transformation Matrix Start ####################################

  # T01
  T01 = transformationMatrixModifiedSymbol(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  print("T01:", T01)
  print()
  T0g = T01

  # T12
  T12 = transformationMatrixModifiedSymbol(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  print("T12:", T12)
  print()
  T0g = T0g * T12

  # T23
  T23 = transformationMatrixModifiedSymbol(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  print("T23:", T23)
  print()
  T0g = T0g * T23

  # T34
  T34 = transformationMatrixModifiedSymbol(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  print("T34:", T34)
  print()
  T0g = T0g * T34

  # T45
  T45 = transformationMatrixModifiedSymbol(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  print("T45:", T45)
  print()
  T0g = T0g * T45

  # T56
  T56 = transformationMatrixModifiedSymbol(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  print("T56:", T56)
  print()
  T0g = T0g * T56

  # T6g (final position and rotation)
  T6g = transformationMatrixModifiedSymbol(0 + robot_theta7, robot_alpha6, robot_a6, robot_d7)
  print("T6g:", T6g)
  print()
  T0g = T0g * T6g
  # print("T0g:", T0g)
  # print()

  ############################## Transformation Matrix End ####################################

  T03 = T01 * T12 * T23
  # print("T03:", simplify(T03))
  # print()

  T03T = T03.T
  # print("T03T:", simplify(T03T))
  # print()

  T13 = T12 * T23
  # print("T13:", simplify(T13))
  # print()

  T13T = T13.T
  # print("T13T:", simplify(T13T))
  # print()

  T36 = T34 * T45 * T56
  # print("T36:", simplify(T36))
  # print()

  T36T = T36.T
  # print("T36T:", simplify(T36T))
  # print()

  T06 = T01 * T12 * T23 * T34 * T45 * T56
  #print("T06: ", Matrix(T06))
  #print()

  T0g = T01 * T12 * T23 * T34 * T45 * T56 * T6g  
  # print("T0g: ", Matrix(T0g))
  # print()


  R1 = eulerAnglesToRotationMatrix(np.array([alpha, beta, gamma]))
  # print("R1")
  # pprint(simplify(R1))

  R0u = Matrix(R1)
  # print("R0u: ", R0u)
  # print()

  # Rotation of urdf_gripper wrt (DH) gripper frame from rotz(pi) * roty(-pi/2) and it's transpose
  # Rgu_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])
  # RguT_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])
  RguT_eval = Rz(math.pi)[:3, :3] * Ry(-math.pi/2)[:3, :3] 

  # Inverse kinematics transformations starts below
  print("Roll: ", roll, ", Pitch: %f", pitch, ", Yaw: %f", yaw)
  #R0u_eval = R0u.evalf(subs = {alpha: pitch, beta: roll, gamma: yaw})
  #R0u_eval = R0u.evalf(subs = {alpha: pitch, beta: yaw, gamma: roll})
  #R0u_eval = R0u.evalf(subs = {alpha: roll, beta: yaw, gamma: pitch})
  R0u_eval = R0u.evalf(subs = {alpha: roll, beta: pitch, gamma: yaw})
  #R0u_eval = R0u.evalf(subs = {alpha: yaw, beta: pitch, gamma: roll})
  #R0u_eval = R0u.evalf(subs = {alpha: yaw, beta: roll, gamma: pitch})
  print("R0u_eval 0 R1:")
  printMatrixE(R0u_eval)

  R0g_eval = R0u_eval
  #R0g_eval = R0u_eval * RguT_eval
  print("R0g_eval:")
  printMatrixE(R0g_eval)









  # wrist_center = get_wrist_center(gripper_point, R0g_eval, robot_d6, robot_d7)
  # print("Wrist center: ", wrist_center)
  # print()

  # j1, j2, j3 = get_first_three_angles(wrist_center)













  px = x - (robot_d6 + robot_d7) * R0g_eval[0, 2]
  py = y - (robot_d6 + robot_d7) * R0g_eval[1, 2]
  pz = z - (robot_d6 + robot_d7) * R0g_eval[2, 2]


  # Variante 1
  # Left arm configuration
  # r = math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(robot_d2, 2))
  r = math.sqrt(math.pow(px, 2) + math.pow(py, 2))
  s = pz - robot_d1

  tmp1 = math.atan2(px, py)
  tmp2 = math.atan2(math.sqrt(math.pow(r, 2) - math.pow(robot_d2, 2)), robot_d2)
  j1 = tmp1 - tmp2
  print("j1_1 OK~: ", j1)

  # Right arm configuration
  j1 = tmp1 - tmp2 + math.pi
  print("j1_2: ", j1)

 
  tmp1 = math.atan2(px, py)
  tmp2 = math.atan2(math.sqrt(math.pow(r, 2) - math.pow(robot_d2, 2)), robot_d2)
  tmp3 = tmp2 + math.pi
  j1 = tmp3 + tmp2
  print("j1_3 OK~: ", j1)

  # 16 S.105
  # tmp1 = math.atan2(px, py)
  # tmp2 = math.atan2(math.sqrt(math.pow(r, 2) - math.pow(robot_d2, 2)), robot_d2)
  # tmp3 = tmp2 + math.pi
  # j1 = tmp1 + tmp3
  # print("j1_4 OK~: ", j1)

  # tmp1 = math.atan2(px, py)
  # tmp2 = math.atan2(-math.sqrt(math.pow(r, 2) - math.pow(robot_d2, 2)), -robot_d2)
  # j1 = tmp1 + tmp2
  # print("j1_5 OK~: ", j1)



  # Ellbow configuration
  r = math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(robot_d2, 2))
  s = pz - robot_d1
  D = (math.pow(r, 2) + math.pow(s, 2) - math.pow(robot_a2, 2) - math.pow(robot_d4, 2)) / (2 * robot_a2 * robot_d4)

  j3 = math.atan2(-math.sqrt(1 - math.pow(D, 2)), D)
  print("j3_2 OK: ", j3)

  j3 = math.atan2(+math.sqrt(1 - math.pow(D, 2)), D)
  print("j3_1 OK: ", j3)
  





  r = math.sqrt(math.pow(px, 2) + math.pow(py, 2))
  r = math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(robot_d2, 2))
  s = pz
  #s = pz - robot_d1

  j2 = math.atan2(r, s) - math.atan2(robot_a2 + robot_d4 * cos(j3), robot_d4 * sin(j3))
  print("j2: ", j2)

  C = (math.pow(r, 2) + math.pow(s, 2) + math.pow(robot_a2, 2) - math.pow(robot_d4, 2)) / (2 * math.sqrt(math.pow(r, 2) + math.pow(s, 2)) * robot_a2)
  tmp1 = math.atan2(s, r)
  tmp2 = math.atan2(math.sqrt(1 - math.pow(C, 2)), C)
  j2 = 0.0
  # if s > 0:
  j2 = tmp2 - tmp1
  print("j2: ", j2)
  # if s < 0:
  j2 = tmp2 + tmp1
  print("j2: ", j2)






  # Variante 2
  # wrist_center = get_wrist_center(gripper_point, R0g_eval, robot_d6, robot_d7)
  # print("Wrist center: ", wrist_center)
  # print()
  xu, yu, zu = gripper_point 

  nx, ny, nz = R0g_eval[0, 2], R0g_eval[1, 2], R0g_eval[2, 2]
  # xw = xu - dg * nx
  # yw = yu - dg * ny
  # zw = zu - dg * nz 

  xw = xu - (robot_d6 + robot_d7) * nx
  yw = yu - (robot_d6 + robot_d7) * ny
  zw = zu - (robot_d6 + robot_d7) * nz


  # j1, j2, j3 = get_first_three_angles(wrist_center)
    
  l = get_hypotenuse(robot_d4, abs(robot_a3))
  phi = atan2(robot_d4, abs(robot_a3))
  
  x_prime = get_hypotenuse(xw, yw)
  mx = x_prime -  robot_a2
  my = zw - robot_d1 
  m = get_hypotenuse(mx, my)
  alpha = atan2(my, mx)
  
  gamma = get_cosine_law_angle(l, robot_a2, m)
  beta = get_cosine_law_angle(m, robot_a2, l)
  
  j1 = atan2(xw, yw)
  j2 = pi/2 - beta - alpha 
  j3 = -(gamma - phi)

  print()
  print()
  print("JJ1: ", j1.evalf())
  print("JJ2: ", j2.evalf())
  print("JJ3: ", j3.evalf())


















  # for testing purpose
  j1 = math.radians(89.29) # 1.558404489105737
  j2 = math.radians(-34.05) # -0.5942846103040692
  j3 = math.radians(125.42) # 2.1889919478512883

  # j1 = math.radians(90.0) # 1.5707963267948966
  # j2 = math.radians(0.0) # 0.0
  # j3 = math.radians(125.0) # 2.181661564992912
    
  R03T = (T03T).evalf(subs={q1: j1, q2: j2, q3: j3})

  Rrpy = np.dot(rot_z(roll), np.dot(rot_y(pitch), rot_z(yaw)))
  R_corr = simplify(rot_z(pi) * rot_y(-pi/2))
  R36 = R03T * Rrpy #* R_corr
  print("Rrpy")
  printMatrixE(Rrpy)
  print("R03T")
  printMatrixE(R03T)
  print("R36")
  printMatrixE(R36)


  # Variant 10  
  print()
  print()
  print("Variant 10")

  j5 = math.atan2(+math.sqrt(math.pow(R36[0,2], 2) + math.pow(R36[1,2], 2)), R36[2,2])
  print("j5_1: ", j5)
  j5 = math.atan2(-math.sqrt(math.pow(R36[0,2], 2) + math.pow(R36[1,2], 2)), R36[2,2])
  print("j5_2: ", j5)

  j4 = math.atan2(R36[2,1], -R36[2,0])
  print("j4_1: ", j4)
  j4 = math.atan2(-R36[2,1], -R36[2,0])
  print("j4_1: ", j4)

  j6 = math.atan2(R36[1,2], R36[0,2])
  print("j6_1: ", j6)
  j6 = math.atan2(-R36[1,2], R36[0,2])
  print("j6_2: ", j6)




  # Variant 11 
  print()
  print()
  print("Variant 11")

  # tf requires a numpy matrix instead of a sympy matrix
  R3_6_np = np.array(R36).astype(np.float64)

  alpha, beta, gamma = transformations.euler_from_matrix(R3_6_np, axes='rzyz')   # xyx, yzx, xyz
  theta4 = alpha
  print(theta4)
  theta5 = beta
  print(theta5)
  theta6 = gamma
  print(theta6)






  # Variant 12 
  print()
  print()
  print("Variant 12")

  tmp1 = math.sin(j1) * R36[0,2] - math.cos(j1) * R36[1,2]
  tmp2 = math.sqrt(1 - math.pow(math.sin(j1) * R36[0,2] - math.cos(j1) * R36[1,2], 2))

  j5 = math.atan2(tmp1, tmp2)
  print(j5)  
  j5 = math.atan2(tmp1, tmp2 * (-1))
  print(j5)
  





  # Variant 13
  print()
  print()
  print("Variant 13")    # j4, j5, j6 = get_last_three_angles(R36_eval)
  sin_q4 = R36[2, 2]
  cos_q4 =  -R36[0, 2]
    
  sin_q5 = sqrt(R36[0, 2]**2 + R36[2, 2]**2) 
  cos_q5 = R36[1, 2]
    
  sin_q6 = -R36[1, 1]
  cos_q6 = R36[1, 0] 
  
  j4 = atan2(sin_q4, cos_q4) * (-1)
  print("q4", j4.evalf())
  j5 = atan2(sin_q5, cos_q5) - math.pi
  print("q5", j5.evalf())
  j6 = atan2(sin_q6, cos_q6) + math.pi
  print("q6", j6.evalf())
  


  # # Variant 14  
  # print()
  # print()
  # print("Variant 14")

  # j5 = +math.atan2(math.sqrt(1 - math.pow(R36[2,2], 2)), R36[2,2])
  # print("j5_1: ", j5)
  # j5 = -math.atan2(math.sqrt(1 - math.pow(R36[2,2], 2)), R36[2,2])
  # print("j5_2: ", j5)

  # j4 = math.atan2(R36[2,1], -R36[2,0])
  # print("j4_1: ", j4)

  # j6 = math.atan2(R36[2,1], R36[2,0])
  # print("j6_1: ", j6)










  print()
  print("Result:")
  print("q1 +/-:", j1, "rad; ", math.degrees(j1), "°")
  print("q2:", j2, "rad; ", math.degrees(j2), "°")
  print("q3 +/-:", j3, "rad; ", math.degrees(j3), "°")
  print("q4:", j4.evalf(), "rad; ", math.degrees(j4.evalf()), "°")
  print("q5 +/-:", j5.evalf(), "rad; ", math.degrees(j5.evalf()), "°")
  print("q6:", j6.evalf(), "rad; ", math.degrees(j6.evalf()), "°")

  print ("89.29, -34.05, 125.42, -21.60, -1.60, 21.63")
  print("Invers kinematic end ...")
  print("########################################")

  return math.degrees(j1), math.degrees(j2), math.degrees(j3), math.degrees(j4), math.degrees(j5), math.degrees(j6)
















# rotation matrices in x axes
def rotx(q):
  sq, cq = sin(q), cos(q)

  r = Matrix([
    [1., 0., 0., 0.],
    [0., cq,-sq, 0.],
    [0., sq, cq, 0.],
    [0., 0. ,0., 1.]])

  r = Matrix([
    [1.,  0.,  0.,  0.],
    [0.,  1.,  0.,  0.],
    [0.,  cq, -sq,  0.],
    [0.,  sq,  cq,  1.]])
    
  return r

# # rotation matrices in y axes
# def roty(q):
#   sq, cq = sin(q), cos(q)

#   r = Matrix([
#     [ cq, 0., sq, 0.],
#     [ 0., 1., 0., 0.],
#     [-sq, 0., cq, 0.],
#     [0., 0. ,0., 1.]])
    
#   return r

# rotation matrices in z axes
def rotz(q):
  sq, cq = sin(q), cos(q)

  r = Matrix([
    [1.,  0.,  0.,  0.],
    [0.,  cq, -sq,  0.],
    [0.,  sq,  cq,  0.],
    [0.,  0.  ,0.,  1.]])
               
  return r


def inversKinematic1(x, y, z, roll, pitch, yaw):
  # input: given position and orientation of the gripper_URDF wrt base frame
  # output: angles q1, q2, q3, q4, q5, q6

  print("Invers kinematic start ...")
  print("Input: ", x, " ", y, " ", z, " ", roll, " ", pitch, " ", yaw)

  roll = math.radians(roll)
  pitch = math.radians(pitch)
  yaw = math.radians(yaw)
    
  gripper_point = x, y, z

  ################################################################################
  # All important symbolic transformations matrices are declared below 
  ################################################################################

  q1, q2, q3, q4, q5, q6 = symbols('q1:7')
  alpha, beta, gamma = symbols('alpha beta gamma', real = True)
  px, py, pz = symbols('px py pz', real = True)

  # X = []
  # Y = []
  # Z = []

  ############################## Transformation Matrix Start ####################################

  # T01
  T01 = transformationMatrixModifiedNewSymbol(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  print("T01:")
  pprint(T01)
  print()
  T0g = T01
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

   # T12
  T12 = transformationMatrixModifiedNewSymbol(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  print("T12:")
  pprint(simplify(T12 * rotx(-math.pi/2) * rotz(-math.pi/2)))
  print()
  T0g = T0g * T12
  # px,py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T23
  T23 = transformationMatrixModifiedNewSymbol(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  print("T23:")
  pprint(T23)
  print()
  T0g = T0g * T23
  # px, py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T34
  T34 = transformationMatrixModifiedNewSymbol(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  print("T34:")
  pprint(T34)
  print()
  T0g = T0g * T34
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T45
  T45 = transformationMatrixModifiedNewSymbol(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  print("T45:")
  pprint(T45)
  print()
  T0g = T0g * T45
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T56
  T56 = transformationMatrixModifiedNewSymbol(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  print("T56:")
  pprint(T56)
  print()
  T0g = T0g * T56
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T6g (final position and rotation)
  T6g = transformationMatrixModifiedNewSymbol(0 + robot_theta7, robot_alpha6, robot_a6, robot_d7)
  print("T6g:")
  pprint(T6g)
  print()
  T0g = T0g * T6g
  # print("T0g:", T0g)
  # print()
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  ############################## Transformation Matrix End ####################################

  T03 = T01 * T12 * T23
  # print("T03:", simplify(T03))
  # print()

  T03T = T03.T
  # print("T03T:", simplify(T03T))
  # print()

  T13 = T12 * T23
  # print("T13:", simplify(T13))
  # print()

  T13T = T13.T
  # print("T13T:", simplify(T13T))
  # print()

  T36 = T34 * T45 * T56
  # print("T36:", simplify(T36))
  # print()

  T36T = T36.T
  # print("T36T:", simplify(T36T))
  # print()

  T06 = T01 * T12 * T23 * T34 * T45 * T56
  # print("T06:")
  # pprint(Matrix(T06))
  # print()

  T0g = T01 * T12 * T23 * T34 * T45 * T56 * T6g  
  # print("T0g: ", Matrix(T0g))
  # print()







  px = T06[1,1]
  print("px:", simplify(px))
  print()

  py = T06[1,2]
  print("py:", simplify(py))
  print()

  pz = T06[1,3]
  print("pz:", simplify(pz))
  print()


  r = math.sqrt(math.pow(px, 2) + math.pow(py, 2))












  R1 = eulerAnglesToRotationMatrix(np.array([gamma, beta, alpha]))
  # print("R1:", simplify(R1))
  # print()














  R0u = Matrix(R1)
  # print("R0u: ", R0u)
  # print()

  R03T = Matrix(T03T[:3, :3])
  # print("R03T: ", R03T)
  # print()

  # Rotation of urdf_gripper wrt (DH) gripper frame from rotz(pi) * roty(-pi/2) and it's transpose
  # Rgu_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])
  # RguT_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])
  RguT_eval = Rz(math.pi)[:3, :3] * Ry(-math.pi/2)[:3, :3] 

  # Inverse kinematics transformations starts below
  #R0u_eval = R0u.evalf(subs = {alpha: pitch, beta: roll, gamma: yaw})
  #R0u_eval = R0u.evalf(subs = {alpha: pitch, beta: yaw, gamma: roll})
  #R0u_eval = R0u.evalf(subs = {alpha: roll, beta: yaw, gamma: pitch})
  R0u_eval = R0u.evalf(subs = {alpha: roll, beta: pitch, gamma: yaw})
  #R0u_eval = R0u.evalf(subs = {alpha: pitch, beta: pitch, gamma: roll})
  #R0u_eval = R0u.evalf(subs = {alpha: pitch, beta: roll, gamma: pitch})
  #R0g_eval = R0u_eval
  R0g_eval = R0u_eval * RguT_eval
  # print("R0g_eval: ", R0u_eval)
  # print()

  wrist_center = get_wrist_center(gripper_point, R0g_eval, robot_d6, robot_d7)
  print("Wrist center: ", wrist_center)
  print()

  j1, j2, j3 = get_first_three_angles(wrist_center)

  R03T_eval = R03T.evalf(subs = {q1: j1.evalf(), q2: j2.evalf(), q3: j3.evalf()})
  R36_eval = R03T_eval * R0g_eval
  print("R03T_eval: ", R03T_eval)
  print()

  j4, j5, j6 = get_last_three_angles(R36_eval)

  j1 = j1.evalf()
  j2 = j2.evalf()
  j3 = j3.evalf()
  j4 = j4.evalf()
  j5 = j5.evalf()
  j6 = j6.evalf()
  
  print("q1 +/-:", j1, "rad; ", math.degrees(j1), "°")
  print("q2:", j2, "rad; ", math.degrees(j2), "°")
  print("q3 +/-:", j3, "rad; ", math.degrees(j3), "°")
  print("q4:", j4, "rad; ", math.degrees(j4), "°")
  print("q5 +/-:", j5, "rad; ", math.degrees(j5), "°")
  print("q6:", j6, "rad; ", math.degrees(j6), "°")

  print ("89.29, -34.05, 125.42, -21.60, -1.60, 21.63")
  print("Invers kinematic end ...")
  print("########################################")

  return j1, j2, j3, j4, j5, j6





def main():
  np.set_printoptions(infstr="(infinity)")
  np.set_printoptions(formatter={'float': "\t{: 0.0f}\t".format})

  #px, py, pz = 0.49792, 1.3673, 3.0000
  #px, py, pz = 0.276668856049955, -0.600333896425126, 0.512767653796484
  #px, py, pz = 0.3253, -0.23290, 0.55470
  #px, py, pz = 0.49792, 1.3673, 2.4988
  #roll, pitch, yaw = 1, -1, -1
  #roll, pitch, yaw = 2.58309, 1.60570, 2.09440
  #roll, pitch, yaw = 0.366, -0.078, 2.561
  #roll, pitch, yaw = 0.17453, 1.57080, 0.00000

  # Calculation 1
  x, y, z = -0.00423, 0.25884, 0.46696
  roll, pitch, yaw = 89.88, 89.89, 0.03
  #
  x = -0.00428667791374644
  y = 0.258835606840153
  z =  0.466983876304506
  roll =  89.87893434791008
  pitch =  89.88231144120364
  yaw =  0.03643662760103062

  # Calculation 1
  # x, y, z = -0.0062, 0.42061, 0.355930
  # roll, pitch, yaw = 90.00, 80.00, 0.00
  # #
  # x = -0.00624999999999999
  # y = 0.420609690412826
  # z =  0.355935300920514
  # roll =  90.0
  # pitch =  80.00000000000001
  # yaw =  2.9438590922394028e-15


  # print("get_angles ...")
  q1, q2, q3, q4, q5, q6 = inversKinematic(x, y, z, roll, pitch, yaw)
  #forwardKinematic(q1, q2, q3, q4, q5, q6)
  #print()
  #print()

  # forwardKinematic(89.29, -34.05, 125.42, -21.60, -1.60, 21.63)        # Result should be: x=-4.23, y=258.84, z=466.96, a=89.88, b=89.89, c=0.03
  # print("----------------------------------------")
  # forwardKinematic(93.76, -54.76, 118.97, -4.89, 30.30, 1.64)          # Result should be: x=-8.60, y=115.87, z=522.94, a=91.28, b=94.42, c=-2.78
  # print("----------------------------------------")
  # forwardKinematic(71.97, -1.45, 106.23, -9.10, -39.35, -23.230)       # Result should be: x=123.430, y=438.550, z=501.370, a=78.28, b=65.92, c=-32.87
  # print("----------------------------------------")
  # forwardKinematic(90.00, 0.0, 125.00, 0.0, -45.00, -0.0)              # Result should be: x=-6.2, y=420.610, z=355.930, a=90.00, b=80.00, c=0.0
  # print("----------------------------------------")
  # forwardKinematic(91.72, -85.60, 85.02, -1.72, 90.59, -0.02)          # Result should be: x=6.2, y=-175.00, z=600.00 a=90.00, b=90.01, c=0.00 (Eneeffector wrong)

  # # float[] _toolCenterPositionKlein = { -52.750f, 0.000f, 118.10f, 0.000f, 0.000f, 0.000f };        // NiederstreiferKlein z, y , y, a, b, c
  # # float[] _toolCenterPositionGross = { -65.500f, 0.000f, 117.500f, 0.000f, 0.000f, 0.000f };       // NiederstreiferGross z, y , y, a, b, c

if __name__=="__main__":
  main()