from calendar import c
import math
from re import A
import rlcompleter
from symtable import Symbol
from sympy import symbols, simplify, pprint,  expand_trig, trigsimp, pi, atan2, cos, sin, tan, sqrt, atan, atanh
from sympy.matrices import Matrix

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
import numpy as np

global robot_d1, robot_d2, robot_d3, robot_d4, robot_d5, robot_d6, robot_dg
global robot_a0, robot_a1, robot_a2, robot_a3, robot_a4, robot_a5, robot_ag
global robot_alpha0, robot_alpha1, robot_alpha2, robot_alpha3, robot_alpha4, robot_alpha5, robot_alphag
global robot_theta1, robot_theta2, robot_theta3, robot_theta4, robot_theta5, robot_theta6, robot_thetag

# Doosan
robot_d1 = 0.135
robot_d2 = 0.00625
robot_d3 = 0.0
robot_d4 = 0.368
robot_d5 = 0.0
robot_d6 = 0.121
robot_dg = 0.0

robot_a0 = 0.0
robot_a1 = 0.0
robot_a2 = 0.411
robot_a3 = 0.0
robot_a4 = 0.0
robot_a5 = 0.0
robot_ag = 0.0

robot_alpha0 = 0.0
robot_alpha1 = -pi/2
robot_alpha2 = 0.0
robot_alpha3 = pi/2
robot_alpha4 = -pi/2
robot_alpha5 = pi/2
robot_alphag = 0.0

robot_theta1 = 0.0
robot_theta2 = -pi/2
robot_theta3 = pi/2
robot_theta4 = 0.0
robot_theta5 = 0.0
robot_theta6 = 0.0
robot_thetag = 0.0

# 6-DOF LinuxCNC
robot_d1 = 0.450
robot_d2 = 0.11
robot_d3 = 0.0
robot_d4 = 0.620
robot_d5 = 0.0
robot_d6 = 0.0
robot_dg = 0.0

robot_a0 = 0.2
robot_a1 = 0.6
robot_a2 = 0.110
robot_a3 = 0.0
robot_a4 = 0.0
robot_a5 = 0.0
robot_ag = 0.0

robot_alpha0 = pi/2
robot_alpha1 = 0.0
robot_alpha2 = pi/2
robot_alpha3 = -pi/2
robot_alpha4 = pi/2
robot_alpha5 = 0.0
robot_alphag = 0.0

robot_theta1 = 0.0
robot_theta2 = 0.0
robot_theta3 = 0.0
robot_theta4 = 0.0
robot_theta5 = 0.0
robot_theta6 = 0.0
robot_thetag = 0.0







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
        print("%6.6f\t" %a[i,j], end = '')
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
  r24 = d * sin(alpha)

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








# Calculates Rotation Matrix given euler angles (https://mathworld.wolfram.com/RotationMatrix.html, https://mathworld.wolfram.com/EulerAngles.html, https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix)
def eulerAnglesToRotationMatrix(angle):
    R_x = np.array([[1,               0,              	  0               ],
                    [0,               cos(angle[0]),      sin(angle[0])   ],
                    [0,               -sin(angle[0]),     cos(angle[0])   ]])

    R_y = np.array([[cos(angle[1]),    0,                 -sin(angle[1])  ],
                    [0,                1,                 0               ],
                    [sin(angle[1]),    0,                 cos(angle[1])    ]])

    R_z = np.array([[cos(angle[2]),    sin(angle[2]),     0               ],
                    [-sin(angle[2]),   cos(angle[2]),     0               ],
                    [0,                0,                 1               ]])

    R = np.dot(R_z, np.dot( R_y, R_x ))

    return R

# Calculates Rotation Matrix given euler angles (https://mathworld.wolfram.com/RotationMatrix.html, https://mathworld.wolfram.com/EulerAngles.html, https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix)
def eulerAnglesToRotationMatrixT(angle):
    R_x = np.array([[1,               0,              	  0               ],
                    [0,               cos(angle[0]),      -sin(angle[0])  ],
                    [0,               sin(angle[0]),     cos(angle[0])    ]])

    R_y = np.array([[cos(angle[1]),    0,                 sin(angle[1])   ],
                    [0,                1,                 0               ],
                    [-sin(angle[1]),   0,                 cos(angle[1])   ]])

    R_z = np.array([[cos(angle[2]),    -sin(angle[2]),    0               ],
                    [sin(angle[2]),    cos(angle[2]),     0               ],
                    [0,                0,                 1               ]])

    R = np.dot(R_z, np.dot( R_y, R_x ))

    return R




# def poseMat(x, y, z, a, b, c):
#   sr =




# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R):
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6

# Calculates rotation matrix to euler angles
# the result is the same as MATLAB except the order
# of the euler angles ( x and z are swapped ).
def rotationMatrixToEulerAngles(R):
    assert(isRotationMatrix(R))

    sy = sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])

    singular = sy < 1e-6

    if  not singular :
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

def get_wrist_center(gripper_point, R0g, dg):
  # get the coordinates of the wrist center wrt to the base frame (xw, yw, zw)
  # given the following info:
  # the coordinates of the gripper (end effector) (x, y, z)
  # the rotation of the gripper in gripper frame wrt to the base frame (R0u)
  # the distance between gripper and wrist center dg which is along common z axis
  # check WRITEUP.pdf for more info
  xu, yu, zu = gripper_point

  nx, ny, nz = R0g[0, 2], R0g[1, 2], R0g[2, 2]
  xw = xu - dg * nx
  yw = yu - dg * ny
  zw = zu - dg * nz

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






def forwardTransformation(q11, q22, q33, q44, q55, q66):
  X = []
  Y = []
  Z = []

  # q1 = q1 * (math.pi / 180)
  q1 = math.radians(q11)
  print("q1: ", q1, "rad; ", q11, "°")
  # q2 = q2 * (math.pi / 180)
  q2 = math.radians(q22)
  print("q1: ", q2, "rad; ", q22, "°")
  #q3 = q3 * (math.pi / 180)
  q3 = math.radians(q33)
  print("q1: ", q3, "rad; ", q33, "°")
  # q4 = q4 * (math.pi / 180)
  q4 = math.radians(q44)
  print("q1: ", q4, "rad; ", q44, "°")
  # q5 = q5 * (math.pi / 180)
  q5 = math.radians(q55)
  print("q1: ", q5, "rad; ", q55, "°")
  # q6 = q6 * (math.pi / 180)
  q6 = math.radians(q66)
  print("q1: ", q6, "rad; ", q66, "°")
  print()

  # T01
  T01 = transformationMatrixModified(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  print("T01:")
  printMatrixE(T01)
  T0g = T01
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

   # T12
  T12 = transformationMatrixModified(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  print("T12:")
  printMatrixE(T12)
  T0g = T0g* T12
  px,py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T23
  T23 = transformationMatrixModified(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  print("T23:")
  printMatrixE(T23)
  T0g = T0g* T23
  px, py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T34
  T34 = transformationMatrixModified(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  print("T34:")
  printMatrixE(T34)
  T0g = T0g* T34
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T45
  T45 = transformationMatrixModified(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  print("T45:")
  printMatrixE(T45)
  T0g = T0g* T45
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T56
  T56 = transformationMatrixModified(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  print("T56:")
  printMatrixE(T56)
  T0g = T0g* T56
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T6g (final position and rotation)
  T6g = transformationMatrixModified(0 + robot_thetag, robot_alphag, robot_ag, robot_dg)
  print("T6g:")
  printMatrixE(T6g)

  T0g = T0g* T6g
  print("T0g:")
  printMatrixE(T0g)
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # print("T0g[0,0]: ", T0g[0,0])
  # print("T0g[1,0]: ", T0g[1,0])
  # print("T0g[2,0]: ", T0g[2,0])

  # print("T0g[2,1]: ", T0g[2,1])
  # print("T0g[2,2]: ", T0g[2,2])

  #fig = plt.figure()
  #ax = fig.add_subplot(111,projection = '3d')
  #ax.set_xlabel('x axis')
  #ax.set_ylabel('y axis')
  #ax.set_zlabel('z axis')

  X = np.reshape(X, (1,7))
  Y = np.reshape(Y, (1,7))
  Z = np.reshape(Z, (1,7))

  #ax.cla()
  #ax.plot_wireframe(X,Y,Z)
  #plt.draw()
  #plt.pause(3)
  #ax.cla()
  #ax.plot_wireframe(Z,Y,X,color='r')

  #plt.show()

  # T13 = T12 * T23
  # print("T13: ", Matrix(T13))

  # W6 = T01 * T12 * T23 * T34 * T45 * T56
  # print("W6: ", Matrix(W6))
  # T0g = T01 * T12 * T23 * T34 * T45 * T56 * T6g
  # print("T0g: ", Matrix(T0g))

  # U6 = T56
  # U5 = T45 * U6
  # U4 = T34 * U5
  # U3 = T23 * U4
  # U2 = T13 * U4
  # U1 = T01 * U2

  # print(T0g[2,0])
  # print(math.pow(T0g[0,0], 2))
  # print(math.pow(T0g[1,0], 2))
  # print(math.sqrt(math.pow(T0g[1,0], 2) + math.pow(T0g[0,0], 2)))

  pitch = math.atan2(-T0g[2,0], math.sqrt(math.pow(T0g[1,0], 2) + math.pow(T0g[0,0], 2)))
  print("Pitch1: ", pitch, " ", pitch * 180 / math.pi)
  print("Pitch1: ", math.pi - pitch, " ", pitch * 180 / math.pi)

  # pitch = math.pi + math.atan(sqrt(pow(T0g[1,0], 2) + pow(T0g[0,0], 2)) / -T0g[2,0])
  # print("Pitch2: ", pitch, " ", pitch * 180 / math.pi)
  # print("Pitch2: ", pitch, " ", (math.pi + pitch) * 180 / math.pi)
  # print("Pitch2: ", pitch, " ", math.degrees(math.pi + pitch))

  roll = 0.0
  if (pitch != math.pi/2 and pitch != -(math.pi/2)):
    roll = math.pi + math.atan(T0g[0,0]/math.cos(pitch) / (T0g[1,0]/math.cos(pitch)))

  yaw = 0.0
  if (pitch == math.pi/2):
    yaw = math.pi + math.atan(T0g[2,2] / T0g[3,2])
  elif (pitch == -math.pi/2):
    yaw = math.pi + -math.atan(T0g[2,2] / T0g[3,2])
  else:
    yaw = math.pi + math.atan(((T0g[2,2]/math.cos(pitch)) / (T0g[2,1]/math.cos(pitch))))

  print("Forward kinematic:")
  print("T0g-X : ", Matrix(X))
  print("T0g-Y : ", Matrix(Y))
  print("T0g-Z : ", Matrix(Z))
  print()
  print("x: ", T0g[0,3])
  print("y: ", T0g[1,3])
  print("z: ", T0g[2,3])
  print("Pitch : ", pitch, "rad; ", math.degrees(pitch), "°")
  print("Roll : ", roll, "rad; ", math.degrees(roll), "°")
  print("Yaw : ", yaw, "rad; ", math.degrees(yaw), "°")
  return X, Y, Z


def inverseTransformation(x, y, z, roll, pitch, yaw):
  # input: given position and orientation of the gripper_URDF wrt base frame
  # output: angles q1, q2, q3, q4, q5, q6

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
  T01 = transformationMatrixModifiedSymbol(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  print("T01:", T01)
  print()
  T0g = T01
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

   # T12
  T12 = transformationMatrixModifiedSymbol(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  print("T12:", T12)
  print()
  T0g = T0g* T12
  # px,py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T23
  T23 = transformationMatrixModifiedSymbol(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  print("T23:", T23)
  print()
  T0g = T0g* T23
  # px, py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T34
  T34 = transformationMatrixModifiedSymbol(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  print("T34:", T34)
  print()
  T0g = T0g* T34
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T45
  T45 = transformationMatrixModifiedSymbol(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  print("T45:", T45)
  print()
  T0g = T0g* T45
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T56
  T56 = transformationMatrixModifiedSymbol(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  print("T56:", T56)
  print()
  T0g = T0g* T56
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T6g (final position and rotation)
  T6g = transformationMatrixModifiedSymbol(0 + robot_thetag, robot_alphag, robot_ag, robot_dg)
  print("T6g:", T6g)
  print()
  T0g = T0g* T6g
  print("T0g:", T0g)
  print()
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  ############################## Transformation Matrix End ####################################

  T03 = T01 * T12 * T23
  print("T03:", simplify(T03))
  print()

  T03T = T03.T
  print("T03T:", simplify(T03T))
  print()

  T13 = T12 * T23
  print("T13:", simplify(T13))
  print()

  T13T = T13.T
  print("T03T:", simplify(T13T))
  print()

  T36 = T34 * T45 * T56
  print("T36:", simplify(T36))
  print()

  T36T = T36.T
  print("T36T:", simplify(T36T))
  print()

  T06 = T01 * T12 * T23 * T34 * T45 * T56
  print("T06: ", Matrix(T06))
  print()

  T0g = T01 * T12 * T23 * T34 * T45 * T56 * T6g
  print("T0g: ", Matrix(T0g))
  print()




  R1 = eulerAnglesToRotationMatrixT(np.array([gamma, beta, alpha]))
  print("R1:", simplify(R1))
  print()




  R0u = Matrix(R1)
  print("R0u: ", R0u)
  print()

  R03T = Matrix(T03T[:3, :3])
  print("R03T: ", R03T)
  print()

  # Rotation of urdf_gripper wrt (DH) gripper frame from rotz(pi) * roty(-pi/2) and it's transpose
  Rgu_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])
  RguT_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])

  # Inverse kinematics transformations starts below
  R0u_eval = R0u.evalf(subs = {alpha: yaw, beta: pitch, gamma: roll })
  R0g_eval = R0u_eval * RguT_eval
  print("R0g_eval: ", R0u_eval)
  print()

  wrist_center = get_wrist_center(gripper_point, R0g_eval, robot_dg)
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

  print("Invers kinematic:")
  print("q1 +/-:", j1, "rad; ", math.degrees(j1), "°")
  print("q2:", j2, "rad; ", math.degrees(j2), "°")
  print("q3 +/-:", j3, "rad; ", math.degrees(j3), "°")
  print("q4:", j4, "rad; ", math.degrees(j4), "°")
  print("q5 +/-:", j5, "rad; ", math.degrees(j5), "°")
  print("q6:", j6, "rad; ", math.degrees(j6), "°")

  return j1, j2, j3, j4, j5, j6














def PoseMat(Xt, Yt, Zt, At, Bt, Ct):
  sr = math.sin(Xt)
  cr = math.cos(Yt)
  sp = math.sin(Zt)
  cp = math.cos(At)
  sy = math.sin(Bt)
  cy = math.cos(Ct)

  t00 = cp*cy
  t10 = cp*sy
  t20 = -sp
  t30 = 0.0

  t01 = sr*sp*cy - cr*sy
  t11 = sr*sp*sy + cr*cy
  t21 = sr*cp
  t31 = 0.0

  t02 = cr*sp*cy + sr*sy
  t12 = cr*sp*sy - sr*cy
  t22 = cr*cp
  t32 = 0.0

  t03 = Xt
  t13 = Yt
  t23 = Zt
  t33 = 1.0

  T = Matrix([[t00, t01, t02, t03],
              [t10, t11, t12, t13],
              [t20, t21, t22, t23],
              [t30, t31, t32, t33]])

  print("T")
  pprint(T)

  return T






def inverseTransformation1(x, y, z, roll, pitch, yaw):
  # input: given position and orientation of the gripper_URDF wrt base frame
  # output: angles q1, q2, q3, q4, q5, q6

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


  # robot_d1, robot_d2, robot_d3, robot_d4, robot_d5, robot_d6, robot_dg = symbols('d1 d2 d3 d4 d5 d6 dg', real = True)
  # robot_a0, robot_a1, robot_a2, robot_a3, robot_a4, robot_a5, robot_ag = symbols('a0 a1 a2 a3 a4 a5 ag', real = True)
  # robot_alpha0, robot_alpha1, robot_alpha2, robot_alpha3, robot_alpha4, robot_alpha5, robot_alphag = symbols('aplha0 aplha1 aplha2 aplha3 aplha4 aplha5 aplhag', real = True)
  # robot_theta1, robot_theta2, robot_theta3, robot_theta4, robot_theta5, robot_theta6, robot_thetag = symbols('theta1 theta2 theta3 theta4 theta5 theta6 thetag', real = True)

  # X = []
  # Y = []
  # Z = []

  ############################## Transformation Matrix Start ####################################

  # T01
  T01 = transformationMatrixClassicSymbol(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  print("T01:")
  pprint(T01)
  print("T01T:")
  pprint(T01.T)
  print()
  T0g = T01
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

   # T12
  T12 = transformationMatrixClassicSymbol(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  print("T12:")
  pprint(T12)
  print()
  T0g = T0g* T12
  # px,py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T23
  T23 = transformationMatrixClassicSymbol(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  print("T23:")
  pprint(T23)
  print()
  T0g = T0g* T23
  # px, py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T34
  T34 = transformationMatrixClassicSymbol(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  print("T34:")
  pprint(T34)
  print()
  T0g = T0g* T34
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T45
  T45 = transformationMatrixClassicSymbol(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  print("T45:")
  pprint(T45)
  print()
  T0g = T0g* T45
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T56
  T56 = transformationMatrixClassicSymbol(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  print("T56:")
  pprint(T56)
  print()
  T0g = T0g* T56
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  # T6g (final position and rotation)
  T6g = transformationMatrixClassicSymbol(0 + robot_thetag, robot_alphag, robot_ag, robot_dg)
  print("T6g:")
  pprint(T6g)
  print()
  T0g = T0g* T6g
  print("T0g:")
  pprint(T0g)
  print()
  # px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  # X.append(px)
  # Y.append(py)
  # Z.append(pz)

  ############################## Transformation Matrix End ####################################

  T03 = T01 * T12 * T23
  print("T03:")
  pprint(simplify(T03))
  print()

  T03T = T03.T
  print("T03T:")
  pprint(simplify(T03T))
  print()

  T13 = T12 * T23
  print("T13:")
  pprint(simplify(T13))
  print()

  T13T = T13.T
  print("T13T:")
  pprint(simplify(T13T))
  print()

  T36 = T34 * T45 * T56
  print("T36:")
  pprint(simplify(T36))
  print()

  T36T = T36.T
  print("T36T:")
  pprint(simplify(T36T))
  print()

  T06 = T01 * T12 * T23 * T34 * T45 * T56
  print("T06:")
  #pprint(simplify(T06))
  print()

  T0g = T01 * T12 * T23 * T34 * T45 * T56 * T6g
  print("T0g:")
  #pprint(simplify(T0g))
  print()










  W6 = T06 # T0g ????
  print("W6:")
  #pprint(simplify(W6))
  print()

  # U6 = simplify(T56)
  # U5 = simplify(T45 * U6)
  # U4 = simplify(T34 * U5)
  # U3 = simplify(T23 * U4)
  # U2 = simplify(T12 * U3)
  # U1 = simplify(T01 * U2)
  U6 = T56
  print("U6:")
  pprint(simplify(U6))
  print()
  U5 = T45 * U6
  print("U5:")
  pprint(simplify(U5))
  print()
  U4 = T34 * U5
  print("U4:")
  pprint(simplify(U4))
  print()
  U3 = T23 * U4
  print("U3:")
  pprint(simplify(U3))
  print()
  U2 = T12 * U3
  print("U2:")
  #pprint(simplify(U2))
  print()
  U1 = T01 * U2
  print("U1:")
  #pprint(simplify(U1))
  print()

  # W6 = V0 = U1
  # print("W6:")
  # #pprint(simplify(W6))
  # print()


  X, Y, Z, M = symbols('X Y Z M', real = True)

  V0j = Matrix([[X], [Y], [Z], [M]])
  W6 = V0j
  V1 = T01.T * W6
  print("V1:")
  #pprint(simplify(V1))
  print()
  V2 = T12.T * V1
  print("V2:")
  #pprint(simplify(V2))
  print()
  V3 = T13.T * V1
  print("V3:")
  #pprint(simplify(V3))
  print()
  V4 = T34.T * V3
  print("V4:")
  #pprint(simplify(V4))
  print()
  V5 = T45.T * V4
  print("V5:")
  #pprint(simplify(V5))
  print()




  Xt = x
  Yt = y
  Zt = z
  At = roll
  Bt = pitch
  Ct = yaw

  n1 = 1 # shoulder on the right
  n2 = 1 # elbow up
  n4 = 1 # wrist not flipped

  elbup = 1
  wrflp = 0
  wsing = 0

  eps = 0.0001



  Tw = PoseMat(Xt, Yt, Zt, At, Bt, Ct)


  ux = Tw[0,0]
  vx = Tw[0,1]
  wx = Tw[0,2]
  qx = Tw[0,3]

  uy = Tw[1,0]
  vy = Tw[1,1]
  wy = Tw[1,2]
  qy = Tw[1,3]

  uz = Tw[2,0]
  vz = Tw[2,1]
  wz = Tw[2,2]
  qz = Tw[2,3]


  # wrist position --------------------------------------------------
  px = qx - robot_d6*wx
  py = qy - robot_d6*wy
  pz = qz - robot_d6*wz

  # solve for th1 ---------------------------------------------------
  r = math.sqrt(px*px + py*py)
  if (r < robot_d2):
    # 'ERROR:--------- point not reachable'
    return 1

  k1 = math.atan2(py, px)
  k2 = math.asin(robot_d2 / r)

  if (n1 == 1):
    th1 = k1 + k2
  else:
    th1 = k1 - k2 + math.pi

  c1 = math.cos(th1)
  s1 = math.sin(th1)

  # solve for th2 ---------------------------------------------------
  print("U2:")
  pprint(U2)
  # v114 = px*c1 + py*s1 - robot_a1
  v114 = U2[0,3]
  print (v114)
  #v124 = pz - robot_d1
  v124 = U2[1,3]
  print (v124)
  r = math.sqrt(v114*v114 + v124*v124)
  k1 = (robot_a1*robot_a1 - robot_d4*robot_d4 - robot_a2*robot_a2 + v114*v114 + v124*v124) / (2*robot_a1*r)

  if (abs(k1) > 1): 
    # 'ERROR:--------- point not reachable'
   return 2
  
  k2 = math.acos(k1)

  if (elbup == 1):
    n2 = 1
  else:
    n2 = -1

  th2 = math.atan2(v124, v114) + n2*k2
  c2 = math.cos(th2)
  s2 = math.sin(th2)

  # solve for th3 -----------------------------------------------
  v214 = c2*v114 + s2*v124 - robot_a1
  v224 = -s2*v114 + c2*v124
  th3 = -math.atan2(robot_a2, robot_d4) + math.atan2(v214, -v224)
  c3 = math.cos(th3)
  s3 = math.sin(th3)

  # solve for th4 -----------------------------------------------
  c23 = math.cos(th2+th3)
  s23 = math.sin(th2+th3)
  v113 = c1*wx + s1*wy
  v313 = c23*v113 + s23*wz
  v323 = s1*wx - c1*wy
  if ((abs(v323) < eps) and (abs(v313) < eps)):
    th4 = 0
  else:
    th4 = math.atan2(n4*v323, n4*v313)

  # take care of singularities and map for continuity
  if ((abs(v323) < eps) and (v313 < eps)):
    th4 = th4old
  if ((v323 > eps) and (v313 < eps)):
    th4 = th4 - 2*math.pi
  if ((abs(v113) < eps) and (abs(v313) < eps) and (abs(v323) < eps)):
    th4 = th4old

  th4old = th4
  
  #
  c4 = math.cos(th4); s4 = math.sin(th4)

  # solve for th5 -------------------------------------------------
  k1 = c4*v313 + s4*v323
  k2 = s23*v113 - c23*wz
  th5 = math.atan2(k1, k2)

  #
  c5 = math.cos(th5)
  s5 = math.sin(th5)

  # solve for th6 --------------------------------------------------
  v111 = c1*ux + s1*uy
  v131 = s1*ux - c1*uy
  v311 = c23*v111 + s23*uz
  v331 = s23*v111 - c23*uz
  v411 = c4*v311 + s4*v131
  v431 = -s4*v311 + c4*v131

  #
  k1 = v431
  k2 = c5*v411 - s5*v331
  th6 = atan2(k1, k2)










  # R1 = eulerAnglesToRotationMatrixT(np.array([gamma, beta, alpha]))
  # print("R1:", simplify(R1))
  # print()




  # R0u = Matrix(R1)
  # print("R0u: ", R0u)
  # print()

  # R03T = Matrix(T03T[:3, :3])
  # print("R03T: ", R03T)
  # print()

  # # Rotation of urdf_gripper wrt (DH) gripper frame from rotz(pi) * roty(-pi/2) and it's transpose
  # Rgu_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])
  # RguT_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])

  # # Inverse kinematics transformations starts below
  # R0u_eval = R0u.evalf(subs = {alpha: yaw, beta: pitch, gamma: roll })
  # R0g_eval = R0u_eval * RguT_eval
  # print("R0g_eval: ", R0u_eval)
  # print()

  # wrist_center = get_wrist_center(gripper_point, R0g_eval, robot_dg)
  # print("Wrist center: ", wrist_center)
  # print()

  # j1, j2, j3 = get_first_three_angles(wrist_center)

  # R03T_eval = R03T.evalf(subs = {q1: j1.evalf(), q2: j2.evalf(), q3: j3.evalf()})
  # R36_eval = R03T_eval * R0g_eval
  # print("R03T_eval: ", R03T_eval)
  # print()


  # j4, j5, j6 = get_last_three_angles(R36_eval)

  # j1 = j1.evalf()
  # j2 = j2.evalf()
  # j3 = j3.evalf()
  # j4 = j4.evalf()
  # j5 = j5.evalf()
  # j6 = j6.evalf()

  print("Invers kinematic:")
  print("q1 +/-:", th1, "rad; ", math.degrees(th1), "°")
  print("q2:", th2, "rad; ", math.degrees(th2), "°")
  print("q3 +/-:", th3, "rad; ", math.degrees(th3), "°")
  print("q4:", th4, "rad; ", math.degrees(th4), "°")
  print("q5 +/-:", th5, "rad; ", math.degrees(th5), "°")
  print("q6:", th6, "rad; ", math.degrees(th6), "°")

  return th1, th2, th3, th4, th5, th6








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

  px, py, pz = -0.00423, 0.25884, 0.46696
  roll, pitch, yaw = 89.88, 89.89, 0.03


  # print("get_angles ...")
  #q1, q2, q3, q4, q5, q6 = inverseTransformation1(px, py, pz, roll, pitch, yaw)
  # forward_kin(q1, q2, q3, q4, q5, q6)
  # print()
  # print()
  # print("########################################")


  forwardTransformation(89.29, -34.05, 125.42, -21.60, -1.60, 21.63)        # Result should be: x=-4.23, y=258.84, z=466.96, a=89.88, b=89.89, c=0.03
  # print("----------------------------------------")
  # forwardTransformation(93.76, -54.76, 118.97, -4.89, 30.30, 1.64)          # Result should be: x=-8.60, y=115.87, z=522.94, a=91.28, b=94.42, c=-2.78

if __name__=="__main__":
  main()