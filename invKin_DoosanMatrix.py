from calendar import c
import math
from re import A
import rlcompleter
from symtable import Symbol
from tkinter import BaseWidget
from tokenize import Double
from sympy import symbols, simplify, pprint,  expand_trig, trigsimp, pi, atan2, cos, sin, tan, sqrt, atan, atanh
from sympy.matrices import Matrix

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


global robot_d1, robot_d2, robot_d3, robot_d4, robot_d5, robot_d6, robot_d7
global robot_a0, robot_a1, robot_a2, robot_a3, robot_a4, robot_a5, robot_a6
global robot_alpha0, robot_alpha1, robot_alpha2, robot_alpha3, robot_alpha4, robot_alpha5, robot_alpha6
global robot_theta1, robot_theta2, robot_theta3, robot_theta4, robot_theta5, robot_theta6, robot_theta7

# Doosan
robot_d1 = 0.135
robot_d2 = 0.0062
#robot_d2 = 0.00625
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
robot_alpha1 = -(pi/2)
robot_alpha2 = 0.0
robot_alpha3 = (pi/2)
robot_alpha4 = -(pi/2)
robot_alpha5 = (pi/2)
robot_alpha6 = 0.0

robot_theta1 = 0.0
robot_theta2 = -(pi/2)
robot_theta3 = (pi/2)
robot_theta4 = 0.0
robot_theta5 = 0.0
robot_theta6 = 0.0
robot_theta7 = 0.0

# float[] _toolCenterPositionKlein = { -52.750f, 0.000f, 118.10f, 0.000f, 0.000f, 0.000f };        // NiederstreiferKlein z, y, x, a, b, c
# float[] _toolCenterPositionGross = { -65.500f, 0.000f, 117.500f, 0.000f, 0.000f, 0.000f };       // NiederstreiferGross z, y, x, a, b, c
l1 = robot_d1
l2 = robot_a2
l3 = robot_d4
l4 = robot_d6
l5 = 0.118100

a2 = robot_a2
#a5 = 0.052750

d1 = robot_d1
d2 = robot_d2
d4 = robot_d4
d6 = robot_d6

# Print matrix in console
def printMatrix(a):
  print("Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]")
  rows = a.shape[0]
  cols = a.shape[1]
  for i in range(0,rows):
    for j in range(0,cols):
        print ("%6.f\t" %a[i,j], end = '')
    print()
  print()

# Print matrix in console
def printMatrixE(a):
  print("Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]")
  rows = a.shape[0]
  cols = a.shape[1]
  for i in range(0,rows):
    for j in range(0,cols):
        print("%6.20f\t" %a[i,j], end = '')
    print()
  print()

# Print matrix in console
def printMatrixEToDegrees(a):
  print("Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]")
  rows = a.shape[0]
  cols = a.shape[1]
  for i in range(0,rows):
    for j in range(0,cols):
        print("%6.18f\t" %math.degrees(a[i,j]), end = '')
    print()
  print()

# Print matrix to file
def printMatrixEtoFile(a, f):
  f.write("Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]\n")
  rows = a.shape[0]
  cols = a.shape[1]
  for i in range(0,rows):
    for j in range(0,cols):
        f.write("%6.18f\t" %(a[i,j]))
    f.write("\n")
  f.write("\n")



# Transformation matrix Classic i-1Ti (https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters))
def transformationMatrixClassic(theta, alpha, a, d):
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
  return T

# Transformation matrix Classic i-1Ti (https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters)
def transformationMatrixClassicSymbol(theta, alpha, a, d):
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
  return T

# Transformation matrix Modified i-1Ti (https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters)
def transformationMatrixModified(theta, alpha, a, d):
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
  return T

# Transformation matrix Modified i-1Ti (https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters)
def transformationMatrixModifiedSymbol(theta, alpha, a, d):
  r11 = cos(theta)
  r12 = -sin(theta)
  r13 = 0.0
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
  return T



def rot_x_Symbolic(q):
    r_x = Matrix([[1,         0,        0,        0],
                  [0,         cos(q),   -sin(q),  0],
                  [0,         sin(q),   cos(q),   0],
                  [0,         0,        0,        1]])
    return r_x

def rot_y_Symbolic(q):
    r_y = Matrix([[cos(q),    0,        sin(q),   0],
                  [0,         1,        0,        0],
                  [-sin(q),   0,        cos(q),   0],
                  [0,         0,        0,        1]])
    return r_y

def rot_z_Symbolic(q):
    r_z = Matrix([[cos(q),  -sin(q),    0,        0],
                  [sin(q),  cos(q),     0,        0],
                  [0,       0,          1,        0],
                  [0,       0,          0,        1]])
    return r_z



def rot_x(q):
    r_x = Matrix([[1,         0,        0,        0],
                  [0,         math.cos(q),   -math.sin(q),  0],
                  [0,         math.sin(q),   math.cos(q),   0],
                  [0,         0,        0,        1]])
    return r_x

def rot_y(q):
    r_y = Matrix([[math.cos(q),    0,        math.sin(q),   0],
                  [0,         1,        0,        0],
                  [-math.sin(q),   0,        math.cos(q),   0],
                  [0,         0,        0,        1]])
    return r_y

def rot_z(q):
    r_z = Matrix([[math.cos(q),  -math.sin(q),    0,        0],
                  [math.sin(q),  math.cos(q),     0,        0],
                  [0,       0,          1,        0],
                  [0,       0,          0,        1]])
    return r_z



# https://github.com/hortovanyi/RoboND-Kinematics-Project/blob/master/FK_debug.py
# https://hive.blog/hive-196387/@juecoree/forward-kinematics-of-puma-560-robot-using-dh-method
def forwardKinematic(q11, q22, q33, q44, q55, q66, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f):
  print("########################################")
  f.writelines("########################################\n")
  print("Forward kinematic start ...")
  f.writelines("Forward kinematic start ...\n")
  print()
  f.writelines("\n")
  print("Input:")

  q1 = math.radians(q11)
  print("q1: ", q1, "rad; ", q11, "°")
  f.writelines("q1: %.20f\n" % q1)

  q2 = math.radians(q22)
  print("q2: ", q2, "rad; ", q22, "°")
  f.writelines("q2: %.20f\n" % q2)

  q3 = math.radians(q33)
  print("q3: ", q3, "rad; ", q33, "°")
  f.writelines("q3: %.20f\n" % q3)

  q4 = math.radians(q44)
  print("q4: ", q4, "rad; ", q44, "°")
  f.writelines("q4: %.20f\n" % q4)

  q5 = math.radians(q55)
  print("q5: ", q5, "rad; ", q55, "°")
  f.writelines("q5: %.20f\n" % q5)

  q6 = math.radians(q66)
  print("q6: ", q6, "rad; ", q66, "°")
  f.writelines("q6: %.20f\n" % q6)

  print()
  f.writelines("\n")
  print("---------------------------------")
  f.writelines("---------------------------------\n")
  print()
  f.writelines("\n")

  # end effector corection (eg.: ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000)
  robot_d7 = ez
  robot_a6 = ex

  # T01
  T01 = transformationMatrixModified(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  # print("T01:")
  # printMatrixE(T01)
  T0g = T01

   # T12
  T12 = transformationMatrixModified(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  # print("T12:")
  # printMatrixE(T12)
  T0g = T0g * T12

  # T23
  T23 = transformationMatrixModified(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  # print("T23:")
  # printMatrixE(T23)
  T0g = T0g * T23

  # T34
  T34 = transformationMatrixModified(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  # print("T34:")
  # printMatrixE(T34))
  T0g = T0g * T34

  # T45
  T45 = transformationMatrixModified(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  # print("T45:")
  # printMatrixE(T45)
  T0g = T0g * T45

  # T56
  T56 = transformationMatrixModified(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  # print("T56:")
  # printMatrixE(T56)
  T0g = T0g * T56

  # T6g (final position and rotation)
  T6g = transformationMatrixModified(0 + robot_theta7, robot_alpha6, robot_a6, robot_d7)
  # print("T6g:")
  # printMatrixE(T6g)

  T0g = T0g * T6g
  print("T0g:")
  printMatrixE(T0g)
  f.writelines("T0g:\n")
  printMatrixEtoFile(T0g, f)

  # Calculate angles (pitch, roll yaw)
  cw = math.atan2(math.sqrt(math.pow(T0g[0,2], 2) + math.pow(T0g[1,2], 2)), T0g[2,2])
  bw = math.atan2(T0g[2,1], -T0g[2,0])
  aw = math.atan2(T0g[1,2], T0g[0,2])

  # Calculates Rotation Matrix given euler angles (https://mathworld.wolfram.com/RotationMatrix.html, https://mathworld.wolfram.com/EulerAngles.html, https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix)
  #R1 = np.dot(rot_z(aw), np.dot(rot_y(cw), rot_z(bw)))[:3, :3]
  R1 = (rot_z_Symbolic(aw) * rot_y_Symbolic(cw) * rot_z_Symbolic(bw))[:3, :3]
  print("R1:")
  printMatrixE(R1)
  f.writelines("R1:")
  printMatrixEtoFile(R1, f)

  so = math.sqrt(1 - math.pow(R1[2,2], 2))
  print("so +: %.20f" % +so)
  print("so -: %.20f" % -so)

  co = R1[2,2]
  print("co: %.20f" % co)

  print("r13: %.20f" % R1[0,2])
  print("r23: %.20f" % R1[1,2])
  print()


  tmp2 = math.atan2(R1[0,2], R1[1,2])
  print("tmp2: %.20f" % math.degrees(tmp2))
  tmp2 = math.atan2(-R1[0,2], -R1[1,2])
  print("tmp2: %.20f" % math.degrees(tmp2))

  tmp1 = math.atan2(R1[2,2], +math.sqrt(1 - math.pow(R1[2,2], 2)))
  print("tmp1: %.20f" % math.degrees(tmp1))
  tmp1 = math.atan2(R1[2,2], -math.sqrt(1 - math.pow(R1[2,2], 2)))
  print("tmp1: %.20f" % math.degrees(tmp1))

  tmp3 = math.atan2(-R1[2,0], R1[2,1])
  print("tmp3: %.20f" % math.degrees(tmp3))
  tmp3 = math.atan2(R1[2,0], -R1[2,1])
  print("tmp3: %.20f" % math.degrees(tmp3))

  # # # end effector corection
  # # # ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000
  # # robot_d7 = ez
  # # #robot_a6 = ex

  # # Matrix for the end effector
  # vector = Matrix([[0, 0, ex]])
  # printMatrixE(vector.T)

  # # RR = np.dot(rot_z(aw)[:3, :3], np.dot(rot_y(cw)[:3, :3], rot_z(bw)[:3, :3]))
  # # printMatrixE(RR)
  # xyz = R1 * vector.T
  # print("xyz:")
  # printMatrixE(xyz)

  # # # Test for the end effector correction
  # # side1 = sin(cw) * ex
  # # print("zz:")
  # # print(side1)

  # # side2 = cos(cw) * ex
  # # print("zz:")
  # # print(side2)

  print("---------------------------------")
  f.writelines("---------------------------------\n")
  print()
  f.writelines("\n")
  print("Summary ...")
  f.writelines("Summary ...\n")

  print("Pose:")
  f.writelines("Pose:\n")
  print("x = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f"                % (T0g[0,3], x, x - T0g[0,3]))
  f.writelines("x = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f\n"       % (T0g[0,3], x, x - T0g[0,3]))

  # print("x = %6.20f\t\toriginal: %d\t\tdifference: %6.20f"                % (T0g[0,3], format(x, '.20f'), x - T0g[0,3]))
  # f.writelines("x = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f\n"       % (T0g[0,3], x, x - T0g[0,3]))

  print("y = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f"                % (T0g[1,3], y, y - T0g[1,3]))
  f.writelines("y = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f\n"       % (T0g[1,3], y, y - T0g[1,3]))
  print("z = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f"                % (T0g[2,3], z, z - T0g[2,3]))
  f.writelines("z = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f\n"       % (T0g[2,3], z, z - T0g[2,3]))
  # print("x = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f"                % (T0g[0,3], x, x - T0g[0,3] + xyz[0,0]))
  # f.writelines("x = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f\n"       % (T0g[0,3], x, x - T0g[0,3] + xyz[0,0]))
  # print("y = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f"                % (T0g[1,3], y, y - T0g[1,3] - xyz[1,0]))
  # f.writelines("y = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f\n"       % (T0g[1,3], y, y - T0g[1,3] - xyz[1,0]))
  # print("z = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f"                % (T0g[2,3], z, z - T0g[2,3] + xyz[2,0]))
  # f.writelines("z = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f\n"       % (T0g[2,3], z, z - T0g[2,3] + xyz[2,0]))
  print("a = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f"                % (math.degrees(aw), a, a - math.degrees(aw)))
  f.writelines("a = %6.20f\toriginal: %6.20f\t\tdifference: %6.20f\n"         % (math.degrees(aw), a, a - math.degrees(aw)))
  print("b = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f"                % (math.degrees(cw), b, b - math.degrees(cw)))
  f.writelines("b = %6.20f\toriginal: %6.20f\t\tdifference: %6.20f\n"         % (math.degrees(cw), b, b - math.degrees(cw)))
  print("c = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f"                % (math.degrees(bw), c, c - math.degrees(bw)))
  f.writelines("c = %6.20f\t\toriginal: %6.20f\t\tdifference: %6.20f\n"       % (math.degrees(bw), c, c - math.degrees(bw)))
  print("Angles:")
  f.writelines("Angles:\n")
  print("q1, q2, q3, q4, q5, q6 = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f"               %(math.degrees(q1), math.degrees(q2), math.degrees(q3), math.degrees(q4), math.degrees(q5), math.degrees(q6)))
  print("q1, q2, q3, q4, q5, q6 = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f"             %(math.degrees(q1), math.degrees(q2), math.degrees(q3), math.degrees(tmp2), math.degrees(tmp1), math.degrees(tmp3)))
  print()
  f.writelines("q1, q2, q3, q4, q5, q6 = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f\n"      %(math.degrees(q1), math.degrees(q2), math.degrees(q3), math.degrees(q4), math.degrees(q5), math.degrees(q6)))
  f.writelines("q1, q2, q3, q4, q5, q6 = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f\n\n"      %(math.degrees(q1), math.degrees(q2), math.degrees(q3), math.degrees(tmp2), math.degrees(tmp1), math.degrees(tmp3)))

  print("---------------------------------")
  f.writelines("---------------------------------\n")
  print()
  f.writelines("\n")

  print("Summary for invers kinematic ...")
  f.writelines("Summary for invers kinematic ...\n")
  print("x, y, z, a, b, c = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f"                     %(T0g[0,3], T0g[1,3], T0g[2,3], math.degrees(aw), math.degrees(cw), math.degrees(bw)))
  f.writelines("x, y, z, a, b, c = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f\n"            %(T0g[0,3], T0g[1,3], T0g[2,3], math.degrees(aw), math.degrees(cw), math.degrees(bw)))
  print("q1, q2, q3, q, q5, q6 = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f"                %(math.degrees(q1), math.degrees(q2), math.degrees(q3), math.degrees(q4), math.degrees(q5), math.degrees(q6)))
  f.writelines("q1, q2, q3, q4, q5, q6 = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f\n"      %(math.degrees(q1), math.degrees(q2), math.degrees(q3), math.degrees(q4), math.degrees(q5), math.degrees(q6)))
  print("ex, ey, ez, ea, eb, ec = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f"               %(ex, ey, ez, ea, eb, ec))
  f.writelines("ex, ey, ez, ea, eb, ec = %6.20f, %6.20f, %6.20f, %6.20f, %6.20f, %6.20f\n\n"      %(ex, ey, ez, ea, eb, ec))

  print()
  f.writelines("\n")
  print("Forward kinematic end ...")
  f.writelines("Forward kinematic end ...\n")
  print("########################################")
  f.writelines("########################################\n")
  return 1


# https://nitishpuri.github.io/posts/robotics/inverse-kinematics-on-kuka-arm-using-ros-and-python/
# https://github.com/hortovanyi/RoboND-Kinematics-Project/blob/master/IK_debug.py
# https://github.com/camisatx/RoboticsND/blob/master/projects/kinematics/kuka_kr210/kuka_ik.py
# https://hive.blog/hive-196387/@juecoree/inverse-kinematics-of-puma-560-robot
# https://github.com/zahid58/Inverse-Kinematics-6-DOF-for-ERC-2019
# https://github.com/mithi/arm-ik
# https://www.hindawi.com/journals/mpe/2021/6647035/
def inversKinematic(x, y, z, roll, pitch, yaw, ex, ey, ez, ea, eb, ec, q1Original, q2Original, q3Original, q4Original, q5Original, q6Original):
  print("Invers kinematic start ...")
  print("########################################")
  print()
  print("Input:")
  print("x: ", x)
  print("y: ", y)
  print("z: ", z)
  print("Roll : ", roll, "°; ", math.radians(roll), "rad")
  print("Pitch : ", pitch, "°; ", math.radians(pitch), "rad")
  print("Yaw : ", yaw, "°; ", math.radians(yaw), "rad")

  print()
  print("---------------------------------")
  print()

  # Convert degree to radians for the angles q1, q2, q3
  roll = math.radians(roll)
  pitch = math.radians(pitch)
  yaw = math.radians(yaw)

  # end effector corection
  # ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000
  robot_d7 = ez
  #robot_a6 = ex

  # Matrix for the end effector
  vector = Matrix([[0, 0, ex]])
  print("vector:")
  printMatrixE(vector.T)

  RR = np.dot(rot_z_Symbolic(roll), np.dot(rot_y_Symbolic(pitch), rot_z_Symbolic(yaw)))[:3, :3]
  print("RR:")
  printMatrixE(RR)

  xyz = RR * vector.T
  print("xyz:")
  printMatrixE(xyz)

  # Test for the end effector correction
  side1 = sin(pitch) * ex
  print("side1: %.20f", side1)

  side2 = cos(pitch) * ex
  print("side2: %.20f", side2)

  print()
  print("---------------------------------")
  print()

  ################################################################################
  # All important symbolic transformations matrices are declared below
  ################################################################################

  q1, q2, q3, q4, q5, q6 = symbols('q1:7')
  alpha, beta, gamma = symbols('alpha beta gamma', real = True)
  px, py, pz = symbols('px py pz', real = True)

  ############################## Transformation Matrix Start ####################################

  # T01
  T01 = transformationMatrixModifiedSymbol(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  # print("T01:", T01)
  # print()
  # T0g = T01

  # T12
  T12 = transformationMatrixModifiedSymbol(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  # print("T12:", T12)
  # print()
  # T0g = T0g * T12

  # T23
  T23 = transformationMatrixModifiedSymbol(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  # print("T23:", T23)
  # print()
  # T0g = T0g * T23

  # T34
  T34 = transformationMatrixModifiedSymbol(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  # print("T34:", T34)
  # print()
  # T0g = T0g * T34

  # T45
  T45 = transformationMatrixModifiedSymbol(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  # print("T45:", T45)
  # print()
  # T0g = T0g * T45

  # T56
  T56 = transformationMatrixModifiedSymbol(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  # print("T56:", T56)
  # print()
  # T0g = T0g * T56

  # T6g (final position and rotation)
  T6g = transformationMatrixModifiedSymbol(0 + robot_theta7, robot_alpha6, robot_a6, robot_d7)
  # print("T6g:", T6g)
  # print()
  # T0g = T0g * T6g
  # print("T0g:", T0g)
  # print()

  ############################## Transformation Matrix End ####################################

  T03 = T01 * T12 * T23
  # print("T03:", simplify(T03))
  # print()

  T03T = T03.T
  # print("T03T:", simplify(T03T))
  # print()

  T04 = T01 * T12 * T23 * T34
  # print("T03:", simplify(T03))
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

  # Calculates Rotation Matrix given euler angles (https://mathworld.wolfram.com/RotationMatrix.html, https://mathworld.wolfram.com/EulerAngles.html, https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix)
  # R1 = np.dot(rot_z(alpha), np.dot(rot_y(beta), rot_z(gamma)))[:3, :3]
  R1 = (rot_z_Symbolic(gamma) * rot_y_Symbolic(beta) *rot_z_Symbolic(alpha))[:3, :3]
  # print("R1")
  # pprint(simplify(R1))
  # print()

  R0u = Matrix(R1)
  # print("R0u: ", R0u)
  # print()

  # Euler Transformation and evaluation
  print("Roll: ", math.degrees(roll),", Pitch: ", math.degrees(pitch),", Yaw: ", math.degrees(yaw))
  print()
  R0u_eval = R0u.evalf(subs = {alpha: roll, beta: pitch, gamma: yaw})
  # print("R0u_eval 0 R1:")
  # printMatrixE(R0u_eval)


  R0g_eval = R0u_eval
  print("R0g_eval:")
  printMatrixE(R0g_eval)
  print()

  print("---------------------------------")
  print()

  # ##################################################### q1, q2, q3 #####################################################





  tmp1 = math.sin(x / math.sqrt(math.pow(x, 2) + math.pow(y, 2)))
  print("tmp1: %.20f" % math.degrees(tmp1))

  tmp2 = math.cos(y / math.sqrt(math.pow(x, 2) + math.pow(y, 2)))
  print("tmp2: %.20f" % math.degrees(tmp2))
  print()



  zetaDefault = math.sin(d2 / math.sqrt(math.pow(x, 2) + math.pow(y, 2)))
  zeta = math.sin(x / math.sqrt(math.pow(x, 2) + math.pow(y, 2)))

  delta1 = zeta - zetaDefault
  delta2 = zeta + zetaDefault




  Rzphi = rot_z(roll)
  print("Rzphi:")
  printMatrixE(Rzphi)
  Rytheta = rot_y(pitch)
  print("Rytheta:")
  printMatrixE(Rytheta)
  Rzpsi = rot_z(yaw)
  print("Rzpsi:")
  printMatrixE(Rzpsi)

  R0g_eval = (rot_z_Symbolic(roll) * rot_y_Symbolic(pitch) * rot_z_Symbolic(yaw))[:3, :3]
  Rrpy = np.dot(rot_z_Symbolic(roll), np.dot(rot_y_Symbolic(pitch), rot_z_Symbolic(yaw)))[:3, :3]
  print("R0g_eval:")
  printMatrixE(R0g_eval)
  print()


  px = x - (robot_d6 + robot_d7) * R0g_eval[0, 2]
  py = y - (robot_d6 + robot_d7) * R0g_eval[1, 2]
  pz = z - (robot_d6 + robot_d7) * R0g_eval[2, 2]

  px = x - (robot_d6 + robot_d7) * R0g_eval[0, 2] + xyz[0,0]
  py = y - (robot_d6 + robot_d7) * R0g_eval[1, 2] + xyz[1,0]
  pz = z - (robot_d6 + robot_d7) * R0g_eval[2, 2] + xyz[2,0]

  print("px: %.20f" % px)
  print("py: %.20f" % py)
  print("pz: %.20f" % pz)
  print()





  # tmp3 = math.atan2(py, px)
  # print("tmp3: %.20f" % math.degrees(tmp3))
  # print()

  # tmp4 = math.asin(robot_d2 / math.sqrt(math.pow(py, 2) + math.pow(px, 2) - math.pow(robot_d2, 2)))
  # print("tmp4: %.20f" % math.degrees(tmp4))
  # print()

  # tmp5 = math.asin(math.sqrt(math.pow(py, 2) + math.pow(px, 2) - math.pow(robot_d2, 2)) / math.sqrt(math.pow(py, 2) + math.pow(px, 2)))
  # print("tmp5: %.20f" % math.degrees(tmp5))
  # print()

  # tmp6 = math.asin(robot_d2 / math.sqrt(math.pow(py, 2) + math.pow(px, 2)))
  # print("tmp6: %.20f" % math.degrees(tmp6))
  # print()

  # px = px - robot_d2 * math.cos(tmp3)
  # py = py - robot_d2 * math.sin(tmp3)
  # pz = pz

  # print("px: %.20f" % px)
  # print("py: %.20f" % py)
  # print("pz: %.20f" % pz)
  # print()




  eps = 0.00001
  ELBOW = 1
  ARM = -1


  r = math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2))
  R = math.sqrt(math.pow(px, 2) + math.pow(py, 2))

  gamma1 = math.asin(py / R)
  gamma2 = math.acos(px / R)

  alpha1 = math.asin(d2 / R)
  alpha2 = math.acos(r / R)





  tmp1 = (py * math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2)) - px * d2) / (math.pow(px, 2) + math.pow(py, 2))
  tmp2 = (px * math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2)) + py * d2) / (math.pow(px, 2) + math.pow(py, 2))

  tmp1 = (py * math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2)) - px * d2)
  tmp2 = (px * math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2)) + py * d2)

  # tmp1 = -px * d2 - ARM * py * math.sqrt(math.pow(px, 2) + math.pow(py, 2) + math.pow(d2, 2))
  # tmp2 = +py * d2 - ARM * px * math.sqrt(math.pow(px, 2) + math.pow(py, 2) + math.pow(d2, 2))

  tmp3 = math.atan2(tmp1, tmp2)
  print("phi1: %.20f" % math.degrees(tmp3))
  print()

  tmp4 = math.atan2(-tmp1, -tmp2)
  print("phi1: %.20f" % math.degrees(tmp4))
  print()






  # px = px - robot_d2 * math.sin(tmp3)
  # py = py - robot_d2 * math.cos(tmp3)
  # pz = pz

  # print("px: %.20f" % px)
  # print("py: %.20f" % py)
  # print("pz: %.20f" % pz)
  # print()







  r = math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2))
  R = math.sqrt(math.pow(px, 2) + math.pow(py, 2) + math.pow(pz, 2) - math.pow(d2, 2))

  sinAlpha = - ((pz) / (R))
  cosAlpha = - ((ARM * r) / (R))

  cosBeta = ((math.pow(robot_a2, 2) + math.pow(R, 2) - (math.pow(robot_d4, 2) + math.pow(robot_a3, 2))) / (2 * robot_a2 * R))
  if((1 - math.pow(cosBeta, 2)) <= eps):
    sinBeta = 0.0
  else:
    sinBeta = math.sqrt(1 - math.pow(cosBeta, 2))

  tmp1 = sinAlpha * cosBeta + (ARM * ELBOW) * cosAlpha * sinBeta
  tmp2 = cosAlpha * cosBeta - (ARM * ELBOW) * sinAlpha * sinBeta

  tmp4 = math.atan2(tmp1, tmp2)
  print("phi2: %.20f" % math.degrees(tmp4))
  print()

















  # Different situations lead to different solutions
  if((math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2)) < 0):
    # There is no solution for theta1.
    th1 = math.nan
    th2 = math.nan
    th3 = math.nan
    th4 = math.nan
    th5 = math.nan
    th6 = math.nan
  else:
    # Do inverse position kinematics and calculate theta1, theta2, theta3.
    #  Calculate theta1.
    if ((math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2)) == 0):
      # Singularity configurations. Only one solution for theta1. Here we
      # assign two identical solutions to theta1, which will be eliminated to one set at the end.
      theta1_1 = math.atan2(py, px)- math.pi/2
      theta1_2 = math.atan2(py, px)- math.pi/2
    else:
      theta1_1 = atan2(py, px) - math.atan2(d2, sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2)))
      theta1_2 = atan2(py, px) + math.atan2(-d2, -sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2)))
    # Calculate theta3.
    r = math.sqrt(math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2))
    s = pz - d1
    cos_temp = (math.pow(d4, 2) + math.pow(a2, 2) - math.pow(r, 2) - math.pow(s, 2)) /2 /d4 /a2
    if (abs(cos_temp) > 1):
      # This is when the desired position is outside of the robot's reachable workspace.
      th1 = math.nan
      th2 = math.nan
      th3 = math.nan
      th4 = math.nan
      th5 = math.nan
      th6 = math.nan
    else:
      temp = math.atan2(math.sqrt(1 - math.pow(cos_temp, 2)), cos_temp)

      theta3_1 = math.pi/2 - temp
      theta3_2 = math.pi/2 + temp

      # Calculate theta2.
      theta2_1 = math.atan2(s, r) - math.atan2(d4 * math.cos(theta3_1), a2 - d4 * math.sin(theta3_1))
      theta2_2 = math.atan2(s, r) - math.atan2(d4 * math.cos(theta3_2), a2 - d4 * math.sin(theta3_2))

      theta3_3 = math.pi - theta3_1
      theta2_3 = math.pi - theta2_1

      theta3_4 = math.pi - theta3_2
      theta2_4 = math.pi - theta2_2

      # Thus the four solutions for the arm are:
      # 1. theta1_1, theta2_1, theta3_1;
      # 2. theta1_1, theta2_2, theta3_2;
      # 3. theta1_2, theta2_3, theta3_3;
      # 4. theta1_2, theta2_4, theta3_4;

      ## Inverse orientation kinematics

      # Calculate R03 using forward kinematics.
      ### 1st pair of solutions
      #A1 = dh_kuchenbe(0,  pi/2,   a, theta1_1);
      A1 = transformationMatrixModifiedSymbol(theta1_1, robot_alpha0, robot_a0, robot_d1)
      #A2 = dh_kuchenbe(c,     0,  -b, theta2_1);
      A2 = transformationMatrixModifiedSymbol(theta2_1, robot_alpha1, robot_a1, robot_d2)
      #A3 = dh_kuchenbe(0, -pi/2,  -d, theta3_1);
      A3 = transformationMatrixModifiedSymbol(theta3_1, robot_alpha2, robot_a2, robot_d3)

      A03 = A1 * A2 * A3
      R03 = A03[:3, :3]

      # Calculate R36.
      #R36 = R03.T * R06
      R36 = R03.T * R0g_eval

      ctheta = R36[2,2]
      stheta_pos = math.sqrt(1 - math.pow(ctheta, 2))
      stheta_neg = -math.sqrt(1 - math.pow(ctheta, 2))
      theta5_1 = -math.atan2(stheta_pos, ctheta)
      theta4_1 = math.atan2(R36[1,2], R36[0,2])
      theta6_1 = math.atan2(R36[2,1], -R36[2,0])

      theta5_2 = -math.atan2(stheta_neg, ctheta)
      theta4_2 = math.atan2(-R36[1,2], -R36[0,2])
      theta6_2 = math.atan2(-R36[2,1], R36[2,0])

      if ((R36[1,2] == 0) and (R36[0,2] == 0)):
        theta4_1 = 0
        theta6_1 = 0
        theta4_2 = math.pi
        theta6_2 = math.pi

      ### 2nd pair of solutions
      #A1 = dh_kuchenbe(0,  pi/2,   a, theta1_1);
      A1 = transformationMatrixModifiedSymbol(theta1_1, robot_alpha0, robot_a0, robot_d1)
      #A2 = dh_kuchenbe(c,     0,  -b, theta2_1);
      A2 = transformationMatrixModifiedSymbol(theta2_1, robot_alpha1, robot_a1, robot_d2)
      #A3 = dh_kuchenbe(0, -pi/2,  -d, theta3_1);
      A3 = transformationMatrixModifiedSymbol(theta3_1, robot_alpha2, robot_a2, robot_d3)

      A03 = A1 * A2 * A3
      R03 = A03[:3, :3]

      # Calculate R36.
      #R36 = R03.T * R06
      R36 = R03.T * R0g_eval

      ctheta = R36[2,2]
      stheta_pos = math.sqrt(1 - math.pow(ctheta, 2))
      stheta_neg = -math.sqrt(1 - math.pow(ctheta, 2))
      theta5_3 = -math.atan2(stheta_pos, ctheta)
      theta4_3 = math.atan2(R36[1,2], R36[0,2])
      theta6_3 = math.atan2(R36[2,1], -R36[2,0])

      theta5_4 = -math.atan2(stheta_neg, ctheta)
      theta4_4 = math.atan2(-R36[1,2], -R36[0,2])
      theta6_4 = math.atan2(-R36[2,1], R36[2,0])

      if ((R36[1,2] == 0) and (R36[0,2] == 0)):
        theta4_3 = 0
        theta6_3 = 0
        theta4_4 = math.pi
        theta6_4 = math.pi

      ### 3rd pair of solutions
      #A1 = dh_kuchenbe(0,  pi/2,   a, theta1_1);
      A1 = transformationMatrixModifiedSymbol(theta1_1, robot_alpha0, robot_a0, robot_d1)
      #A2 = dh_kuchenbe(c,     0,  -b, theta2_1);
      A2 = transformationMatrixModifiedSymbol(theta2_1, robot_alpha1, robot_a1, robot_d2)
      #A3 = dh_kuchenbe(0, -pi/2,  -d, theta3_1);
      A3 = transformationMatrixModifiedSymbol(theta3_1, robot_alpha2, robot_a2, robot_d3)

      A03 = A1 * A2* A3
      R03 = A03[:3, :3]

      # Calculate R36.
      #R36 = R03.T * R06
      R36 = R03.T * R0g_eval

      ctheta = R36[2,2]
      stheta_pos = math.sqrt(1 - math.pow(ctheta, 2))
      stheta_neg = -math.sqrt(1 - math.pow(ctheta, 2))
      theta5_5 = -math.atan2(stheta_pos, ctheta)
      theta4_5 = math.atan2(R36[1,2], R36[0,2])
      theta6_5 = math.atan2(R36[2,1], -R36[2,0])

      theta5_6 = -atan2(stheta_neg, ctheta)
      theta4_6 = math.atan2(-R36[1,2], -R36[0,2])
      theta6_6 = math.atan2(-R36[2,1], R36[2,0])

      if (R36[1,2] == 0 and R36[0,2] == 0):
        theta4_5 = 0
        theta6_5 = 0
        theta4_6 = math.pi
        theta6_6 = math.pi

      ### 4th pair of solutions
      #A1 = dh_kuchenbe(0,  pi/2,   a, theta1_1);
      A1 = transformationMatrixModifiedSymbol(theta1_1, robot_alpha0, robot_a0, robot_d1)
      #A2 = dh_kuchenbe(c,     0,  -b, theta2_1);
      A2 = transformationMatrixModifiedSymbol(theta2_1, robot_alpha1, robot_a1, robot_d2)
      #A3 = dh_kuchenbe(0, -pi/2,  -d, theta3_1);
      A3 = transformationMatrixModifiedSymbol(theta3_1, robot_alpha2, robot_a2, robot_d3)

      A03 = A1 * A2* A3
      R03 = A03[:3, :3]

      # Calculate R36.
      #R36 = R03.T * R06
      R36 = R03.T * R0g_eval

      ctheta = R36[2,2]
      stheta_pos = math.sqrt(1 - math.pow(ctheta, 2))
      stheta_neg = -math.sqrt(1 - math.pow(ctheta, 2))
      theta5_7 = -math.atan2(stheta_pos, ctheta)
      theta4_7 = math.atan2(R36[1,2], R36[0,2])
      theta6_7 = math.atan2(R36[2,1], -R36[2,0])

      theta5_8 = -math.atan2(stheta_neg, ctheta)
      theta4_8 = math.atan2(-R36[1,2], -R36[0,2])
      theta6_8 = math.atan2(-R36[2,1], R36[2,0])

      if (R36[1,2] == 0 and R36[0,2] == 0):
        theta4_7 = 0
        theta6_7 = 0
        theta4_8 = math.pi
        theta6_8 = math.pi

      # Put all the solutions of theta together
      if ((math.pow(px, 2) + math.pow(py, 2) - math.pow(d2, 2)) == 0):
        # Under singularity configurations, there are four solutions.
        th1 = [theta1_1, theta1_1, theta1_1, theta1_1]
        th2 = [theta2_1, theta2_2, theta2_1, theta2_2]
        th3 = [theta3_1, theta3_2, theta3_1, theta3_2]
        th4 = [theta4_1, theta4_3, theta4_2, theta4_4]
        th5 = [theta5_1, theta5_3, theta5_2, theta5_4]
        th6 = [theta6_1, theta6_3, theta6_2, theta6_4]
      else:
        # Normal configurations, there are eight solutions.
        th1 = [theta1_1, theta1_1, theta1_2, theta1_2, theta1_1, theta1_1, theta1_2, theta1_2]
        th2 = [theta2_1, theta2_2, theta2_3, theta2_4, theta2_1, theta2_2, theta2_3, theta2_4]
        th3 = [theta3_1, theta3_2, theta3_3, theta3_4, theta3_1, theta3_2, theta3_3, theta3_4]
        th4 = [theta4_1, theta4_3, theta4_5, theta4_7, theta4_2, theta4_4, theta4_6, theta4_8]
        th5 = [theta5_1, theta5_3, theta5_5, theta5_7, theta5_2, theta5_4, theta5_6, theta5_8]
        th6 = [theta6_1, theta6_3, theta6_5, theta6_7, theta6_2, theta6_4, theta6_6, theta6_8]

      thetas =  Matrix([th1, th2, th3, th4, th5, th6]).T
      print("thetas:")
      printMatrixEToDegrees(thetas)
      print()






















  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! J1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Left arm configuration ... "Robot Modeling and Control (Mark W. Spong, Seth Hutchinson, and M. Vidyasagar)" S. 92 (3.47, 3.48)
  eps = 0.00001

  tmp1 = math.pow(px, 2) + math.pow(py, 2) - math.pow(robot_d2, 2)
  if(tmp1 <= eps):
    r = 0.000
  else :
    r = math.sqrt(tmp1)

  s = pz - robot_d1

  tmp1 = math.atan2(px, py)
  tmp2 = math.atan2(r, robot_d2)
  tmp3 = math.atan2(math.sqrt(math.pow(px, 2) + math.pow(py, 2)), robot_d2)
  j1 = tmp1 + tmp2
  print("j1_1: %.20f" % math.degrees(j1))
  #jj1 = j1

  j1 = tmp1 - tmp2
  print("j1_2: %.20f" % math.degrees(j1))
  jj1 = j1

  # Right arm configuration ... "Robot Modeling and Control (Mark W. Spong, Seth Hutchinson, and M. Vidyasagar)" S. 93 (3.49)
  j1 = tmp1 + tmp2 + math.pi
  print("j1_3: %.20f" % math.degrees(j1))
  #jj1 = j1

  # OR ...

  tmp1 = math.pow(px, 2) + math.pow(py, 2) - math.pow(robot_d2, 2)
  if(tmp1 <= eps):
    r = 0.000
  else :
    r = math.sqrt(tmp1)

  s = pz - robot_d1

  tmp1 = math.atan2(px, py)
  tmp2 = math.atan2(-r, -robot_d2)
  j1 = tmp1 + tmp2
  print("j1_4: %.20f" % math.degrees(j1))

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! J3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  eps = 0.00001

  # Ellbow up & down configuration ... "Robot Modeling and Control (Mark W. Spong, Seth Hutchinson, and M. Vidyasagar)" S. 94 (3.55)
  tmp1 = math.pow(px, 2) + math.pow(py, 2)
  tmp2 = math.pow(robot_d2, 2)
  if((tmp1 - tmp2) >= eps):
    r = math.sqrt(tmp1 - tmp2)
  else :
    r = 0.000

  s = pz - robot_d1
  D = (math.pow(r, 2) + math.pow(s, 2) - math.pow(robot_a2, 2) - math.pow(robot_d4, 2)) / (2 * robot_a2 * robot_d4)

  if s >= 0: # up
    # OLD: j3 = math.atan2(-math.sqrt(1 - math.pow(D, 2)), D)
    if((1 - D) <= eps):
      j3 = math.atan2(0, D)
      #j3 = j3 - math.pi/2 # Hack
    else:
      j3 = math.atan2(+math.sqrt(1 - math.pow(D, 2)), D)
      #j3 = j3 - math.pi/2 # Hack
    print("j3_1 up: %.20f" % math.degrees(j3))
    jj3 = j3

  # j3 = math.atan2(D, 0)
  # #j3 = j3 - math.pi/2 # Hack
  # print("j3_1 up: %.20f" % math.degrees(j3))
  # ###############
  # tmp1 = math.pow(D, 2)
  # if((tmp1 <= (1.000 + eps)) and (tmp1 >= (1.000 - eps))):
  #   tmp1 = 0.000
  # else :
  #   tmp1 = math.sqrt(1 - tmp1)
  # j3 = math.atan2(D, +tmp1)
  # # j3 = math.atan2(D, +math.sqrt(1 - math.pow(D, 2)))
  # #j3 = j3 - math.pi/2 # Hack
  # print("j3_1 up: %.20f" % math.degrees(j3))



  if s < 0: # down
    # OLD: # j3 = math.atan2(+math.sqrt(1 - math.pow(D, 2)), D)
    if((1 - D) <= eps):
      j3 = math.atan2(0, D) - math.pi/2
      #j3 = j3 - math.pi/2 # Hack
    else:
      j3 = math.atan2(-math.sqrt(1 - math.pow(D, 2)), D)
      #j3 = j3 - math.pi/2 # Hack
    print("j3_2 down: %.20f" % math.degrees(j3))
    jj3 = j3

  # j3 = math.atan2(D, 0) - math.pi/2
  # #j3 = j3 - math.pi/2 # Hack
  # print("j3_2 down: %.20f" % math.degrees(j3))
  # ###############
  # tmp1 = math.pow(D, 2)
  # if((tmp1 <= (1.000 + eps)) and (tmp1 >= (1.000 - eps))):
  #   tmp1 = 0.000
  # else :
  #   tmp1 = math.sqrt(1 - tmp1)
  # j3 = math.atan2(D, -tmp1)
  # # j3 = math.atan2(D, -math.sqrt(1 - math.pow(D, 2)))
  # #j3 = j3 - math.pi/2 # Hack
  # print("j3_2 down: %.20f" % math.degrees(j3))
  # jj3 = j3

  #jj3 = j3 = 0.0
  #jj3 = j3 = math.pi/2

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! J2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  eps = 0.00001

  tmp1 = math.pow(px, 2) + math.pow(py, 2)
  tmp2 = math.pow(robot_d2, 2)

  if((tmp1 - tmp2) >= eps):
    r = math.sqrt(tmp1 - tmp2)
  else :
    r = 0.000

  s = pz - robot_d1

  test1 = math.atan2(r, s)
  test2 = math.atan2(robot_a2 + robot_d4 * cos(j3), robot_d4 * sin(j3))
  test2 = math.atan2(robot_d4 * sin(j3), robot_a2 + robot_d4 * cos(j3))

  j2 = math.atan2(r, s) - math.atan2(robot_a2 + robot_d4 * cos(j3), robot_d4 * sin(j3))
  print("j2_1: %.20f" % math.degrees(j2))
  j2 = math.atan2(r, s) - math.atan2(robot_d4 * sin(j3), robot_a2 + robot_d4 * cos(j3))
  print("j2_1: %.20f" % math.degrees(j2))
  jj2 = j2
  j2 = math.atan2(r, s) + math.atan2(robot_a2 + robot_d4 * cos(j3), robot_d4 * sin(j3))
  print("j2_1: %.20f" % math.degrees(j2))
  j2 = math.atan2(r, s) + math.atan2(robot_d4 * sin(j3), robot_a2 + robot_d4 * cos(j3))
  print("j2_1: %.20f" % math.degrees(j2))
  #j2 = j2 + math.pi/2 # Hack



  # Getting the right j1, j2, j3
  j1 = jj1
  j2 = jj2
  j3 = jj3

  # ##################################################### q4, q5, q6 #####################################################

  # Only for testing
  #j1, j2, j3 = math.radians(q1Original),  math.radians(q2Original), math.radians(q3Original)

  R03T = (T03T).evalf(subs={q1: j1, q2: j2, q3: j3})[:3, :3]

  Rrpy = np.dot(rot_z_Symbolic(roll), np.dot(rot_y_Symbolic(pitch), rot_z_Symbolic(yaw)))[:3, :3]
  Rrpy = (rot_z_Symbolic(roll) * rot_y_Symbolic(pitch) * rot_z_Symbolic(yaw))[:3, :3]
  Rrpy = (rot_x_Symbolic(roll) * rot_y_Symbolic(pitch) * rot_z_Symbolic(yaw))[:3, :3]
  #R_corr = simplify(rot_z(pi))
  R36 = R03T * Rrpy #* R_corr
  print("Rrpy")
  printMatrixE(Rrpy)
  print("R03T")
  printMatrixE(R03T)
  print("R36")
  printMatrixE(R36)

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! J5 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  eps = 0.00001

  tmp1 = sin(j1) * R36[0,2] - cos(j1) * R36[1,2]
  if ((1 - math.pow(tmp1, 2)) <= eps):
    tmp2 = 0.00
  else:
    tmp2 = math.sqrt(1 - math.pow(tmp1, 2))

  j5 = atan2(tmp1, +tmp2)
  #j5 = j5 - math.pi/2 # Hack
  print("q5: %.20f" %  math.degrees(j5))

  j5 = atan2(tmp1, -tmp2)
  j5 = j5 - math.pi/2 # Hack
  print("q5: %.20f" %  math.degrees(j5))

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! J4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  eps = 0.00001

  tmp1 = (cos(j1) * cos(j2 + j3) * R36[0,2] + sin(j1) * cos(j2 + j3) * R36[1,2] + sin(j2 + j3) * R36[2,2]).evalf()
  tmp2 = (-cos(j1) * sin(j2 * j3) * R36[0,2] - sin(j1) * sin(j2 + j3) * R36[1,2] + cos(j2 + j3) * R36[2,2]).evalf()
  # if ((1 - math.pow((sin(j1) * R36[0,1] - cos(j1) * R36[1,2]), 2)) <= eps):
  #   tmp2 = 0.00
  # else:
  #   tmp2 = math.sqrt(1 - math.pow((sin(j1) * R36[0,1] - cos(j1) * R36[1,2]), 2))

  if((tmp1 <= eps) and (tmp2 <= eps)):
    j4 = 0.0
  else:
    j4 = atan2(tmp1, tmp2)
  j4 = j4 - math.pi/2 # Hack
  print("q4: %.20f" %  math.degrees(j4))

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! J6 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  eps = 0.00001

  tmp1 = (-sin(j1) * R36[0,0] + cos(j1) * R36[1,0]).evalf()
  tmp2 = (sin(j1) * R36[0,1] - cos(j1) * R36[1,1]).evalf()
  # if ((1 - math.pow((sin(j1) * R36[0,1] - cos(j1) * R36[1,2]), 2)) <= eps):
  #   tmp2 = 0.00
  # else:
  #   tmp2 = math.sqrt(1 - math.pow((sin(j1) * R36[0,1] - cos(j1) * R36[1,2]), 2))

  if((tmp1 <= eps) and (tmp2 <= eps)):
    j6 = 0.0
  else:
    j6 = atan2(tmp1, tmp2)
  j6 = j6 - math.pi/2 # Hack
  print("q6: %.20f" %  math.degrees(j6))











  # # Variant 13
  # print()
  # print()
  # print("Variant 13")    # j4, j5, j6 = get_last_three_angles(R36_eval)
  # sin_q4 = R36[2, 2]
  # cos_q4 =  -R36[0, 2]

  # sin_q5 = sqrt(R36[0, 2]**2 + R36[2, 2]**2)
  # cos_q5 = R36[1, 2]

  # sin_q6 = -R36[1, 1]
  # cos_q6 = R36[1, 0]


  # # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! J4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # # j4 = atan2(cos_q4, sin_q4) + (math.pi / 2)
  # # print("q4", math.degrees(j4.evalf()))
  # j4 = atan2(cos_q4, sin_q4)
  # print("q4", math.degrees(j4.evalf()))

  # j5 = atan2(sin_q5, cos_q5) - math.pi
  # # print("q5", math.degrees(j5.evalf()))
  # # j5 = atan2(sin_q5, cos_q5)
  # # print("q5", math.degrees(j5.evalf()))
  # # j5 = atan2(sin_q5, cos_q5) + math.pi
  # # print("q5", math.degrees(j5.evalf()))
  # # j5 = atan2(sin_q5, cos_q5) + math.pi / 2
  # # print("q5", math.degrees(j5.evalf()))
  # # j5 = atan2(sin_q5, cos_q5) - math.pi / 2
  # # print("q5", math.degrees(j5.evalf()))

  # j6 = atan2(sin_q6, cos_q6)
  # print("q6", math.degrees(j6.evalf()))
  # j6 = (atan2(sin_q6, cos_q6) + math.pi) - math.pi
  # print("q6", math.degrees(j6.evalf()))
  # j6 = (atan2(sin_q6, cos_q6) + math.pi)
  # print("q6", math.degrees(j6.evalf()))
  # j6 = (atan2(sin_q6, cos_q6) + math.pi / 2)
  # print("q6", math.degrees(j6.evalf()))






  print()
  print("Result:")
  print("q1 +/-: %6.20frad;\t %6.20f°\t\toriginal: %6.20f°\t\tdifference: %6.20f° (left arm / right arm)" %(j1, math.degrees(j1), q1Original, (abs(q1Original) - abs(math.degrees(j1)))))
  print("q2: %6.20frad;\t\t %6.20f°\t\toriginal: %6.20f°\t\tdifference: %6.20f° (elbow up / elbow down)°" %(j2, math.degrees(j2), q2Original, (abs(q2Original) - abs(math.degrees(j2)))))
  print("q3 +/-: %6.20frad;\t %6.20f°\t\toriginal: %6.20f°\t\tdifference: %6.20f°" %(j3, math.degrees(j3), q3Original, (abs(q3Original) - abs(math.degrees(j3)))))
  print("q4: %6.20frad;\t\t %6.20f°\t\toriginal: %6.20f°\t\tdifference: %6.20f°" %(j4, math.degrees(j4), q4Original, abs(q4Original) - abs(math.degrees(j4))))
  print("q5 +/-: %6.20frad\t %6.20f°\t\toriginal: %6.20f°\t\tdifference: %6.20f°" %(j5, math.degrees(j5), q5Original, abs(q5Original) - abs(math.degrees(j5))))
  print("q6: %6.20frad\t\t %6.20f°\t\toriginal: %6.20f°\t\tdifference: %6.20f°" %(j6, math.degrees(j6), q6Original, abs(q6Original) - abs(math.degrees(j6))))

  print("Invers kinematic end ...")
  print("########################################")

  return math.degrees(j1), math.degrees(j2), math.degrees(j3), math.degrees(j4), math.degrees(j5), math.degrees(j6)
  #return j1, j2, j3, j4, j5, j6

def main():
  np.set_printoptions(infstr="(infinity)")
  np.set_printoptions(formatter={'float': "\t{: 0.0f}\t".format})

  # Calculation 1
  x, y, z, a, b, c = 0.00622573099065294126, -0.29285994375053248095, 0.54743996268304373487, 90.00280364221974593875, 90.03040263869252157747, -0.00041279850247499621
  q1, q2, q3, q4, q5, q6 = 91.72786000000000683485, -83.77836000000000638011, 82.95619999999999549800, -1.72524699999999997502, 90.85218999999999311967, -0.02516447000000000128
  ex, ey, ez, ea, eb, ec = ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # # Calculation 2
  # x, y, z, a, b, c = 0.00621957061202105688, -0.17473197000900078657, 0.60012728831046313616, 90.00280364221974593875, 90.03040263869252157747, -0.00041279850247499621
  # q1, q2, q3, q4, q5, q6 = 91.72786000000000683485, -83.77836000000000638011, 82.95619999999999549800, -1.72524699999999997502, 90.85218999999999311967, -0.02516447000000000128
  # ex, ey, ez, ea, eb, ec = -0.052750, 0.000000, 0.118100, 0.000000, 0.000000, 0.000000

  # # Calculation 3
  # x, y, z, a, b, c = -0.00624342358044325044, -0.25208865719984852038, 0.56246459269437065753, 89.99841969328601010147, 95.00647999671636512176, -0.00032959675925155326
  # q1, q2, q3, q4, q5, q6 = 90.00050000000000238742, -79.99165000000000702585, 84.99872999999999478860, -0.00207237000000000003, 89.99939999999999429292, -0.00014803000000000000
  # ex, ey, ez, ea, eb, ec = ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # # Calculation 4
  # x, y, z, a, b, c = -0.00624035510397932033, -0.12983581968671117912, 0.60470694360635102171, 89.99841969328601010147, 95.00647999671636512176, -0.00032959675925155326
  # q1, q2, q3, q4, q5, q6 = 90.00050000000000238742, -79.99165000000000702585, 84.99872999999999478860, -0.00207237000000000003, 89.99939999999999429292, -0.00014803000000000000
  # ex, ey, ez, ea, eb, ec = -0.052750, 0.000000, 0.118100, 0.000000, 0.000000, 0.000000

  # # Calculation 5
  x, y, z, a, b, c = -0.00597660090787822430, -0.23297860973622386282, 0.54091510737302239686, 89.99237127958376447623, 63.22621510689209856082, -0.01213066769093150421
  q1, q2, q3, q4, q5, q6 = 90.04352000000000089130, -91.37971000000000287855, 102.32599999999999340616, -0.05773026000000000546, 52.27993000000000023419, 0.00014803000000000000
  ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # # Calculation 6
  # x, y, z, a, b, c = -0.00597703117717071110, -0.15130123242096643610, 0.64121055007635763268, 89.99237132462658905752, 63.22631510685228306556, -0.01227879768204425889
  # q1, q2, q3, q4, q5, q6 = 90.04352000000000089130, -91.37960999999999955890, 102.32599999999999340616, -0.05773026000000000546, 52.27993000000000023419, 0.00000000000000000000
  # ex, ey, ez, ea, eb, ec = -0.052750, 0.000000, 0.118100, 0.000000, 0.000000, 0.000000
  # #ex, ey, ez, ea, eb, ec = -0.065500, 0.000000, 0.117500, 0.000000, 0.000000, 0.000000

  # # Roboter über der mitte

  # # Claculation 7
  # x, y, z, a, b, c = 0.36712896686322232798, -0.38852550549334541330, 0.51482615071429593812, 129.18533184969555804855, 39.81645892012360121726, 2.42908311224520678806
  # q1, q2, q3, q4, q5, q6 = 133.43289999999998940439, -92.94849000000000671662, 59.77634000000000469299, -2.84373399999999998400, 72.93094000000000676209, -0.00059210999999999997
  # ex, ey, ez, ea, eb, ec = 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # Claculation 8
  # x, y, z, a, b, c = 0.34665795811081900890, -0.35987401010406627755, 0.63928579121119410988, 129.18533184969555804855, 39.81645892012360121726, 2.42908311224520678806
  # q1, q2, q3, q4, q5, q6 = 133.43289999999998940439, -92.94849000000000671662, 59.77634000000000469299, -2.84373399999999998400, 72.93094000000000676209, -0.00059210999999999997
  # ex, ey, ez, ea, eb, ec = -0.05274999999999999828, 0.00000000000000000000, 0.11809999999999999665, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # Claculation 9
  # x, y, z, a, b, c = -0.00916987627183628133, 0.10843739655188605608, 0.53243025984828662445, 91.75991050690176109583, 90.74974599004613651232, 12.54158341021471834154
  # q1, q2, q3, q4, q5, q6 = 93.70431000000000665295, -55.74374999999999857891, 118.63970000000000482032, -4.15554299999999976478, 27.91657999999999972829, 16.24041000000000067871
  # ex, ey, ez, ea, eb, ec = 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # Claculation 10
  # x, y, z, a, b, c = -0.00136811476125670191, 0.22749682185081582242, 0.58237180833069346342, 91.75991050690176109583, 90.74974599004613651232, 12.54158341021471834154
  # q1, q2, q3, q4, q5, q6 = 93.70431000000000665295, -55.74374999999999857891, 118.63970000000000482032, -4.15554299999999976478, 27.91657999999999972829, 16.24041000000000067871
  # ex, ey, ez, ea, eb, ec = -0.05274999999999999828, 0.00000000000000000000, 0.11809999999999999665, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # Claculation 11 (Nice angles)
  # x, y, z, a, b, c = -0.00626517223557047139, -0.28993992890120567907, 0.50305293011553597893, 89.99679196999855435024, 90.00811999999997681243, -0.00029607402757792005
  # q1, q2, q3, q4, q5, q6 = 89.99693999999999505235, -89.99022999999999683496, 89.99952999999999292413, -0.00014803000000000000, 89.99881999999999493411, -0.00029605000000000001
  # ex, ey, ez, ea, eb, ec = 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # Claculation 12 (Nice angles)
  # x, y, z, a, b, c = -0.00625883190091027983, -0.17183245448832001889, 0.55578619236631765510, 89.99679196999855435024, 90.00811999999997681243, -0.00029607402757792005
  # q1, q2, q3, q4, q5, q6 = 89.99693999999999505235, -89.99022999999999683496, 89.99952999999999292413, -0.00014803000000000000, 89.99881999999999493411, -0.00029605000000000001
  # ex, ey, ez, ea, eb, ec = -0.05274999999999999828, 0.00000000000000000000, 0.11809999999999999665, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # Claculation 13 (Singularity)
  # x, y, z, a, b, c = -0.00619999143738556718, 0.10256245199084459974, 0.56531478715374294808, 89.99995118947381911312, 84.99698999887598915848, 0.00125486003076745369
  # q1, q2, q3, q4, q5, q6 = 89.99970000000000425189, -54.99857999999999691454, 114.99600000000000932232, 0.00059210999999999997, 24.99956999999999851525, 0.00074012999999999998
  # ex, ey, ez, ea, eb, ec = 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # Claculation 14 (Singularity)
  # x, y, z, a, b, c = -0.00619873982882447692, 0.21561227886185954650, 0.62816308977892942877, 89.99995118947381911312, 84.99698999887598915848, 0.00125486003076745369
  # q1, q2, q3, q4, q5, q6 = 89.99970000000000425189, -54.99857999999999691454, 114.99600000000000932232, 0.00059210999999999997, 24.99956999999999851525, 0.00074012999999999998
  # ex, ey, ez, ea, eb, ec = -0.05274999999999999828, 0.00000000000000000000, 0.11809999999999999665, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # # # Claculation 15 (Singularity)
  # x, y, z, a, b, c = -0.00620039839162776860, 0.10302216280906202617, 0.55477226839635329636, 90.00014407275156713695, 89.99647999701934963923, 0.00150927602607089885
  # q1, q2, q3, q4, q5, q6 = 89.99970000000000425189, -54.99857999999999691454, 114.99559999999999604370, 0.00088816000000000003, 29.99945999999999912689, 0.00074012999999999998
  # ex, ey, ez, ea, eb, ec = 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # # Claculation 16 (Singularity)
  # x, y, z, a, b, c = -0.00619930581926256448, 0.22111892185921194764, 0.60752952382778302454, 90.00014407275156713695, 89.99647999701934963923, 0.00150927602607089885
  # q1, q2, q3, q4, q5, q6 = 89.99970000000000425189, -54.99857999999999691454, 114.99559999999999604370, 0.00088816000000000003, 29.99945999999999912689, 0.00074012999999999998
  # ex, ey, ez, ea, eb, ec = -0.05274999999999999828, 0.00000000000000000000, 0.11809999999999999665, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000
  # #ex, ey, ez, ea, eb, ec = 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000










  # # Claculation 25
  x, y, z, a, b, c = 0.25717679234691975809, 0.35268800117509269132, 0.49040289872098646873, 68.01716188013186581429, 44.35993344881290312287, -3.82589720563539481546
  q1, q2, q3, q4, q5, q6 = 49.56877999999999673264, 2.15042799999999978411, 110.51659999999999683951, -13.58867000000000047066, -70.33765999999999962711, 14.23864999999999980673
  ex, ey, ez, ea, eb, ec = 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # # Claculation 26
  # x, y, z, a, b, c = 0.27073567778784979332, 0.39567917927792539690, 0.61163856892959400646, 68.01716188013186581429, 44.35993344881290312287, -3.82589720563539481546
  # q1, q2, q3, q4, q5, q6 = 49.56877999999999673264, 2.15042799999999978411, 110.51659999999999683951, -13.58867000000000047066, -70.33765999999999962711, 14.23864999999999980673
  # ex, ey, ez, ea, eb, ec = -0.05274999999999999828, 0.00000000000000000000, 0.11809999999999999665, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

















  # Robot movement

  # # Claculation 50
  # x, y, z, a, b, c = -0.00520236283541147728, 0.37612751728306470511, 0.51699737845966109440, 89.87802321880285205680, 90.01168052804928265687, 0.02578639325435625948
  # x, y, z, a, b, c = 0.00625000, 0.37600000, 0.51758940, 90.00000000, 90.00000000, 0.00000000
  # q1, q2, q3, q4, q5, q6 = 89.80559999999999831743, -34.17128000000000298542, 125.84520000000000550244, -2.49513100000000020984, -1.66381599999999996164, 2.51985200000000020282
  # ex, ey, ez, ea, eb, ec = -0.05274999999999999828, 0.00000000000000000000, 0.11809999999999999665, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # Claculation 51
  # x, y, z, a, b, c = 0.00619859793296640565, 0.38850605099275370424, 0.63636969881697891260, 90.00013014400718702746, 44.99590540961028750644, -0.00268843385874517355
  # x, y, z, a, b, c = 0.00620000, 0.38850000, 0.63638680, 90.00000000, 45.00000000, 0.00000000
  # q1, q2, q3, q4, q5, q6 = 87.23309000000000423825, -13.51667999999999913996, 120.03409999999999513420, -2.22439099999999978508, -61.57301999999999253532, 3.01411199999999990240
  # ex, ey, ez, ea, eb, ec = -0.05274999999999999828, 0.00000000000000000000, 0.11809999999999999665, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000

  # # Claculation 52
  # x, y, z, a, b, c = -0.00592167572153881258, 0.36386196706491724662, 0.53870674552562813897, 89.79285355604409346597, 89.99367089646374040512, 0.05603143098136382366
  # x, y, z, a, b, c = 0.00620000, 0.36350000, 0.55203630, 90.00000000, 90.00000000, 0.00000000
  # q1, q2, q3, q4, q5, q6 = 90.29278999999999655302, -36.19223999999999819011, 123.19580000000000552518, -9.49633200000000066154, 3.03157899999999980167, 9.53926000000000051671
  # ex, ey, ez, ea, eb, ec = -0.05274999999999999828, 0.00000000000000000000, 0.11809999999999999665, 0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000







  # Claculation 80
  x, y, z, a, b, c = 0.000000, (d2), (l1 + l2 + l3 + l4) - 0.0000001, 0.000000, 0.000000, 0.000000
  q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
  ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # Claculation 81
  x, y, z, a, b, c = (l4), (d2), (l1 + l2 + l3) - 0.0000001, 0.000000, 90.000000, 0.000000
  q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 0.000000, 0.000000, 90.000000, 0.000000
  ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # # Claculation 82
  # x, y, z, a, b, c = (l3 + l4), (d2), (l1 + l2), 0.000000, 0.000000, 0.000000
  # q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 90.000000, 0.000000, 90.000000, 0.000000
  # ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # # Claculation 83
  # x, y, z, a, b, c = (l2 + l3 + l4), (d2), (l1), 0.000000, 90.000000, 0.000000
  # q1, q2, q3, q4, q5, q6 = 0.000000, 90.000000, 0.000000, 0.000000, 90.000000, 0.000000
  # ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000




  # # Claculation 80_
  # x, y, z, a, b, c = (d2), 0.000000, (l1 + l2 + l3 + l4) - 0.0000001, 0.000000, 0.000000, 0.000000
  # q1, q2, q3, q4, q5, q6 = 90.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
  # ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # # Claculation 81_
  # x, y, z, a, b, c = (d2), (l4), (l1 + l2 + l3) - 0.0000001, 0.000000, -90.000000, 0.000000
  # q1, q2, q3, q4, q5, q6 = 90.000000, 0.000000, 0.000000, 0.000000, 90.000000, 90.000000
  # ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # # Claculation 82_
  # x, y, z, a, b, c = (d2), (l3 + l4), (l1 + l2), 0.000000, 90.000000, 0.000000
  # q1, q2, q3, q4, q5, q6 = 90.000000, 0.000000, 90.000000, 0.000000, 90.000000, 0.000000
  # ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # # Claculation 83_
  # x, y, z, a, b, c = (d2), (l2 + l3 + l4), (l1), 0.000000, 90.000000, 0.000000
  # q1, q2, q3, q4, q5, q6 = 90.000000, 90.000000, 0.000000, 0.000000, 90.000000, 0.000000
  # ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000






  # # # Claculation 84
  # x, y, z, a, b, c = 0.36799999999999999378, 0.00620000000000001193, 0.42500000000000004441, 45.00000000000000000000, 180.00000000000000000000, 45.00000000000000000000
  # q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 90.000000, 0.000000, 90.000000, 0.000000
  # ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

  # # # Claculation 85
  # x, y, z, a, b, c = 0.00620000000000003448, -0.36799999999999999378, 0.42500000000000004441, -44.99999999999999289457, 180.00000000000000000000, 45.00000000000000000000
  # q1, q2, q3, q4, q5, q6 = -90.000000, 0.000000, 90.000000, 0.000000, 90.00000000000000000000, 0.000000
  # ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000















  # # Calculation 1 (Test)
  # x = -0.00428667791374643704
  # y = 0.25883560684015277653
  # z = 0.46698387630450588492
  # roll = 89.87893434791007507556
  # pitch = 89.88231144120364035643
  # yaw = 0.03643662760103062032
  # q1, q2, q3, q4, q5, q6 = math.degrees(1.55840448910573692309), math.degrees(-0.59428461030406920518), math.degrees(2.18899194785128825558), math.degrees(-0.37699111843077520723), math.degrees(-0.02792526803190927345), math.degrees(0.37751471720637347351)






  q1, q2, q3, q4, q5, q6 = inversKinematic(x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, q1, q2, q3, q4, q5, q6)
  # q5 = q5 * (-1)
  # f = open(r"01_Calculation\Calculation.txt", "w")
  # forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
  # f.close()
  # print()
  # print()

  term = 80


  match term:
      case 1:
        f = open(r"01_Calculation\Calculation 1.txt", "w")
        print("Forward kinematic #1")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 91.72786000, -83.77836000, 82.95620000, -1.72524700, 90.85219000, -0.02516447
        # Result should be:
        x, y, z, a, b, c = 0.00627572, -0.29285840, 0.54744000, 90.00281000, 90.03040000, -0.00041287
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000, 0.000, 0.000, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 2:
        f = open(r"01_Calculation\Calculation 2.txt", "w")
        print("Forward kinematic #2")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 91.72786000, -83.77836000, 82.95620000, -1.72524700, 90.85219000, -0.02516447
        # Result should be:
        x, y, z, a, b, c = 0.00626955, -0.17473040, 0.60012730, 90.00281000, 90.03040000, -0.00041287
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 3:
        f = open(r"01_Calculation\Calculation 3.txt", "w")
        print("Forward kinematic #3")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 90.00050000, -79.99165000, 84.99873000, -0.00207237, 89.99940000, -0.00014803
        # Result should be:
        x, y, z, a, b, c = -0.00619345, -0.25208870, 0.56246460, 89.99841000, 95.00648000, -0.00032959
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000, 0.000, 0.000, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 4:
        f = open(r"01_Calculation\Calculation 4.txt", "w")
        print("Forward kinematic #4")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 90.00050000, -79.99165000, 84.99873000, -0.00207237, 89.99940000, -0.00014803
        # Result should be:
        x, y, z, a, b, c = -0.00619037, -0.12983580, 0.60470690, 89.99841000, 95.00648000, -0.00032959
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 5:
        f = open(r"01_Calculation\Calculation 5.txt", "w")
        print("Forward kinematic #5")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 90.04352000, -91.37971000, 102.32600000, -0.05773026, 52.27993000, 0.00014803
        # Result should be:
        x, y, z, a, b, c = -0.00592660, -0.23297870, 0.54091530, 89.99237000, 63.22620000, -0.01213065
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000, 0.000, 0.000, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 6:
        f = open(r"01_Calculation\Calculation 6.txt", "w")
        print("Forward kinematic #6")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 90.04352000, -91.37961000, 102.32600000, -0.05773026, 52.27993000, 0.00000000
        # Result should be:
        x, y, z, a, b, c = -0.00592689, -0.15130220, 0.64121040, 89.99237000, 63.22620000, -0.01213065
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 7:
        f = open(r"01_Calculation\Calculation 7.txt", "w")
        print("Forward kinematic #8")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 133.43290000, -92.94849000, 59.77634000, -2.84373400, 72.93094000, -0.00059211
        # Result should be:
        x, y, z, a, b, c = 0.36716570, -0.38849160, 0.51482530, 129.18530000, 39.81636000, 2.42909400
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000, 0.000, 0.000, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 8:
        f = open(r"01_Calculation\Calculation 8.txt", "w")
        print("Forward kinematic #8")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 133.43290000, -92.94849000, 59.77634000, -2.84373400, 72.93094000, -0.00059211
        # Result should be:
        x, y, z, a, b, c = 0.34669420, -0.35983960, 0.63928580, 129.18530000, 39.81646000, 2.42908200
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 9:
        f = open(r"01_Calculation\Calculation 9.txt", "w")
        print("Forward kinematic #9")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 93.70431000, -55.74375000, 118.63970000, -4.15554300, 27.91658000, 16.24041000
        # Result should be:
        x, y, z, a, b, c = -0.00911997, 0.10844050, 0.53243060, 91.75991000, 90.74970000, 12.54159000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000, 0.000, 0.000, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      # ...
      case 10:
        f = open(r"01_Calculation\Calculation 10.txt", "w")
        print("Forward kinematic #10")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 93.70431000, -55.74375000, 118.63970000, -4.15554300, 27.91658000, 16.24041000
        # Result should be:
        x, y, z, a, b, c = -0.00131821, 0.22749990, 0.58237220, 91.75991000, 90.74970000, 12.54159000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 11:
        f = open(r"01_Calculation\Calculation 11.txt", "w")
        print("Forward kinematic #11")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 89.99694000, -89.99023000, 89.99953000, -0.00014803, 89.99882000, -0.00029605
        # Result should be:
        x, y, z, a, b, c = -0.00621517, -0.28994000, 0.50305290, 89.99680000, 90.00811000, -0.00029608
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000, 0.000, 0.000, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 12:
        f = open(r"01_Calculation\Calculation 12.txt", "w")
        print("Forward kinematic #12")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 89.99694000, -89.99023000, 89.99953000, -0.00014803, 89.99882000, -0.00029605
        # Result should be:
        x, y, z, a, b, c = -0.00620883, -0.17183250, 0.55578630, 89.99680000, 90.00811000, -0.00029608
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 13:
        f = open(r"01_Calculation\Calculation 13.txt", "w")
        print("Forward kinematic #13")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 89.99970000, -54.99858000, 114.99600000, 0.00059211, 24.99957000, 0.00074013
        # Result should be:
        x, y, z, a, b, c = -0.00620000, 0.10256230, 0.56531510, 89.99995000, 84.99695000, 0.00125486
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000, 0.000, 0.000, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 14:
        f = open(r"01_Calculation\Calculation 14.txt", "w")
        print("Forward kinematic #14")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 89.99970000, -54.99858000, 114.99600000, 0.00059211, 24.99957000, 0.00074013
        # Result should be:
        x, y, z, a, b, c = -0.00619848, 0.21390240, 0.64081260, 89.99995000, 84.99695000, 0.00125486
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 15:
        f = open(r"01_Calculation\Calculation 15.txt", "w")
        print("Forward kinematic #15")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 89.99970000, -54.99858000, 114.99560000, 0.00088816, 29.99946000, 0.00074013
        # Result should be:
        x, y, z, a, b, c = -0.00620041, 0.10302210, 0.55477220, 90.00014000, 89.99648000, 0.00150928
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000, 0.000, 0.000, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 16:
        f = open(r"01_Calculation\Calculation 16.txt", "w")
        print("Forward kinematic #16")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 89.99970000, -54.99858000, 114.99560000, 0.00088816, 29.99946000, 0.00074013
        # Result should be:
        x, y, z, a, b, c = -0.00619898, 0.22051810, 0.62027950, 90.00014000, 89.99648000, 0.00150928
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()



      case 25:
        f = open(r"01_Calculation\Calculation 25.txt", "w")
        print("Forward kinematic #25")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 49.56878000, 2.15042800, 110.51660000, -13.58867000, -70.33766000, 14.23865000
        # Result should be:
        x, y, z, a, b, c = 0.25717680, 0.35268800, 0.49040320, 68.01717000, 44.35990000, -3.82591000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000, 0.000, 0.000, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 26:
        f = open(r"01_Calculation\Calculation 26.txt", "w")
        print("Forward kinematic #26")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 49.56878000, 2.15042800, 110.51660000, -13.58867000, -70.33766000, 14.23865000
        # Result should be:
        x, y, z, a, b, c = 0.27073560, 0.39567910, 0.61163880, 68.01717000, 44.35990000, -3.82591000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()




      # Robot moving V2
      case 50:
        f = open(r"01_Calculation\Calculation 50.txt", "w")
        print("Forward kinematic #50")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 89.80560000, -34.17128000, 125.84520000, -2.49513100, -1.66381600, 2.51985200
        # Result should be:
        x, y, z, a, b, c = 0.00625000, 0.37600000, 0.51758940, 90.00000000, 90.00000000, 0.00000000
        #x, y, z, a, b, c = -0.00520603, 0.37612730, 0.51699840, 89.87802000, 90.01141000, 0.02178904
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 51:
        f = open(r"01_Calculation\Calculation 51.txt", "w")
        print("Forward kinematic #51")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 87.23309000, -13.51668000, 120.03410000, -2.22439100, -61.57302000, 3.01411200
        # Result should be:
        x, y, z, a, b, c = 0.00620000, 0.38850000, 0.63638680, 90.00000000, 45.00000000, 0.00000000
        #x, y, z, a, b, c = 0.00619823, 0.38851090, 0.63637290, 90.00014000, 44.99681000, -0.00265480
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 52:
        f = open(r"01_Calculation\Calculation 52.txt", "w")
        print("Forward kinematic #52")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 90.29279000, -36.19224000, 123.19580000, -9.49633200, 3.03157900, 9.53926000
        # Result should be:
        x, y, z, a, b, c = 0.00620000, 0.36350000, 0.55203630, 90.00000000, 90.00000000, 0.00000000
        #x, y, z, a, b, c = -0.00518637, 0.37616280, 0.51703570, 89.87478000, 90.00785000, 0.01956835
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -0.052750, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()






      case 80:
        f = open(r"01_Calculation\Calculation 80.txt", "w")
        print("Forward kinematic #80")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 90.000000
        # Result should be:
        x, y, z, a, b, c = 0.000000, d2, (l1 + l2 + l3 + l4), 0.000000, 0.000000, 0.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 81:
        f = open(r"01_Calculation\Calculation 81.txt", "w")
        print("Forward kinematic #81")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 0.000000, 0.000000, 90.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = (l4), d2, (l1 + l2 + l3), 0.000000, 90.000000, 0.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 82:
        f = open(r"01_Calculation\Calculation 82.txt", "w")
        print("Forward kinematic #82")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 90.000000, 0.000000, 0.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = (l3 + l4), d2, (l1 + l2), 0.000000, 90.000000, 0.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 83:
        f = open(r"01_Calculation\Calculation 83.txt", "w")
        print("Forward kinematic #83")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 90.000000, 0.000000, 90.000000, 0.000000, 0.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = -d2, (l3 + l4), (l1 + l2), 90.000000, 90.000000, 0.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 84:
        f = open(r"01_Calculation\Calculation 84.txt", "w")
        print("Forward kinematic #84")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 90.000000, 0.000000, 0.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = (l3 + l4), d2, (l1 + l2), 0.000000, 90.000000, 0.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 85:
        f = open(r"01_Calculation\Calculation 85.txt", "w")
        print("Forward kinematic #85")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 90.000000, 0.000000, 90.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = l3, d2, (l1 + l2) - l4, 45.000000, 180.000000, 45.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 86:
        f = open(r"01_Calculation\Calculation 86.txt", "w")
        print("Forward kinematic #86")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = -90.000000, 0.000000, 90.000000, 0.000000, 90.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = d2, -l3, (l1 + l2) - l4, -45.000000, 180.000000, 45.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 87:
        f = open(r"01_Calculation\Calculation 87.txt", "w")
        print("Forward kinematic #87")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = -90.000000, 0.000000, 90.000000, 0.000000, 90.000000, 0.000000
        # Result should be:
        z = l1 + (sin(90 - q2)) * l2 + (cos(q3 - q2) * l3) - l4
        x, y, z, a, b, c = d2, -l3, z, -45.000000, 180.000000, 45.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()




      case 90:
        f = open(r"01_Calculation\Calculation 80.txt", "w")
        print("Forward kinematic #80")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = -a5, d2, (l1 + l2 + l3 + l4) + l5, 0.000000, 0.000000, 0.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -a5, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 91:
        f = open(r"01_Calculation\Calculation 91.txt", "w")
        print("Forward kinematic #91")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 90.000000, 0.000000, 0.000000, 90.000000, 0.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = -d2 + a5, 0.000000, (l1 + l2 + l3 + l4) + l5, 0.000000, 0.000000, 0.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -a5, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()
      case 92:
        f = open(r"01_Calculation\Calculation 92.txt", "w")
        print("Forward kinematic #92")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = 180.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = a5, -d2, (l1 + l2 + l3 + l4) + l5, 0.000000, 0.000000, 0.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -a5, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()



      case 95:
        f = open(r"01_Calculation\Calculation 95.txt", "w")
        print("Forward kinematic #95")
        # Angels are:
        q1, q2, q3, q4, q5, q6 = -90.000000, 0.000000, 90.000000, 0.000000, 90.000000, 0.000000
        # Result should be:
        x, y, z, a, b, c = d2, -l3 - a5, (l1 + l2) - l4 - l5, -45.000000, 180.000000, 45.000000
        # Endeffector are:
        ex, ey, ez, ea, eb, ec = -a5, 0.000, 0.11810, 0.000, 0.000, 0.000

        forwardKinematic(q1, q2, q3, q4, q5, q6, x, y, z, a, b, c, ex, ey, ez, ea, eb, ec, f)
        f.close()



if __name__=="__main__":
  main()
