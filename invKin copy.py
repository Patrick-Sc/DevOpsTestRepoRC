from distutils.command.clean import clean
from sympy import symbols, cos, sin, pi, simplify, pprint, tan, expand_trig, sqrt, trigsimp, atan2
from sympy.matrices import Matrix

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
import numpy as np

global robot_d1, robot_d2, robot_d3, robot_d4, robot_d5, robot_d6, robot_d7
global robot_a1, robot_a2, robot_a3, robot_a4, robot_a5, robot_a6, robot_a7
global robot_alpha1, robot_alpha2, robot_alpha3, robot_alpha4, robot_alpha5, robot_alpha6, robot_alpha7
global robot_theta1, robot_theta2, robot_theta3, robot_theta4, robot_theta5, robot_theta6, robot_theta7

# Original
robot_d1 = 0.75
robot_d2 = 0.0
robot_d3 = 0.0
robot_d4 = 1.5
robot_d5 = 0.0
robot_d6 = 0.0
robot_d7 = 0.303

robot_a1 = 0.0
robot_a2 = 0.35
robot_a3 = 1.25
robot_a4 = -0.054
robot_a5 = 0.0
robot_a6 = 0.0
robot_a7 = 0.0

robot_alpha1 = 0.0
robot_alpha2 = -pi/2  # -90.0
robot_alpha3 = 0.0
robot_alpha4 = -pi/2  # -90.0
robot_alpha5 = pi/2   # 90.0
robot_alpha6 = -pi/2  # -90.0
robot_alpha7 = 0.0

robot_theta1 = 0.0 
robot_theta2 = -pi/2 
robot_theta3 = 0.0 
robot_theta4 = 0.0
robot_theta5 = 0.0
robot_theta6 = 0.0
robot_theta7 = 0.0

# Doosan
# robot_d1 = 0.0
# robot_d2 = 0.135
# robot_d3 = 0.00625
# robot_d4 = 0.0
# robot_d5 = 0.368
# robot_d6 = 0.0
# robot_d7 = 0.121

# robot_a1 = 0.0
# robot_a2 = 0.0
# robot_a3 = 0.441
# robot_a4 = 0.0
# robot_a5 = 0.0
# robot_a6 = 0.0
# robot_a7 = 0.0

# robot_alpha1 = 0.0
# robot_alpha2 = -pi/2
# robot_alpha3 = 0.0
# robot_alpha4 = pi/2
# robot_alpha5 = -pi/2
# robot_alpha6 = pi/2
# robot_alpha7 = 0.0

# robot_theta1 = 0.0 
# robot_theta2 = -pi/2 
# robot_theta3 = pi/2 
# robot_theta4 = 0.0
# robot_theta5 = 0.0
# robot_theta6 = 0.0
# robot_theta7 = 0.0

# PUMA 560
# robot_d1 = 0.67183
# robot_d2 = 0.13970
# robot_d3 = 0.0
# robot_d4 = 0.43180
# robot_d5 = 0.0
# robot_d6 = 0.0565
# robot_d7 = 0.0

# robot_a1 = 0.0
# robot_a2 = 0.4318
# robot_a3 = -0.2032
# robot_a4 = 0.0
# robot_a5 = 0.0
# robot_a6 = 0.0
# robot_a7 = 0.0

# robot_alpha1 = -pi/2 
# robot_alpha2 = 0.0
# robot_alpha3 = pi/2 
# robot_alpha4 = -pi/2
# robot_alpha5 = pi/2
# robot_alpha6 = 0.0
# robot_alpha7 = 0.0

# robot_theta1 = 0.0 
# robot_theta2 = -pi/2 
# robot_theta3 = pi/2 
# robot_theta4 = 0.0
# robot_theta5 = 0.0
# robot_theta6 = 0.0
# robot_theta7 = 0.0

# # Excel (UR5)
# robot_d1 = 0.1625
# robot_d2 = 0.0
# robot_d3 = 0.0
# robot_d4 = 0.1333
# robot_d5 = 0.0997
# robot_d6 = 0.0996
# robot_d7 = 0.0

# robot_a1 = 0.0
# robot_a2 = -0.4250
# robot_a3 = -0.3922
# robot_a4 = 0.0
# robot_a5 = 0.0
# robot_a6 = 0.0
# robot_a7 = 0.0

# robot_alpha1 = pi/2 
# robot_alpha2 = 0.0
# robot_alpha3 = 0.0
# robot_alpha4 = pi/2
# robot_alpha5 = -pi/2
# robot_alpha6 = 0.0
# robot_alpha7 = 0.0

# robot_theta1 = 0.0 
# robot_theta2 = 0.0
# robot_theta3 = 0.0 
# robot_theta4 = 0.0
# robot_theta5 = 0.0
# robot_theta6 = 0.0
# robot_theta7 = 0.0

# Pose
def pose(theta, alpha, a, d):
  # returns the pose T of one joint frame i with respect to the previous joint frame (i - 1)
  # given the parameters:
  # theta: theta[i]
  # alpha: alpha[i-1]
  # a: a[i-1]
  # d: d[i]

  # r11 = cos(theta)
  # r12 = -sin(theta)
  # r21 = sin(theta) * cos(alpha)
  # r22 = cos(theta) * cos(alpha)
  # r23 = -sin(alpha)
  # r31 = sin(theta) * sin(alpha)
  # r32 = cos(theta) * sin(alpha)
  # r33 = cos(alpha)
  # y = -d * sin(alpha)
  # z = d * cos(alpha)
    
  # T = Matrix([
  #   [r11, r12, 0.0, a],
  #   [r21, r22, r23, y],
  #   [r31, r32, r33, z],
  #   [0.0, 0.0, 0.0, 1.0]
  # ])
  
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
    
  T = Matrix([
    [r11, r12, r13, r14],
    [r21, r22, r23, r24],
    [r31, r32, r33, r34],
    [0.0, 0.0, 0.0, 1.0]
  ])

  T = simplify(T)

  # print("Transformation for each link")
  # print("T : ", T)

  return T

# Forward kinematic
def forward_kin(q1, q2, q3, q4, q5, q6,):
  X = []
  Y = []
  Z = []

  # T01
  T01 = pose(q1 + robot_theta1, robot_alpha1, robot_a1, robot_d1)   # pose(theta[i], alpha[i-1], a[i-1], d[i]) 
  T0g = T01
  px,py,pz = T0g[0,3], T0g[1,3], T0g[2,3]  
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T12
  T12 = pose(q2 + robot_theta2, robot_alpha2, robot_a2, robot_d2)   # pose(theta[i], alpha[i-1], a[i-1], d[i])
  T0g = T0g* T12
  px,py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T23
  T23 = pose(q3 + robot_theta3, robot_alpha3, robot_a3, robot_d3)   # pose(theta[i], alpha[i-1], a[i-1], d[i])
  T0g = T0g* T23
  px,py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T34
  T34 = pose(q4 + robot_theta4, robot_alpha4, robot_a4, robot_d4)   # pose(theta[i], alpha[i-1], a[i-1], d[i])
  T0g = T0g* T34
  px,py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T45
  T45 = pose(q5 + robot_theta5, robot_alpha5, robot_a5, robot_d5)   # pose(theta[i], alpha[i-1], a[i-1], d[i])
  T0g = T0g* T45
  px,py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T56
  T56 = pose(q6 + robot_theta6, robot_alpha6, robot_a6, robot_d6)   # pose(theta[i], alpha[i-1], a[i-1], d[i])
  T0g = T0g* T56
  px,py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T6g (final position and rotation)
  T6g = pose(0 + robot_theta7, robot_alpha7, robot_a7, robot_d7)   # pose(theta[i], alpha[i-1], a[i-1], d[i])
  T0g = T0g* T6g
  px,py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # fig = plt.figure()
  # ax = fig.add_subplot(111,projection = '3d')
  # ax.set_xlabel('x axis')
  # ax.set_ylabel('y axis')
  # ax.set_zlabel('z axis')
  
  X = np.reshape(X, (1,7))
  Y = np.reshape(Y, (1,7))
  Z = np.reshape(Z, (1,7))
  
  # ax.cla()
  # ax.plot_wireframe(X,Y,Z)
  # plt.draw()
  # plt.pause(3)
  # ax.cla()
  # ax.plot_wireframe(Z,Y,X,color='r')
  
  # plt.show()

  print("Forward kinematic")
  print("X : ", Matrix(X))
  print("Y : ", Matrix(Y))
  print("Z : ", Matrix(Z))
  return X,Y,Z


def create_plot():
  fig = plt.figure()
  ax = fig.add_subplot(111,projection = '3d')
  ax.set_xlabel('x axis')
  ax.set_ylabel('y axis')
  ax.set_zlabel('z axis')
  ax.set_xlim3d([0, 2])
  ax.set_ylim3d([0, 3])
  ax.set_zlim3d([0, 3])
  ax.set_autoscale_on(False)
  return fig,ax

def update_plot(X,Y,Z,fig,ax):
  X = np.reshape(X,(1,7))
  Y = np.reshape(Y,(1,7))
  Z = np.reshape(Z,(1,7))
  ax.cla()
  ax.plot_wireframe(X,Y,Z)
  #plt.draw()
  ax.set_xlabel('x axis')
  ax.set_ylabel('y axis')
  ax.set_zlabel('z axis')
  ax.set_xlim3d([0, 2])
  ax.set_ylim3d([0, 3])
  ax.set_zlim3d([0, 3])
  ax.set_autoscale_on(False)
  fig.canvas.draw()
  fig.canvas.flush_events()
  #plt.pause(3)
  #ax.cla()
  #ax.plot_wireframe(Z,Y,X,color='r')

# Calculates Rotation Matrix given euler angles (Rotation Matrix and Euler Angles)
def eulerAnglesToRotationMatrix(theta) :     
    R_x = np.array([[1,         0,             0              ],
                    [0,         cos(theta[0]), -sin(theta[0]) ],
                    [0,         sin(theta[0]), cos(theta[0])  ]
                    ])   

    R_y = np.array([[cos(theta[1]),    0,      sin(theta[1])  ],
                    [0,                1,      0              ],
                    [-sin(theta[1]),   0,      cos(theta[1])  ]
                    ])
                 
    R_z = np.array([[cos(theta[2]),    -sin(theta[2]),    0],
                    [sin(theta[2]),    cos(theta[2]),     0],
                    [0,                0,                 1]
                    ])                     
                     
    R = np.dot(R_z, np.dot( R_y, R_x ))
 
    return R

# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6

# Calculates rotation matrix to euler angles
# The result is the same as MATLAB except the order
# of the euler angles ( x and z are swapped ).
def rotationMatrixToEulerAngles(R) :
 
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
#---------------------------------------


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
    
  #a1, a2, a3 = 0.35, 1.25, -0.054
  #d1, d4 = 0.75, 1.5
  #l = 1.50097168527591 
  l = get_hypotenuse(robot_d4, abs(robot_a3))
  print("l : ",l)
  #phi = 1.53481186671284 
  phi = atan2(robot_d4, abs(robot_a3))
  print("phi : ",phi)
  
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
  #Matrix([
  #[-sin(q4)*sin(q6) + cos(q4)*cos(q5)*cos(q6), -sin(q4)*cos(q6) - sin(q6)*cos(q4)*cos(q5), -sin(q5)*cos(q4)],
  #[                           sin(q5)*cos(q6),                           -sin(q5)*sin(q6),          cos(q5)],
  #[-sin(q4)*cos(q5)*cos(q6) - sin(q6)*cos(q4),  sin(q4)*sin(q6)*cos(q5) - cos(q4)*cos(q6),  sin(q4)*sin(q5)]])
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

def get_angles1(x, y, z, roll, pitch, yaw):
  # input: given position and orientation of the gripper_URDF wrt base frame
  # output: angles q1, q2, q3, q4, q5, q6
    
  gripper_point = x, y, z

  ################################################################################
  # All important symbolic transformations matrices are declared below 
  ################################################################################

  q1, q2, q3, q4, q5, q6 = symbols('q1:7')
  alpha, beta, gamma = symbols('alpha beta gamma', real = True)
  px, py, pz = symbols('px py pz', real = True)

  # Rotation of joint 3 wrt to the base frame interms the first three angles q1, q2, q3
  R03 = Matrix([
    [sin(q2 + q3)*cos(q1), cos(q1)*cos(q2 + q3), -sin(q1)],
    [sin(q1)*sin(q2 + q3), sin(q1)*cos(q2 + q3),  cos(q1)],
    [        cos(q2 + q3),        -sin(q2 + q3),        0]])

  # Transpose of R03 
  R03T = Matrix([
    [sin(q2 + q3)*cos(q1), sin(q1)*sin(q2 + q3),  cos(q2 + q3)],
    [cos(q1)*cos(q2 + q3), sin(q1)*cos(q2 + q3), -sin(q2 + q3)],
    [            -sin(q1),              cos(q1),             0]])

  # Rotation of joint 6 wrt to frame of joint 3 interms of the last three angles q4, q5, q6
  R36 = Matrix([
    [-sin(q4)*sin(q6) + cos(q4)*cos(q5)*cos(q6), -sin(q4)*cos(q6) - sin(q6)*cos(q4)*cos(q5), -sin(q5)*cos(q4)],
    [                           sin(q5)*cos(q6),                           -sin(q5)*sin(q6),          cos(q5)],
    [-sin(q4)*cos(q5)*cos(q6) - sin(q6)*cos(q4),  sin(q4)*sin(q6)*cos(q5) - cos(q4)*cos(q6),  sin(q4)*sin(q5)]])

  # Rotation of urdf_gripper with respect to the base frame interms of alpha = yaw, beta = pitch, gamma = roll
  R0u = Matrix([
    [1.0*cos(alpha)*cos(beta), -1.0*sin(alpha)*cos(gamma) + sin(beta)*sin(gamma)*cos(alpha), 1.0*sin(alpha)*sin(gamma) + sin(beta)*cos(alpha)*cos(gamma)],
    [1.0*sin(alpha)*cos(beta),  sin(alpha)*sin(beta)*sin(gamma) + 1.0*cos(alpha)*cos(gamma), sin(alpha)*sin(beta)*cos(gamma) - 1.0*sin(gamma)*cos(alpha)],
    [          -1.0*sin(beta),                                     1.0*sin(gamma)*cos(beta),                                    1.0*cos(beta)*cos(gamma)]])

  # Total transform of gripper wrt to base frame given orientation yaw (alpha), pitch (beta), roll (gamma) and position px, py, pz
  T0g_b = Matrix([
    [1.0*sin(alpha)*sin(gamma) + sin(beta)*cos(alpha)*cos(gamma),  1.0*sin(alpha)*cos(gamma) - 1.0*sin(beta)*sin(gamma)*cos(alpha), 1.0*cos(alpha)*cos(beta), px],
    [sin(alpha)*sin(beta)*cos(gamma) - 1.0*sin(gamma)*cos(alpha), -1.0*sin(alpha)*sin(beta)*sin(gamma) - 1.0*cos(alpha)*cos(gamma), 1.0*sin(alpha)*cos(beta), py],
    [                                   1.0*cos(beta)*cos(gamma),                                        -1.0*sin(gamma)*cos(beta),           -1.0*sin(beta), pz],
    [                                                          0,                                                                0,                        0,  1]])

  # Total transform of gripper wrt to base frame given angles q1, q2, q3, q4, q5, q6
#  T0g_a = Matrix([
#    [((sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*cos(q5) + sin(q5)*cos(q1)*cos(q2 + q3))*cos(q6) - (-sin(q1)*cos(q4) + sin(q4)*sin(q2 + q3)*cos(q1))*sin(q6), -((sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*cos(q5) + sin(q5)*cos(q1)*cos(q2 + q3))*sin(q6) + (sin(q1)*cos(q4) - sin(q4)*sin(q2 + q3)*cos(q1))*cos(q6), -(sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*sin(q5) + cos(q1)*cos(q5)*cos(q2 + q3), -0.303*sin(q1)*sin(q4)*sin(q5) + 1.25*sin(q2)*cos(q1) - 0.303*sin(q5)*sin(q2 + q3)*cos(q1)*cos(q4) - 0.054*sin(q2 + q3)*cos(q1) + 0.303*cos(q1)*cos(q5)*cos(q2 + q3) + 1.5*cos(q1)*cos(q2 + q3) + 0.35*cos(q1)],
#    [ ((sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*cos(q5) + sin(q1)*sin(q5)*cos(q2 + q3))*cos(q6) - (sin(q1)*sin(q4)*sin(q2 + q3) + cos(q1)*cos(q4))*sin(q6), -((sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*cos(q5) + sin(q1)*sin(q5)*cos(q2 + q3))*sin(q6) - (sin(q1)*sin(q4)*sin(q2 + q3) + cos(q1)*cos(q4))*cos(q6), -(sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*sin(q5) + sin(q1)*cos(q5)*cos(q2 + q3),  1.25*sin(q1)*sin(q2) - 0.303*sin(q1)*sin(q5)*sin(q2 + q3)*cos(q4) - 0.054*sin(q1)*sin(q2 + q3) + 0.303*sin(q1)*cos(q5)*cos(q2 + q3) + 1.5*sin(q1)*cos(q2 + q3) + 0.35*sin(q1) + 0.303*sin(q4)*sin(q5)*cos(q1)],
#    [                                                                -(sin(q5)*sin(q2 + q3) - cos(q4)*cos(q5)*cos(q2 + q3))*cos(q6) - sin(q4)*sin(q6)*cos(q2 + q3),                                                                  (sin(q5)*sin(q2 + q3) - cos(q4)*cos(q5)*cos(q2 + q3))*sin(q6) - sin(q4)*cos(q6)*cos(q2 + q3),                                     -sin(q5)*cos(q4)*cos(q2 + q3) - sin(q2 + q3)*cos(q5),                                                                                 -0.303*sin(q5)*cos(q4)*cos(q2 + q3) - 0.303*sin(q2 + q3)*cos(q5) - 1.5*sin(q2 + q3) + 1.25*cos(q2) - 0.054*cos(q2 + q3) + 0.75],
#    [                                                                                                                                                            0,                                                                                                                                                             0,                                                                                        0,                                                                                                                                                                                                              1]])

  # Rotation of urdf_gripper wrt (DH) gripper frame from rotz(pi) * roty(-pi/2) and it's transpose
  Rgu_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])
  RguT_eval = Matrix([[0, 0, 1], [0, -1.00000000000000, 0], [1.00000000000000, 0, 0]])

  # Inverse kinematics transformations starts below
  R0u_eval = R0u.evalf(subs = {alpha: yaw, beta: pitch, gamma: roll})
  R0g_eval = R0u_eval * RguT_eval

  wrist_center = get_wrist_center(gripper_point, R0g_eval, robot_d7)

  j1, j2, j3 = get_first_three_angles(wrist_center)

  R03T_eval = R03T.evalf(subs = {q1: j1.evalf(), q2: j2.evalf(), q3: j3.evalf()})
  R36_eval = R03T_eval * R0g_eval

  j4, j5, j6 = get_last_three_angles(R36_eval)

  j1 = j1.evalf()
  j2 = j2.evalf()
  j3 = j3.evalf()
  j4 = j4.evalf()
  j5 = j5.evalf()
  j6 = j6.evalf()

  print("Inverse kinematic")
  print("q1 +/-: ", j1)
  print("q2 : ", j2)
  print("q3 +/-: ", j3)
  print("q4 : ", j4)
  print("q5 +/-: ", j5)
  print("q6 : ", j6)

  return j1, j2, j3, j4, j5, j6














# rotation matrices in x axes
def rotx(q):
  sq, cq = sin(q), cos(q)

  r = Matrix([
    [1., 0., 0.],
    [0., cq,-sq],
    [0., sq, cq]
  ])
    
  return r

# rotation matrices in y axes
def roty(q):
  sq, cq = sin(q), cos(q)

  r = Matrix([
    [ cq, 0., sq],
    [ 0., 1., 0.],
    [-sq, 0., cq]
  ])
    
  return r

# rotation matrices in  z axes
def rotz(q):
  sq, cq = sin(q), cos(q)

  r = Matrix([
    [cq,-sq, 0.],
    [sq, cq, 0.],
    [0., 0., 1.]
  ])
    
  return r

def get_angles2(x, y, z, roll, pitch, yaw):
  # input: given position and orientation of the gripper_URDF wrt base frame
  # output: angles q1, q2, q3, q4, q5, q6

  

    # This is given position and orientation of the gripper wrt to URDFrame

  #roll, pitch, yaw = 0.366, -0.078, 2.561

  gripper_point = x, y, z
  
  
  
  
  #P5 = [ x - d6 * nx, y - d6 * ny, z - d6 * nz ]











    
  

  ################################################################################
  # All important symbolic transformations matrices are declared below 
  ################################################################################

  q1, q2, q3, q4, q5, q6 = symbols('q1:7')
  alpha, beta, gamma = symbols('alpha beta gamma', real = True)
  px, py, pz = symbols('px py pz', real = True)

  # R01 = Matrix([
  #   [cos(q1),   -sin(q1) * cos(a1),   sin(q1) * sin(a1)],
  #   [sin(q1),   cos(q1) * cos(a1),    -cos(q1) * sin(a1)],
  #   [0,         sin(a1),              cos(a1)]
  # ])



  d1 = robot_d1
  d2 = robot_d2
  d3 = robot_d3
  d4 = robot_d4
  d5 = robot_d5
  d6 = robot_d6
  d7 = robot_d7
  
  a1 = robot_a1
  a2 = robot_a2
  a3 = robot_a3
  a4 = robot_a4
  a5 = robot_a5
  a6 = robot_a6
  a7 = robot_a7
  
  alpha1 = robot_alpha1
  alpha2 = robot_alpha2
  alpha3 = robot_alpha3
  alpha4 = robot_alpha4
  alpha5 = robot_alpha5
  alpha6 = robot_alpha6
  alpha7 = robot_alpha7

  theta1 = robot_theta1
  theta2 = robot_theta2
  theta3 = robot_theta3
  theta4 = robot_theta4
  theta5 = robot_theta5
  theta6 = robot_theta6
  theta7 = robot_theta7


  # get the pose (homogenous transforms) of each joint wrt to previous joint
  T01 = pose(q1 + theta1, alpha1, a1, d1)
  T12 = pose(q2 + theta2, alpha2, a2, d2)
  T23 = pose(q3 + theta3, alpha3, a3, d3)
  T34 = pose(q4 + theta4, alpha4, a4, d4)
  T45 = pose(q5 + theta5, alpha5, a5, d5)
  T56 = pose(q6 + theta6, alpha6, a6, d6)
  T6g = pose(0 + theta7,  alpha7, a7, d7)

  # From the poses, get the rotation of joint 3 wrt to the base frame and the transpose 
  # We will need this later
  T03 = simplify(T01 * T12 * T23)
  R03 = T03[:3, :3]
  R03T = R03.T
  
  print("R03 = ")
  print(R03)
  print()

  print("R03.T = ")
  print(R03T)
  print()

  # From the poses, get the rotation of joint 6 wrt to the joint 3
  # We will need this later 
  T36 = simplify(T34 * T45 * T56)
  R36 = T36[:3, :3]

  print("R36 = ")
  print(R36)
  print()



  # the yaw, pitch roll is given wrt to the URDF frame 
  # We must convert this to gripper frame by performing
  # a rotation of 180 degrees ccw about the z axis and then 
  # a rotation of 90 degrees cw about the new y axis

  # This is the transpose of the rotation of the urdf frame wrt to gripper frame and its transpose
  # ( which is strangely the same) which is important later

  Rgu = (rotz(pi) * roty(-pi/2)).T
  RguT = Rgu.T
  print("RguT = ")
  print(RguT)
  print()
  print("Rgu == RguT = ")
  print(Rgu == RguT)
  print()

  #  euler_R is the composite rotation matrix of the following
  # a rotation of alpha in the z axis
  # a rotation of beta in the new y axis
  # a rotation of gamma in the new x axis 
  # this will be useful later 

  #alpha, beta, gamma = symbols('alpha beta gamma', real = True)
  euler_R = simplify(rotz(alpha) * roty(beta) * rotx(gamma))
  print("Euler R = ")
  print(euler_R)
  print()


  # This is the rotation of the gripper in URDF wrt to base frame 
  R0u_eval = euler_R.evalf(subs = {alpha: yaw, beta: pitch, gamma: roll})

  # R0g * Rgu = R0u 
  R0g_eval = R0u_eval * RguT

  # calculate wrist center
  wrist_center = get_wrist_center(gripper_point, R0g_eval, dg = d7)
  print("wrist_center =")
  print(wrist_center)
  print()

  # evaluated R0g
  print("evaluated R0g =")
  print(R0g_eval)
  print()

  j1, j2, j3 = get_first_three_angles(wrist_center)
  #j1, j2, j3 = get_first_three_angles(wrist_center, robot_d1, robot_d2, robot_d3, robot_d4, robot_d5, robot_d6, robot_d7, robot_a1, robot_a2, robot_a3, robot_a4, robot_a5, robot_a6, robot_a7, robot_alpha1, robot_alpha2, robot_alpha3, robot_alpha4, robot_alpha5, robot_alpha6, robot_alpha7)

  R03T_eval = R03T.evalf(subs = {q1: j1.evalf(), q2: j2.evalf(), q3: j3.evalf()})
  R36_eval = R03T_eval * R0g_eval

  j4, j5, j6 = get_last_three_angles(R36_eval)

  j1 = j1.evalf()
  j2 = j2.evalf()
  j3 = j3.evalf()
  j4 = j4.evalf()
  j5 = j5.evalf()
  j6 = j6.evalf()

  return j1, j2, j3, j4, j5, j6

def main():
  px, py, pz = 0.49792, 1.3673, 3.0000
  #px, py, pz = 0.276668856049955, -0.600333896425126, 0.512767653796484
  #px, py, pz = 0.49792, 1.3673, 2.4988
  roll, pitch, yaw = 1, -1, -1
  #roll, pitch, yaw =  1.60570, 2.58309, 2.09440
  #roll, pitch, yaw = 0.366, -0.078, 2.561

  q1, q2, q3, q4, q5, q6 = get_angles1(px, py, pz, roll, pitch, yaw)
  forward_kin(q1, q2, q3, q4, q5, q6)
  print()

  # q1, q2, q3, q4, q5, q6 = get_angles2(px, py, pz, roll, pitch, yaw)
  # forward_kin(q1, q2, q3, q4, q5, q6)
  # print()

  # now plot the arm with updated angles
  # forward_kin(q1, q2, q3, q4, q5, q6)
  # print()
  # forward_kin(0, 0, 0, 0, 0, 0)
  # print()
  # forward_kin(5.270894341, 3.316125579, 1.029744259, 3.473205211, 2.094395102, 1.570796327)

if __name__=="__main__":
  main()
