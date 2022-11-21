from re import A
import rlcompleter
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

# Kuka KR210 (Original)
# robot_d1 = 0.75
# robot_d2 = 0.0
# robot_d3 = 0.0
# robot_d4 = 1.5
# robot_d5 = 0.0
# robot_d6 = 0.0
# robot_d7 = 0.303

# robot_a0 = 0.0
# robot_a1 = 0.35
# robot_a2 = 1.25
# robot_a3 = -0.054
# robot_a4 = 0.0
# robot_a5 = 0.0
# robot_a6 = 0.0

# robot_alpha0 = 0.0
# robot_alpha1 = -pi/2  # -90.0
# robot_alpha2 = 0.0
# robot_alpha3 = -pi/2  # -90.0
# robot_alpha4 = pi/2   # 90.0
# robot_alpha5 = -pi/2  # -90.0
# robot_alpha6 = 0.0

# robot_theta1 = 0.0 
# robot_theta2 = -pi/2 
# robot_theta3 = 0.0 
# robot_theta4 = 0.0
# robot_theta5 = 0.0
# robot_theta6 = 0.0
# robot_theta7 = 0.0

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

# # Excel (UR5e)
# robot_d1 = 0.1625
# robot_d2 = 0.0
# robot_d3 = 0.0
# robot_d4 = 0.1333
# robot_d5 = 0.0997
# robot_d6 = 0.0996
# robot_d7 = 0.0

# robot_a0 = 0.0
# robot_a1 = -0.4250
# robot_a2 = -0.3922
# robot_a3 = 0.0
# robot_a4 = 0.0
# robot_a5 = 0.0
# robot_a6 = 0.0

# robot_alpha0 = pi/2 
# robot_alpha1 = 0.0
# robot_alpha2 = 0.0
# robot_alpha3 = pi/2
# robot_alpha4 = -pi/2
# robot_alpha5 = 0.0
# robot_alpha6 = 0.0

# robot_theta1 = 0.0 
# robot_theta2 = 0.0
# robot_theta3 = 0.0 
# robot_theta4 = 0.0
# robot_theta5 = 0.0
# robot_theta6 = 0.0
# robot_theta7 = 0.0

# ...
robot_d1 = 0.1625
robot_d2 = 0.0
robot_d3 = 0.0
robot_d4 = 0.1333
robot_d5 = 0.0997
robot_d6 = 0.0996
robot_d7 = 0.0

robot_a1 = 0.0
robot_a2 = -0.4250
robot_a3 = -0.3922
robot_a4 = 0.0
robot_a5 = 0.0
robot_a6 = 0.0
robot_a7 = 0.0

robot_alpha1 = pi/2 
robot_alpha2 = 0.0
robot_alpha3 = 0.0
robot_alpha4 = pi/2
robot_alpha5 = -pi/2
robot_alpha6 = 0.0
robot_alpha7 = 0.0

robot_theta1 = 0.0 
robot_theta2 = 0.0
robot_theta3 = 0.0 
robot_theta4 = 0.0
robot_theta5 = 0.0
robot_theta6 = 0.0
robot_theta7 = 0.0





def pose3(theta, alpha, a, d):
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
    
  T = Matrix([
    [r11, r12, r13, r14],
    [r21, r22, r23, r24],
    [r31, r32, r33, r34],
    [0.0, 0.0, 0.0, 1.0]
  ])
  
  #T = simplify(T)

  return T

def forward_kin3(q1, q2, q3, q4, q5, q6):
  X = []
  Y = []
  Z = []

  # T01
  T01 = pose3(q1 + robot_theta1, robot_alpha1, robot_a1, robot_d1)
  print("T01: ", Matrix(T01))
  T0g = T01
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

   # T12
  T12 = pose3(q2 + robot_theta2, robot_alpha2, robot_a2, robot_d2)
  print("T12: ", Matrix(T12))
  T0g = T0g* T12
  px,py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T23
  T23 = pose3(q3 + robot_theta3, robot_alpha3, robot_a3, robot_d3)
  print("T23: ", Matrix(T23))
  T0g = T0g* T23
  px, py,pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T34
  T34 = pose3(q4 + robot_theta4, robot_alpha4, robot_a4, robot_d4)
  print("T34: ", Matrix(T34))
  T0g = T0g* T34
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T45
  T45 = pose3(q5 + robot_theta5, robot_alpha5, robot_a5, robot_d5)
  print("T45: ", Matrix(T45))
  T0g = T0g* T45
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

  # T56
  T56 = pose3(q6 + robot_theta6, robot_alpha6, robot_a6, robot_d6)
  print("T56: ", Matrix(T56))
  T0g = T0g* T56
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)
  
  # T6g (final position and rotation)
  T6g = pose3(0 + robot_theta7, robot_alpha7, robot_a7, robot_d7)
  print("T6g: ", Matrix(T6g))
  T0g = T0g* T6g
  print("T0g: ", Matrix(T0g))
  px, py, pz = T0g[0,3], T0g[1,3], T0g[2,3]
  X.append(px)
  Y.append(py)
  Z.append(pz)

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

  T13 = T12 * T23
  print("T13: ", Matrix(T13))

  W6 = T01 * T12 * T23 * T34 * T45 * T56
  print("W6: ", Matrix(W6))
  T0g = T01 * T12 * T23 * T34 * T45 * T56 * T6g  
  print("T0g: ", Matrix(T0g))

  U6 = T56
  U5 = T45 * U6
  U4 = T34 * U5
  U3 = T23 * U4
  U2 = T13 * U4
  U1 = T01 * U2

  pitch = atan2(T0g[1,3], sqrt(T0g[2,3] * T0g[2,3] + T0g[3,3] * T0g[3,3] ))

  roll = 0.0
  if (pitch != pi/2 and pitch != -(pi/2)):
    roll = atan2(-(T0g[2,3]/cos(pitch)), T0g[3,3]/cos(pitch))

  yaw = 0.0
  if (pitch == pi/2):
    yaw = atan2(T0g[3,2], T0g[2,2])
  elif (pitch == -pi/2):
    yaw = -atan2(T0g[3,2], T0g[2,2])
  else:
    yaw = atan2(-(T0g[1,2]/cos(pitch)), T0g[1,1]/cos(pitch))

  print("Forward kinematic 2:")
  print(Matrix(T0g))
  print("X : ", Matrix(X))
  print("Y : ", Matrix(Y))
  print("Z : ", Matrix(Z))
  print("Pitch : ", pitch, " ", pitch * 57.29577951307854999853275233034)
  print("Roll : ", roll, " ", roll * 57.29577951307854999853275233034)
  print("Yaw : ", yaw, " ", yaw * 57.29577951307854999853275233034)
  return X, Y, Z









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


#------------  Rotation Matrix and Euler Angles -----------
# Calculates Rotation Matrix given euler angles.
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
    
  # a1, a2, a3 = 0.35, 1.25, -0.054
  # d1, d4 = 0.75, 1.5
  # l = 1.50097168527591
  l = get_hypotenuse(robot_d4, abs(robot_a3))
  # phi = 1.53481186671284
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

  # Total transform of gripper wrt to base frame given orientation yaw (alpha), pitch (beta), roll (beta) and position px, py, pz
  T0g_b = Matrix([
    [1.0*sin(alpha)*sin(gamma) + sin(beta)*cos(alpha)*cos(gamma),  1.0*sin(alpha)*cos(gamma) - 1.0*sin(beta)*sin(gamma)*cos(alpha), 1.0*cos(alpha)*cos(beta), px],
    [sin(alpha)*sin(beta)*cos(gamma) - 1.0*sin(gamma)*cos(alpha), -1.0*sin(alpha)*sin(beta)*sin(gamma) - 1.0*cos(alpha)*cos(gamma), 1.0*sin(alpha)*cos(beta), py],
    [                                   1.0*cos(beta)*cos(gamma),                                        -1.0*sin(gamma)*cos(beta),           -1.0*sin(beta), pz],
    [                                                          0,                                                                0,                        0,  1]])

  # Total transform of gripper wrt to base frame given angles q1, q2, q3, q4, q5, q6
  T0g_a = Matrix([
    [((sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*cos(q5) + sin(q5)*cos(q1)*cos(q2 + q3))*cos(q6) - (-sin(q1)*cos(q4) + sin(q4)*sin(q2 + q3)*cos(q1))*sin(q6), -((sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*cos(q5) + sin(q5)*cos(q1)*cos(q2 + q3))*sin(q6) + (sin(q1)*cos(q4) - sin(q4)*sin(q2 + q3)*cos(q1))*cos(q6), -(sin(q1)*sin(q4) + sin(q2 + q3)*cos(q1)*cos(q4))*sin(q5) + cos(q1)*cos(q5)*cos(q2 + q3), -0.303*sin(q1)*sin(q4)*sin(q5) + 1.25*sin(q2)*cos(q1) - 0.303*sin(q5)*sin(q2 + q3)*cos(q1)*cos(q4) - 0.054*sin(q2 + q3)*cos(q1) + 0.303*cos(q1)*cos(q5)*cos(q2 + q3) + 1.5*cos(q1)*cos(q2 + q3) + 0.35*cos(q1)],
    [ ((sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*cos(q5) + sin(q1)*sin(q5)*cos(q2 + q3))*cos(q6) - (sin(q1)*sin(q4)*sin(q2 + q3) + cos(q1)*cos(q4))*sin(q6), -((sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*cos(q5) + sin(q1)*sin(q5)*cos(q2 + q3))*sin(q6) - (sin(q1)*sin(q4)*sin(q2 + q3) + cos(q1)*cos(q4))*cos(q6), -(sin(q1)*sin(q2 + q3)*cos(q4) - sin(q4)*cos(q1))*sin(q5) + sin(q1)*cos(q5)*cos(q2 + q3),  1.25*sin(q1)*sin(q2) - 0.303*sin(q1)*sin(q5)*sin(q2 + q3)*cos(q4) - 0.054*sin(q1)*sin(q2 + q3) + 0.303*sin(q1)*cos(q5)*cos(q2 + q3) + 1.5*sin(q1)*cos(q2 + q3) + 0.35*sin(q1) + 0.303*sin(q4)*sin(q5)*cos(q1)],
    [                                                                -(sin(q5)*sin(q2 + q3) - cos(q4)*cos(q5)*cos(q2 + q3))*cos(q6) - sin(q4)*sin(q6)*cos(q2 + q3),                                                                  (sin(q5)*sin(q2 + q3) - cos(q4)*cos(q5)*cos(q2 + q3))*sin(q6) - sin(q4)*cos(q6)*cos(q2 + q3),                                     -sin(q5)*cos(q4)*cos(q2 + q3) - sin(q2 + q3)*cos(q5),                                                                                 -0.303*sin(q5)*cos(q4)*cos(q2 + q3) - 0.303*sin(q2 + q3)*cos(q5) - 1.5*sin(q2 + q3) + 1.25*cos(q2) - 0.054*cos(q2 + q3) + 0.75],
    [                                                                                                                                                            0,                                                                                                                                                             0,                                                                                        0,                                                                                                                                                                                                              1]])

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
  
  print("Invers kinematic 1:")
  print("q1 +/-:", j1)
  print("q2:", j2)
  print("q3 +/-:", j3)
  print("q4:", j4)
  print("q5 +/-:", j5)
  print("q6:", j6)

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

# rotation matrices in z axes
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

  gripper_point = x, y, z

  ################################################################################
  # All important symbolic transformations matrices are declared below 
  ################################################################################

  q1, q2, q3, q4, q5, q6 = symbols('q1:7')
  alpha, beta, gamma = symbols('alpha beta gamma', real = True)
  px, py, pz = symbols('px py pz', real = True)

  # get the pose (homogenous transforms) of each joint wrt to previous joint
  # T01 = pose1(q1 + robot_theta1, robot_alpha1, robot_a1, robot_d1)
  # T12 = pose1(q2 + robot_theta2, robot_alpha2, robot_a2, robot_d2)
  # T23 = pose1(q3 + robot_theta3, robot_alpha3, robot_a3, robot_d3)
  # T34 = pose1(q4 + robot_theta4, robot_alpha4, robot_a4, robot_d4)
  # T45 = pose1(q5 + robot_theta5, robot_alpha5, robot_a5, robot_d5)
  # T56 = pose1(q6 + robot_theta6, robot_alpha6, robot_a6, robot_d6)
  # T6g = pose1(0 + robot_theta7,  robot_alpha7, robot_a7, robot_d7)
  T01 = pose2(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  T12 = pose2(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  T23 = pose2(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  T34 = pose2(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  T45 = pose2(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  T56 = pose2(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  T6g = pose2(0 + robot_theta7,  robot_alpha6, robot_a6, robot_d7)

  # From the poses, get the rotation of joint 3 wrt to the base frame and the transpose 
  # We will need this later
  T03 = simplify(T01 * T12 * T23)
  R03 = T03[:3, :3]
  R03T = R03.T
  
  # print("R03 = ")
  # print(R03)
  # print()

  # print("R03.T = ")
  # print(R03T)
  # print()

  # From the poses, get the rotation of joint 6 wrt to the joint 3
  # We will need this later 
  T36 = simplify(T34 * T45 * T56)
  R36 = T36[:3, :3]
  R36T = R36.T
  # print("R36 = ")
  # print(R36)
  # print()

  # the yaw, pitch roll is given wrt to the URDF frame 
  # We must convert this to gripper frame by performing
  # a rotation of 180 degrees ccw about the z axis and then 
  # a rotation of 90 degrees cw about the new y axis

  # This is the transpose of the rotation of the urdf frame wrt to gripper frame and its transpose
  # ( which is strangely the same) which is important later
  Rgu = (rotz(pi) * roty(-pi/2)).T
  RguT = Rgu.T
  # print(RguT)
  # print(Rgu == RguT)

  #  euler_R is the composite rotation matrix of the following
  # a rotation of alpha in the z axis
  # a rotation of beta in the new y axis
  # a rotation of gamma in the new x axis 
  # this will be useful later 
  euler_R = simplify(rotz(alpha) * roty(beta) * rotx(gamma))
  # print(euler_R)

  # This is the rotation of the gripper in URDF wrt to base frame 
  R0u_eval = euler_R.evalf(subs = {alpha: yaw, beta: pitch, gamma: roll})

  # R0g * Rgu = R0u 
  R0g_eval = R0u_eval * RguT
  # print(R0g_eval)

  # calculate wrist center
  wrist_center = get_wrist_center(gripper_point, R0g_eval, robot_d7)
  # print("wrist_center", wrist_center)

  # evaluated R0g
  # print("evaluated R0g:")
  # print(R0g_eval)

  j1, j2, j3 = get_first_three_angles(wrist_center)

  j1 = j1.evalf()
  j2 = j2.evalf()
  j3 = j3.evalf()

  R03T_eval = R03T.evalf(subs = {q1: j1, q2: j2, q3: j3})
  # print(R03T_eval)
  R36_eval = R03T_eval * R0g_eval

  j4, j5, j6 = get_last_three_angles(R36_eval)

  j4 = j4.evalf()
  j5 = j5.evalf()
  j6 = j6.evalf()

  print("Invers kinematic 2:")
  print("q1 +/-:", j1)
  print("q2:", j2)
  print("q3 +/-:", j3)
  print("q4:", j4)
  print("q5 +/-:", j5)
  print("q6:", j6)

  return j1, j2, j3, j4, j5, j6

def get_angles3(x, y, z, roll, pitch, yaw):
  # input: given position and orientation of the gripper_URDF wrt base frame
  # output: angles q1, q2, q3, q4, q5, q6

  gripper_point = x, y, z

  ################################################################################
  # All important symbolic transformations matrices are declared below 
  ################################################################################

  q1, q2, q3, q4, q5, q6 = symbols('q1:7')
  alpha, beta, gamma = symbols('alpha beta gamma', real = True)
  px, py, pz = symbols('px py pz', real = True)

  # get the pose (homogenous transforms) of each joint wrt to previous joint
  # T01 = pose1(q1 + robot_theta1, robot_alpha1, robot_a1, robot_d1)
  # T12 = pose1(q2 + robot_theta2, robot_alpha2, robot_a2, robot_d2)
  # T23 = pose1(q3 + robot_theta3, robot_alpha3, robot_a3, robot_d3)
  # T34 = pose1(q4 + robot_theta4, robot_alpha4, robot_a4, robot_d4)
  # T45 = pose1(q5 + robot_theta5, robot_alpha5, robot_a5, robot_d5)
  # T56 = pose1(q6 + robot_theta6, robot_alpha6, robot_a6, robot_d6)
  # T6g = pose1(0 + robot_theta7,  robot_alpha7, robot_a7, robot_d7)
  T01 = pose2(q1 + robot_theta1, robot_alpha0, robot_a0, robot_d1)
  T12 = pose2(q2 + robot_theta2, robot_alpha1, robot_a1, robot_d2)
  T23 = pose2(q3 + robot_theta3, robot_alpha2, robot_a2, robot_d3)
  T34 = pose2(q4 + robot_theta4, robot_alpha3, robot_a3, robot_d4)
  T45 = pose2(q5 + robot_theta5, robot_alpha4, robot_a4, robot_d5)
  T56 = pose2(q6 + robot_theta6, robot_alpha5, robot_a5, robot_d6)
  T6g = pose2(0 + robot_theta7,  robot_alpha6, robot_a6, robot_d7)

  T06 = T01 * T12 * T23 * T34 * T45 * T56
  print(Matrix(T06))

  P05 = T06 * [0, 0, -robot_a6, 1].T
  print(Matrix(P05))






  # # From the poses, get the rotation of joint 3 wrt to the base frame and the transpose 
  # # We will need this later
  # T03 = simplify(T01 * T12 * T23)
  # R03 = T03[:3, :3]
  # R03T = R03.T
  
  # # print("R03 = ")
  # # print(R03)
  # # print()

  # # print("R03.T = ")
  # # print(R03T)
  # # print()

  # # From the poses, get the rotation of joint 6 wrt to the joint 3
  # # We will need this later 
  # T36 = simplify(T34 * T45 * T56)
  # R36 = T36[:3, :3]
  # R36T = R36.T
  # # print("R36 = ")
  # # print(R36)
  # # print()

  # # the yaw, pitch roll is given wrt to the URDF frame 
  # # We must convert this to gripper frame by performing
  # # a rotation of 180 degrees ccw about the z axis and then 
  # # a rotation of 90 degrees cw about the new y axis

  # # This is the transpose of the rotation of the urdf frame wrt to gripper frame and its transpose
  # # ( which is strangely the same) which is important later
  # Rgu = (rotz(pi) * roty(-pi/2)).T
  # RguT = Rgu.T
  # # print(RguT)
  # # print(Rgu == RguT)

  # #  euler_R is the composite rotation matrix of the following
  # # a rotation of alpha in the z axis
  # # a rotation of beta in the new y axis
  # # a rotation of gamma in the new x axis 
  # # this will be useful later 
  # euler_R = simplify(rotz(alpha) * roty(beta) * rotx(gamma))
  # # print(euler_R)

  # # This is the rotation of the gripper in URDF wrt to base frame 
  # R0u_eval = euler_R.evalf(subs = {alpha: yaw, beta: pitch, gamma: roll})

  # # R0g * Rgu = R0u 
  # R0g_eval = R0u_eval * RguT
  # # print(R0g_eval)

  # # calculate wrist center
  # wrist_center = get_wrist_center(gripper_point, R0g_eval, robot_d7)
  # # print("wrist_center", wrist_center)

  # # evaluated R0g
  # # print("evaluated R0g:")
  # # print(R0g_eval)

  # j1, j2, j3 = get_first_three_angles(wrist_center)

  # j1 = j1.evalf()
  # j2 = j2.evalf()
  # j3 = j3.evalf()

  # R03T_eval = R03T.evalf(subs = {q1: j1, q2: j2, q3: j3})
  # # print(R03T_eval)
  # R36_eval = R03T_eval * R0g_eval

  # j4, j5, j6 = get_last_three_angles(R36_eval)

  # j4 = j4.evalf()
  # j5 = j5.evalf()
  # j6 = j6.evalf()

  # print("Invers kinematic 2:")
  # print("q1 +/-:", j1)
  # print("q2:", j2)
  # print("q3 +/-:", j3)
  # print("q4:", j4)
  # print("q5 +/-:", j5)
  # print("q6:", j6)

  return j1, j2, j3, j4, j5, j6




def pmMatQuatConvert(m, q):
  # from Stephe's "space" book e1 = (c32 - c23) / 4*e4 e2 = (c13 - c31) /
  # 4*e4 e3 = (c21 - c12) / 4*e4 e4 = sqrt(1 + c11 + c22 + c33) / 2

  # if e4 == 0 e1 = sqrt(1 + c11 - c33 - c22) / 2 e2 = sqrt(1 + c22 - c33
  # - c11) / 2 e3 = sqrt(1 + c33 - c11 - c22) / 2 to determine whether to
  # take the positive or negative sqrt value since e4 == 0 indicates a
  # 180* rotation then (0 x y z) = (0 -x -y -z). Thus some generallities
  # can be used: 1) find which of e1, e2, or e3 has the largest magnitude
  # and leave it pos. 2) if e1 is largest then if c21 < 0 then take the
  # negative for e2 if c31 < 0 then take the negative for e3 3) else if e2 
  # is largest then if c21 < 0 then take the negative for e1 if c32 < 0
  # then take the negative for e3 4) else if e3 is larger then if c31 < 0
  # then take the negative for e1 if c32 < 0 then take the negative for e2

  # Note: c21 in the space book is m->x.y in this C code

  a = 0.0

  # #ifdef PM_DEBUG
  #   if (!pmMatIsNorm(m)) {
  # #ifdef PM_PRINT_ERROR
  # 	pmPrintError("Bad matrix in pmMatQuatConvert\n");
  # #endif
  # 	return pmErrno = PM_NORM_ERR;
  #   }
  # #endif

    q->s = 0.5 * pmSqrt(1.0 + m->x.x + m->y.y + m->z.z);

    if (fabs(q->s) > QS_FUZZ) {
	q->x = (m->y.z - m->z.y) / (a = 4 * q->s);
	q->y = (m->z.x - m->x.z) / a;
	q->z = (m->x.y - m->y.x) / a;
    } else {
	q->s = 0;
	q->x = pmSqrt(1.0 + m->x.x - m->y.y - m->z.z) / 2.0;
	q->y = pmSqrt(1.0 + m->y.y - m->x.x - m->z.z) / 2.0;
	q->z = pmSqrt(1.0 + m->z.z - m->y.y - m->x.x) / 2.0;

	if (q->x > q->y && q->x > q->z) {
	    if (m->x.y < 0.0) {
		q->y *= -1;
	    }
	    if (m->x.z < 0.0) {
		q->z *= -1;
	    }
	} else if (q->y > q->z) {
	    if (m->x.y < 0.0) {
		q->x *= -1;
	    }
	    if (m->y.z < 0.0) {
		q->z *= -1;
	    }
	} else {
	    if (m->x.z < 0.0) {
		q->x *= -1;
	    }
	    if (m->y.z < 0.0) {
		q->y *= -1;
	    }
	}
    }

    pmErrno = PM_OK;
    return pmQuatNorm(q, q);

def pmHomPoseConvert(PmHomogeneous const * const h, PmPose * const p):
  r1 = 0.0

  p->tran = h->tran;
  r1 = pmMatQuatConvert(&h->rot, &p->rot);

  return r1




def pmQuatRpyConvert(PmQuaternion const * const q, PmRpy * const rpy):
    PmRotationMatrix m;
    r1, r2 = 0.0, 0.0

    #! \todo FIXME-- need direct equations */
    r1 = pmQuatMatConvert(q, &m);
    r2 = pmMatRpyConvert(&m, rpy);

    return pmErrno = (r1 || r2) ? PM_NORM_ERR : PM_OK;



def pumaKinematicsForward(joint):
  PUMA_A2 = 300.0
  PUMA_A3 = 50.0
  PUMA_D3 = 70.0
  PUMA_D4 = 400.0
  PUMA_D6 = 70.0

  SINGULAR_FUZZ = 0.000001
  FLAG_FUZZ = 0.000001

  # flags for inverse kinematics
  PUMA_SHOULDER_RIGHT = 1
  PUMA_ELBOW_DOWN = 2
  PUMA_WRIST_FLIP = 4
  PUMA_SINGULAR = 8  # joints at a singularity

  s1, s2, s3, s4, s5, s6 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  c1, c2, c3, c4, c5, c6 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  s23 = 0.0
  c23 = 0.0
  t1, t2, t3, t4, t5= 0.0, 0.0, 0.0, 0.0, 0.0
  # sumSq, k = 0.0
  # PmHomogeneous hom;
  # PmPose worldPose;
  # PmRpy rpy;

  # Calculate sin of joints for future use
  s1 = sin(joint[0] * pi/180)
  s2 = sin(joint[1] * pi/180)
  s3 = sin(joint[2] * pi/180)
  s4 = sin(joint[3] * pi/180)
  s5 = sin(joint[4] * pi/180)
  s6 = sin(joint[5] * pi/180)

  # Calculate cos of joints for future use
  c1 = cos(joint[0] * pi/180)
  c2 = cos(joint[1] * pi/180)
  c3 = cos(joint[2] * pi/180)
  c4 = cos(joint[3] * pi/180)
  c5 = cos(joint[4] * pi/180)
  c6 = cos(joint[5] * pi/180)

  s23 = c2 * s3 + s2 * c3
  c23 = c2 * c3 - s2 * s3

  # Calculate terms to be used in definition of... 
  # first column of rotation matrix.
  t1 = c4 * c5 * c6 - s4 * s6
  t2 = s23 * s5 * c6
  t3 = s4 * c5 * c6 + c4 * s6
  t4 = c23 * t1 - t2
  t5 = c23 * s5 * c6

  # Define first column of rotation matrix
  Rxx = c1 * t4 + s1 * t3
  Rxy = s1 * t4 - c1 * t3
  Rxz = -s23 * t1 - t5

  # Calculate terms to be used in definition of... 
  # second column of rotation matrix.
  t1 = -c4 * c5 * s6 - s4 * c6
  t2 = s23 * s5 * s6
  t3 = c4 * c6 - s4 * c5 * s6
  t4 = c23 * t1 + t2
  t5 = c23 * s5 * s6

  # Define second column of rotation matrix */
  Ryx= c1 * t4 + s1 * t3
  Ryy = s1 * t4 - c1 * t3
  Ryz = -s23 * t1 + t5

  # Calculate term to be used in definition of...
  # third column of rotation matrix.
  t1 = c23 * c4 * s5 + s23 * c5

  # Define third column of rotation matrix
  Rzx = -c1 * t1 - s1 * s4 * s5
  Rzy = -s1 * t1 + c1 * s4 * s5
  Rzz = s23 * c4 * s5 - c23 * c5

  # Calculate term to be used in definition of...
  # position vector.
  t1 = PUMA_A2 * c2 + PUMA_A3 * c23 - PUMA_D4 * s23

  # Define position vector
  Rtx= c1 * t1 - PUMA_D3 * s1
  Rty = s1 * t1 + PUMA_D3 * c1
  Rtz = -PUMA_A3 * s23 - PUMA_A2 * s2 - PUMA_D4 * c23

  # Calculate terms to be used to...
  # determine flags.
  sumSq = Rtx * Rtx + Rty * Rty - PUMA_D3 * PUMA_D3
  k = (sumSq + Rtz * Rtz - PUMA_A2 * PUMA_A2 - PUMA_A3 * PUMA_A3 - PUMA_D4 * PUMA_D4) / (2.0 * PUMA_A2)

  # reset flags
  iflags = 0

  # Set shoulder-up flag if necessary
  if (abs(joint[0]*pi/180 - atan2(Rty, Rtx) + atan2(PUMA_D3, -sqrt(sumSq))) < FLAG_FUZZ):
    iflags |= PUMA_SHOULDER_RIGHT

  # Set elbow down flag if necessary
  if (abs(joint[2]*pi/180 - atan2(PUMA_A3, PUMA_D4) + atan2(k, -sqrt(PUMA_A3 * PUMA_A3 + PUMA_D4 * PUMA_D4 - k * k))) < FLAG_FUZZ):
    iflags |= PUMA_ELBOW_DOWN

  # set singular flag if necessary
  t1 = -Rzx * s1 + Rzy * c1
  t2 = -Rzx * c1 * c23 - Rzy * s1 * c23 + Rzz * s23
  if (abs(t1) < SINGULAR_FUZZ and abs(t2) < SINGULAR_FUZZ):
    iflags |= PUMA_SINGULAR
   # if not singular set wrist flip flag if necessary 
  else:
    if (! (abs(joint[3]*PM_PI/180 - atan2(t1, t2)) < FLAG_FUZZ)):
      iflags |= PUMA_WRIST_FLIP

  #  add effect of d6 parameter
  Rtx = Rtx + Rzx*PUMA_D6
  Rty = Rty + Rzy*PUMA_D6
  Rtz = Rty + Rzz*PUMA_D6

  # convert hom.rot to world->quat
  pmHomPoseConvert(&hom, &worldPose);
  pmQuatRpyConvert(&worldPose.rot,&rpy);
  world->tran = worldPose.tran;
  world->a = rpy.r * 180.0/pi;
  world->b = rpy.p * 180.0/pi;
  world->c = rpy.y * 180.0/pi;




  return 








def main():
  #px, py, pz = 0.49792, 1.3673, 3.0000
  #px, py, pz = 0.276668856049955, -0.600333896425126, 0.512767653796484
  px, py, pz = 0.3253, -0.23290, 0.55470
  #px, py, pz = 0.49792, 1.3673, 2.4988
  #roll, pitch, yaw = 1, -1, -1
  #roll, pitch, yaw = 2.58309, 1.60570, 2.09440
  #roll, pitch, yaw = 0.366, -0.078, 2.561
  roll, pitch, yaw = 0.17453, 1.57080, 0.00000

  print("Input: ", px, " ", py, " ", pz, " ", roll, " ", pitch, " ", yaw)
  print()
  print()


  # print("get_angles1 ...")
  # q1, q2, q3, q4, q5, q6 = get_angles1(px, py, pz, roll, pitch, yaw)
  # forward_kin1(q1, q2, q3, q4, q5, q6)
  # forward_kin2(q1, q2, q3, q4, q5, q6)
  # print()
  # print()


  # print("get_angles2 ...")
  # q1, q2, q3, q4, q5, q6 = get_angles2(px, py, pz, roll, pitch, yaw)
  # forward_kin1(q1, q2, q3, q4, q5, q6)
  # forward_kin2(q1, q2, q3, q4, q5, q6)
  # print()
  # print()

  
  # print("get_angles3 ...")
  # q1, q2, q3, q4, q5, q6 = get_angles3(px, py, pz, roll, pitch, yaw)
  # forward_kin1(q1, q2, q3, q4, q5, q6)
  # forward_kin2(q1, q2, q3, q4, q5, q6)
  # print()
  # print()



  # forward_kin1(5.270894341, 3.316125579, 1.029744259, 3.473205211, 2.094395102, 1.570796327)
  # forward_kin2(5.270894341, 3.316125579, 1.029744259, 3.473205211, 2.094395102, 1.570796327)


  forward_kin3(0.0, 3.141592654, 1.570796327, 0.0, 0.0, 0.174532925)

if __name__=="__main__":
  main()