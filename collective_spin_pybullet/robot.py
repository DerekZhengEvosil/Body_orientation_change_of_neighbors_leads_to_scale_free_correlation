from collections import deque
import numpy as np
import pybullet as p
import itertools
import math

class Robot():
    def __init__(self, init_pos, init_vel, robot_id, dt):
        self.id = robot_id
        self.dt = dt

        startPos = init_pos
        startOrientation = p.getQuaternionFromEuler(init_vel)
        self.pybullet_id = p.loadURDF("../models/swmrob_v26.urdf", startPos, startOrientation)
        self.wheelRadius = 0.0175  # m
        self.leftWheelJointIdx = 2
        self.rightWheelJointIdx = 1

        self.joint_ids = list(range(p.getNumJoints(self.pybullet_id)))
        self.initial_position = startPos
        self.initial_Orientation = startOrientation
        self.reset()
        self.tail = []
        self.past_posvel = deque()
        self.desSpeed = None
        self.rotateRate = None
        self.multiple_TR = None
        self.leftSpeed = 0
        self.rightSpeed = 0
        p.changeDynamics(self.pybullet_id, -1, mass=0.18)
        # No friction between bbody and surface.
        p.changeDynamics(self.pybullet_id, -1, lateralFriction=0, rollingFriction=0.)
        # Friction between joint links and surface.
        for i in range(p.getNumJoints(self.pybullet_id)):
            p.changeDynamics(self.pybullet_id, i, lateralFriction=1., rollingFriction=0.)

    def reset(self):
        p.resetBasePositionAndOrientation(self.pybullet_id, self.initial_position, self.initial_Orientation)
            
    def set_wheel_velocity(self, leftVelocity, rightVelocity):
        maxForce = 1000
        leftVelocity = leftVelocity / 1000   # swarmbang运动速度为毫米，在bulletEngine中需要转化为米进行计算
        rightVelocity = rightVelocity / 1000
        leftWheel_expRotateSpeed = round((leftVelocity / (2 * math.pi * self.wheelRadius)) * math.pi * 2, 2)
        rightWheel_expRotateSpeed = round((rightVelocity / (2 * math.pi * self.wheelRadius)) * math.pi * 2, 2)

        p.setJointMotorControl2(self.pybullet_id, self.leftWheelJointIdx, p.VELOCITY_CONTROL, targetVelocity=leftWheel_expRotateSpeed, force=maxForce)
        p.setJointMotorControl2(self.pybullet_id, self.rightWheelJointIdx, p.VELOCITY_CONTROL, targetVelocity=rightWheel_expRotateSpeed, force=maxForce)

    def get_pos_and_orientation(self):
        pos, rot = p.getBasePositionAndOrientation(self.pybullet_id)
        euler = p.getEulerFromQuaternion(rot)
        return np.array(pos), euler[2]

    def compute_controller(self, leftVelocity, rightVelocity):
        self.set_wheel_velocity(leftVelocity, rightVelocity)


    

    
       
