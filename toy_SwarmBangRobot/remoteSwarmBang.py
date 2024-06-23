import pybullet as p
import time
import pybullet_data
import math

cid = p.connect(p.SHARED_MEMORY)
if (cid < 0):
    p.connect(p.GUI)
p.setAdditionalSearchPath(pybullet_data.getDataPath())
p.resetSimulation()
p.setGravity(0, 0, -10)
useRealTimeSim = 1
p.setRealTimeSimulation(useRealTimeSim)  # either this
p.loadURDF("plane.urdf")
car = p.loadURDF("../models/swmrob_v26.urdf", [0, 0, .3])
for i in range(p.getNumJoints(car)):
    print(p.getJointInfo(car, i))
for wheel in range(p.getNumJoints(car)):
    p.setJointMotorControl2(car, wheel, p.VELOCITY_CONTROL, targetVelocity=0, force=0)
    p.getJointInfo(car, wheel)

p.changeDynamics(car, -1, mass=0.18)
p.changeDynamics(car, -1, rollingFriction=0)
p.changeDynamics(car, -1, lateralFriction=0)

p.changeDynamics(car, 1, rollingFriction=0)
p.changeDynamics(car, 1, lateralFriction=1)

p.changeDynamics(car, 2, rollingFriction=0)
p.changeDynamics(car, 2, lateralFriction=1)

wheels = [1, 2]
print("----------------")

wheelRadius = 0.0175  # m
secDis = 2 * math.pi * wheelRadius
linearVelocitybar = p.addUserDebugParameter("linearvelocity", -25, 25, 0)
angularVelocitybar = p.addUserDebugParameter("angularvelocity", -0.83, 0.83, 0)
while (True):
    linearVelocity = p.readUserDebugParameter(linearVelocitybar)
    angularVelocity = p.readUserDebugParameter(angularVelocitybar)
    maxForce = 1000
    leftVelocity = (linearVelocity - math.floor(
        math.degrees(angularVelocity) / 1.91)) / maxForce
    rightVelocity = (linearVelocity + math.floor(math.degrees(angularVelocity) / 1.91)) / maxForce
    leftWheel_expRotateSpeed = round((leftVelocity / (2 * math.pi * wheelRadius)) * math.pi * 2, 2)
    rightWheel_expRotateSpeed = round((rightVelocity / (2 * math.pi * wheelRadius)) * math.pi * 2, 2)
    p.setJointMotorControl2(car, 2, p.VELOCITY_CONTROL,
                            targetVelocity=leftWheel_expRotateSpeed, force=maxForce)
    p.setJointMotorControl2(car, 1, p.VELOCITY_CONTROL,
                            targetVelocity=rightWheel_expRotateSpeed, force=maxForce)
    vel = p.getBaseVelocity(car)
    dynInfo = p.getDynamicsInfo(car, -1)
    print("vx:{}-vy:{}-vz:{}".format(vel[0][0] * 1000, vel[0][1] * 1000, vel[0][2] * 1000))
    if (useRealTimeSim == 0):
        p.stepSimulation()
    time.sleep(0.01)
