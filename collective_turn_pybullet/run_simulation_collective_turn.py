import time
import configparser
from collective_turn_simulation import World
plot_ini_setting = 0
animation = 1
pybullet_tag = 'DIRECT'  # GUI DIRECT, GUI is for showing the simulation in pybullet enviroment
align_type = 'selective'
select_type = 'BOC'   # BOC random
tau = 10
turning_angle = 180 # 90 degree
expNum = 100
sim_index = 1
maxSimStep = 200
dist_thresh = 10
activateTime = 50
parameter_config = configparser.ConfigParser()
parameter_version = 1
config_name = "Inner_parameter_v{}.ini".format(parameter_version)
parameter_config.read(config_name)
time_tag = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
world = World(pybullet_tag, sim_index, time_tag, maxSimStep, expNum, align_type, select_type, tau,
                    turning_angle, plot_ini_setting, animation, dist_thresh, parameter_config, activateTime)
start_time = time.time()
while True:
    runSec = world.stepSimulation()
    # print(runSec)
    if runSec > maxSimStep:
        break

end_time = time.time()
execution_time = end_time - start_time
print("step：{} Elapsed {} s".format(runSec, execution_time))