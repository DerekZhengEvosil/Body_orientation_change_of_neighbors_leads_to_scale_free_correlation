import time
import configparser
from collective_spin_sim import World
plot_ini_setting = 0
animation = 1
pybullet_tag = 'DIRECT'  # GUI/DIRECT GUI for showing simulation in pybullet environment
align_type = 'selective'
select_type = 'BOC'  # BOC random
sim_index = 1
maxSimStep = 60 #  time that the spin initiator just spun 2pi angle.
activateTime = 25
expNum = 100
tau = 5
dist_thresh = 10
file_index = 1
parameter_config = configparser.ConfigParser()
parameter_version = 1
config_name = "Inner_parameter_v{}.ini".format(parameter_version)
parameter_config.read(config_name)
time_tag = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
world = World(pybullet_tag, sim_index, time_tag, maxSimStep, expNum, align_type, select_type, tau, plot_ini_setting, animation, dist_thresh, parameter_config, activateTime)
start_time = time.time()
while True:
    runSec = world.stepSimulation()
    if runSec > maxSimStep:
        break

end_time = time.time()
execution_time = end_time - start_time
print("step：{} Elapsed {} s".format(runSec, execution_time))