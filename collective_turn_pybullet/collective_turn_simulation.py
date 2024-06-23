import pybullet as p
from robot import Robot
from matplotlib.patches import Ellipse, Polygon, FancyArrowPatch
from matplotlib.collections import PatchCollection
from sawrmAlgorithm_collective_turn import swarmAlgorithm
from get_ellipse_visual_area import _cast_to_pm_pi, _subjpol_to_objcart, _line_intersect
import cv2
import h5py
from OpenGL.GL import *
from pylab import *
import os


class World():
    def __init__(self, pybullet_tag, sim_index, time_tag, maxSimStep, expNum, align_type, select_type,
                 tau, turning_angle, plot_ini_setting, animation, dist_thresh, parameter_config,
                 activateTime):
        # create the physics simulator
        self.informed_dir = None
        self.informed_coords = None
        if pybullet_tag == 'DIRECT':
            self.physicsClient = p.connect(p.DIRECT)  # p.GUI DIRECT
        elif pybullet_tag == 'GUI':
            self.physicsClient = p.connect(p.GUI)  # p.GUI DIRECT
        p.setGravity(0, 0, -9.81)
        self.dt = 1. / 50.
        p.setPhysicsEngineParameter(self.dt, numSubSteps=1)
        p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
        self.pybullet_tag = pybullet_tag
        self.time = 0.0
        self.runSec = 0
        self.break_flag = False
        self.turning_angle = turning_angle
        self.real_wall = False
        self.r_agent = 0.042
        self.informed_id = None
        self.informed_id_list = []
        self.BL = 0.06
        self.w = 0.3
        self.activateTime = activateTime
        self.scaleDiv = 1000
        self.expNum = expNum
        self.cycTime = 1
        self.maxSimStep = maxSimStep
        self.tau = tau
        self.dist_thresh = dist_thresh
        self.vision_memory = []
        self.turning_direction = None
        self.expected_2agent_dist = 3 * 0.05 + 0.01
        self.informed_sin = None
        self.x_coords = []
        self.y_coords = []
        self.rotated_vector_1 = []
        self.rotated_vector_2 = []
        self.turn_cnt = 0
        self.first_dir = None
        self.swarmAlgorothm = swarmAlgorithm(self.r_agent, align_type, select_type, \
                                             tau, [], [], \
                                             self.cycTime, self.BL, self.w, self.dist_thresh, parameter_config)
        self.info_flag = 1
        self.animation = animation
        current_directory = os.path.dirname(__file__)
        randPos_dir = current_directory + "/randomPosition_N{}/".format(self.expNum)
        self.file_index = [1]
        init_position_filename = randPos_dir + "ranPos_{}".format(self.file_index[0]) + ".mat"
        all_pos_mat = h5py.File(init_position_filename, 'r')
        all_pos = np.array(all_pos_mat['pos'])
        all_pos = all_pos / 1000
        all_pos = all_pos.T
        # all_pos[0, :] = all_pos[0, :] - all_pos[0, -1] + self.ini_pos_x
        all_v = np.zeros((2, self.expNum))
        all_v[0, :] = 1
        plt.rcParams['toolbar'] = 'None'

        all_pos_BL = all_pos / self.BL
        # draw the wall in particle
        if plot_ini_setting:
            fig, ax = plt.subplots(figsize=(10, 4))
            radii = 1 * np.ones(self.expNum)
            patches1 = []
            for x1, y1, r in zip(all_pos_BL[0, :], all_pos_BL[1, :], radii):
                circle = Circle((x1, y1), r)
                patches1.append(circle)
            pp1 = PatchCollection(patches1, alpha=1)
            pp1.set_color("orange")
            ax.add_collection(pp1)  # 注意参数的赋值
            quiver(all_pos_BL[0, :], all_pos_BL[1, :], 1 * all_v[0, :], 1 * all_v[1, :], \
                   angles='xy', scale_units='xy', scale=1, color='black')
            # for i in range(self.expNum):
            #     text(all_pos_BL[0, i], all_pos_BL[1, i], str(i), horizontalalignment='center', verticalalignment='center', fontsize=8)
            # ax.set_xlim(-2, 6)
            # ax.set_ylim(-1, 1)
            axis('equal')
            plt.show()

        if self.animation:
            self.animation_fig = plt.figure(figsize=(6, 6.5), facecolor='w')
            self.animation_fig.patch.set_facecolor('white')
            self.animation_fig.patch.set_alpha(1)
            mngr = plt.get_current_fig_manager()
            # mngr.window.wm_geometry("+300+300")
            # plt.grid(True)
            plt.ion()
            plt.show()

        p.resetDebugVisualizerCamera(5.1, 0, -89.99, (0.8, -0.3, -2.27))
        # p.removeUserDebugItem(self.physicsClient)
        self.planeId = p.loadURDF("../models/plane.urdf")
        p.changeDynamics(self.planeId, -1, lateralFriction=5., rollingFriction=0)

        if align_type == 'selective':
            self.dirName = "U_turn_Data_selective/{}/{}-N={}-tau={}-distTH={}-freq={}/Sim_{}".format(time_tag,
                                                                                                           select_type,
                                                                                                           self.expNum,
                                                                                                           self.tau,
                                                                                                           dist_thresh,
                                                                                                           self.turning_angle,
                                                                                                           sim_index)
            self.res_dir = "U_turn_Data_selective/{}/{}-N={}-tau={}-freq={}/time_res".format(time_tag,
                                                                                                   select_type,
                                                                                                   self.expNum,
                                                                                                   self.tau,
                                                                                                   self.turning_angle)
        self.picDir = self.dirName + "/snapShot"
        if not os.path.exists(self.dirName):
            os.makedirs(self.dirName)
            os.makedirs(self.picDir)

        self.robots = []
        for idCnt in range(self.expNum):
            x_pos = all_pos[0][idCnt]
            y_pos = all_pos[1][idCnt]
            x_vel = all_v[0][idCnt]
            y_vel = all_v[1][idCnt]
            robotPos = [x_pos, y_pos, 0]
            robotHeading = [0, 0, np.pi / 2]
            self.robots.append(Robot(robotPos, robotHeading, idCnt, self.cycTime))
            p.stepSimulation()
            p.stepSimulation()

        self.stepSimulation()
        self.stepSimulation()

    def getFrameFromGame(self, w, h, curStep):
        data = []
        glReadBuffer(GL_FRONT)
        # 从缓冲区中的读出的数据是字节数组
        data = glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE)
        arr = np.zeros((h * w * 3), dtype=np.uint8)
        for i in range(0, len(data), 3):
            arr[i] = data[i + 2]
            arr[i + 1] = data[i + 1]
            arr[i + 2] = data[i]
        arr = np.reshape(arr, (h, w, 3))
        cv2.flip(arr, 0, arr)
        # cv2.imshow('scene', arr)
        picName = self.picDir + "/snapshot_{}.jpeg".format(curStep)
        cv2.imwrite(picName, arr)

    def unzip(self, items):
        return ([item[i] for item in items] for i in range(len(items[0])))

    def reset(self):
        """
        Resets the position of all the robots
        """
        for r in self.robots:
            r.reset()
        p.stepSimulation()

    def rotate_vector(self, vector, angle_degrees):
        angle_radians = np.radians(angle_degrees)
        rotation_matrix = np.array([[np.cos(angle_radians), -np.sin(angle_radians)],
                                    [np.sin(angle_radians), np.cos(angle_radians)]])
        rotated_vector = np.dot(rotation_matrix, vector)
        return rotated_vector

    def stepSimulation(self):
        """
        Simulates one step simulation
        """
        if self.time > self.cycTime:
            # Initialize the pos and velocity
            pos = np.zeros([2, len(self.robots)])
            v = np.zeros([2, len(self.robots)])
            all_past_posvel = np.zeros([4, len(self.robots)])
            # get position and velocity from the bullet engine under the unit of m
            # position and velocity memory
            for r in self.robots:
                posTmp, oriTmp = r.get_pos_and_orientation()
                pos[0, r.id] = posTmp[0]
                pos[1, r.id] = posTmp[1]
                v[0, r.id] = math.sin(oriTmp)
                v[1, r.id] = -1 * math.cos(oriTmp)
                robotInfo = "{},{},{}".format(posTmp[0], posTmp[1], oriTmp)
                savefile = open(self.dirName + "/simData_{}.txt".format(r.id), "a")
                savefile.write(robotInfo + '\n')
                if self.runSec >= self.tau:
                    r.past_posvel.popleft()
                    r.past_posvel.append(np.append(pos[:, r.id], v[:, r.id], axis=0))
                else:
                    r.past_posvel.append(np.append(pos[:, r.id], v[:, r.id], axis=0))
                all_past_posvel[:, r.id] = np.array(r.past_posvel[0])

            pos_BL = pos / self.BL
            angular_area, adjacency_matrix, visual_field, polygon_param_cell, polygon_param_belonging_cell, \
            tangent_pt_obj_cart, metric_distance_center, adjacency_dist_thresh = \
                self.swarmAlgorothm.calculateVision(self.robots, pos_BL, v)

            ellipse_vis_param, vis_angle = self.swarmAlgorothm.calculateProjection(self.robots, angular_area,
                                                                                   polygon_param_cell,
                                                                                   polygon_param_belonging_cell)
            self.vision_memory.append(ellipse_vis_param)
            if self.runSec >= self.tau:
                past_vision = self.vision_memory[self.runSec - self.tau:self.runSec]
            else:
                past_vision = self.vision_memory[0:self.runSec]

            # update the controllers
            self.robots = self.swarmAlgorothm.calculateAction(self.robots, pos_BL, v, all_past_posvel, angular_area,
                                                              past_vision, adjacency_dist_thresh, self.informed_id)

            # update the informed robot
            # if self.runSec % self.activateTime == 0 and self.runSec > 0 and self.info_flag:
            if self.runSec == self.activateTime:
                self.turn_cnt = self.turn_cnt + 1
                self.info_flag = 0
                group_velocity = np.mean(v, axis=1)
                self.rotated_vector_1 = self.rotate_vector(group_velocity, -self.turning_angle)
                self.rotated_vector_2 = self.rotate_vector(group_velocity, self.turning_angle)
                projections = np.dot(pos_BL.T, group_velocity) / np.linalg.norm(group_velocity)
                self.informed_id = np.argmax(projections)  # sorted_indices[len(sorted_indices) // 2]  #
                self.informed_id_list.append(self.informed_id)
                self.informed_dir = self.rotated_vector_1

            if self.runSec > self.activateTime:
                # the rule for the informed individual
                informed_pos, informed_ori = self.robots[self.informed_id].get_pos_and_orientation()
                informed_vel = [math.sin(informed_ori), -1 * math.cos(informed_ori)]
                self.turning_direction = self.informed_dir
                # print("informed_dir:{}".format(self.informed_dir))
                informed_rotateDir = np.arcsin(
                    informed_vel[0] * self.turning_direction[1] - informed_vel[1] * self.turning_direction[
                        0])  # x1y2 – x2y1
                informed_dot_temp = informed_vel[0] * self.turning_direction[0] + informed_vel[1] * \
                                    self.turning_direction[1]
                if informed_dot_temp > 1.0:
                    informed_dot_temp = 1.0
                elif informed_dot_temp < -1.0:
                    informed_dot_temp = -1.0
                informed_rotateAnlge = np.arccos(informed_dot_temp) * 180 / np.pi
                informed_rotRate = informed_rotateAnlge / self.cycTime
                multiple_TR = 0
                informed_maxRotRate = 25 * 1.91
                if informed_rotRate > informed_maxRotRate:
                    rotRate = informed_maxRotRate
                    multiple_TR = math.floor(rotRate / 1.91)
                else:
                    rotRate = informed_rotRate
                    multiple_TR = math.floor(rotRate / 1.91)

                if informed_rotateDir > 0:
                    multiple_TR = multiple_TR
                else:
                    multiple_TR = -1 * multiple_TR
                informed_linear_vel = 20
                self.robots[self.informed_id].leftSpeed = informed_linear_vel + multiple_TR
                self.robots[self.informed_id].rightSpeed = informed_linear_vel - multiple_TR

            for r in self.robots:
                r.compute_controller(r.leftSpeed, r.rightSpeed)

            self.runSec = self.runSec + 1
            self.time = 0.0
            if self.runSec == self.maxSimStep:
                savefile = open(self.dirName + "/informed_id.txt", "a")
                id_tmp = "{}".format(self.informed_id)
                savefile.write(id_tmp + '\n')
                savefile = open(self.dirName + "/informed_{}_velocity.txt".format(self.informed_id), "a")
                idx_record = open(self.dirName + "/idx_record.txt", "a")
                idx_record.write("{}".format(self.file_index[0]) + '\n')
                vel_tmp = "{},{}".format(self.informed_dir[0], self.informed_dir[1])
                savefile.write(vel_tmp + '\n')

            if self.animation:
                phi = np.arctan2(v[1, :], v[0, :])
                plt.clf()
                ax = self.animation_fig.add_subplot(111)
                alpha = 0.5
                edgecolor = 'none'
                edgewidth = 1
                ellipse_color = "r"
                zorder = 100
                ellipses_list = []
                for i in range(self.expNum):
                    if self.runSec > self.activateTime and i == self.informed_id:
                        ellipses = Ellipse(pos_BL[:, i], 1.0, self.w,
                                           _cast_to_pm_pi(phi[i]) * 180.0 / np.pi)
                        ax.add_artist(ellipses)
                        ellipses.set_clip_box(ax.bbox)
                        informed_ellipse_color = "b"
                        ellipses.set_facecolor(informed_ellipse_color)
                        ellipses.set_alpha(alpha)
                        ellipses.set_edgecolor(edgecolor)
                        ellipses.set_linewidth(edgewidth)
                        ellipses.set_zorder(zorder)
                        ellipses_list.append(ellipses)
                    else:
                        ellipses = Ellipse(pos_BL[:, i], 1.0, self.w,
                                           _cast_to_pm_pi(phi[i]) * 180.0 / np.pi)
                        ax.add_artist(ellipses)
                        ellipses.set_clip_box(ax.bbox)
                        ellipses.set_facecolor(ellipse_color)
                        ellipses.set_alpha(alpha)
                        ellipses.set_edgecolor(edgecolor)
                        ellipses.set_linewidth(edgewidth)
                        ellipses.set_zorder(zorder)
                        ellipses_list.append(ellipses)

                viewer_id = 0
                vf_color = 'darkseagreen'
                polygon_param_cell_tmp = polygon_param_cell[viewer_id + 1]
                polygon_param_belonging_cell_tmp = polygon_param_belonging_cell[viewer_id + 1]
                for k in range(len(polygon_param_cell_tmp)):
                    polygon_tmp = polygon_param_cell_tmp[k + 1]
                    p1 = polygon_tmp[0]
                    p2 = polygon_tmp[1]
                    focal = polygon_tmp[2]
                    visual_area = plt.Polygon([p1, p2, focal])
                    ax.add_artist(visual_area)
                    visual_area.set_facecolor(vf_color)
                    visual_area.set_alpha(alpha)
                    ellipses_list.append(visual_area)
                # ax.set_xlim(-2, 6)
                # ax.set_ylim(-1, 1)
                quiver(pos_BL[0, :], pos_BL[1, :], 1 * v[0, :], 1 * v[1, :], \
                       angles='xy', scale_units='xy', scale=2, color='black')  # 注意参数的赋值
                axis('equal')
                # plt.plot(self.x_coords, self.y_coords)
                # plt.axis('off')
                title("step={}".format(self.runSec))
                figName = self.picDir + "/traj_{}.pdf".format(self.runSec)
                self.animation_fig.savefig(figName, dpi=600, format='pdf')
                plt.pause(0.001)
        # do one simulation step
        p.stepSimulation()
        self.time += self.dt
        return self.runSec
