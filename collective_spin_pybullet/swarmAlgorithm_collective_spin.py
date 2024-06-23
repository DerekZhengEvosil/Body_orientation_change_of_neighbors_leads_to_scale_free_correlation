import copy
from scipy.spatial.distance import pdist, squareform
from get_ellipse_visual_area import get_ellipse_visual_area
from get_ellipse_visual_area import _cast_to_pm_pi, _subjpol_to_objcart, _line_intersect
from pylab import *
from copy import deepcopy


class swarmAlgorithm():
    def __init__(self, r_agent, align_type, select_type, tau, wall_agent_vel, cycTime,
                 fish_BL, w, dist_thresh, config):
        self.align_type = align_type
        self.dist_thresh = dist_thresh
        self.tau = tau
        self.select_type = select_type
        self.wall_agent_vel = wall_agent_vel
        self.fish_BL = fish_BL
        self.w = w
        self.r_agent = r_agent
        self.cycTime = cycTime
        self.v0 = 0
        self.v0_min = 0
        self.rigid_rep_fov = np.pi
        self.maxRotRate = 10 * 1.91
        self.maxRotRate_low = 10 * 1.91
        self.nearest_alignmates = None

    def calculateVision(self, robotsList, pos_BL, v):
        heading = np.arctan2(v[1, :], v[0, :])
        angular_area, adjacency_matrix, visual_field, polygon_param_cell, polygon_param_belonging_cell, tangent_pt_obj_cart, metric_distance_center = \
            get_ellipse_visual_area(self.w, len(robotsList), self.fish_BL / self.fish_BL, heading, pos_BL, l=1.,
                                    threshold=0.0, vis_thresh=0.0, dist_thresh=self.dist_thresh,
                                    non_overlap_angles=None)
        adjacency_dist_thresh = {}
        for robot_id in range(1, len(robotsList) + 1):
            belong_tmp = polygon_param_belonging_cell[robot_id]
            neighbor_id_list = [id_tmp for id_tmp in belong_tmp.values()]
            if any(neighbor_id_list):
                adjacency_dist_thresh[robot_id] = np.unique(neighbor_id_list)
            else:
                adjacency_dist_thresh[robot_id] = []

        return angular_area, adjacency_matrix, visual_field, polygon_param_cell, \
               polygon_param_belonging_cell, tangent_pt_obj_cart, metric_distance_center, adjacency_dist_thresh

    def calculateProjection(self, robotsList, angular_area, polygon_param_cell, polygon_param_belonging_cell):
        ellipse_vis_param = {}
        vis_angle = {}
        if np.any(angular_area):
            for r in robotsList:
                polygon_param = polygon_param_cell[r.id + 1]
                polygon_param_belonging = polygon_param_belonging_cell[r.id + 1]
                tangent_point_dis = np.zeros(len(robotsList))
                vis_angle_tmp = zeros(len(robotsList))
                belonging_id_tmp = {}
                polygon_info_tmp = {}
                for pp in range(len(polygon_param)):
                    pp = pp + 1  # index starts from 1
                    belonging_id = polygon_param_belonging[pp]
                    belonging_id_tmp[pp] = belonging_id
                    polygon_param_tmp = polygon_param[pp]
                    polygon_param_p1 = polygon_param_tmp[0]
                    polygon_param_p2 = polygon_param_tmp[1]
                    polygon_param_focal = polygon_param_tmp[2]
                    polygon_info_tmp[pp] = [polygon_param_p1, polygon_param_p2, polygon_param_focal]
                    focal_2_p1 = polygon_param_p1 - polygon_param_focal
                    focal_2_p1 = focal_2_p1 / np.linalg.norm(focal_2_p1)
                    focal_2_p2 = polygon_param_p2 - polygon_param_focal
                    focal_2_p2 = focal_2_p2 / np.linalg.norm(focal_2_p2)
                    tangent_point_dis[belonging_id] = tangent_point_dis[belonging_id] + \
                                                      np.linalg.norm(polygon_param_p1 - polygon_param_p2)
                    dot_temp = np.dot(focal_2_p1, focal_2_p2)
                    if dot_temp > 1.0:
                        dot_temp = 1.0
                    elif dot_temp < -1.0:
                        dot_temp = -1.0
                    vis_angle_tmp[belonging_id] = vis_angle_tmp[belonging_id] + math.acos(dot_temp)
                ellipse_vis_param[r.id] = tangent_point_dis
                vis_angle[r.id] = vis_angle_tmp
        else:
            ellipse_vis_param = {}
            vis_angle = {}
        return ellipse_vis_param, vis_angle

    def calculateAction(self, robotsList, pos, v, all_past_posvel, angular_area, past_vision, adjacency_dist_thresh,
                        informed_id, runSec):
        # the unit of position has already been changed to BL
        heading = np.arctan2(v[1, :], v[0, :])
        dist_xy = squareform(pdist(pos.T))
        align_mate = []
        old_heading = heading
        if self.align_type == "selective":
            if self.select_type == "BOC":
                for r in robotsList:
                    # neighbor = np.where(angular_area[:, r.id] > 0)[0]
                    # Nneig = neighbor.shape[0]
                    neighbor = adjacency_dist_thresh[r.id + 1]
                    if any(neighbor):
                        Nneig = neighbor.shape[0]
                        projection_rate = np.zeros(len(robotsList))
                        # accumulate the projection rate within the period
                        for qq in range(1, len(past_vision)):
                            projection_rate = projection_rate + (
                                    abs(past_vision[qq][r.id] - past_vision[qq - 1][r.id]) / self.cycTime)
                        cij = np.ones(Nneig)
                        iDir = v[:, r.id]
                        rij = pos[:, neighbor] - pos[:, r.id].reshape(2, 1)
                        dist_ij = np.linalg.norm(rij, axis=0)
                        nij = rij / dist_ij
                        front_preference = np.power(((1 + np.dot(iDir, nij)) / 2), 0)
                        cij = projection_rate[neighbor] * front_preference
                        salience = np.vstack((neighbor[np.argsort(cij)[::-1]], cij[np.argsort(cij)[::-1]])).T
                        salience[:, 1] = salience[:, 1] / sum(salience[:, 1])
                        align_mate.append(salience)
                    else:
                        align_mate.append(np.array([r.id]))

            elif self.select_type == "random":
                for r in robotsList:
                    neighbor = adjacency_dist_thresh[r.id + 1]
                    random_select = np.random.randint(0, len(neighbor))
                    align_mate.append(np.array([neighbor[random_select], 1]))

            elif self.select_type == "nearest":
                if runSec == 0:
                    for r in robotsList:
                        dist_tmp = dist_xy[r.id, :]
                        dist_tmp[dist_tmp == 0] = np.inf
                        nearest_neighbor = np.argmin(dist_tmp)
                        align_mate.append(np.array([nearest_neighbor, 1]))
                    self.nearest_alignmates = deepcopy(align_mate)
                else:
                    align_mate = self.nearest_alignmates

        desSpeed = self.v0 * np.ones([1, len(robotsList)])
        mate_list = []
        for r in robotsList:
            if r.id == informed_id:
                mate_list.append(r.id)
                continue
            dist_agent = dist_xy[r.id, :]
            dist_agent[r.id] = np.nan
            iDir = v[:, r.id]
            rij = pos[:, :] - pos[:, r.id].reshape(2, 1)
            dist_ij = np.linalg.norm(rij, axis=0)
            dist_ij[r.id] = nan
            nij = rij / dist_ij
            relative_heading = np.arccos(np.dot(iDir, nij))
            r.desSpeed = desSpeed[0, r.id]
            try:
                mate_id = align_mate[r.id][0, 0]
            except:
                mate_id = align_mate[r.id][0]
            mate_list.append(mate_id)
            r.leftSpeed = robotsList[int(mate_id)].leftSpeed
            r.rightSpeed = robotsList[int(mate_id)].rightSpeed
            #
            # # if r.leftSpeed == 0 and r.rightSpeed == 0:
            # #     r.leftSpeed = robotsList[int(mate_id)].leftSpeed
            # #     r.rightSpeed = robotsList[int(mate_id)].rightSpeed
            # # else:
            # #     r.leftSpeed = r.leftSpeed
            # #     r.rightSpeed = r.rightSpeed
        return robotsList, mate_list
