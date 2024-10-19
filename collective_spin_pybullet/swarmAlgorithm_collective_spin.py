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
